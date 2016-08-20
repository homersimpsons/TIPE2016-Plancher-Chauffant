import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import pickle as pick
import copy as c

def Simulation(chemin):     #Chemin de sauvegarde
    '''Constantes'''
    rho_fluide = 1034       #Masse volumique du fluide (kg/m^3)
    C_fluide = 3600         #Capacité thermique du fluide (J/kg.K)
    lambda_tuyau = 0.46     #Conductvié thermique du tuyau (W/m.K)
    lambda_beton = 0.92     #Conductvié thermique du béton (W/m.K)
    h = 20                  #Coefficient de conducto-convection sol-air (W/k.m^2)
    
    try:
        NewPlancher,Params = read(chemin)
        emur,L1,L2,Z,z,l,v,R,ModelTuyau,Theta,Tair,Text,ThetaTuyau,epsilon,dl,ScaleDiff,TypeTuyau = Params
        Tuyau = Model(L1,L2,l+R,emur,ModelTuyau)
        print('z=',z*dl,' l=',l*dl*ScaleDiff,' v=',v,' Thetatuyau=',ThetaTuyau,' TypeTuyau=',TypeTuyau,' ModelTuyau=',ModelTuyau,' LongueurTuyau=',(len(Tuyau)-1)*dl*ScaleDiff)
        ExtendTuyau(Tuyau,R,z,ScaleDiff)
        
    except:
        '''Paramètres'''
        emur = 0.2      #Épaisseur du mur (m)
        L1 = 5+2*emur   #Longueur (m)
        L2 = 3+emur     #Largeur (m)
        Z = 0.08        #Hauteur (m)
        
        z = 0.02        #Position du tuyau (m)
        l = 0.1         #Espacement entre les tuyaux (m)
        v = 0.4         #Vitesse du fluide dans le tuyau (m/s)
        R = 0.008       #Rayon Tuyau (m)
        ModelTuyau = 1  #Disposition des tuyaux
        
        Theta = 283     #Température du plancher au début (°K)
        Tair = 291      #Température de l'air (°K)
        Text = 278      #Température extérieure (°K)
        ThetaTuyau = 318#Température du Tuyau (°K)
        
        epsilon = 0.01  #condition de fin
        dl = 0.004      #Pas d'intégration (m)
        ScaleDiff = 3   #Facteur d'échelle entre les axes
        TypeTuyau = 2   #Type de refroidissement du tuyau
        
        '''Calculs'''
        Dm = v*np.pi*(R**2)*rho_fluide                      #(kg/s)
        alpha = (lambda_tuyau*2*np.pi*R)/(Dm*C_fluide*2*R)  #(1/m)
        
        '''Discrétisation'''
        L1 = int(L1/(ScaleDiff*dl))
        L2 = int(L2/(ScaleDiff*dl))
        Z = int(Z/dl)
        emur = int(emur/(dl*ScaleDiff))
        l = int(l/(ScaleDiff*dl))
        R = max(int(R/dl),1)
        z = int(z/dl)
        
        Params = [emur,L1,L2,Z,z,l,v,R,ModelTuyau,Theta,Tair,Text,ThetaTuyau,epsilon,dl,ScaleDiff,TypeTuyau]
        
        Tuyau = Model(L1,L2,l+2*R,emur,ModelTuyau)
        Tuyau = ExtendTuyau(Tuyau,R,z,ScaleDiff)
        
        NewPlancher = TuyauLim(Limites(Creation(L1,L2,Z,Theta),Tair,Text,h,lambda_beton,dl,ScaleDiff),Tuyau,ThetaTuyau,alpha,ScaleDiff,dl,TypeTuyau)
        i = 0
        while i==0 or Norme(Plancher,NewPlancher,epsilon):  #Simulation
            i += 1
            Plancher = c.deepcopy(NewPlancher)
            NewPlancher = Laplace(NewPlancher,Plancher,ScaleDiff**2)
            NewPlancher = TuyauLim(Limites(NewPlancher,Tair,Text,h,lambda_beton,dl,ScaleDiff),Tuyau,ThetaTuyau,alpha,ScaleDiff,dl,TypeTuyau)
        save(chemin,NewPlancher,Params)
    Affichage(NewPlancher,dl,ScaleDiff)
    return (Flux(NewPlancher,dl,ScaleDiff,lambda_beton,emur),FluxTuyau(NewPlancher,Tuyau,v,R,dl,C_fluide,rho_fluide))

'''     Initialisation  '''
    
def Creation(L1,L2,Z,Theta):    #Créer la matrice des points du plancher
    return [[[Theta for i in range (L1)] for j in range (L2)] for k in range (Z)]

def ExtendTuyau(Tuyau,R,z,ScaleDiff):   #Créer les portions de tuyau
    INT = [[0,0]]
    EXT = [[0,R],[0,-R],[max(R//ScaleDiff,1),0],[-max(R//ScaleDiff,1),0]]
    for j in range(1,R):
        INT.append([0,j])
        INT.append([0,-j])
        maxi = max(((R**2-j**2)**0.5)//ScaleDiff,1)
        for i in range(1,maxi):
            INT.append([i,j])
            INT.append([-i,j])
            INT.append([i,-j])
            INT.append([-i,-j])
        EXT.append([j,maxi])
        EXT.append([-j,maxi])
        EXT.append([j,-maxi])
        EXT.append([-j,maxi])
    Tuyau[0] = [ExtP(Tuyau[0],z,INT,True),ExtP(Tuyau[0],z,EXT,True)]
    for i in range (1,len(Tuyau)):
        direction = (Tuyau[i][0]==Tuyau[i-1][0][0][2])
        Tuyau[i] = [ExtP(Tuyau[i],z,INT,direction),ExtP(Tuyau[i],z,EXT,direction)]
    return INTEXT(Tuyau)

def ExtP(point,z,SET,direction):
    if direction:
        return [[x[1]+z,point[1],x[0]+point[0]] for x in SET]
    return [[x[1]+z,x[0]+point[1],point[0]] for x in SET]

def INTEXT(Tuyau):  #Retire les points extérieurs qui sont aussi à l'intérieur
    INT = [Tuyau[i][0] for i in range (len(Tuyau))]
    for i in range (len(Tuyau)):
        Tuyau[i][1] = [ext for ext in Tuyau[i][1] if ext not in INT]
    return Tuyau

'''     Calculs         '''

def Laplace(NewPlancher,Plancher,fact):     #Applique l'équation de la chaleur
    for i in range (1,len(Plancher[0][0])-1):
        for j in range (1,len(Plancher[0])-1):
            for k in range (1,len(Plancher)-1):
                NewPlancher[k][j][i] = (Plancher[k][j][i-1]+Plancher[k][j][i+1]+Plancher[k][j-1][i]+Plancher[k][j+1][i]+fact*(Plancher[k-1][j][i]+Plancher[k+1][j][i]))/(fact*2+4)
    return NewPlancher

'''     Limites         '''

def Limites(Plancher,Tair,Text,h,lamb,dl,ScaleDiff):    #Limites sur les bords
    for k in range (len(Plancher)):
        for i in range (len(Plancher[0][0])):
            Plancher[k][0][i] = Tair
            Plancher[k][-1][i] = Text
        for j in range (len(Plancher[0])):
            Plancher[k][j][0] = Text
            Plancher[k][j][-1] = Text
    for i in range (len(Plancher[0][0])):
        for j in range (len(Plancher[0])):
            Plancher[-1][j][i] = (lamb*Plancher[-2][j][i]+h*dl*Tair)/(lamb+h*dl)
            Plancher[0][j][i] = Plancher[1][j][i]
    return Plancher

def TuyauLim(Plancher,Tuyau,ThetaTuyau,alpha,ScaleDiff,dl,Type):    #Limite du tuyau
    precision=5     #precision pour la méthode d'Euler
    if Type==1:     #Tuyau à température constante
        TEMP = [ThetaTuyau for i in range (len(Tuyau))]
    else:           #Résolution avec Euler
        TMoyP = []
        for i in range (len(Tuyau)):
            Tmoy = 0
            for x in Tuyau[i][1]:   #Température moyenne autour de chaque portion
                Tmoy += Plancher[x[0]][x[1]][x[2]]
            TMoyP.append(Tmoy/len(Tuyau[i][1]))
        ftemp = interp.interp1d([i*ScaleDiff*dl for i in range (len(Tuyau))] , TMoyP)
        def deriv(x,T,alpha):
            return alpha*(ftemp(x)-T)
        TEMP = Euler(deriv,ThetaTuyau,[i*ScaleDiff*dl/precision for i in range ((len(Tuyau)-1)*precision+1)],alpha)[::precision]
    for i in range (len(Tuyau)):
        for x in Tuyau[i][0]:
            Plancher[x[0]][x[1]][x[2]] = TEMP[i]
    return Plancher

'''     Résultats       '''

def Flux(Plancher,dl,ScaleDiff,lamb,emur):  #Flux sortant du plancher
    flux_haut = 0
    flux_total = 0
    for i in range (emur,len(Plancher[0][0])-emur):
        for j in range (1,len(Plancher[0])-emur):
            flux_haut += (Plancher[-2][j][i]-Plancher[-1][j][i])*dl
            flux_total += (Plancher[1][j][i]-Plancher[0][j][i])*dl
    flux_total += flux_haut
    for j in range (len(Plancher[0])-emur,len(Plancher[0])-1):
        for i in range (1,emur):
            flux_total += (Plancher[-2][j][i]-Plancher[-1][j][i])*dl
            flux_total += (Plancher[1][j][i]-Plancher[0][j][i])*dl
        for i in range (len(Plancher[0][0])-emur,len(Plancher[0][0])-1):
            flux_total += (Plancher[-2][j][i]-Plancher[-1][j][i])*dl
            flux_total += (Plancher[1][j][i]-Plancher[0][j][i])*dl
    for k in range (1,len(Plancher)-1):
        for i in range (1,len(Plancher[0][0])-1):
            flux_total += (Plancher[k][-2][i]-Plancher[k][-1][i])*dl*ScaleDiff
            flux_total += (Plancher[k][1][i]-Plancher[k][0][i])*dl*ScaleDiff
        for j in range (1,len(Plancher[0])-1):
            flux_total += (Plancher[k][j][-2]-Plancher[k][j][-1])*dl*ScaleDiff
            flux_total += (Plancher[k][j][1]-Plancher[k][j][0])*dl*ScaleDiff
    return (lamb*flux_haut,lamb*flux_total)

def FluxTuyau(Plancher,Tuyau,v,R,dl,C_fluide,rho_fluide):   #Flux sortant du tuyau
    x = Tuyau[0][0][0]
    y = Tuyau[-1][0][0]
    TMoyP = [Plancher[x[0]][x[1]][x[2]],Plancher[y[0]][y[1]][y[2]]]
    return v*np.pi*((R*dl)**2)*C_fluide*rho_fluide*(TMoyP[0]-TMoyP[-1])

def Affichage(Plancher,dl,ScaleDiff):
    fig = plt.figure()
    im = plt.imshow(Plancher[-1],cmap='jet',interpolation='nearest')
    fig.colorbar(im)
    plt.show()

'''     Outils          '''

def Norme(P1,P2,value=-1):  #Norme infinie (distance)
    lim = not(value==-1)
    max = 0
    for i in range (len(P1[0][0])):
        for j in range (len(P1[0])):
            for k in range (len(P1)):
                if abs(P1[k][j][i]-P2[k][j][i])>max:
                    if lim and abs(P1[k][j][i]-P2[k][j][i])>value:
                        return True
                    max = abs(P1[k][j][i]-P2[k][j][i])
    if lim:
        return False
    return max

def Euler(deriv,T0,X,alpha):
    TEMP = [T0]
    T = T0
    for i in range(1,len(X)) :
        T = T+(X[i]-X[i-1])*deriv(X[i-1],T,alpha)
        TEMP.append(T)
    return TEMP

def save(chemin,Plancher,Params):   #Sauvegarde
    file = open(chemin+'Plancher.bin', 'wb')
    pick.dump(Plancher,file)
    file.close()
    file = open(chemin+'Params.bin', 'wb')
    pick.dump(Params,file)
    file.close()
    return

def read(chemin):   #Restauration
    file = open(chemin+'Plancher.bin', 'rb')
    Plancher = pick.load(file)
    file.close()
    file = open(chemin+'Params.bin', 'rb')
    Params = pick.load(file)
    file.close()
    return Plancher,Params