def Model(L1,L2,l,emur,ModelTuyau):
    if ModelTuyau==1:
        return Model1(L1,L2,l,emur)
    if ModelTuyau==2:
        return Model2(L1,L2,l,emur)

def Model1(L1,L2,l,emur):   #Serpentin
    x=l+emur
    y=-1
    Tuyau=[]
    for i in range (L2-l-emur):
        y+=1
        Tuyau.append((x,y))
    while x<L1-3*l-emur:
        for i in range (l):
            x+=1
            Tuyau.append((x,y))
        for i in range (L2-3*l-emur):
            y-=1
            Tuyau.append((x,y))
        for i in range (l):
            x+=1
            Tuyau.append((x,y))
        for i in range (L2-3*l-emur):
            y+=1
            Tuyau.append((x,y))
    while x<L1-l-emur:
        x+=1
        Tuyau.append((x,y))
    for i in range (L2-2*l-emur):
        y-=1
        Tuyau.append((x,y))
    for i in range (L1-3*l-2*emur):
        x-=1
        Tuyau.append((x,y))
    while y>0:
        y-=1
        Tuyau.append((x,y))
    return Tuyau

def nearturn(x,y,dep,l,T):
    for i in range (3*l):
        x,y = av(x,y,dep)
        if [x,y] in T:
            return i+1
    return False

def av(x,y,dep):
    x+=((dep==0)-(dep==2))
    y+=((dep==3)-(dep==1))
    return x,y

def Model2(L1,L2,l,emur):   #Escargot
    x=l+emur
    y=-1
    sens=1
    dep=3
    Tuyau=[]
    for i in range (L2-emur-l):
        y+=1
        Tuyau.append([x,y])
    for i in range (L1-2*(l+emur)):
        x+=1
        Tuyau.append([x,y])
    for i in range (L2-emur-2*l):
        y-=1
        Tuyau.append([x,y])
    for i in range (L1-2*(l+emur)-2*l):
        x-=1
        Tuyau.append([x,y])
    
    while y:
        while y and ([x+l*(1+(sens==1))*((dep==0)-(dep==2)),y+l*(1+(sens==1))*((dep==3)-(dep==1))] not in Tuyau) :
            x,y = av(x,y,dep)
            Tuyau.append([x,y])
        dep = (dep+sens)%4
        dist = nearturn(x,y,dep,l,Tuyau)
        if dist :
            '''demi-tour'''
            for i in range (dist//3):
                x,y = av(x,y,dep)
                Tuyau.append([x,y])
            dep = (dep+sens)%4
            while ([x+l*((dep==0)-(dep==2)),y+l*((dep==3)-(dep==1))] not in Tuyau):
                x,y = av(x,y,dep)
                Tuyau.append([x,y])
            sens = -1
            dep = (dep+sens)%4
            for i in range (dist//3):
                x,y = av(x,y,dep)
                Tuyau.append([x,y])
            dep = (dep+sens)%4
            while ([x+l*((dep==0)-(dep==2)),y+l*((dep==3)-(dep==1))] not in Tuyau):
                x,y = av(x,y,dep)
                Tuyau.append([x,y])
            dep = (dep-1)%4
            while dep!=1:
                while ([x+l*((dep==0)-(dep==2)),y+l*((dep==3)-(dep==1))] not in Tuyau):
                    Tuyau.append([x,y])
                    x,y = av(x,y,dep)
                dep = (dep+sens)%4
            sens=-1
    return Tuyau
