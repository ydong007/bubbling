import random
import copy
import numpy as np
import pandas as pd
import math
import energy
import tools
from param import sigma

def grow(membrane,system,log,th0):
    
    chosen_angle = chooseangle(membrane)
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
    print(chosen_angle, len(membrane.triangle))
    newmembrane = nextmove(membrane, chosen_angle, system, log, th0)
    return newmembrane

def chooseangle(membrane):

    open_angle = openangle(membrane)
    mintheta = 2 * np.pi
    for i, row in open_angle.iterrows():
        mintheta = min(mintheta, open_angle['theta'].loc[i])
    angleweight = np.exp(-.5 * ((mintheta - open_angle['theta']) / sigma ) ** 2)
    chosen_angle=open_angle.sample(weights=angleweight, random_state=1) #random.seed(9001)
        
    return chosen_angle

def openangle(membrane):

    lvectors=tools.linevec(membrane)
    nhat=tools.normalvec(membrane)

    open_angle = pd.DataFrame(columns=['oindex', 'theta', 'angle_type' ]) #angle_type:0 is vertex angle and 1 is for line angle
    for m in membrane.vertex[membrane.vertex.edge==1].index:
        if (membrane.vertex.loc[m,'pnt']>=4): 
            for i in membrane.line[membrane.line.lnt==1].index:
                if (membrane.line.loc[i,'lv1']==m):
                    c=i
                elif (membrane.line.loc[i,'lv2']==m):
                    d=i
            angle=tools.vertex_angle(d,c,0,m,membrane,lvectors,nhat)
            if (membrane.vertex.loc[m,'pnt']!=4):
                n_o=len(open_angle)
                open_angle.loc[n_o]=[m,angle,0]
            elif(angle<= (np.pi/3)):
                n_o=len(open_angle)
                open_angle.loc[n_o]=[m,angle,0]

    for l in membrane.line[membrane.line.lnt==1].index:
        j=membrane.line['ll'].loc[l]
        k=membrane.line['lr'].loc[l]
        lv=membrane.line['lv2'].loc[l]
        rv=membrane.line['lv1'].loc[l]
        if(membrane.vertex['pnt'].loc[lv]<=4 and membrane.vertex['pnt'].loc[rv]<=4):
            theta1 = tools.vertex_angle(l, j, 0, lv, membrane, lvectors, nhat)
            theta2 = tools.vertex_angle(l, k, 1, rv, membrane, lvectors, nhat)
            angle = min(theta1, theta2)
            n_o=len(open_angle)
            open_angle.loc[n_o]=[l,angle,1]

    return open_angle

def nextmove(membrane, chosen_angle, system, log,th0):

    print("*****angle type*****")
    print(chosen_angle['angle_type'].iloc[0])
    print(chosen_angle['oindex'].iloc[0])
    print(chosen_angle['theta'].iloc[0])
    
    if(chosen_angle['angle_type'].iloc[0]==0):
        v = chosen_angle['oindex'].iloc[0]
        if (membrane.vertex['pnt'].loc[v]==5):
                newmembrane = insert_tri(v,membrane)
                return newmembrane
        
    else:
        l = chosen_angle['oindex'].iloc[0]
        newmembrane = nextvertex(membrane,l,th0)
        return newmembrane

def insert_tri(v,membrane):

    for i in membrane.line[membrane.line.lnt==1].index:
        if (membrane.line.loc[i,'lv1']==v):
            c=i
        elif (membrane.line.loc[i,'lv2']==v):
            d=i
        else:
            print("check the next line")

    # d=dd # so stupid! but I had to put this line otherwise it gives error
    a=membrane.line.loc[d,'lv1']
    b=membrane.line.loc[c,'lv2']
    rline=membrane.line.loc[d,'lr']
    lline=membrane.line.loc[c,'ll']
    n_l=len(membrane.line)
    n_t=len(membrane.triangle)
    membrane.line.loc[1+n_l]=[a,b,1+n_t,1000,1,rline,lline,0]
    membrane.line.loc[d,'lnt']=2
    membrane.line.loc[c,'lnt']=2
    membrane.line.loc[d,'lt2']=1+n_t
    membrane.line.loc[c,'lt2']=1+n_t
    membrane.vertex.loc[v,'pnt']=1+membrane.vertex.loc[v,'pnt']
    membrane.vertex.loc[v,'edge']=0
    membrane.vertex.loc[a,'pnt']=1+membrane.vertex.loc[a,'pnt']
    membrane.vertex.loc[b,'pnt']=1+membrane.vertex.loc[b,'pnt']
    membrane.line.loc[rline,'ll']=1+n_l
    membrane.line.loc[lline,'lr']=1+n_l
    membrane.triangle.loc[1+n_t]=[v,a,b,3+n_l,1+n_l,2+n_l,0]
    membrane.line.loc[2+n_l]=[b,membrane.line.loc[c,'lv1'],1+n_t,membrane.line.loc[c,'lt1'],2,rline,lline,c]
    membrane.line.loc[3+n_l]=[membrane.line.loc[d,'lv2'],a,1+n_t,membrane.line.loc[d,'lt1'],2,rline,lline,d]

        
    print("insert tri")
    return membrane

def nextvertex(membrane,l,th0):

    lvectors=tools.linevec(membrane)
    nhat=tools.normalvec(membrane)
    tn=membrane.line['lt1'].loc[l]
    a=membrane.line['lv1'].loc[l]
    b=membrane.line['lv2'].loc[l]
    [c,d]=np.array(membrane.vertex.loc[[a,b],list('xyz')])
    llp=abs(1-(.25*(lvectors['mag'].loc[l]**2)))
    ll=math.sqrt(llp)
    lv=(np.array(lvectors.loc[[l],list('xyz')]))/lvectors['mag'].loc[l]
    nh=nhat.loc[tn].as_matrix(columns=None)
    nh=nh[:3]
    bvec=np.cross(lv,nh)
    new=((c+d)*0.5)+ll*(np.cos(th0)*bvec-np.sin(th0)*nh)
    lt2=membrane.line['lt1'].loc[l]

    [x,y,z]=new[0,:]
    n_v=len(membrane.vertex)
    membrane.vertex.loc[1+n_v]= [x,y,z,1,1]

    n_v=len(membrane.vertex)
    n_t=len(membrane.triangle)
    n_l=len(membrane.line)
    rline=membrane.line['lr'].loc[l]
    lline=membrane.line['ll'].loc[l]
    membrane.line['ll'].loc[rline]=n_l+1
    membrane.line['lr'].loc[lline]=n_l+2
    membrane.line['lt2'].loc[l]=1+n_t
    membrane.line['lnt'].loc[l]=2
    membrane.vertex['pnt'].loc[a]=1+membrane.vertex['pnt'].loc[a]
    membrane.vertex['pnt'].loc[b]=1+membrane.vertex['pnt'].loc[b]
    membrane.line.loc[1+n_l]=[a,n_v,1+n_t,1000,1,rline,2+n_l,0]
    membrane.line.loc[2+n_l]=[n_v,b,1+n_t,1000,1,1+n_l,lline,0]
    membrane.line.loc[3+n_l]=[b,a,1+n_t,lt2,2,2+n_l,1+n_l,l]

    membrane.triangle.loc[1+n_t]=[a,n_v,b,3+n_l,1+n_l,2+n_l,0]
    
    # print(membrane.vertex)
    # print(membrane.line)
    # print(membrane.triangle)

    return membrane


def backup(membarne):

    backupmembrane=copy.deepcopy(membarne)
    print("backup")
    return backupmembrane

def restore(membrane,backupmembrane):
    
    membrane=copy.deepcopy(backupmembrane)
    print("restore")
    return membrane


