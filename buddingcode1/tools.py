import pandas as pd
import numpy as np
from numpy import linalg as LA
import math


def linevec(membrane):
	lvectors = pd.DataFrame(columns=['x', 'y', 'z', 'mag'])

	for i, row in membrane.line.iterrows():
		w=membrane.line['lv1'].loc[i]
		o=membrane.line['lv2'].loc[i]
		vx=membrane.vertex['x'].loc[o]-membrane.vertex['x'].loc[w]
		vy=membrane.vertex['y'].loc[o]-membrane.vertex['y'].loc[w]
		vz=membrane.vertex['z'].loc[o]-membrane.vertex['z'].loc[w]
		a=[vx,vy,vz]
		mag=LA.norm(a)
		lvectors.loc[i] = [vx, vy, vz, mag]
	return lvectors

def normalvec(membrane):
	nhat = pd.DataFrame(columns=['nx', 'ny', 'nz', 'mag' ])
	for i, row in membrane.triangle.iterrows():
		v1=membrane.triangle['tv1'].loc[i]
		v2=membrane.triangle['tv2'].loc[i]
		v3=membrane.triangle['tv3'].loc[i]
		[p1,p2,p3]=np.array(membrane.vertex.loc[[v1,v2,v3],list('xyz')])

		[nx,ny,nz]=np.cross(p1,p2)+np.cross(p2,p3)+np.cross(p3,p1)
		a=[nx,ny,nz]
		mag=LA.norm(a)
		nhat.loc[i] = [nx/mag, ny/mag, nz/mag, mag]
	return nhat

def vertex_angle(d,c,tag,m,membrane,lvectors,nhat):

	

	n1=membrane.line['lt1'].loc[d]
	n2=membrane.line['lt1'].loc[c]
	nhat_add = nhat.loc[n1].add(nhat.loc[n2], fill_value=0).as_matrix(columns=None)
	nhat_add = nhat_add[:3]
	[l1,l2]=np.array(lvectors.loc[[d,c],list('xyz')])
	dot_pro=np.dot(l1,l2)
	angle = -dot_pro/((lvectors['mag'].loc[c])*(lvectors['mag'].loc[d]))
	angle = np.arccos(max(min(angle, 1.0), -1.0))

	cross_pro=-np.cross(l1,l2)/((lvectors['mag'].loc[c])*(lvectors['mag'].loc[d]))
	if(tag==1):
		cross_pro=-cross_pro

	tdot=np.dot(nhat_add,cross_pro)
	# print("tdot", tdot)
	if (tdot<0):
		if(membrane.vertex['pnt'].loc[m]>4):
			return 0.001
		else:
			angle = (2*np.pi)-angle
			# print("if", angle)
			return angle
	else:
		# print("else", angle)
		return angle	
	
def calRave(membrane):
	rsq=0
	xc=0
	yc=0
	zc=0
	for i, row in membrane.vertex.iterrows():
		[xi,yi,zi]=np.array(membrane.vertex.loc[i,list('xyz')])
		xc=xc+xi
		yc=yc+yi
		zc=zc+zi
	xc=xc/len(membrane.vertex)
	yc=yc/len(membrane.vertex)
	zc=zc/len(membrane.vertex)

	for i, row in membrane.vertex.iterrows():
		[xi,yi,zi]=np.array(membrane.vertex.loc[i,list('xyz')])
		print("cal Rave", xi,yi,zi)
		rsqi=((xi-xc)**2)+((yi-yc)**2)+((zi-zc)**2)
		print(rsqi)
		rsq=rsq+np.sqrt(rsqi)
	# print(rsq)
	rave=rsq/len(membrane.vertex)
	return rave

def calPrimeter(membrane):
	lvectors=linevec(membrane)
	peri=0
	for i in membrane.line[membrane.line.lnt==1].index:
		l=lvectors.loc[i,'mag']
		peri=l+peri

	return peri

