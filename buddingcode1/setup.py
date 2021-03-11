import param
from param import  R3_2, Rcore, Rvertex, l0, stretchingp, bendingp, ljc, Rcut, epslj
import hoomd.deprecated as hoomdold
from hoomd.integrate import _integrator, _integration_method
import hoomd
from hoomd import md
import pandas as pd
import numpy as np
import energy

class sphere():
	"""docstring for polymer"""
	def __init__(self,center):
	
		self.n=1
		self.rp=Rcore
		self.particles=center
		self.pids=0
		self.typeid=[0]
		self.bodyid=-1 
		

class mesh():
   
    def __init__(self,vertex,line,triangle):
    	self.rv=Rvertex
    	self.l0=l0
    	self.bodyid=-1
    	self.bondid=0
    	self.dihedralid=0
    	self.vertex=vertex
    	self.line=line
    	self.triangle=triangle
    	self.n=len(self.vertex)
    	self.particles=self.vertex[list('xyz')].values
    	self.pids=[1]*self.n
    	self.typeid=[1]
    	self.bonds=np.int_(self.line[['lv1','lv2']])
    	self.dihedrals=[]
    	self.dihedralid=0

    def update_mesh_info(self):
    	self.n=len(self.vertex)
    	self.pids=self.typeid*len(self.vertex)
    	self.particles=np.array(self.vertex[list('xyz')])
    	self.bonds=np.int_(self.line[['lv1','lv2']])
    	dihedralgroup=[]
    	for l, row in self.line.iterrows():
    		if l == 0: continue
    		if (self.line['lnt'].loc[l])==2 and (self.line['double'].loc[l])==0:
    			v0 = self.line['lv1'].loc[l]
    			v1 = self.line['lv2'].loc[l]
    			t0 = self.line['lt1'].loc[l]
    			t1 = self.line['lt2'].loc[l]
    			[tv1,tv2,tv3]=self.triangle.loc[t0,['tv1','tv2','tv3']]
    			v2 = tv1+tv2+tv3-v0-v1
    			[tv1,tv2,tv3]=self.triangle.loc[t1,['tv1','tv2','tv3']]
    			v3 = tv1+tv2+tv3-v0-v1
    			dihedralgroup=dihedralgroup+[[v2,v1,v0,v3]]
    			self.dihedrals=np.int_(dihedralgroup)
		#-----------------------------------

def init_snap(th0,ks,kb):
	hoomd.context.initialize("")
	#define and assign snapshot
	snap=hoomd.data.make_snapshot(N=0,
		box=hoomd.data.boxdim(L=200*l0),
		particle_types=['C','M'], # M:mesh->1, C:core->0
		bond_types=['mesh'],
		dihedral_types=['capsomer'],
		angle_types=['ds'],
		dtype='double')
	core, membrane, system=structure_assignments(snap,th0,ks,kb)
	print("sanaz 2", th0)
	return snap, core, membrane, system

def structure_assignments(snap,th0,ks,kb):
	# print("before assignment",snap.particles.N)
	core=assign_core(snap)
	membrane=assign_mesh(snap)
	# print("after assignment",snap.particles.N)
	# print_snap(snap)
	system=assign_hoomd(snap,np.pi-th0,ks,kb)
	return core, membrane, system

def assign_hoomd(snap, t0,ks,kb):
	system=hoomd.init.read_snapshot(snap)
	if stretchingp:
		bond_harmonic = hoomd.md.bond.harmonic()
		bond_harmonic.bond_coeff.set('mesh', k=ks, r0=l0)
		
	if bendingp:
		dtable = hoomd.md.dihedral.table(width=1000)
		dtable.dihedral_coeff.set('capsomer', func=energy.dihedral_harmonic, coeff=dict(kappa= kb, theta0=t0))
		print("in hoomd", t0)
	
	nl=hoomd.md.nlist.cell()
	if ljc:
		lj_force = hoomd.md.pair.lj(r_cut=Rcut, nlist=nl)
		lj_force.set_params(mode="shift")
		lj_force.pair_coeff.set('C', 'M', epsilon=epslj, sigma=Rcore+Rvertex, alpha=2) #sigma=Rcore+Rvertex
		lj_force.pair_coeff.set('C', 'C', epsilon=0, sigma=0)
		lj_force.pair_coeff.set('M', 'M', epsilon=0, sigma=0)

	return system
def assign_core(snap):
	center=np.array([[0.,0.,-Rcore]])
	core=sphere(center)
	pstart=snap.particles.N
	pend=pstart+core.n
	snap.particles.resize(pend)
	snap.particles.diameter[0] = 2*Rcore
	snap.particles.position[0]=core.particles
	snap.particles.typeid[0]= core.pids
	snap.particles.body[0]= core.bodyid
	# print("sanaz24",pstart,pend)
	return core
	

def assign_mesh(snap):
	vertex=pd.DataFrame(columns=['x','y','z','pnt','edge'])
	line=pd.DataFrame(columns=['lv1', 'lv2', 'lt1', 'lt2', 'lnt', 'lr', 'll', 'double'])
	triangle=pd.DataFrame(columns=['tv1', 'tv2', 'tv3', 'tl1', 'tl2', 'tl3', 'sign'])
	# print("before")
	membrane=mesh(vertex,line,triangle)
	# print("after")
	#-----------------

	membrane.vertex.loc[1] = [0., 0., 0., 1, 1]
	membrane.vertex.loc[2] = [0.5, R3_2, 0.0, 1, 1]
	membrane.vertex.loc[3] = [1.0, 0., 0., 1, 1]

	membrane.line.loc[1] = [1, 2, 1, 1000, 1, 3, 2, 0]
	membrane.line.loc[2] = [2, 3, 1, 1000, 1, 1, 3, 0]
	membrane.line.loc[3] = [3, 1, 1, 1000, 1, 2, 1, 0]
	membrane.triangle.loc[1] = [1, 2, 3, 1, 2, 3, 0]
	#-------------------
	membrane.update_mesh_info()
	#-------------------
	pstart=snap.particles.N
	pend=pstart+membrane.n
	snap.particles.resize(pend)
	
	for i, coordinates in enumerate(membrane.particles, start=pstart):
		snap.particles.diameter[i] = 2*membrane.rv
		snap.particles.position[i]=coordinates
		snap.particles.body[i] = membrane.bodyid
	for i, pid in enumerate(membrane.pids, start=pstart):
		snap.particles.typeid[i] = pid
	#-----------------
	bstart=snap.bonds.N
	snap.bonds.resize(bstart+len(membrane.bonds))
	for i, bond in enumerate(membrane.bonds,start=bstart):
		snap.bonds.group[i]=bond 
		snap.bonds.typeid[i]=membrane.bondid
	#------------------
	dstart=snap.dihedrals.N

	snap.dihedrals.resize(dstart+len(membrane.dihedrals))
	for i, dihedral in enumerate(membrane.dihedrals,start=dstart):
		snap.dihedrals.group[i]=dihedral
		snap.dihedrals.typeid[i]=membrane.dihedralid
	return membrane

def snap_assignments(snap,membrane):
	
	pstart=1
	pend=pstart+membrane.n
	snap.particles.resize(pend)
	
	for i, coordinates in enumerate(membrane.particles, start=pstart):
		snap.particles.diameter[i] = 2*membrane.rv
		snap.particles.position[i]=coordinates
		snap.particles.body[i] = membrane.bodyid
	for i, pid in enumerate(membrane.pids, start=pstart):
		snap.particles.typeid[i] = pid
	
	#---------bonds
	bstart=0
	snap.bonds.resize(bstart+len(membrane.bonds))
	for i, bond in enumerate(membrane.bonds,start=bstart):
		snap.bonds.group[i]=bond
		snap.bonds.typeid[i]=membrane.bondid
	#---------dihedrals
	dstart=0
	snap.dihedrals.resize(dstart+len(membrane.dihedrals))
	for i, dihedral in enumerate(membrane.dihedrals,start=dstart):
		snap.dihedrals.group[i]=dihedral
		snap.dihedrals.typeid[i]=membrane.dihedralid

def print_snap(snap):
	print("number of particles",snap.particles.N)
	print("position",snap.particles.position)
	print("typeid",snap.particles.typeid)
	print("number of bonds", snap.bonds.group)
	print("number of dihedrals", snap.dihedrals.N)
	print("number of dihedrals", snap.dihedrals.group)