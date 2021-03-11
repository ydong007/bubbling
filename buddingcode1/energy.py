
import hoomd.deprecated as hoomdold
from hoomd.integrate import _integrator, _integration_method
import hoomd
from hoomd import md
import pandas as pd
import numpy as np
import math
import setup
from numpy import linalg as LA
from param import ljc, epslj, Rcore, Rvertex


def energy(membrane, log):
	
	total_potential_energy=log.query('potential_energy')
	stretching=log.query('bond_harmonic_energy')
	bending=log.query('dihedral_table_energy')
	if ljc:
		lj=log.query('pair_lj_energy')
	else:
		lj=0

	Energy = pd.DataFrame(columns=['potential_energy', 'bond_harmonic_energy', 'dihedral_harmonic_energy','pair_lj_energy'])
	Energy.loc[0, 'potential_energy']=log.query('potential_energy')
	Energy.loc[0, 'bond_harmonic_energy']=log.query('bond_harmonic_energy')
	Energy.loc[0, 'dihedral_harmonic_energy']=bending
	
	if ljc:
		Energy.loc[0, 'pair_lj_energy']=log.query('pair_lj_energy')
	
	return stretching, bending, stretching+bending
def dihedral_harmonic(theta, kappa, theta0):
	V =  kappa * (1 - np.cos(theta-theta0))
	F = - kappa*np.sin(theta-theta0)
	return (V, F)

def relax(membrane,system):
	snapshot=system.take_snapshot(all=True)
	# print_snap(snapshot)
	membrane.update_mesh_info()
	# print("check # in snapshot E")
	# print_snap(snapshot)
	setup.snap_assignments(snapshot,membrane)
	print("check # in snapshot F")
	# setup.print_snap(snapshot)
	system.restore_snapshot(snapshot)
	print("check # in snapshot G")
	# setup.print_snap(snapshot)
	# print("after update")
	# #---------------
	
	print("check # in snapshot H")
	# setup.print_snap(snapshot)
	# print("sys pA: ",system.particles)
	# print(membrane.vertex)
	# errorcode=relaxation()

	relaxation()
	for p in system.particles: 
		if p.typeid in membrane.typeid:
				# print("p.tag  ", p.tag)
				# print("p.position ", p.position)
				membrane.vertex.loc[p.tag,list('xyz')]=p.position
	#I checked the relaxation by running the linevec function here and comparing the linevec magnitude before and after relaxation. It is really relaxed
	print("check # in snapshot I")
	# setup.print_snap(snapshot)
	# print("sys pB: ",system.particles)

	print("after relaxation")
	# print(membrane.vertex)
	# print(membrane.line)
	# print(membrane.triangle)
	return membrane

def relaxation():
	
	fire=hoomd.md.integrate.mode_minimize_fire(dt=0.01, ftol=1e-4, Etol=1e-4)
	nve=hoomd.md.integrate.nve(group=hoomd.group.nonrigid())
	hoomd.context.current.thermos[-1].disable() 
	print("check the nve")

	# i=0
	count=0
	while not(fire.has_converged()) and count<1:
		hoomd.run(1000)
		count+=1
		
	nve.disable()
	return 1
