import param 
# from param import file_name, r0, l0, gamma, ks, kb
from param import Rcore,R0,l0,gamma,th0
import hoomd
from hoomd import md
import math
import hoomd.deprecated as hoomdold
from hoomd.integrate import _integrator, _integration_method
import setup
import movements
import energy
import output
import tools

def start(ks,kb,file_name,paramfile_name,lastframe_name):

	print("hi")
	output.output_params(paramfile_name, R0, l0, gamma, ks, kb, th0)
	snap, core, membrane, system=setup.init_snap(th0,ks,kb)
	log = hoomd.analyze.log(filename='None', quantities=['N','potential_energy','pair_lj_energy','bond_harmonic_energy','dihedral_table_energy'], period=1, header_prefix='#', overwrite=True)

	hoomd.dump.gsd(filename=file_name,group=hoomd.group.all(),period=None,static=[])
	#--------------
	
	n_t = len(membrane.triangle)
	while n_t <10000:
		n_t = len(membrane.triangle)
		n_v= len(membrane.vertex)
		n_l = len(membrane.line)
		membrane = movements.grow(membrane,system,log,th0)
		membrane=energy.relax(membrane,system)
		hoomd.dump.gsd(filename=file_name,group=hoomd.group.all(),period=None,static=[])
		
		print("end of one growth step")
		
	output.lastframe(lastframe_name,gamma,R0,membrane)
	