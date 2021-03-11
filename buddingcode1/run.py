
import numpy as np
import os
from param import gamma, R0
import start
import math
import output

def run():
	

	new="/Users/sanazpanahandeh/Desktop/buddingcode1/flatmesh"+"/R0="+str(R0)+",gamma="+str(gamma)
	os.mkdir(new)
	lastframe_name=new+"/lastf"+"R0="+str(R0)+",gamma="+str(gamma)+".csv"
	file_name=new+"/membarne"+".gsd"
	# frame_name=new+"/frames"+".csv"
	# energy_file=new+"/energy"+".csv"
	paramfile_name=new+"/paramfile"+".csv"
	# th0 =2 * np.arcsin(l0 / np.sqrt((12 * R0 * R0) - 3*l0*l0))
	ks=np.sqrt(gamma)
	kb=1/ks
	start.start(ks,kb,file_name,paramfile_name,lastframe_name)

run()