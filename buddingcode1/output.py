import pandas as pd
import csv

def output_params(paramfile_name, r0, l0, gamma, ks, kb, th0):
	
	with open(paramfile_name,'a') as timing:
		
		timing.write("R0:%f\n"%r0)
		timing.write("l0:%f\n"%l0)
		timing.write("ks:%f\n"%ks)
		timing.write("kb:%f\n"%kb)
		timing.write("gamma:%f\n"%gamma)
		timing.write("th0:%f\n"%th0)
		
def lastframe(lastframe_name,gamma,r0,membrane):

	with open(lastframe_name,'a') as lastf:
		writer = csv.writer(lastf, delimiter='\t')
		if membrane=="null":
			writer.writerow([gamma,r0,0,0])
		else:	
			writer.writerow([gamma,r0,len(membrane.triangle),len(membrane.vertex)])
			for j, row in membrane.triangle.iterrows():
 				[v1,v2,v3]=membrane.triangle.loc[j,['tv1','tv2','tv3']].astype(int)
 				writer.writerow([v1, v2, v3])
 				
 			