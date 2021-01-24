import numpy as np
import pandas as pd
import os

def validate_fluxes(fileName, baseDir):
	data = pd.read_csv(os.path.join(baseDir, fileName))
	col = list(data.columns)[-3]

	X = sorted(list(set([ round(x,5) for x in data["X1"] ])))
	x1_in, x1_out = X[0], X[-2]
	x2_in, x2_out = X[1], X[-1]

	delta = (x2_in-x1_in)/4

	s1 = sum( [ qx*(x2_in/2) for qx, x1, x2 in zip(data[ col ], data["X1"], data["X2"]) if abs(x1-x1_in)<delta and abs(x2-x2_in)<delta ] )
	s2 = sum( [ qx*(x2_in/2) for qx, x1, x2 in zip(data[ col ], data["X1"], data["X2"]) if abs(x1-x1_out)<delta and abs(x2-x2_out)<delta] )

	k1=float( fileName.split('=')[1].split(',')[0] )
	k2=float( fileName.split('=')[2].split('-')[0] )
	size=int( fileName.split('-')[2].split('x')[0] )

	dfData["k1"].append(k1);dfData["k2"].append(k2);dfData["size"].append(size);dfData["sum_in"].append(s1);dfData["sum_out"].append(s2)

	# print(s1, s2)
	# print(f"{ fileName } - Σq(x={(x1_in+x2_in)/2})={s1} - Σq(x={(x1_out+x2_out)/2})={s2} - (Σq(x2)-Σq(x1))/Σq(x2)={(s2-s1):.5e}")

if __name__ == "__main__":
	baseDir = "fluxos - vec\\results"
	dfData = {"k1":[], "k2":[], "size":[], "sum_in":[], "sum_out":[]}
	for fileName in os.listdir(baseDir):
		validate_fluxes(os.path.join(baseDir, fileName).replace(baseDir,'')[1:], baseDir)

	pd.DataFrame(dfData).to_csv( "soma dos fluxos.csv", index=False )
