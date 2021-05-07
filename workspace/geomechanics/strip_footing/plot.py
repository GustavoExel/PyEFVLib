import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

Cv = 2.268843885e-08

kPa = 1000
mm = 1e-3

width = 16.0
height = 8.0
load_width = 1.0

load = 1e+05

dts = [4407.5,8815.1,22038,44075,88151,220380,440750,881510,2203800,4407500,8815100,22038000,44075000]
Ns = [10,5,6,5,5,6,5,5,6,5,5,6,5]
timeSteps = [t for dt,N in zip(dts,Ns) for t in N*[dt]]
times = np.array([sum(timeSteps[:i]) for i in range(len(timeSteps))])

fig1, ax1 = plt.subplots(2,2)
fig1.subplots_adjust(hspace=0.4)

fig2, ax2 = plt.subplots(2,2)
fig2.subplots_adjust(hspace=0.4)

fig3, ax3 = plt.subplots(2,2)
fig3.subplots_adjust(hspace=0.4)


dirname = os.path.realpath(os.path.dirname(__file__))
def plot(ax1,ax2,ax3, name,marker,color,label):
	filePath = f"{name}/{name}.csv"
	df = pd.read_csv(os.path.join(dirname,filePath))
	X = df["X"]
	Y = df["Y"]
	numberOfTimeSteps = int(df.keys()[-1].split(" - ")[0].replace("TimeStep",""))

	time = times[1:1+numberOfTimeSteps]
	dst2 = lambda x,y,x0,y0:(x-x0)**2+(y-y0)**2
	idxA = min([(idx,dst2(x,y,0.0,height)) for idx,(x,y) in enumerate(zip(X,Y))],key=lambda t:t[1])[0]
	idxB = min([(idx,dst2(x,y,0.0,height-load_width)) for idx,(x,y) in enumerate(zip(X,Y))],key=lambda t:t[1])[0]
	idxC = min([(idx,dst2(x,y,load_width,height-load_width)) for idx,(x,y) in enumerate(zip(X,Y)) if x<=load_width],key=lambda t:t[1])[0]
	idxD = min([(idx,dst2(x,y,load_width,height)) for idx,(x,y) in enumerate(zip(X,Y)) if x<=load_width],key=lambda t:t[1])[0]

	# Pressure
	pA = np.array([df[f"TimeStep{s} - p"][idxA] for s in range(1,numberOfTimeSteps+1)])
	ax1[0][0].scatter(time*Cv/(load_width**2), pA/load, marker=marker,color=color, label=label)
	ax1[0][0].set_title("A")
	ax1[0][0].set_ylabel(r"Normalized Pressure ($p/p_0$)")
	ax1[0][0].set_xlabel("Time Factor (T)")
	ax1[0][0].semilogx()
	ax1[0][0].legend()

	pB = np.array([df[f"TimeStep{s} - p"][idxB] for s in range(1,numberOfTimeSteps+1)])
	ax1[1][0].scatter(time*Cv/(load_width**2), pB/load, marker=marker,color=color)
	ax1[1][0].set_title("B")
	ax1[1][0].set_ylabel(r"Normalized Pressure ($p/p_0$)")
	ax1[1][0].set_xlabel("Time Factor (T)")
	ax1[1][0].semilogx()

	pC = np.array([df[f"TimeStep{s} - p"][idxC] for s in range(1,numberOfTimeSteps+1)])
	ax1[1][1].scatter(time*Cv/(load_width**2), pC/load, marker=marker,color=color)
	ax1[1][1].set_title("C")
	ax1[1][1].set_ylabel(r"Normalized Pressure ($p/p_0$)")
	ax1[1][1].set_xlabel("Time Factor (T)")
	ax1[1][1].semilogx()

	pD = np.array([df[f"TimeStep{s} - p"][idxD] for s in range(1,numberOfTimeSteps+1)])
	ax1[0][1].scatter(time*Cv/(load_width**2), pD/load, marker=marker,color=color)
	ax1[0][1].set_title("D")
	ax1[0][1].set_ylabel(r"Normalized Pressure ($p/p_0$)")
	ax1[0][1].set_xlabel("Time Factor (T)")
	ax1[0][1].semilogx()


	# X Displacements
	uA = np.array([df[f"TimeStep{s} - u_x"][idxA] for s in range(1,numberOfTimeSteps+1)])
	ax2[0][0].scatter(time, uA/mm, marker=marker,color=color,label=label)
	ax2[0][0].set_title("A")
	ax2[0][0].set_ylabel("X Displacement (mm)")
	ax2[0][0].set_xlabel("Time (s)")
	ax2[0][0].set_ylim((-14e-4,14e-4))	
	ax2[0][0].semilogx()
	ax2[0][0].legend()

	uB = np.array([df[f"TimeStep{s} - u_x"][idxB] for s in range(1,numberOfTimeSteps+1)])
	ax2[1][0].scatter(time, uB/mm, marker=marker,color=color)
	ax2[1][0].set_title("B")
	ax2[1][0].set_ylabel("X Displacement (mm)")
	ax2[1][0].set_xlabel("Time (s)")
	ax2[1][0].set_ylim((-14e-4,14e-4))	
	ax2[1][0].semilogx()

	uC = np.array([df[f"TimeStep{s} - u_x"][idxC] for s in range(1,numberOfTimeSteps+1)])
	ax2[1][1].scatter(time, uC/mm, marker=marker,color=color)
	ax2[1][1].set_title("C")
	ax2[1][1].set_ylabel("X Displacement (mm)")
	ax2[1][1].set_xlabel("Time (s)")
	ax2[1][1].semilogx()

	uD = np.array([df[f"TimeStep{s} - u_x"][idxD] for s in range(1,numberOfTimeSteps+1)])
	ax2[0][1].scatter(time, uD/mm, marker=marker,color=color)
	ax2[0][1].set_title("D")
	ax2[0][1].set_ylabel("X Displacement (mm)")
	ax2[0][1].set_xlabel("Time (s)")
	ax2[0][1].semilogx()

	# Y Displacements
	vA = np.array([df[f"TimeStep{s} - u_y"][idxA] for s in range(1,numberOfTimeSteps+1)])
	ax3[0][0].scatter(time, vA/mm, marker=marker,color=color, label=label)
	ax3[0][0].set_title("A")
	ax3[0][0].set_ylabel("Y Displacement (mm)")
	ax3[0][0].set_xlabel("Time (s)")
	ax3[0][0].semilogx()
	ax3[0][0].legend()

	vB = np.array([df[f"TimeStep{s} - u_y"][idxB] for s in range(1,numberOfTimeSteps+1)])
	ax3[1][0].scatter(time, vB/mm, marker=marker,color=color)
	ax3[1][0].set_title("B")
	ax3[1][0].set_ylabel("Y Displacement (mm)")
	ax3[1][0].set_xlabel("Time (s)")
	ax3[1][0].semilogx()

	vC = np.array([df[f"TimeStep{s} - u_y"][idxC] for s in range(1,numberOfTimeSteps+1)])
	ax3[1][1].scatter(time, vC/mm, marker=marker,color=color)
	ax3[1][1].set_title("C")
	ax3[1][1].set_ylabel("Y Displacement (mm)")
	ax3[1][1].set_xlabel("Time (s)")
	ax3[1][1].semilogx()

	vD = np.array([df[f"TimeStep{s} - u_y"][idxD] for s in range(1,numberOfTimeSteps+1)])
	ax3[0][1].scatter(time, vD/mm, marker=marker,color=color)
	ax3[0][1].set_title("D")
	ax3[0][1].set_ylabel("Y Displacement (mm)")
	ax3[0][1].set_xlabel("Time (s)")
	ax3[0][1].semilogx()

	df = pd.DataFrame({"time":time,
						"pA":pA,"pB":pB,"pC":pC,"pD":pD,
						"uA":uA,"uB":uB,"uC":uC,"uD":uD,
						"vA":vA,"vB":vB,"vC":vC,"vD":vD })
	df.to_csv(f"{name}.csv",index=False)



plot(ax1,ax2,ax3,"strip_footing_uneq_spacing_xy", 's', 'b',"Unequal Spacing XY")
plot(ax1,ax2,ax3,"strip_footing_uneq_spacing_x", 'v', 'k',"Unequal Spacing X")
plot(ax1,ax2,ax3,"strip_footing_eq_spacing", 'x', 'r',"Equal Spacing")
plot(ax1,ax2,ax3,"strip_footing_super_fine", 'x', 'g',"More refined")

TAa = [1.0000e-04,1.9971e-04,2.9806e-04,3.9883e-04,4.9619e-04,5.9526e-04,7.0121e-04,8.1113e-04,9.0474e-04,9.9094e-04,1.9790e-03,2.9536e-03,3.9522e-03,5.0073e-03,6.0070e-03,6.9486e-03,7.8928e-03,8.9654e-03,9.6425e-03,1.9611e-02,2.9269e-02,3.9164e-02,4.9619e-02,5.9526e-02,6.8856e-02,7.9650e-02,8.8842e-02,9.9094e-02,1.9790e-01,2.9536e-01,3.9522e-01,4.9170e-01,5.7922e-01,6.8233e-01,7.8928e-01,8.9654e-01,9.8196e-01,1.9611e+00,2.9269e+00,3.9164e+00,4.7846e+00,5.8452e+00,6.6395e+00]
pAa = [9.5618e-03,1.7530e-02,2.5498e-02,3.1474e-02,4.1434e-02,4.7410e-02,5.7371e-02,6.5339e-02,7.3307e-02,8.1275e-02,1.6295e-01,2.4462e-01,3.2829e-01,4.1394e-01,4.9960e-01,5.8526e-01,6.6892e-01,7.5458e-01,8.3825e-01,8.5020e-01,8.4622e-01,8.3426e-01,8.1833e-01,8.0040e-01,7.8247e-01,7.6056e-01,7.4263e-01,7.2271e-01,5.6733e-01,4.5976e-01,3.8606e-01,3.3028e-01,2.8645e-01,2.5259e-01,2.2470e-01,2.0080e-01,1.8287e-01,9.3227e-02,5.7371e-02,4.1434e-02,3.5458e-02,2.9482e-02,2.5498e-02]
ax1[0][0].plot(TAa, pAa, marker='.',color='y', label="Manoharan")
ax1[0][0].legend()

TBa = [1.0000e-04,2.0338e-04,3.0354e-04,4.0616e-04,5.0531e-04,6.0619e-04,7.1409e-04,8.1113e-04,9.0474e-04,1.0091e-03,2.0153e-03,3.0079e-03,4.0248e-03,5.0073e-03,6.0070e-03,7.0762e-03,8.0378e-03,8.9654e-03,9.8196e-03,1.9611e-02,2.9269e-02,3.9883e-02,4.9619e-02,5.9526e-02,6.8856e-02,7.8213e-02,9.9094e-02,1.9790e-01,2.9536e-01,3.9522e-01,4.9170e-01,5.7922e-01,6.8233e-01,7.7505e-01,8.8037e-01,9.6425e-01,1.9611e+00,2.9269e+00,3.9164e+00,4.8724e+00,5.8452e+00,6.7614e+00,7.8213e+00,9.5551e+00,1.9790e+01,2.9536e+01,3.9522e+01,5.0073e+01,5.8986e+01,6.9486e+01,8.0378e+01,9.1301e+01,1.0184e+02]
pBa = [5.5215e-03,9.2025e-03,1.2883e-02,1.7791e-02,2.2699e-02,2.7607e-02,3.2515e-02,3.7423e-02,4.1104e-02,4.7239e-02,9.5092e-02,1.4049e-01,1.8834e-01,2.3620e-01,2.8405e-01,3.3190e-01,3.7853e-01,4.2638e-01,4.7423e-01,4.7669e-01,4.7914e-01,4.8405e-01,4.8773e-01,4.9141e-01,4.9387e-01,4.9632e-01,4.9509e-01,4.5583e-01,4.0675e-01,3.6135e-01,3.2454e-01,2.9387e-01,2.6564e-01,2.4233e-01,2.2393e-01,2.0675e-01,1.1595e-01,7.9141e-02,6.1963e-02,5.0920e-02,4.2331e-02,3.8650e-02,3.3742e-02,2.8834e-02,1.2883e-02,5.5215e-03,3.0675e-03,3.0675e-03,1.8405e-03,1.8405e-03,6.1350e-04,6.1350e-04,6.1350e-04]
ax1[1][0].plot(TBa, pBa, marker='.',color='y')


TCa = [1.0000e-04,1.9916e-04,3.0221e-04,4.0392e-04,5.0210e-04,6.0190e-04,6.9586e-04,8.0447e-04,9.1333e-04,1.0183e-03,1.9916e-03,3.0221e-03,4.0392e-03,5.0210e-03,6.0190e-03,6.9586e-03,7.7583e-03,8.8081e-03,9.8203e-03,1.9559e-02,2.9145e-02,3.8254e-02,4.8422e-02,5.7004e-02,6.7108e-02,7.6189e-02,8.6498e-02,9.6439e-02,1.9559e-01,2.9145e-01,3.8954e-01,4.8422e-01,5.8047e-01,6.8335e-01,7.9002e-01,8.8081e-01,9.6439e-01,1.9559e+00,2.9145e+00,3.9667e+00,4.9308e+00,5.8047e+00,6.8335e+00,7.9002e+00,9.8203e+00,1.9916e+01,2.9678e+01,3.9667e+01,4.8422e+01,5.8047e+01,6.8335e+01]
pCa = [2.7833e-03,5.1690e-03,8.7475e-03,1.2326e-02,1.5905e-02,1.8290e-02,2.0676e-02,2.3062e-02,2.6640e-02,3.0219e-02,5.8847e-02,8.9861e-02,1.2207e-01,1.5308e-01,1.8529e-01,2.1511e-01,2.4732e-01,2.7952e-01,3.1292e-01,3.3320e-01,3.4155e-01,3.4513e-01,3.4632e-01,3.4632e-01,3.4513e-01,3.4274e-01,3.4036e-01,3.3797e-01,3.0934e-01,2.8429e-01,2.5924e-01,2.3897e-01,2.1988e-01,2.0318e-01,1.8767e-01,1.7575e-01,1.6382e-01,9.8211e-02,6.9583e-02,5.4076e-02,4.5726e-02,3.9761e-02,3.4990e-02,3.1412e-02,2.5447e-02,1.4712e-02,6.3618e-03,3.9761e-03,2.7833e-03,2.7833e-03,2.7833e-03]
ax1[1][1].plot(TCa, pCa, marker='.',color='y')


plt.show()