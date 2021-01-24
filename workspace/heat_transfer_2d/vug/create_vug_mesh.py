import os, numpy as np

def writeVug(width, height, nxT, nyT, nxI, nyI, outputPath):
	# nxT, nxT stands for the total number of divisions
	# nxI, nyI stands for the number of divisions in the inner domain

	physicalNames = { "Outer":2, "Inner":2, "North":1, "South":1, "West":1, "East":1 }
	nodes = [(x,y,0.0) for y in np.linspace(0,height,(nyT+1)) for x in np.linspace(0,width,(nxT+1))]

	def tellRegion(idx):
		# Tells whether element belongs to outer or inner region
		# Outer = 1, Inner = 2
		row = idx // (nxT+1)
		col = idx % (nxT+1)
		firstInnerX = (nxT - nxI + 1)//2
		firstInnerY = (nyT - nyI + 1)//2

		lastInnerX = firstInnerX + nxI
		lastInnerY = firstInnerY + nyI

		if firstInnerX <= col < lastInnerX and firstInnerY <= row < lastInnerY:
			return 2
		else:
			return 1

	elements = []
	elements += [[len(elements)+i+1, *(1, 2), *(3, 3), *(vIdx+1,vIdx+2)] for i, vIdx in enumerate( range(nyT*(nxT+1), (nxT+1)*(nyT+1)-1) )]
	elements += [[len(elements)+i+1, *(1, 2), *(4, 4), *(vIdx+1,vIdx+2)] for i, vIdx in enumerate( range(0, nxT) )]
	elements += [[len(elements)+i+1, *(1, 2), *(5, 5), *(vIdx+1,vIdx+nxT+2)] for i, vIdx in enumerate( range(0, nyT*(nxT+1), (nxT+1)) )]
	elements += [[len(elements)+i+1, *(1, 2), *(6, 6), *(vIdx+1,vIdx+nxT+2)] for i, vIdx in enumerate( range(nxT, (nxT+1)*(nyT+1)-1, (nxT+1)) )]
	elements += [[len(elements)+i+1, *(3, 2), *(tellRegion(vIdx), tellRegion(vIdx)), *(vIdx+1, vIdx+2, vIdx+nxT+3, vIdx+nxT+2)] for i, vIdx in enumerate( [ vIdx for vIdx in range((nxT+1)*(nyT+1)-(nxT+2)) if (vIdx+1)%(nxT+1)!=0 ] )]

	text = ""
	text += "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n"
	text += str( len(physicalNames) ) + "\n"
	text += "".join([ "{} {} \"{}\"\n".format(physicalNames[name], idx+1, name) for idx, name in enumerate(physicalNames.keys()) ])
	text += "$EndPhysicalNames\n$Nodes\n"
	text += str( len(nodes) ) + "\n"
	text += "".join([ "{} {} {} {}\n".format(idx+1, x, y, z) for idx,(x,y,z) in enumerate(nodes) ])
	text += "$EndNodes\n$Elements\n"
	text += str( len(elements) ) + "\n"
	text += "\n".join([" ".join([str(ev) for ev in e]) for e in elements])
	text += "\n$EndElements\n"

	with open(outputPath, "w") as f:
		f.write(text)

if __name__ == "__main__":
	writeVug(1.0, 1.0, 3, 3, 1, 1, os.path.join( os.path.dirname(__file__), "vug.msh" ))