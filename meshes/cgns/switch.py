import os
if 'Blank.cgns' in list(os.walk('.'))[0][2]:
	os.system('mv blank.cgns b.cgns && mv Blank.cgns blank.cgns')
else:
	os.system('mv blank.cgns Blank.cgns && mv b.cgns blank.cgns')

