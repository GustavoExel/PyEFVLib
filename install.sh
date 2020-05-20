rm -r PyEFVLib/simulation/CGNS/build/*
cmake -B PyEFVLib/simulation/CGNS/build/ -S PyEFVLib/simulation/CGNS
make -C PyEFVLib/simulation/CGNS/build/ -s