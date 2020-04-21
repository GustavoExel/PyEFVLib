#include "include/CgnsOpener.hpp"
#include <iostream>
#include <boost/filesystem.hpp>
// #include <boost/shared_ptr.hpp>
#include <cgnslib.h>
#include <cgns_io.h>


int main()
{
	CgnsOpener cgnsOpener("/home/gustavoe/Documents/Sinmec/GTRelated/GridReader/results/Results.cgns", "Modify");
	std::string p="home/gustavoe/Documents/Sinmec/GTRelated/GridReader/results";
	boost::filesystem::path file(p);
	std::cout << file.root_name() << std::endl;



	int fileIndex;
	cg_open(boost::filesystem::absolute(file).c_str(), 0, &fileIndex);

	cg_close(fileIndex);

	return 0;
}
