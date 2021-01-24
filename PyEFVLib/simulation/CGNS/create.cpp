#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <cgnslib.h>
#include <boost/filesystem.hpp>


void readData(std::string &filePath, int *sizes, std::string basePath)
{
    std::ifstream data(basePath + "/data.txt");
    data >> filePath;
    data >> *sizes;
    data >> *(sizes+1);
    data >> *(sizes+2);
}

void writeCoords(int &fileIndex, int &baseIndex, int &zoneIndex, int &coordinateIndex, std::string basePath)
{
    std::vector< std::vector<double> > coords(3);
    std::ifstream file(basePath + "/coords.txt");
    std::string line;
    double coordVal;
    int i=0;
    while( getline(file, line) )
    {
        std::stringstream ss(line);
        int j=0;
        while(ss >> coordVal ) { j++;
            coords[i].emplace_back(coordVal);
        }
        i++;
    }

    cg_coord_write(fileIndex, baseIndex, zoneIndex, RealDouble, "CoordinateX", &coords[0][0], &coordinateIndex);
    cg_coord_write(fileIndex, baseIndex, zoneIndex, RealDouble, "CoordinateY", &coords[1][0], &coordinateIndex);    
    cg_coord_write(fileIndex, baseIndex, zoneIndex, RealDouble, "CoordinateZ", &coords[2][0], &coordinateIndex);    
    file.close();
}

struct SectionData {
    std::string name;
    int elementType;
    int begin;
    int end;
};

void writeSections(int &fileIndex, int &baseIndex, int &zoneIndex, int &sectionIndex, std::string basePath) {
    std::ifstream connectivityFile(basePath + "/connectivity.txt");
    std::ifstream sectionsFile(basePath + "/sections.txt");
    std::vector< std::vector<int> > connectivities;
    std::vector< SectionData > sections;

    std::stringstream ss;
    std::string line; int elementIndex; std::vector<int> e;
    while (getline(connectivityFile, line)) {
        ss.clear(); ss << line;
        e.clear();
        while( ss >> elementIndex ) e.emplace_back( elementIndex + 1 );
        connectivities.emplace_back(e);
    }

    std::string wd;
    while (getline(sectionsFile, line)) {
        ss.clear(); ss << line;
        SectionData section;
        ss >> section.name;        
        ss >> section.elementType;        
        ss >> section.begin;        
        ss >> section.end;        
        sections.emplace_back( section );
    }

    for (auto section : sections)
    {
        std::vector<int> sectionConnectivities;
        for (auto e = connectivities.begin() + section.begin; e != connectivities.begin() + section.end+1; ++e) {
            std::vector<int> myvec = *e;
            for (auto v : *e) {
                sectionConnectivities.emplace_back( v );
            }
        }

        std::unordered_map<int,int> sizeType = {{2,BAR_2},{3,TRI_3},{4,TETRA_4},{5,PYRA_5},{6,PENTA_6},{8,HEXA_8}};
        cg_section_write(fileIndex, baseIndex, zoneIndex, section.name.c_str(), ElementType_t(sizeType[section.elementType]), section.begin+1,
                         section.end+1, 0, &sectionConnectivities[0],  &sectionIndex);
    }
}

int main(int argc, char* argv[])
{
    std::string fileName;
    int baseIndex = 1, zoneIndex, fileIndex, coordinateIndex, sectionIndex;
    int physicalDimension = 3;
    int cellDimension = 3;
    int sizes[3];

    std::string baseName = "BASE";
    std::string zoneName = "ZONE";
    std::string basePath(static_cast<boost::filesystem::path>(argv[0]).parent_path().c_str());

    readData(fileName, sizes, basePath);
    cg_open(fileName.c_str(), CG_MODE_WRITE, &fileIndex);
    cg_base_write(fileIndex, baseName.c_str(), cellDimension, physicalDimension, &baseIndex);
    cg_zone_write(fileIndex, baseIndex, zoneName.c_str(), sizes, Unstructured, &zoneIndex);
    writeCoords(fileIndex, baseIndex, zoneIndex, coordinateIndex, basePath);
    writeSections(fileIndex, baseIndex, zoneIndex, sectionIndex, basePath);
}