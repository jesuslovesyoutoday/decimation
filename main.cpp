#include <iostream>
#include <string>
#include "decimation.h"

int main()
{
    std::string filename;
    std::cout << "Enter .vtk filename ..." << std::endl;
    std::cin >> filename;
    
    std::string filename_out;
    std::cout << "Enter .vtk output filename ..." << std::endl;
    std::cin >> filename_out;

    double t = 0.001;
    double param = 0.1; // want to remove param*100% of all vertices
    
    Decimation decimation;
    
    decimation.ReadFromVtk(filename);
    decimation.SelectValidPairs(t);
    decimation.ComputeQMatrices();
    decimation.SortValidPairs();
    decimation.RemoveVertices(param);
    decimation.WriteNewVtk(filename, filename_out);
    
    return 0;
}
