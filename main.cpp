#include <iostream>
#include <string>
#include "decimation.h"

int main()
{
    std::string filename;
    std::cout << "Enter .vtk filename ..." << std::endl;
    std::cin >> filename;

    double t = 1e-10;
    double param = 0.2; // want to remove 20% of all vertices
    
    Decimation decimation;
    
    decimation.ReadFromVtk(filename);
    decimation.SelectValidPairs(t);
    decimation.ComputeQMatrices();
    decimation.SortValidPairs();
    std::cout << "hi" << std::endl;
    decimation.RemoveVertices(param);
    
    return 0;
}
