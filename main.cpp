#include <iostream>
#include <string>
#include "decimation.h"

int main()
{
    std::string filename;
    std::cout << "Enter .vtk filename ..." << std::endl;
    std::cin >> filename;
    
    /*int* vertices = ReadFromVtk(filename);
    
    double* qmatrices = QMatrices(vertices);
    
    int* vertices_sort = QMatricesSort(vertices, qmatrices);
    
    int* vertices_new = RemoveVertices(vertices_sort);
    
    
    delete vertices;
    delete qmatrics;
    delete vertices_sort;
    delete vertices_new;*/
    double t = 1e-10;
    Decimation decimation;
    decimation.ReadFromVtk(filename);
    decimation.SelectValidPairs(t);
    decimation.ComputeQMatrices();
    
    return 0;
}
