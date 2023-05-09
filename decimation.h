#ifndef DECIMATION
#define DECIMATION

#include <string>
#include <vector>
#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>

class Decimation
{
    private:
        double t;
        int numVertices;
        
        std::vector<std::vector<double>> vertices;
        std::vector<int> verticesId;
        std::vector<std::vector<std::vector<double>>> triangles;

        std::vector<std::vector<int>> validPairs;
        std::vector<int> validPairs1D;
        
        std::vector<double> qmatrices;
        
    public:
        void ReadFromVtk(std::string filename);
        void SelectValidPairs(double treshold);
        void ComputeQMatrices();
};

/*int* LinearSysSolve();

int* ReadFromVtk(string filename);

double* QMatrices(int* vertices);

int* QMatricesSort(int* vertices, double* qmatrices);

int* Remove vertices(int* vertices_sort);*/



#endif // DECIMATION
