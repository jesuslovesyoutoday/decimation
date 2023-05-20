#ifndef DECIMATION
#define DECIMATION

#include <string>
#include <vector>
#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>

class Decimation
{
    private:
        double t;
        int numVertices;
        int newNumVertices;
        
        std::vector<std::vector<double>> vertices;
        std::vector<std::vector<double>> newVertices;
        std::vector<int> verticesId;
        std::vector<int> newVerticesId;
        std::vector<std::vector<std::vector<double>>> triangles;
        std::vector<std::vector<int>> trId;

        std::vector<std::vector<int>> validPairs;
        
        std::vector<std::pair<int, double>> qmatrices;
        
    public:
        void ReadFromVtk(std::string filename);
        void SelectValidPairs(double treshold);
        void ComputeQMatrices();
        void SortValidPairs();
        void RemoveVertices(double param);
        void WriteNewVtk(std::string filename, std::string filename_out);
};

#endif // DECIMATION
