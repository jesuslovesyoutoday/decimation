#include "decimation.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <xtensor/xtensor.hpp>
#include <xtensor/xio.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <Eigen/Dense>

double DistanceBPoints(std::vector<double> p1, std::vector<double> p2)
{   
    double dx = p1[0] - p2[0];
    double dy = p1[1] - p2[1];
    double dz = p1[2] - p2[2];
    
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

void Decimation::ReadFromVtk(std::string filename)
{
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = reader->GetOutput();

    vtkSmartPointer<vtkPoints> points = unstructuredGrid->GetPoints();
    
    this->numVertices = points->GetNumberOfPoints();
    for (int i = 0; i < this->numVertices * 3; i++)
    {
        double vertex[3];
        points->GetPoint(i, vertex);
        std::vector<double> point;
        for (int j = 0; j < 3; j++)
        {
            point.push_back(vertex[i]);
        }
        
        this->vertices.push_back(point);
    }
    
    
    vtkSmartPointer<vtkCellArray> cells = unstructuredGrid->GetCells();
    int ncells = cells->GetNumberOfCells();
    int n = 0;
    vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
    for (int i = 0; i < ncells; i++) 
    {
        cells->GetCell(i, idList);
        int nverts = idList->GetNumberOfIds();
        for (int j=0; j<nverts; j++) 
        {
            this->verticesId.push_back(idList->GetId(j));
        }
    }
    
    for (int i = 0; i < this->verticesId.size(); i=i+3)
    {
        std::vector<std::vector<double>> triang;
        triang.push_back(this->vertices[this->verticesId[i]]);
        triang.push_back(this->vertices[this->verticesId[i+1]]);
        triang.push_back(this->vertices[this->verticesId[i+2]]);
        this->triangles.push_back(triang);
    }  
    /*for (int i = 0; i < this->numVertices; i++)
    {
        std::cout << this->vertices[i][0] << this->vertices[i][1] << this->vertices[i][2] << std::endl;
    }*/
};

void Decimation::SelectValidPairs(double treshold)
{
    this->t = treshold;
    for (int i = 0; i < this->numVertices; i++)
    {
        if (i != 0)
        {
            int point1 = this->verticesId[i-1];
            int point2 = this->verticesId[i];
            if(DistanceBPoints(this->vertices[point1], this->vertices[point2]) < this->t)
            {
                std::vector<int> pair {point1, point2};
                this->validPairs.push_back(pair);
                
                this->validPairs1D.push_back(point1);
                this->validPairs1D.push_back(point2);
            }
        }
        int point1 = this->verticesId[i];
        int point2 = this->verticesId[i+1];
        int point3 = this->verticesId[i+2];
        
        std::vector<int> pair1 {point1, point2};
        std::vector<int> pair2 {point2, point3};
        std::vector<int> pair3 {point1, point3};
        
        this->validPairs.push_back(pair1);
        this->validPairs.push_back(pair2);
        this->validPairs.push_back(pair3);
        
        this->validPairs1D.push_back(point1);
        this->validPairs1D.push_back(point2);
        this->validPairs1D.push_back(point2);
        this->validPairs1D.push_back(point3);
        this->validPairs1D.push_back(point1);
        this->validPairs1D.push_back(point3);
    }
    /*for (int i = 0; i < this->numVertices; i++)
    {
        std::cout << this->validPairs[i][0] << "  " << this->validPairs[i][1] << std::endl;
    }*/
}

void Decimation::ComputeQMatrices()
{
    for (int i = 0; i < this->validPairs.size(); i++)
    {
        std::vector<int> pair = validPairs[i];
        
        std::vector<std::vector<std::vector<double>>> tr1;
        std::vector<std::vector<std::vector<double>>> tr2;
        
        for (int j = 0; j < this->triangles.size(); j++)
        {   
            for (int k = 0; k < 3; k++)
            {
                if (DistanceBPoints(this->vertices[pair[0]], this->triangles[j][k]) == 0)
                {
                    tr1.push_back(this->triangles[j]);
                }
            }
        }
        
        for (int j = 0; j < this->triangles.size(); j++)
        {   
            for (int k = 0; k < 3; k++)
            {
                if (DistanceBPoints(this->vertices[pair[1]], this->triangles[j][k]) == 0)
                {
                    tr2.push_back(this->triangles[j]);
                }
            }
        }
        
        Eigen::Matrix4d K1(4, 4);           
        
        for (int j = 0; j < tr1.size(); j++)
        {
            Eigen::Matrix3d A(3, 3);
            Eigen::Vector3d B(3), x(3);
            
            A << tr1[j][0][0], tr1[j][0][1], tr1[j][0][2], tr1[j][1][0], tr1[j][1][1], tr1[j][1][2], tr1[j][2][0], tr1[j][2][1], tr1[j][2][2];
            B << -1, -1, -1;
            x = A.colPivHouseholderQr().solve(B);
            double a = x[0];
            double b = x[1];
            double c = x[2];
            double d = 1;
            
            Eigen::Matrix4d k1(4, 4);  
            k1 << a*a, a*b, a*c, a*d, a*b, b*b, b*c, b*d, a*c, b*c, c*c, c*d, d*a, d*b, d*c, d*d;
            K1 = K1+ k1;       
        }

        Eigen::Matrix4d K2(4, 4);           
        
        for (int j = 0; j < tr2.size(); j++)
        {
            Eigen::Matrix3d A(3, 3);
            Eigen::Vector3d B(3), x(3);
            
            A << tr2[j][0][0], tr2[j][0][1], tr2[j][0][2], tr2[j][1][0], tr2[j][1][1], tr2[j][1][2], tr2[j][2][0], tr2[j][2][1], tr2[j][2][2];
            B << -1, -1, -1;
            x = A.colPivHouseholderQr().solve(B);
            double a = x[0];
            double b = x[1];
            double c = x[2];
            double d = 1;
            
            Eigen::Matrix4d k2(4, 4);  
            k2 << a*a, a*b, a*c, a*d, a*b, b*b, b*c, b*d, a*c, b*c, c*c, c*d, d*a, d*b, d*c, d*d;
            K2= K2+ k2;  
        }
        
        Eigen::Vector4d v1, v2;
        v1 << this->vertices[pair[0]][0], this->vertices[pair[0]][1], this->vertices[pair[0]][2], 1;
        v1 << this->vertices[pair[1]][0], this->vertices[pair[1]][1], this->vertices[pair[1]][2], 1;
        double q1 = v1.transpose() * K1 * v1;
        double q2 = v2.transpose() * K2 * v2;
        double q = q1 + q2;
        this->qmatrices.push_back(q);
        
    }
}
