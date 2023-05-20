#include "decimation.h"
#include <iostream>
#include <cmath>
#include <utility>
#include <algorithm>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkUnstructuredGridWriter.h>

#include <xtensor/xtensor.hpp>
#include <xtensor/xio.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <Eigen/Dense>


bool Sorting(std::pair<int, double> &a, std::pair<int, double> &b)
{   
    return (a.second < b.second);
}


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
    for (int i = 0; i < this->numVertices; i++)
    {
        double vertex[3];
        points->GetPoint(i, vertex);
        std::vector<double> point;
        for (int j = 0; j < 3; j++)
        {
            point.push_back(vertex[j]);
        }
        
        this->vertices.push_back(point);
    }
    
    /*for (int i = 0; i < this->numVertices; i++)
    {
        std::vector<double> point = this->vertices[i];
        std::cout << point[0] << " " << point[1] << " " << point[2] << std::endl;
    }*/
    
    vtkSmartPointer<vtkCellArray> cells = unstructuredGrid->GetCells();
    int ncells = cells->GetNumberOfCells();
    int n = 0;
    vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
    for (int i = 0; i < ncells; i++) 
    {
        cells->GetCellAtId(i, idList);
        int nverts = idList->GetNumberOfIds();
        for (int j=0; j<nverts; j++) 
        {
            this->verticesId.push_back(idList->GetId(j));
        }
    }
    
    for (int i = 0; i < this->verticesId.size(); i=i+3)
    {
        std::vector<std::vector<double>> triang;
        std::vector<int> id;
        triang.push_back(this->vertices[this->verticesId[i]]);
        id.push_back(this->verticesId[i]);
        triang.push_back(this->vertices[this->verticesId[i+1]]);
        id.push_back(this->verticesId[i+1]);
        triang.push_back(this->vertices[this->verticesId[i+2]]);
        id.push_back(this->verticesId[i+2]);
        this->triangles.push_back(triang);
        this->trId.push_back(id);
    }  
};

void Decimation::SelectValidPairs(double treshold)
{
    this->t = treshold;
    for (int i = 0; i < this->verticesId.size(); i = i + 3)
    {
        if (i != 0)
        {
            int point1 = this->verticesId[i-1];
            int point2 = this->verticesId[i];
            if(DistanceBPoints(this->vertices[point1], this->vertices[point2]) < this->t)
            {
                std::vector<int> pair {point1, point2};
                this->validPairs.push_back(pair);
                
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
        
    }
    /*for (int i = 0; i < this->validPairs.size(); i++)
    {
        std::vector<double> p1 = this->vertices[this->validPairs[i][0]];
        std::vector<double> p2 = this->vertices[this->validPairs[i][1]];
        std::cout << DistanceBPoints(p1, p2) << std::endl;
        //std::cout << this->validPairs[i][0] << "  " << this->validPairs[i][1] << std::endl;
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
            
            A << tr1[j][0][0], tr1[j][0][1], tr1[j][0][2], 
                 tr1[j][1][0], tr1[j][1][1], tr1[j][1][2], 
                 tr1[j][2][0], tr1[j][2][1], tr1[j][2][2];
            B << -1, -1, -1;
            x = A.colPivHouseholderQr().solve(B);
            double a = x[0];
            double b = x[1];
            double c = x[2];
            double d = 1;
            
            Eigen::Matrix4d k1(4, 4);  
            k1 << a*a, a*b, a*c, a*d, 
                  a*b, b*b, b*c, b*d, 
                  a*c, b*c, c*c, c*d, 
                  d*a, d*b, d*c, d*d;
            K1 = K1+ k1;       
        }

        Eigen::Matrix4d K2(4, 4);           
        
        for (int j = 0; j < tr2.size(); j++)
        {
            Eigen::Matrix3d A(3, 3);
            Eigen::Vector3d B(3), x(3);
            
            A << tr2[j][0][0], tr2[j][0][1], tr2[j][0][2], 
                 tr2[j][1][0], tr2[j][1][1], tr2[j][1][2], 
                 tr2[j][2][0], tr2[j][2][1], tr2[j][2][2];
            B << -1, -1, -1;
            x = A.colPivHouseholderQr().solve(B);
            double a = x[0];
            double b = x[1];
            double c = x[2];
            double d = 1;
            
            Eigen::Matrix4d k2(4, 4);  
            k2 << a*a, a*b, a*c, a*d,
                  a*b, b*b, b*c, b*d, 
                  a*c, b*c, c*c, c*d, 
                  d*a, d*b, d*c, d*d;
            K2= K2+ k2;  
        }
        
        Eigen::Vector4d v1, v2;
        v1 << this->vertices[pair[0]][0], 
              this->vertices[pair[0]][1], 
              this->vertices[pair[0]][2], 1;
        v1 << this->vertices[pair[1]][0], 
              this->vertices[pair[1]][1], 
              this->vertices[pair[1]][2], 1;
        double q1 = v1.transpose() * K1 * v1;
        double q2 = v2.transpose() * K2 * v2;
        double q = q1 + q2;
        std::pair <int, double> p = std::make_pair(i, q);
        this->qmatrices.push_back(p);
        
    }
}

void Decimation::SortValidPairs()
{
    std::sort(this->qmatrices.begin(), this->qmatrices.end(), Sorting);
}

void Decimation::RemoveVertices(double param)
{
    std::vector<std::vector<double>> vert;
    for (int i = 0; i < this->numVertices; i ++)
    {
        vert.push_back(this->vertices[i]);
    }
    
    /*for (int i = 0; i < this->verticesId.size(); i++)
    {
        this->newVerticesId.push_back(this->verticesId[i]);
    }*/

    this->newNumVertices = int(this->numVertices * (1 - param));
    int removeVert = this->numVertices - this->newNumVertices;
    
    if (this->validPairs.size() < removeVert)
    {
        removeVert = this->validPairs.size();
    }
    
    removeVert = 1;
    
    std::vector<int> removedId;
    
    for (int i = 0; i < removeVert; i++)
    {
        int ind1 = this->validPairs[std::get<0>(this->qmatrices[i])][0];
        int ind2 = this->validPairs[std::get<0>(this->qmatrices[i])][1];
        
        vert[ind2] = vert[ind1];
        removedId.push_back(ind2);
        
        for (int j = 0; j < this->trId.size(); j++)
        {
            std::vector<int> tr = this->trId[j];
            if ((std::find(tr.begin(), tr.end(), ind1) != tr.end()) && (std::find(tr.begin(), tr.end(), ind2) != tr.end()))
            {
                this->trId.erase(this->trId.begin() + j);
            }
        }
        for (int j = 0; j < this->trId.size(); j++)
        {
            for (int k = 0; k < 3; k++)
            {
                if (this->trId[j][k] == ind2)
                {
                    trId[j][k] = ind1;
                }
            }
        }
    }
    
    for (int i = 0; i < this->trId.size(); i++)
    {
        this->newVerticesId.push_back(this->trId[i][0]);
        this->newVerticesId.push_back(this->trId[i][1]);
        this->newVerticesId.push_back(this->trId[i][0]);
    }

    /*for (int i = 0; i < this->newVerticesId.size(); i ++)
    {
        if (std::find(removedId.begin(), removedId.end(), this->newVerticesId[i]) != removedId.end())
        {
            this->newVerticesId.erase(this->newVerticesId.begin() + i);
            i = i - 1;
        }
    }

    for (int i = 0; i < this->newVerticesId.size(); i ++)
    {
        std::cout << this->newVerticesId[i] << std::endl;
    }
    
    for (int i = 0; i < vert.size(); i++)
    {
        if (std::find(this->newVerticesId.begin(), this->newVerticesId.end(), i) != this->newVerticesId.end())
        {
            //this->newVertices.push_back(vert[i]);
        }
        else
        {
            vert.erase(vert.begin() + i);
            for (int j = 0; j < this->newVerticesId.size(); j++)
            {
                //std::cout << "hi" << std::endl;
                if (this->newVerticesId[j] > i)
                {
                    this->newVerticesId[j] -= 1;
                }
            }
        }
    }*/
    
    for (int i = 0; i < vert.size(); i++)
    {
        this->newVertices.push_back(vert[i]);
    }

    this->newNumVertices = this->newVertices.size();
    
    /*for (int i = 0; i < this->newVerticesId.size(); i++)
    {
        while (this->newVerticesId[i] == this->newVerticesId[i+1])
        {
            this->newVerticesId.erase(this->newVerticesId.begin() + i + 1);
            if(i+1 == this->newVerticesId.size())
            {
                break;
            }
        }
    }*/
}

void Decimation::WriteNewVtk(std::string filename, std::string filename_out)
{
    vtkSmartPointer<vtkPoints> newVtkVertices = vtkSmartPointer<vtkPoints>::New();
    newVtkVertices->SetNumberOfPoints(this->newNumVertices);
    for (int i = 0; i < this->newNumVertices; i++)
    {   
        double coords[3];
        for (int j = 0; j < 3; j++)
        {
            coords[j] = this->newVertices[i][j];
        }
        newVtkVertices->SetPoint(i, coords);
        
    }
    
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    for (int i = 0; i < this->newVerticesId.size(); i = i + 3)
    {
        vtkIdType pointIds[3];
        pointIds[0] = this->newVerticesId[i];
        pointIds[1] = this->newVerticesId[i+1];
        pointIds[2] = this->newVerticesId[i+2];
        cells->InsertNextCell(3, pointIds);
    }
    
    /*for (int i = 0; i < this->newVerticesId.size(); i++)
    {
        std::cout << this->newVerticesId[i] << std:: endl;
    }
    std::cout << this->newVerticesId.size() << std:: endl;*/
    
   vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    
    vtkSmartPointer<vtkUnstructuredGrid> uGrid = reader->GetOutput();
    
    uGrid->SetPoints(newVtkVertices);
    uGrid->SetCells(VTK_TRIANGLE, cells);


    vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetInputData(uGrid);
    writer->SetFileName(filename_out.c_str());
    writer->Write();
}
