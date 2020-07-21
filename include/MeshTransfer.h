#ifndef MESH_TRANSFER_H
#define MESH_TRANSFER_H

#include "Common.h"
#include "TriangleGradient.h"

#include "Utils.h"

class MeshTransfer{
    
public:
    MeshTransfer() {}
    ~MeshTransfer(){}
    
    void setSource(trimesh::TriMesh _S0);
    void setTarget(trimesh::TriMesh _T0);
    void setStationaryVertices(const std::vector<int>& sv);

    // using a target shape
    trimesh::TriMesh transfer(trimesh::TriMesh _S1);
    
    // using a per-face deformation gradient
    trimesh::TriMesh transfer(const std::vector<Eigen::Matrix3d>& S1grad,  bool sameFlag = false);
    
private:
    void computeS0grad();
    void computeT0grad();
    
private:
    trimesh::TriMesh             m_S0, m_T0;
    std::vector<Eigen::Matrix3d> m_S0grad, m_T0grad;
    std::vector<double>          m_Ds;
    std::vector<int>             m_stationary_vertices;
};

#endif
