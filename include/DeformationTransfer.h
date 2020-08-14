// Deformation Transfer
#ifndef DEFORMATION_TRANSFER
#define DEFORMATION_TRANSFER

#include "Common.h"

typedef Eigen::Matrix<double, -1,  1> Vec;
typedef Eigen::Matrix<double, -1, -1> Mat;
typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMat;
typedef Eigen::Vector3f Float3;
typedef Eigen::Matrix3f Mat3f;
typedef Eigen::Vector3i Int3;
typedef Eigen::Vector4i Int4;

// Input:  A0, A1, B0
// Output: B1
// Where Ai, Bi are triangle meshes with the same topology
class DeformationTransfer{
public:
    DeformationTransfer(){}
    ~DeformationTransfer(){}
    
    void init(const trimesh::TriMesh& _A0,
              const trimesh::TriMesh& _B0);
    
    // given A1, to solve B1
    bool transfer(const trimesh::TriMesh& _A1,
                  trimesh::TriMesh& _B1);
    
    
private:
    void setupAnchorMat();
    void setupAnchorRhs();
    
    void setupRegularizationMat();
    void setupRegularizationRhs();
    
    void setupE1Mat();
    void setupE1Rhs(const trimesh::TriMesh& A1);
    
    void saveResultMesh(const Eigen::Matrix<double, -1, 1>& m_x,
                                  trimesh::TriMesh& _B1);
    
    Eigen::Matrix3d triangleGradient(const trimesh::TriMesh& mesh,
                                    const int& fidx);

    
private:
    bool m_initFlag;
    trimesh::TriMesh m_A0, m_B0; // A0 and B0 share the same topology
    int m_numVerts, m_numFaces, m_numTotal;
    std::vector<int> m_anchors;  //index of anchor points
    
    
    //weigths
    double w_anchor;
    double w_reg;
    double w_1;
    
    // for anchor
    SpMat m_anchorMat;
    SpMat m_anchorMatT;
    Eigen::Matrix<double, -1, 1> m_anchorRhs;
    
    // regularization
    SpMat m_regAtA;
    Eigen::Matrix<double, -1, 1> m_regAtb;
    Eigen::Matrix<double, -1, 1> m_anchorRegSumAtb;
    
    // energy for source-target triangle correspondences
    SpMat m_E1Mat, m_E1MatT;
    Eigen::Matrix<double, -1, 1> m_E1Rhs;
    
    // total energy 
    SpMat m_AtA;
    Eigen::Matrix<double, -1, 1> m_Atb;
    Eigen::Matrix<double, -1, 1> m_x;
    
    // solver
    Eigen::SimplicialCholesky<SpMat> m_solver;
};



#endif
