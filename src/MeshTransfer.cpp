#include "MeshTransfer.h"
#include <TriMesh.h>

void MeshTransfer::setStationaryVertices(const std::vector<int>& sv){
    m_stationary_vertices = sv;
}

void MeshTransfer::setSource(trimesh::TriMesh _S0){
    m_S0 = _S0;
    computeS0grad();
}

void MeshTransfer::setTarget(trimesh::TriMesh _T0){
    m_T0 = _T0;
    computeT0grad();
}


void MeshTransfer::computeS0grad(){
    int numFaces = m_S0.faces.size();
    m_S0grad.resize(numFaces);
    m_Ds.resize(numFaces);
    
    for(int i=0; i<numFaces; i++){
        pair<Eigen::Matrix3d, double> tmp = triangleGradient2(m_S0, i);
        m_S0grad[i] = tmp.first;
        m_Ds[i]     = tmp.second;
    }
}

void MeshTransfer::computeT0grad(){
    int numFaces = m_T0.faces.size();
    m_T0grad.resize(numFaces);
    
    for(int i=0; i<numFaces; i++){
        Eigen::Matrix3d TV = triangleGradient(m_T0, i);
        Eigen::Matrix3d TVinv = TV.inverse();
        
        Eigen::Vector3d s;
        s[0] = -TVinv(0, 0) - TVinv(1, 0);
        s[1] = -TVinv(0, 1) - TVinv(1, 1);
        s[2] = -TVinv(0, 2) - TVinv(1, 2);
        
        m_T0grad[i] << s[0], TVinv(0,0), TVinv(1,0),
                       s[1], TVinv(0,1), TVinv(1,1),
                       s[2], TVinv(0,2), TVinv(1,2);
    }
}


trimesh::TriMesh MeshTransfer::transfer(trimesh::TriMesh _S1){  
    int numFaces = _S1.faces.size();
    
    // compute the deformation gradients of S1
    // G_{s->t} : M_t * M_s^-1
    std::vector<Eigen::Matrix3d> S1grad(numFaces);
    for(int i=0; i<numFaces; i++){
        Eigen::Matrix3d V = m_S0grad[i];
        Eigen::Matrix3d Vt = triangleGradient(_S1, i);
        S1grad[i] = (Vt * V.inverse()).transpose();
    }
    
    // delegate to transfer by deformation gradients
    return transfer(S1grad);
}


/* reconstruct the blendshape vertices based on the local frame
 * @input:
 *      S1grad = G_{s0->s1}
 */
// G_{t0->t1} = G_{s0->s1}  
trimesh::TriMesh MeshTransfer::transfer(const std::vector<Eigen::Matrix3d>& S1grad, bool sameFlag)
{    
    
    int numFaces = m_S0.faces.size();
    int numVerts = m_S0.vertices.size();
    
    // assemble sparse matrix A
    int nrowsA = numFaces * 3;
    int nsv = m_stationary_vertices.size();
    int nrowsC = nsv;
    int nrows = nrowsA + nrowsC;
    int ncols = numVerts;
    
//     debug("numFaces", numFaces);
//     debug("numVerts", numVerts);
//     debug("nsv", nsv);
//     debug("nrows", nrows);
//     debug("ncols", ncols);
    
    using Tripletd = Eigen::Triplet<double>;
    using SparseMatrixd = Eigen::SparseMatrix<double>;  //Eigen::SparseMatrix<double, Eigen::RowMajor>;//
    
    std::vector<Tripletd> A_coeffs;
    A_coeffs.reserve(9 * numFaces + nrowsC);
    

    
    SparseMatrixd A(nrows, ncols);
    for(int i=0, ioffset=0; i<numFaces; i++){
        trimesh::TriMesh::Face f = m_S0.faces[i];
        Eigen::Matrix3d Ti = m_T0grad[i]; 
        
        A_coeffs.push_back(Tripletd(ioffset, f[0], Ti(0,0)));
        A_coeffs.push_back(Tripletd(ioffset, f[1], Ti(0,1)));
        A_coeffs.push_back(Tripletd(ioffset, f[2], Ti(0,2)));
        ioffset++;
        
        A_coeffs.push_back(Tripletd(ioffset, f[0], Ti(1,0)));
        A_coeffs.push_back(Tripletd(ioffset, f[1], Ti(1,1)));
        A_coeffs.push_back(Tripletd(ioffset, f[2], Ti(1,2)));
        ioffset++;
        
        A_coeffs.push_back(Tripletd(ioffset, f[0], Ti(2,0)));
        A_coeffs.push_back(Tripletd(ioffset, f[1], Ti(2,1)));
        A_coeffs.push_back(Tripletd(ioffset, f[2], Ti(2,2)));
        ioffset++;
    }
    
    const double w_stationary = 0.1;
    for(int i=0; i<nsv; i++){
        A_coeffs.push_back(Tripletd(nrowsA + i, m_stationary_vertices[i], w_stationary));
    }
    A.setFromTriplets(A_coeffs.begin(), A_coeffs.end());
    

    // fill C
    Eigen::MatrixXd C(nrows, 3);
    for(int i=0; i<3; i++){
        for(int j=0, joffset=0; j<numFaces; j++){
            Eigen::Matrix3d Sj = S1grad[j];
            C(joffset, i) = Sj(0, i); joffset++;
            C(joffset, i) = Sj(1, i); joffset++;
            C(joffset, i) = Sj(2, i); joffset++;
        }
    }
    for(int i=0; i<3; i++){
        for(int j=0; j<nsv; j++){
            trimesh::vec3 vj = m_T0.vertices[j];
            C(nrowsA+j, i) = vj[i] * w_stationary;
        }
    }
    
    // fill D
    std::vector<Tripletd> D_coeffs;
    D_coeffs.reserve(nrows);
    for(int i=0; i<nrowsA; i++){
        int fidx = i / 3;
        D_coeffs.push_back(Tripletd(i, i, m_Ds[fidx]));
    }
    for(int i=nrowsA; i<nrows; i++){
        D_coeffs.push_back(Tripletd(i,i,1));
    }
    SparseMatrixd D(nrows, nrows);
    D.setFromTriplets(D_coeffs.begin(), D_coeffs.end());
    

    SparseMatrixd At = A.transpose();    
    SparseMatrixd AtD = At * D;
    SparseMatrixd AtA = (At * A).pruned();
    SparseMatrixd AtDA = (AtD * A).pruned();
    Eigen::MatrixXd AtDC = AtD * C;
    Eigen::MatrixXd AtC = At * C;
    
        
        
    const double epsilon = 1e-6;
    Eigen::SparseMatrix<double> eye(ncols, ncols);
    for(int j=0; j<ncols; j++) eye.insert(j,j) = epsilon;
    AtA += eye;
    

    //solve for GtDG \ GtDC
    Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> solver;
    //cout<<sameFlag<<endl;
    if(sameFlag) solver.compute(AtA);
    else solver.compute(AtDA);
    
    if(solver.info() != Eigen::Success) {exit(-1);}
    
    
    Eigen::MatrixXd x;
    if(sameFlag) x = solver.solve(AtC);
    else x = solver.solve(AtDC);
    
    if(solver.info() != Eigen::Success) {exit(-1);}
    
    
    trimesh::TriMesh ret = m_T0;
    for(int i=0; i<numVerts; i++){
        ret.vertices[i] = trimesh::vec3(x(i,0),x(i,1), x(i,2));
    }
    
    //ret.write("ret.ply");
    return ret;
}



