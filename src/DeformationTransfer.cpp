#include "DeformationTransfer.h"

void DeformationTransfer::init(const trimesh::TriMesh& _A0,
                               const trimesh::TriMesh& _B0)
{
    m_A0 = _A0;
    m_B0 = _B0;      
    
    m_numVerts = m_A0.vertices.size();
    m_numFaces = m_A0.faces.size();
    m_numTotal  = m_numVerts + m_numFaces;
    
    cout<<"numVerts:"<<m_numVerts<<"  numFaces:"<<m_numFaces<<endl;
    
    // simply use point 0 as anchor points for global position
    m_anchors.clear();
    //m_anchors.push_back(0);
    
    for(int i=0; i<m_A0.vertices.size(); i++){
        if(m_A0.is_bdy(i)) m_anchors.push_back(i);
    }
    cout<<"# anchor: "<<m_anchors.size()<<endl;
    
    
    //pre-compute
    setupAnchorMat();
    setupAnchorRhs();
    
    setupRegularizationMat();
    setupRegularizationRhs();
    
    setupE1Mat();
    
    /**********************************************
     * parameters
    ***********************************************/
//  const double Transfer_Weight_Smoothness = 1e-3;
// 	const double Transfer_Weight_Identity   = 1e-6;
// 	const double Transfer_Weight_Correspond = 1.0;
// 	const double Transfer_Weight_Anchor     = 1e8;
// 	const double Transfer_Weight_Regularization = 1e-8;

//     w_anchor = Transfer_Weight_Anchor         / (1e-3f + m_anchorMat.rows());
//     w_reg    = Transfer_Weight_Regularization / (1e-3f + m_regAtA.rows());
//     w_corres = Transfer_Weight_Correspond     / (1e-3f + m_E1Mat.rows());

    w_anchor = 0; 
    w_reg    = 0; 
    w_corres = 1; 
    
    cout<<"Parameters = [anchor: " << w_anchor <<  "]  [reglar: " << w_reg<< "]  [corres: " << w_corres <<"]"<<endl;
    

    
    m_anchorMatT = m_anchorMat.transpose();  
    m_E1MatT = m_E1Mat.transpose();
    
    m_AtA = w_corres * m_E1MatT * m_E1Mat + w_anchor * m_anchorMatT * m_anchorMat  + w_reg * m_regAtA;
    m_solver.compute(m_AtA);

    
    m_anchorRegSumAtb = w_anchor * m_anchorMatT * m_anchorRhs + w_reg * m_regAtb;
    
    m_initFlag = true;
    
    return;
}


bool DeformationTransfer::transfer(const trimesh::TriMesh& _A1, 
                                   trimesh::TriMesh& _B1)
{
    
    if(!m_initFlag) return false;
    
    // computing all energy matrices
    setupE1Rhs(_A1);
    
    //// sum all the energy terms
    m_Atb = m_anchorRegSumAtb + w_corres * m_E1MatT * m_E1Rhs;
    
    //solve
    m_x = m_solver.solve(m_Atb);
    
    
    _B1 = m_B0;
    //save results
    saveResultMesh(m_x, _B1);
    
    return true;
}

void DeformationTransfer::saveResultMesh(const Eigen::Matrix<double, -1, 1>& m_x,
                                  trimesh::TriMesh& _B1)
{
    //cout<<m_x.size()<<endl;
    for(int i=0; i<m_numVerts; i++){
        for(int k=0; k<3; k++){
            _B1.vertices[i][k] = m_x[k * (m_x.size()/3) + i];
        }
    }
    _B1.write("ret_deformTransfer.ply");
}


void DeformationTransfer::setupAnchorMat(){
    m_anchorMat.resize(3*m_anchors.size(), 3*m_numTotal);
    
    for(int i=0; i<m_anchors.size(); i++){
        int idx = m_anchors[i];
        for(int k=0; k<3; k++){
            m_anchorMat.insert(i*3+k, idx + k*m_numTotal) = 1;
        }
    }
    m_anchorMat.finalize();
}


void DeformationTransfer::setupAnchorRhs(){
     m_anchorRhs.resize(3*m_anchors.size());
     m_anchorRhs.setZero();
     
     for(int i=0; i<m_anchors.size(); i++){
         int idx = m_anchors[i];
         for(int k=0; k<3; k++){
             m_anchorRhs[i*3+k] = m_B0.vertices[idx][k];
        }
    }
}

void DeformationTransfer::setupRegularizationMat(){
    m_regAtA.resize(3*m_numTotal, 3*m_numTotal);
    m_regAtA.reserve(3*m_numTotal);
    
    for(int i=0; i<3*m_numTotal; i++){
        m_regAtA.insert(i,i) = 1;
    }
    m_regAtA.finalize();
}


Eigen::Matrix<double, 3, 4> triangleLocalFrame(
    const trimesh::TriMesh& mesh, const int& fidx)
{
	trimesh::TriMesh::Face f = mesh.faces[fidx];
	
	trimesh::point v0 = mesh.vertices[f[0]];
	trimesh::point v1 = mesh.vertices[f[1]];
	trimesh::point v2 = mesh.vertices[f[2]];
	
	Eigen::Vector3d v1v0 = Eigen::Vector3d(v1.x-v0.x, v1.y-v0.y, v1.z-v0.z);
	Eigen::Vector3d v2v0 = Eigen::Vector3d(v2.x-v0.x, v2.y-v0.y, v2.z-v0.z);
	
	Eigen::Vector3d n = v1v0.cross(v2v0);
	Eigen::Vector3d nn = n.normalized();
	
	Eigen::Matrix<double, 3, 4> G;

	G << v0[0], v1[0], v2[0], v0[0] + nn[0],
		 v0[1], v1[1], v2[1], v0[1] + nn[1],
		 v0[2], v1[2], v2[2], v0[2] + nn[2];

	return G;
}


void DeformationTransfer::setupRegularizationRhs(){
    m_regAtb.resize(3*m_numTotal);
    m_regAtb.setZero();
    
    for(int i=0; i<m_numVerts; i++){
        for(int k=0; k<3; k++){ 
            m_regAtb[k * m_numTotal + i] = m_B0.vertices[i][k];
        }
    }
    
    for(int i=0; i<m_numFaces; i++){
        Eigen::Matrix<double, 3, 4> localFrame = triangleLocalFrame(m_B0, i);
        for(int k=0; k<3; k++){
            m_regAtb[k * m_numTotal + m_numVerts + i] = localFrame(k, 3);
        }
    }
}



Eigen::Matrix3d DeformationTransfer::triangleGradient(
    const trimesh::TriMesh& mesh, const int& fidx)
{
	trimesh::TriMesh::Face f = mesh.faces[fidx];
	
	trimesh::point v0 = mesh.vertices[f[0]];
	trimesh::point v1 = mesh.vertices[f[1]];
	trimesh::point v2 = mesh.vertices[f[2]];
	
	Eigen::Vector3d v1v0 = Eigen::Vector3d(v1.x-v0.x, v1.y-v0.y, v1.z-v0.z);
	Eigen::Vector3d v2v0 = Eigen::Vector3d(v2.x-v0.x, v2.y-v0.y, v2.z-v0.z);
	
	Eigen::Vector3d n = v1v0.cross(v2v0);
	Eigen::Vector3d nn = n.normalized();
	
	Eigen::Matrix3d G;

	G << v1v0[0], v2v0[0], nn[0],
		 v1v0[1], v2v0[1], nn[1],
		 v1v0[2], v2v0[2], nn[2];

	return G;
}



// G
// |v1v0[0], v2v0[0], nn[0]|
// |v1v0[1], v2v0[1], nn[1]|
// |v1v0[2], v2v0[2], nn[2]|
inline void getMatrixT(const Eigen::Matrix3d& G,
                       Eigen::Matrix<double, 3, 4>& A)
{
    // G inverse()
    Eigen::Matrix3d V = G.inverse();
    //|-V(0, 0) - V(1, 0) - V(2, 0),   V(0, 0),   V(1, 0),   V(2, 0)|
    //|-V(0, 1) - V(1, 1) - V(2, 1),   V(0, 1),   V(1, 1),   V(2, 1)|
    //|-V(0, 2) - V(1, 2) - V(2, 2),   V(0, 2),   V(1, 2),   V(2, 2)|
    
    A(0, 0) = -V(0, 0) - V(1, 0) - V(2, 0);
	A(0, 1) =  V(0, 0);
	A(0, 2) =  V(1, 0);
	A(0, 3) =  V(2, 0);
    
	A(1, 0) = -V(0, 1) - V(1, 1) - V(2, 1);
	A(1, 1) =  V(0, 1);
	A(1, 2) =  V(1, 1);
	A(1, 3) =  V(2, 1);
    
	A(2, 0) = -V(0, 2) - V(1, 2) - V(2, 2);
	A(2, 1) =  V(0, 2);
	A(2, 2) =  V(1, 2);
	A(2, 3) =  V(2, 2);
}


//  Eigen::Matrix3d TVinv = TV.inverse();
//         Eigen::Vector3d s;
//         s[0] = -TVinv(0, 0) - TVinv(1, 0);
//         s[1] = -TVinv(0, 1) - TVinv(1, 1);
//         s[2] = -TVinv(0, 2) - TVinv(1, 2);
//         
//         m_T0grad[i] << s[0], TVinv(0,0), TVinv(1,0),
//                        s[1], TVinv(0,1), TVinv(1,1),
//                        s[2], TVinv(0,2), TVinv(1,2);



// The matrix T is in block diag style:
// | A 0 0 |
// | 0 A 0 |
// | 0 0 A |
// where each A is a 3x4 matrix  
inline void fillCooSys(std::vector<Eigen::Triplet<double> >& cooSys,
                       const Eigen::Matrix<double, 3, 4>& T,
                       const int& row, const int& m_numTotal,
                       const vector<int>& id)
{
    int pos = row * 4;
    for(int k = 0; k < 3; k++) // 3
    {  
        int yb = k * 3;        
        for(int y = 0; y < 3; y++)   // 3
        {
            for(int x = 0; x<4; x++) // 4
            {
                int col = m_numTotal * k + id[x];  
                cooSys[pos++] = Eigen::Triplet<double>(row+yb+y, col, T(y, x));
            }
        }
    }
    return;
}


// gradient transfer for mesh B0
// Ax=b
// A(m_numFaces * 9, m_numTotal * 3)
void DeformationTransfer::setupE1Mat(){
    std::vector<Eigen::Triplet<double> > cooSys;
    cooSys.resize(m_numFaces * 9 * 4);   
    
    Eigen::Matrix<double, 3, 4> Ti;
    Ti.setZero();
    
    for(int fidx=0; fidx<m_numFaces; fidx++){
        Eigen::Matrix3d G = triangleGradient(m_B0, fidx);
        getMatrixT(G, Ti); //Ti (3, 4)
        
        trimesh::TriMesh::Face f = m_B0.faces[fidx];
        vector<int> pidx = {f[0], f[1], f[2], m_numVerts + fidx};
        
        const int row = fidx * 9;
        fillCooSys(cooSys, Ti, row, m_numTotal, pidx);
    }
    
    
    //fill matrix A
    m_E1Mat.resize(m_numFaces * 9, m_numTotal * 3);
    if(cooSys.size()>0){  
        m_E1Mat.setFromTriplets(cooSys.begin(), cooSys.end());  
    }
}



void DeformationTransfer::setupE1Rhs(const trimesh::TriMesh& A1){
    m_E1Rhs.resize(m_numFaces * 9);
    
    Eigen::Matrix<double, 3, 4> Si_A;
    Si_A.setZero();
    Eigen::Matrix<double, 4, 1> Si_x[3];
    Eigen::Matrix<double, 3, 1> Si_b[3];
    
    for(int fidx = 0; fidx < m_numFaces; fidx++){
        Eigen::Matrix3d G0 = triangleGradient(m_A0, fidx);        
        getMatrixT(G0, Si_A);
        
        Eigen::Matrix<double, 3, 4> localFrame = triangleLocalFrame(A1, fidx);
        
        for(int k=0; k<4; k++){
            Si_x[0][k] = localFrame(0, k);
            Si_x[1][k] = localFrame(1, k);
            Si_x[2][k] = localFrame(2, k);
        }
        
        Si_b[0] = Si_A * Si_x[0];
		Si_b[1] = Si_A * Si_x[1];
		Si_b[2] = Si_A * Si_x[2];
        
        // push matrix
		const int row = 9 * fidx;
		m_E1Rhs[row + 0] = Si_b[0][0];
		m_E1Rhs[row + 1] = Si_b[0][1];
		m_E1Rhs[row + 2] = Si_b[0][2];

		m_E1Rhs[row + 3] = Si_b[1][0];
		m_E1Rhs[row + 4] = Si_b[1][1];
		m_E1Rhs[row + 5] = Si_b[1][2];

		m_E1Rhs[row + 6] = Si_b[2][0];
		m_E1Rhs[row + 7] = Si_b[2][1];
		m_E1Rhs[row + 8] = Si_b[2][2];
    }
}




