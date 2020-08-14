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
    
//     for(int i=0; i<m_A0.vertices.size(); i++){
//         if(m_A0.is_bdy(i)) m_anchors.push_back(i);
//     }
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
//  const static double Transfer_Weight_Smoothness = 1e-3;
// 	const static double Transfer_Weight_Identity   = 1e-6;
// 	const static double Transfer_Weight_Correspond = 1.0;
// 	const static double Transfer_Weight_Anchor     = 1e8;
// 	const static double Transfer_Weight_Regularization = 1e-8;

    w_anchor = 1;//Transfer_Weight_Anchor         / (1e-3f + m_anchorMat.rows());
    w_reg    = 0;//Transfer_Weight_Regularization / (1e-3f + m_regAtA.rows());
    w_1      = 1;//Transfer_Weight_Correspond     / (1e-3f + m_E1Mat.rows());
    
    cout<<"Parameters = [anchor: " << w_anchor <<  "]  [reglar: " << w_reg<< "]  [corres: " << w_1 <<"]"<<endl;
    

    
    m_anchorMatT = m_anchorMat.transpose();  
    m_E1MatT = m_E1Mat.transpose();
    
    m_AtA = w_1 * m_E1MatT * m_E1Mat + w_anchor * m_anchorMatT * m_anchorMat  + w_reg * m_regAtA;
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
    m_Atb = m_anchorRegSumAtb + w_1 * m_E1MatT * m_E1Rhs;
    
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
    _B1.write("output.ply");
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



void DeformationTransfer::setupRegularizationRhs(){
    m_regAtb.resize(3*m_numTotal);
    m_regAtb.setZero();
    
    for(int i=0; i<m_numVerts; i++){
        for(int k=0; k<3; k++){ 
            m_regAtb[k * m_numTotal + i] = m_B0.vertices[i][k];
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
    
    //cout<<G<<endl;

	return G;
}



// G
// |v1v0[0], v2v0[0], nn[0]|
// |v1v0[1], v2v0[1], nn[1]|
// |v1v0[2], v2v0[2], nn[2]|
inline void getMatrixT(const Eigen::Matrix3d& G,
                       Eigen::Matrix<double, 3, 4>& A)
{
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



inline void fillCooSys(std::vector<Eigen::Triplet<double> >& cooSys,
                       const Eigen::Matrix<double, 3, 4>& T,
                       const int& row, const int& m_numTotal,
                       const vector<int>& id)
{
    // The matrix T is in block diag style:
	// | A 0 0 |
	// | 0 A 0 |
	// | 0 0 A |
	// where each A is a 3x4 matrix  
    const static int nBlocks = 3;
    const static int nPoints = 4;
    const static int nCoords = 3;
    
    int pos = row * 4;
    for(int iBlock=0; iBlock<nBlocks; iBlock++) // 3
    {  
        int yb = iBlock * nCoords;        
        for(int y = 0; y < nCoords; y++)   // 3
        {
            for(int x = 0; x<nPoints; x++) // 4
            {
                int col = m_numTotal * iBlock + id[x];  
                cooSys[pos++] = Eigen::Triplet<double>(row+yb+y, col, T(y, x));
            }
        }
    }
    return;
}


// gradient transfer for mesh B0
void DeformationTransfer::setupE1Mat(){
    std::vector<Eigen::Triplet<double> > cooSys;
    
    cooSys.resize(m_numFaces * 9 * 4);   
    
    Eigen::Matrix<double, 3, 4> Ti;
    Ti.setZero();
    
    for(int fidx=0; fidx<m_numFaces; fidx++){
        Eigen::Matrix3d G = triangleGradient(m_B0, fidx);
        getMatrixT(G, Ti);
        
        trimesh::TriMesh::Face f = m_B0.faces[fidx];
        vector<int> pidx = {f[0], f[1], f[2], m_numVerts + fidx};
        
        const int row = fidx * 9;
        fillCooSys(cooSys, Ti, row, m_numTotal, pidx);
    }
    
    //equation size
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
        
        Eigen::Matrix3d G1 = triangleGradient(A1, fidx);
        
        int idx0 = A1.faces[fidx][0];
        trimesh::point p0 =A1.vertices[idx0];  //mesh.vertices[f[0]];
        
        for(int k=0; k<4; k++){
            if(k==0){
                Si_x[0][k] = p0.x;
                Si_x[1][k] = p0.y;
                Si_x[2][k] = p0.z;
            }else{
                Si_x[0][k] = G1(0, k-1) + p0.x;
                Si_x[1][k] = G1(1, k-1) + p0.y;
                Si_x[2][k] = G1(2, k-1) + p0.z;
            }
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




