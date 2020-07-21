#include "BlendshapeGenerator.h"

const int NUM_POSES  = 20;   // number of poses
const int NUM_SHAPES = 47;   // number of blendshapes


struct PointResidual{
    PointResidual(double x, double y, double z, int idx, const std::vector<trimesh::TriMesh> &dB)
        : m_x(x), m_y(y), m_z(z), m_idx(idx), m_dB(dB){}
    
    template <typename T>
    bool operator()(const T* const alpha, T* residual) const{
        T p[3];  
        p[0] = T(0); p[1] = T(0); p[2] = T(0);
        
        for(int i=0; i<NUM_SHAPES-1; i++){
            trimesh::vec3 pt = m_dB[i].vertices[m_idx];
            p[0] += T(pt[0]) * alpha[i];
            p[1] += T(pt[1]) * alpha[i];
            p[2] += T(pt[2]) * alpha[i];
        }
        
        residual[0] = T(m_x) - p[0];
        residual[1] = T(m_y) - p[1];
        residual[2] = T(m_z) - p[2];
        return true;
    }
    
private:
    const double m_x, m_y, m_z;   // S-B0
    const int    m_idx;           // vertex index
    const std::vector<trimesh::TriMesh> &m_dB;
};


struct PriorResidual{
    PriorResidual(double* prior):m_prior(prior){}
    
    template <typename T>
    bool operator()(const T* const alpha, T* residual) const{
        for(int i=0; i<NUM_SHAPES-1; i++){
            residual[i] = T(m_prior[i]) - alpha[i];
        }
        return true;
    }
private:
    const double* m_prior;
};



std::vector<double> BlendshapeGenerator::estimateWeights(const trimesh::TriMesh& S, 
                                                         const trimesh::TriMesh& B0,
                                                         const std::vector<trimesh::TriMesh>& dB,
                                                         const std::vector<double>& w0_vec,   //46 init weight
                                                         const std::vector<double>& wp_vec,   //46 prior weight
                                                         bool  prior_weight_flag)
{    
    
    ColorStream(ColorOutput::Blue)<<"estimate weight ...";
    
    double w[NUM_SHAPES-1];
    for(int i=0; i<NUM_SHAPES-1; i++) w[i] = w0_vec[i];   
    
    double wp[NUM_SHAPES-1];
    for(int i=0; i<NUM_SHAPES-1; i++) wp[i] = wp_vec[i];
    
    ceres::Problem problem;
    
    int numVerts = S.vertices.size();
    for(int i=0; i<numVerts; i++){
        trimesh::vec3 vS  = S.vertices[i];    // S
        trimesh::vec3 vB0 = B0.vertices[i];   // B0
        double dx = vS[0] - vB0[0];           // S- B0
        double dy = vS[1] - vB0[1];
        double dz = vS[2] - vB0[2];
        
        ceres::CostFunction *costFunc = new AutoDiffCostFunction<PointResidual, 3, 46>(new PointResidual(dx, dy, dz, i, dB));
        problem.AddResidualBlock(costFunc, NULL, w);
    }
    
    if(prior_weight_flag){
        ceres::CostFunction *costFunc = new AutoDiffCostFunction<PriorResidual, 46, 46>(new PriorResidual(wp));
        problem.AddResidualBlock(costFunc, NULL, w);
    }
    for(int i=0; i<NUM_SHAPES-1; i++){
        problem.SetParameterLowerBound(w, i, 0.0);
        problem.SetParameterUpperBound(w, i, 1.0);
    }
    
    // solve 
    ceres::Solver::Options options;
    options.max_num_iterations  = 10;
	options.linear_solver_type = ceres::DENSE_QR;
	options.minimizer_progress_to_stdout = false;
    options.num_threads = 8;
    options.num_linear_solver_threads = 8;
    options.initial_trust_region_radius = 1.0;
    options.min_lm_diagonal = 1.0;
    options.max_lm_diagonal = 1.0;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
    //std::cout<<summary.BriefReport()<<std::endl;
    
    std::vector<double> ret;
    std::cout<<"    done";
    for(int i=0; i<NUM_SHAPES-1; i++){
        //std::cout<<w[i]<<" ";
        ret.push_back(w[i]);
    }
    //std::cout<<std::endl;
    return ret;
}



trimesh::TriMesh BlendshapeGenerator::reconstructMesh(const std::vector<double>& weight, 
                                                      const trimesh::TriMesh& B0,
                                                      const std::vector<trimesh::TriMesh>& dB)
{
    trimesh::TriMesh ret = B0;
    for(int i=0; i<ret.vertices.size(); i++)
    {
        for(int j=0; j<NUM_SHAPES-1; j++){
            ret.vertices[i].x += weight[j] * dB[j].vertices[i].x; 
            ret.vertices[i].y += weight[j] * dB[j].vertices[i].y;
            ret.vertices[i].z += weight[j] * dB[j].vertices[i].z;
        }
    }
    return ret;
}


/*
 * @input:
 * 
 */
std::vector<trimesh::TriMesh> BlendshapeGenerator::optimizeBlendshapes(
        const std::vector<trimesh::TriMesh>& S,                   // 20
        const std::vector<std::vector<Eigen::Matrix3d> >& Sgrad,  // 20 * 22800 * mat33
        const trimesh::TriMesh& B0,                               // B0
        const std::vector<std::vector<double> >& alpha,           // 20 * 46
        double beta, double gamma,
        const std::vector<std::vector<Eigen::Matrix3d> >& prior,  // 46 * 22800 * mat33
        const std::vector<std::vector<double> >& w_prior)         // 46 * 22800
{

    cout<<"Optimize blendshapes ..."<<endl;
    
    int nfaces = B0.faces.size();
    int nverts = B0.vertices.size();
    
    cout<<"nfaces: "<<nfaces<<"     nverts: "<<nverts<<endl;
    
    int nposes  = S.size();     // 20
    int nshapes = prior.size(); // 46
    
    
    // the triangle gradient of B0 
    std::vector<Eigen::Matrix3d> B0grad(nfaces);
    //std::vector<double> B0d(nfaces);
    for(int i=0; i<nfaces; i++){
//         std::pair<Eigen::Matrix3d, double> GD = triangleGradient2(B0, i);
//         B0grad[i] = GD.first;
//         B0d[i] = GD.second;
        B0grad[i] = triangleGradient(B0, i);
    }
    
    
    using Tripletd = Eigen::Triplet<double>;
    using SparseMatrixd = Eigen::SparseMatrix<double, Eigen::RowMajor>;
    std::vector<Tripletd> Adata_coeffs;

    // upper part of A
    // sigma (alpha_ij * M^B_i) = M^S_j - M^B_0 
    // [npose * 9, nshapes * 9]
    for(int i=0; i<nposes; i++){
        int rowoffset = 9 * i;
        for(int j=0; j<nshapes; j++){   // 46
            int coloffset = 9 * j;
            for(int k=0; k<9; k++){
                Adata_coeffs.push_back(Tripletd(rowoffset+k, coloffset+k, alpha[i][j]));
            }
        }
    }
    
    
    int nrows = (nposes + nshapes) * 9;
    int ncols = nshapes * 9;
    
    std::vector<std::vector<Eigen::Matrix3d> > M(nfaces);
    //#pragma omp parallel for
    for(int j=0; j<nfaces; j++){
        if(j%1000 == 0) std::cout<<j<<std::endl<<std::flush;
        //cout<<"face "<<j<<endl;    
        
        Eigen::VectorXd b(nrows);
        
        Eigen::Matrix3d B0j = B0grad[j];
        for(int i=0; i<nposes; i++){
            Eigen::Matrix3d Sgrad_ij = Sgrad[i][j]; //20 * 22800 * mat33
            for(int k=0; k<9; k++){
                b(i * 9 + k) = Sgrad_ij(k/3, k%3) - B0j(k/3, k%3);
            }
        }
        
        // lower part of A
        std::vector<Tripletd> A_coeffs = Adata_coeffs;
        for(int i=0; i<nshapes; i++){
            int rowoffset = (nposes + i) * 9;
            int coloffset = i * 9;
            for(int k=0; k<9; k++){
                A_coeffs.push_back(Tripletd(rowoffset+k, coloffset+k, beta * w_prior[i][j]));
            }
        }
        
        // lower part of B
        for(int i=0; i<nshapes; i++){
            for(int k=0; k<9; k++){
                b( (nposes + i) * 9 + k) = beta * w_prior[i][j] * prior[i][j](k/3, k%3);
            }
        }
        
        //cout<<"Constructing linear system done..."<<endl;
        
        Eigen::SparseMatrix<double> A(nrows, ncols);
        A.setFromTriplets(A_coeffs.begin(), A_coeffs.end());
        A.makeCompressed();
        Eigen::SparseMatrix<double> AtA = (A.transpose() * A).pruned();
        
        // epsilon
        const double epsilon = 1e-6;
        Eigen::SparseMatrix<double> eye(ncols, ncols);
        for(int j=0; j<ncols; j++) eye.insert(j,j) = epsilon;
        AtA += eye;
        
        // solve Ax = b
        Eigen::CholmodSupernodalLLT< Eigen::SparseMatrix<double> > solver;
        solver.compute(AtA);
        
        Eigen::VectorXd x = solver.solve(A.transpose() * b);  //vectorXd //46 * 9
        
        M[j] = std::vector<Eigen::Matrix3d>(nshapes);  
        for(int i=0; i<nshapes; i++){
            for(int k=0; k<9; k++){
                M[j][i](k/3, k%3) = x(9*i + k, 0); 
            }
        }
    }
    
    // reconstruct the blendshape
    ColorStream(ColorOutput::Blue)<<"reconstruct the blendshape from local gradient frame ...";
    MeshTransfer transfer;
    transfer.setSource(B0);
    transfer.setTarget(B0);
    
    std::vector<trimesh::TriMesh> B_new(nshapes);
    for(int i=0; i<nshapes; i++){        
        std::vector<Eigen::Matrix3d> Bgrad_i(nfaces);
        for(int j=0; j<nfaces; j++){
            Eigen::Matrix3d M0j = B0grad[j];
            Eigen::Matrix3d Mij = M[j][i];
            
            Bgrad_i[j] = ((Mij + M0j) * M0j.inverse()).transpose();       
        }
        B_new[i] = transfer.transfer(Bgrad_i, true);
        B_new[i].write("B/B_" + std::to_string(i) + ".obj");
    }
    return B_new;
}


void BlendshapeGenerator::blendshapeGeneration(){
    std::string A_path = "../data/Tester_1/Blendshape/";     // 47 template blenshape
    std::string T_path = "../data/Tester_1/TrainingPose/";   // 20 template poses
    std::string S_path = "../data/Tester_101/TrainingPose/"; // 20 training poses
    std::string B_path = "../data/Tester_101/Blendshape/";   // 47 blendshape to be solved (unknown)
    
    const int nshapes = NUM_SHAPES;  //47
    const int nposes  = NUM_POSES;   //20

    // given A0-A46, T0-T19, S0-S19
    // -> B0-B46
    std::vector<trimesh::TriMesh> A(nshapes); // 47 example blendshapes
    std::vector<trimesh::TriMesh> T(nposes);  // 20 example poses
    std::vector<trimesh::TriMesh> S(nposes);  // 20 training poses
    std::vector<trimesh::TriMesh> B(nshapes); // 47 unknown variables
    
    std::vector<trimesh::TriMesh> B_new(nshapes-1);
    
    // load the template blendshapes and groundtruth blendshapes
    for(int i=0; i<nshapes; i++){
        A[i]    = *trimesh::TriMesh::read(A_path + "shape_" + std::to_string(i) + ".obj");
    }
    
    for(int i=0; i<nshapes-1; i++){
        
        //B_new[i] = *trimesh::TriMesh::read(B_path + "shape_" + std::to_string(i+1) + ".obj");
    }
    
    // load training poses
    for(int i=0; i<nposes; i++){
        T[i] = *trimesh::TriMesh::read(T_path + "pose_" + std::to_string(i) + ".obj");
        S[i] = *trimesh::TriMesh::read(S_path + "pose_" + std::to_string(i) + ".obj");
    }
    
    trimesh::TriMesh A0 = A[0];
    trimesh::TriMesh B0 = S[0];
    
    // step 1.1: initialize B blendshape
    std::cout<<"1.1 Intialize Blendshape B ... "<<std::endl;
    MeshTransfer transfer;
    transfer.setSource(A0);
    transfer.setTarget(B0);
    std::vector<trimesh::TriMesh> B_init(nshapes);
    B_init[0] = B0;
    for(int i=1; i<nshapes; i++){
        std::cout << i <<" ";
        B_init[i] = transfer.transfer(A[i]);
        B_init[i].write("initBlendshape/"+std::to_string(i)+".obj");
    }
    std::cout<<"    Done"<<std::endl;
    
    
    // step 1.2: the delta shapes
    std::vector<trimesh::TriMesh> dB(nshapes-1);  // 46
    for(int i=0; i<nshapes-1; i++){
        dB[i] = B_init[i+1];
        for(int j=0; j<dB[i].vertices.size(); j++){
            dB[i].vertices[j] -= B0.vertices[j];
        }
    }
    
    
    // step 2.1: estimate the initial weights
    std::cout<<"2.1 Estimate initial weight ..."<<std::endl;
    std::vector<double> w0_vec(nshapes-1, 0.0), wp_vec(nshapes-1, 0.0);
    w0_vec[0] = 1.0;
    std::vector<std::vector<double> > alpha_refs(nposes, std::vector<double>(nshapes-1, 0.0));  // 20 * 46
    for(int i=0; i<nposes; i++)
    {
        std::cout << i <<" ";
//        alpha_refs[i] = estimateWeights(S[i], B0, dB, w0_vec, wp_vec, true);
//         
//         trimesh::TriMesh ret = reconstructMesh(alpha_refs[i], B0, dB);
//         ret.write("reconstruct/"+std::to_string(i)+".ply");
        
    }
    std::cout<<"    Done"<<std::endl;
    
    
    // step 2.2: triangle gradient
    std::cout<<"2.2 TriangleGradient for A0 and B0 ...";
    int nfaces = A0.faces.size();
    std::vector<Eigen::Matrix3d> MA0, MB0;  // 22800 * mat33 // local frame for each triangle
    for(int j=0; j<nfaces; j++){
        Eigen::Matrix3d MA0j = triangleGradient(A0, j);
        Eigen::Matrix3d MB0j = triangleGradient(B0, j);
        MA0.push_back(MA0j);
        MB0.push_back(MB0j);
    }
    std::cout<<"    Done"<<std::endl;
    
    
    // step 3.1: prior parameters
    std::cout<<"3.1 w_prior and prior ..."<<std::endl;
    double kappa = 0.1, theta = 2;
    std::vector<std::vector<Eigen::Matrix3d> > prior(nshapes-1, std::vector<Eigen::Matrix3d>(nfaces)); // 46 * 22800 * mat33
    std::vector<std::vector<double> > w_prior(nshapes-1, std::vector<double>(nfaces, 0.0));            // 46 * 22800
    for(int i=0; i<nshapes-1; i++){   
        std::cout<<i<<" ";
        trimesh::TriMesh Ai = A[i+1];   // for A1 - A46
        
        for(int j=0; j<nfaces; j++){
            Eigen::Matrix3d MA0j = MA0[j];
            Eigen::Matrix3d MAij = triangleGradient(Ai, j);
            Eigen::Matrix3d GA0Ai = MAij * MA0j.inverse();
            Eigen::Matrix3d MB0j = MB0[j];
            prior[i][j] = GA0Ai * MB0j - MB0j;   // localframe for B
            
            // GA0Ai * MA0j - MA0j = MAij * MA0j.inverse() * MA0j - MA0j =  MAij * MA0j - MA0j
            double MAij_norm = (MAij - MA0j).norm();   
            w_prior[i][j] = pow( (1 + MAij_norm) / (kappa + MAij_norm), theta);
        }
    }
    std::cout<<"    Done"<<std::endl;
    
    
    // triangleGradient for S
    std::vector<std::vector<Eigen::Matrix3d> > Sgrad(nposes, std::vector<Eigen::Matrix3d>(nfaces));   // 20 * 22800 * mat33
    for(int i=0; i<nposes; i++){
        for(int j=0; j<nfaces; j++){
            Eigen::Matrix3d Sij = triangleGradient(S[i], j); 
            Sgrad[i][j] = Sij;
        }
    }
    
    
    /***********************************************
     *************** main loop *********************
     * ********************************************/
    bool   converged = false;
    int    maxIters = 10;//10;
    double beta_max  = 0.5,  beta_min  = 0.1;
    double gamma_max = 1000, gamma_min = 100;
    int numVerts = B[0].vertices.size();    // 11510
    
    std::vector<std::vector<double> > alpha = alpha_refs;   // 20*46
    int iters = 0;
    while(iters < maxIters && !converged){
        cout<<"Iteration: "<<iters<<endl;
        converged = false;
        
        double beta  = beta_max  - 1.0*iters/maxIters * (beta_max  - beta_min);
        double gamma = gamma_max - 1.0*iters/maxIters * (gamma_max - gamma_min);
        
        // refine blendshapes
        B_new = optimizeBlendshapes(S, Sgrad, B0, alpha, beta, gamma, prior, w_prior);
        for(int i=0; i<nshapes-1; i++) B[i+1] = B_new[i];
        
        
        // update delta shapes
        for(int i=0; i<nshapes-1; i++){
            for(int k=0; k<numVerts; k++){
                dB[i].vertices[k] = B[i+1].vertices[k] - B[0].vertices[k];
            }
        }
        
        // update weights
        for(int i=0; i<nposes; i++){
            alpha[i] = estimateWeights(S[i], B0, dB, alpha[i], wp_vec, true);
        }
        
        // update reconstruction vertices
        for(int i=0; i<nposes; i++){
            trimesh::TriMesh tmp = B0;
            for(int j=0; j<nshapes-1; j++){
                for(int k=0; k<numVerts; k++){
                    tmp.vertices[k] += alpha[i][j] * dB[j].vertices[k];
                }
            }
            S[i] = tmp;
        }
    
        
        // compute deformation gradient for each triangle of S
        for(int i=0; i<nposes; i++){
            for(int j=0; j<nfaces; j++){
                Sgrad[i][j] = triangleGradient(S[i], j);
            }
        }
    
        iters++;
        if(converged || iters == maxIters) break;
        
    }
    
    
    // save output
    for(int i=0; i<nshapes; i++){
        B[i].write("B/B_" + std::to_string(i) + ".obj");
    }
    for(int i=0; i<nposes; i++){
        S[i].write("S/S_" + std::to_string(i) + ".obj");
    }
    
    return;
    
}
