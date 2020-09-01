#include "TestCases.h"


void TestCases::testCeres()
{
	double init_x = 5.0;
	double x = init_x;
	
	//build the problem
	ceres::Problem problem;
	
	// set up the cost function 
	// use the autodifferentation to obtain the derivate
	ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<CostFunctor, 1, 1>(new CostFunctor);
	problem.AddResidualBlock(cost_function, NULL, &x);
	
	// run the solver
	ceres::Solver::Options options;
	options.linear_solver_type = ceres::DENSE_QR;
	options.minimizer_progress_to_stdout = false;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	
	std::cout<<summary.BriefReport()<<std::endl;
	std::cout<<"init x:"<<init_x<<"   ->  "<<x<<std::endl;
	
	return;
}


void TestCases::testEigenMatrix()
{
	/*
	2 -1 0 0 0
	-1 2 -1 0 0
	0 -1 2 -1 0
	0 0 -1 2 -1
	0 0 0 -1 2
	*/
	Eigen::SparseMatrix<double> A(5, 5);
	
	std::vector<Eigen::Triplet<double> > _A;
	_A.push_back(Eigen::Triplet<double>(0, 0, 2));
	_A.push_back(Eigen::Triplet<double>(0, 1,-1));
	_A.push_back(Eigen::Triplet<double>(1, 0,-1));
	_A.push_back(Eigen::Triplet<double>(1, 1, 2));
	_A.push_back(Eigen::Triplet<double>(1, 2,-1));
	_A.push_back(Eigen::Triplet<double>(2, 1,-1));
	_A.push_back(Eigen::Triplet<double>(2, 2, 2));
	_A.push_back(Eigen::Triplet<double>(2, 3,-1));
	_A.push_back(Eigen::Triplet<double>(3, 2,-1));
	_A.push_back(Eigen::Triplet<double>(3, 3, 2));
	_A.push_back(Eigen::Triplet<double>(3, 4,-1));
	_A.push_back(Eigen::Triplet<double>(4, 3,-1));
	_A.push_back(Eigen::Triplet<double>(4, 4, 2));
	_A.push_back(Eigen::Triplet<double>(2, 1, 2));
	
	A.setFromTriplets(_A.begin(), _A.end());

	
	auto AtA = A.transpose() * A;
	
	Eigen::Matrix<double,5,1> b;
	for(int i=0; i<5; i++) b(i,0)=1;
	
	Eigen::CholmodSupernodalLLT< Eigen::SparseMatrix<double> > solver;
	solver.compute(AtA);
	Eigen::VectorXd x = solver.solve(A.transpose()*b);
	
	for(int i=0; i<5; i++) std::cout<<x(i,0)<<" ";
	std::cout<<std::endl;
	
	Eigen::VectorXd b1 = A * x;
	for(int i=0; i<5; i++) std::cout<<b1(i,0)<<" ";
	std::cout<<std::endl;
}

void TestCases::testTriangleGradient()
{
	trimesh::TriMesh* mesh = trimesh::TriMesh::read("../emily-23686.obj");
	Eigen::Matrix3d tmp = triangleGradient(*mesh, 0);
    cout<<tmp * tmp.inverse()<<endl;
	
	pair<Eigen::Matrix3d, double> ret = triangleGradient2(*mesh, 0);
	cout<<ret.first<<endl;
	cout<<ret.second<<endl;
	
	return;
}


void TestCases::testMeshTransfer(){
    cout<<"Test Mesh Transfer"<<endl;
    
    MeshTransfer transfer;
    trimesh::TriMesh S0 = *trimesh::TriMesh::read("../data/Tester_1/Blendshape/shape_0.obj");
    trimesh::TriMesh S1 = *trimesh::TriMesh::read("../data/Tester_1/Blendshape/shape_1.obj");
    trimesh::TriMesh T0 = *trimesh::TriMesh::read("../data/Tester_101/Blendshape/shape_0.obj");
    transfer.setSource(S0);
    transfer.setTarget(T0);
    transfer.transfer(S1);
}   

void TestCases::testBlendshapeGeneration(){
    BlendshapeGenerator bg;
    bg.blendshapeGeneration();
}

void TestCases::testDeformationTransfer(){  
    cout<<"Test Deformation Transfer"<<endl;
    
    trimesh::TriMesh A0 = *trimesh::TriMesh::read("../data/Tester_1/Blendshape/shape_0.obj");
    trimesh::TriMesh A1 = *trimesh::TriMesh::read("../data/Tester_1/Blendshape/shape_1.obj");
    trimesh::TriMesh B0 = *trimesh::TriMesh::read("../data/Tester_101/Blendshape/shape_0.obj");
    trimesh::TriMesh B1;
    
    DeformationTransfer dt;
    dt.init(A0, B0);
    dt.transfer(A1, B1);
    
    /*
    #pragma omp parallel for num_threads(8)
	for (int iMesh = 0; iMesh < 47; iMesh++)
	{
		const int tid = omp_get_thread_num();
	}
     */
    
    return;
}

int main(int argc, char** argv)
{
//     if (argc != 3 && argc != 4)
// 	{
// 		printf("Usage: TestCase [src_folder] [target0.obj] [result_folder]");
// 		return -1;
// 	}
    
	TestCases test;
    
    //test.testEigenMatrix();
    test.testMeshTransfer();
    
    //test.testBlendshapeGeneration();
    test.testDeformationTransfer();
    
    return 0;
}
