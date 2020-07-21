#include "TriangleGradient.h"

Eigen::Matrix3d triangleGradient(const trimesh::TriMesh& mesh, const int& fidx)
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


pair<Eigen::Matrix3d, double> triangleGradient2(const trimesh::TriMesh& mesh, const int& fidx)
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
		 
	double d = 0.5 * n.norm();

	return make_pair(G,d);
}
