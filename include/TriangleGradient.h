#ifndef TRIANGLE_GRADIENT_H
#define TRIANGLE_GRADIENT_H

#include "Common.h"


Eigen::Matrix3d triangleGradient(const trimesh::TriMesh& mesh, const int& fidx);

pair<Eigen::Matrix3d, double> triangleGradient2(const trimesh::TriMesh& mesh, const int& fidx);

#endif
