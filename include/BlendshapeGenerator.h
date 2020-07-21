#ifndef BLENDSHAPE_GENERATOR_H
#define BLENDSHAPE_GENERATOR_H

#include "Common.h"
#include "MeshTransfer.h"

class BlendshapeGenerator{
public:
    void blendshapeGeneration();
    
    
protected:
    std::vector<double> estimateWeights(const trimesh::TriMesh& S, 
                                        const trimesh::TriMesh& B0,
                                        const std::vector<trimesh::TriMesh>& dB,
                                        const std::vector<double>& w0,
                                        const std::vector<double>& wp,
                                        bool  prior_weight_flag = true);
    
    trimesh::TriMesh reconstructMesh(const std::vector<double>& weight, 
                                     const trimesh::TriMesh& B0,
                                     const std::vector<trimesh::TriMesh>& dB);
    
    std::vector<trimesh::TriMesh> optimizeBlendshapes(const std::vector<trimesh::TriMesh>& S,
                                                      const std::vector<std::vector<Eigen::Matrix3d> >& Sgrad,
                                                      const trimesh::TriMesh& B0,
                                                      const std::vector<std::vector<double> >& alpha,
                                                      double beta, double gamma,
                                                      const std::vector<std::vector<Eigen::Matrix3d> >& prior,
                                                      const std::vector<std::vector<double> >& w_prior);

private:
    std::vector<trimesh::TriMesh> m_S;
    
};


#endif
