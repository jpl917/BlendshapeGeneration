#ifndef TEST_CASES_H
#define TEST_CASES_H

#include "Common.h"
#include "TriangleGradient.h"
#include "MeshTransfer.h"
#include "BlendshapeGenerator.h"


struct CostFunctor {
	template <typename T>
	bool operator()(const T* const x, T* residual) const 
	{
		residual[0] = T(10.0) - x[0];
		return true;
	}
};


class TestCases
{
public:
	static void testCeres();
	static void testEigenMatrix();
	static void testTriangleGradient();
    static void testMeshTransfer();
    static void testBlendshapeGeneration();
};



#endif


