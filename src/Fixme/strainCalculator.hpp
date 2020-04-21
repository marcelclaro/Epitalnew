/*
 * strainCalculator.hpp
 *
 *  Created on: Dec 18, 2014
 *      Author: marcelclaro
 */

#include "Grid.hpp"
#include "CustomTypes.hpp"
#include "CMaterial.hpp"
#include "Materials.hpp"
#include "Units.hpp"
#include "DFunction3D.hpp"
#include "Heterostructure3D.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SuperLUSupport>
#include <Eigen/UmfPackSupport>
#include <Eigen/CholmodSupport>
#include <unsupported/Eigen/src/SparseExtra/MarketIO.h>

// Must be set prior to any ViennaCL includes if you want to use ViennaCL algorithms on Eigen objects
#define VIENNACL_WITH_EIGEN 1
// ViennaCL headers
#include <viennacl/linalg/ilu.hpp>
#include <viennacl/linalg/cg.hpp>
#include <viennacl/linalg/bicgstab.hpp>
#include <viennacl/linalg/gmres.hpp>
#include <viennacl/io/matrix_market.hpp>

#include <paralution/paralution.hpp>


#include <vector>

#ifndef STRAINCALCULATOR_HPP_
#define STRAINCALCULATOR_HPP_

namespace epital {

template <typename TYPE>
class strainCalculator {
public:
	strainCalculator();
	virtual ~strainCalculator();

	static void solveDisplace(Heterostructure3D<TYPE>& structure, DiscreteFunction3D<TYPE,TYPE>& ux,DiscreteFunction3D<TYPE,TYPE>& uy,DiscreteFunction3D<TYPE,TYPE>& uz);

};

} /* namespace epital */

#include "strainCalculator.cpp"

#endif /* STRAINCALCULATOR_HPP_ */
