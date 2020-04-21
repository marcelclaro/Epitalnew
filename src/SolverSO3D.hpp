/*
 * SolverSO3D.hpp
 *
 *  Created on: Nov 11, 2014
 *      Author: marcel
 */

#ifndef SOLVERSO3D_HPP_
#define SOLVERSO3D_HPP_

namespace epital {

template<typename complextype,typename gridtype>
class SolverSO3D {
public:
	SolverSO3D();
	virtual ~SolverSO3D();

};

} /* namespace epital */

#include "SolverSO3D.cpp"

#endif /* SOLVERSO3D_HPP_ */
