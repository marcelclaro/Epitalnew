/*
 * QCD.hpp
 *
 *  Created on: Dec 11, 2012
 *      Author: marcel
 */

#include <thread>
#include <iostream>
#include <limits>
#include <string>
#include <cstring>
#include <cstdio>
#include <list>
#include <complex>
#include <chrono>
#include <ratio>
#include <vector>

#include "CustomTypes.hpp"
#include "CMaterial.hpp"
#include "Materials.hpp"
#include "Units.hpp"
#include "Alloy.hpp"
#include "Grid.hpp"
#include "DFunction.hpp"
#include "Graphics.hpp"
#include "Heterostructure.hpp"
#include "WaveTimeEvolution.hpp"
#include "Solution.hpp"
#include "SolverSO.hpp"

#ifndef QCD_HPP_
#define QCD_HPP_

namespace epital {

template<typename complextype,typename xType>
class QCD {
public:
	QCD(vector<typename Heterostructure<complextype>::Epilayer>& activeregion_, vector<typename Heterostructure<complextype>::Epilayer>& contacts_, shared_ptr<Material> substrate_ = make_shared<GaAs>(77));

	virtual ~QCD();

	complextype ResistanceArea(Temperature T);

	Solution<complextype,xType> Solve();

	Solution<complextype,xType> getSolution();

	void setSolveParams(Band selectedband, long int levels, complextype timedifusion = 100.0_fs, complextype timestep = 0.01_fs, long int gridpoints=2000, std::function<std::complex<complextype>(xType)> inifuncpointer = nullptr);

private:
	Heterostructure<complextype>* sample;
	Heterostructure<complextype>* doubledsample;
	Solution<complextype,xType> simplesolution;
	Solution<complextype,xType> doubledsolution;
	bool solved;
	bool dsolved;

	bool solveinitializated;
	Band selectedband;
	long int levels;
	complextype timedifusion;
	complextype timestep;
	long int gridpoints;
	std::function<std::complex<complextype>(xType)> inifuncpointer;

	static std::complex<complextype> localizefunction(xType x,xType begin, xType end);


};

} /* namespace epital */

#include "QCD.cpp"

#endif /* QCD_HPP_ */
