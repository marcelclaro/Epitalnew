/*
 * SolverTM.hpp
 *
 *  Created on: Apr 22, 2013
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
#include <memory>

#include "CustomTypes.hpp"
#include "CMaterial.hpp"
#include "Materials.hpp"
#include "Units.hpp"
#include "Alloy.hpp"
#include "Grid.hpp"
#include "DFunction.hpp"
#include "Heterostructure.hpp"
#include "WaveTimeEvolution.hpp"
#include "Solution.hpp"
#include "SolutionPW.hpp"
#include "SolutionWannier.hpp"



#ifndef SOLVERTM_HPP_
#define SOLVERTM_HPP_


template<typename complextype, typename xType = double>
class PlaneWavesFunction;

#include "PWFunction.hpp"

using namespace std;

namespace epital {

template<typename complextype, typename xType>
class SolverTM {
public:
	SolverTM(Heterostructure<xType>& sample, Band selectedband);
	virtual ~SolverTM();


	/**
	 * Solve a periodic system
	 * @param sample pointer to sample
	 * @param selectedband band to solve
	 * @param levels number of levels to solve
	 * @return The solution
	 */
	template<typename yType>
	static SolutionPW<complextype,xType> SolvePeriodic(Heterostructure<yType>& sample, Band selectedband, long int levels, complextype levelseparation, complextype initialenergy, complextype finalenergy);

	template<typename yType>
	static vector<complextype> getminibandwidth(Heterostructure<yType>& sample, Band selectedband,SolutionPW<complextype,xType> solution, long int boundlevels);

	template<typename yType>
	static SolutionWannier<complextype,xType> getWannierFunctions(Heterostructure<yType>& sample, Band selectedband,SolutionPW<complextype,xType> solution, long int brillouinpoints=1000);

	//complextype Jp(complextype energy);

	//complextype Jm(complextype energy);
private:
	shared_ptr<PlaneWavesFunction<complextype,xType>> getWavefunction(double energy, complextype bloch);

	double Determinantz(double energy);
	double Determinantnonzero(double energy,double& value);

	vector<xType> potential;
	vector<xType> effectivemasses;
	vector<xType> widths;
	vector<pair<xType,xType>> limitslst;
	complextype carriersignal;

	//inline complextype Lpq(long int p, long int q,vector<std::complex<complextype>>& fi, vector<std::complex<complextype>>& hi , vector<std::complex<complextype>>& gi);

	//inline complextype Spq(long int p, long int q);


};

} /* namespace epital */

#include "SolverTM.cpp"
#endif /* SOLVERTM_HPP_ */
