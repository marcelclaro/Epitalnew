/*
 * SolverSO.hpp
 *
 *  Created on: Nov 6, 2012
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

#ifndef SOLVERSO_HPP_
#define SOLVERSO_HPP_

namespace epital{

/**
 * \brief Solver using Split-Operator formalism
 */
template<typename complextype,typename xType>
class SolverSO {
public:
	/**
	 * Specify the band to solve;
	 */

	SolverSO();
	virtual ~SolverSO();

	/**
	 * Solve a finite (aperiodic) system
	 * @param sample pointer to sample
	 * @param selectedband band to solve
	 * @param levels number of levels to solve
	 * @param timedifusion time to get solution stabilization
	 * @param timestep time step of time operator
	 * @param gridpoints grid points of potential and wavefunction
	 * @param inifuncpointer pointer to initial wavefunction
	 * @return The solution
	 */
	template<typename yType>
	static Solution<complextype,xType> Solve(Heterostructure<yType>& sample, Band selectedband, long int levels, complextype timedifusion = 100.0_fs, complextype timestep = 0.01_fs, long int gridpoints=2000, std::function<std::complex<complextype>(xType)> inifuncpointer = nullptr);

	/**
	 * Solve a periodic system
	 * @param sample pointer to sample
	 * @param selectedband band to solve
	 * @param levels number of levels to solve
	 * @param timedifusion time to get solution stabilization
	 * @param timestep time step of time operator
	 * @param gridpoints grid points of potential and wavefunction
	 * @param inifuncpointer pointer to initial wavefunction
	 * @return The solution
	 */
	template<typename yType>
	static Solution<complextype,xType> SolvePeriodic(Heterostructure<yType>& sample, Band selectedband, long int levels, complextype timedifusion = 100.0_fs, complextype timestep = 0.01_fs, long int gridpoints=2000, std::function<std::complex<complextype>(xType)> inifuncpointer = nullptr);

	/**
	 * Solve a periodic system
	 * @param sample pointer to sample
	 * @param selectedband band to solve
	 * @param levels number of levels to solve
	 * @param timedifusion time to get solution stabilization
	 * @param timestep time step of time operator
	 * @param gridpoints grid points of potential and wavefunction
	 * @param inifuncpointer pointer to initial wavefunction
	 * @return The solution
	 */
	template<typename yType>
	static Solution<complextype,xType> SolvePeriodicFFT(Heterostructure<yType>& sample, Band selectedband, long int levels, complextype timedifusion = 100.0_fs, complextype timestep = 0.01_fs, long int gridpoints=2000, std::function<std::complex<complextype>(xType)> inifuncpointer = nullptr);



private:
	/**
	 * For use as initial wavefunction ( A gaussian)
	 * @param x variable
	 * @param sigma standard error
	 * @param mean mean of the function
	 * @param norm norm of the function
	 * @return function value at x
	 */
	static std::complex<double> commongaussian(double x,double sigma, double mean, double norm);

	/**
		 * For use as initial wavefunction ( A even and odd gaussian)
		 * @param x variable
		 * @param sigma standard error
		 * @param mean mean of the function
		 * @param norm norm of the function
		 * @return function value at x
		 */
	static std::complex<double> oddevengaussian(double x,double sigma, double mean, double norm);

};

#include "SolverSO.cpp"

}

#endif /* SOLVERSO_HPP_ */

