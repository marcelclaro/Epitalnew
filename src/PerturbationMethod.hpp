/*
 * PerturbationMethod.hpp
 *
 *  Created on: Apr 15, 2013
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

#ifndef PERTURBATIONMETHOD_HPP_
#define PERTURBATIONMETHOD_HPP_

namespace epital {

template <typename timetype, typename complextype>
class PerturbationMethod {
public:
	/**
		 * Type of carrier of wavefunction
		 */
		enum Carrier{
			HOLE=-1,  //!< HOLE
			ELECTRON=1//!< ELECTRON
		} carriertype;


	PerturbationMethod(PerturbationMethod::Carrier carrier = Carrier::ELECTRON);
	virtual ~PerturbationMethod();


	template<typename yType, typename xType>
	std::complex<complextype> EnergyBrilluoin(long int level, Solution<complextype,xType> solution, DiscreteFunction<yType,xType>& effectivemass,yType blochvector = 0.0);


private:
	template<typename yType, typename xType>
	std::complex<complextype> Vab(Solution<complextype,xType>& solution, DiscreteFunction<yType,xType>& effectivemass,yType blochvector,long int a, long int b);

	template<typename xType>
	std::complex<complextype> Eab(Solution<complextype,xType>& solution,long int a, long int b);

};

} /* namespace epital */

#include "PerturbationMethod.cpp"
#endif /* PERTURBATIONMETHOD_HPP_ */
