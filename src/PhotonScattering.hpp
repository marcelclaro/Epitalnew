/*
 * PhotonScattering.hpp
 *
 *  Created on: Aug 20, 2013
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


#ifndef PHOTONSCATTERING_HPP_
#define PHOTONSCATTERING_HPP_

namespace epital {

template<typename complextype>
class PhotonScattering {
public:
	PhotonScattering();
	virtual ~PhotonScattering();

	template<typename xType>
	static complextype couplingSquare(shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial);

	template<typename xType>
	static complextype MeanTransRate(Temperature T, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,
			complextype energyfinallevel, complextype energyinitiallevel, complextype finaldensity,complextype initialdensity, complextype power, shared_ptr<Material> material, Band band);



};

} /* namespace epital */

#include "PhotonScattering.cpp"
#endif /* PHOTONSCATTERING_HPP_ */
