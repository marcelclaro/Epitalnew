/*
 * PhotonScattering.cpp
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

#include <CustomTypes.hpp>
#include <CMaterial.hpp>
#include <Materials.hpp>
#include <Units.hpp>
#include <Alloy.hpp>
#include <Grid.hpp>
#include <DFunction.hpp>
#include <Graphics.hpp>
#include <Heterostructure.hpp>


#include <PhotonScattering.hpp>

#ifndef PHOTONSCATTERING_CPP_
#define PHOTONSCATTERING_CPP_

namespace epital {

template<typename complextype>
PhotonScattering<complextype>::PhotonScattering() {


}

template<typename complextype>
PhotonScattering<complextype>::~PhotonScattering() {
	// TODO Auto-generated destructor stub
}

template<typename complextype>
template<typename xType>
complextype PhotonScattering<complextype>::couplingSquare(shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial){
	//Assert invalid wavefunctions
	if(wavefunctionfinal->getGrid()!=(wavefunctioninitial->getGrid())){
			throw std::runtime_error("Incompatible grid at FormFactor on PhononScattering");
	};

	//Temporary for multipication
	DiscreteFunction<std::complex<complextype>,xType> temp = wavefunctioninitial->getDifferential();


	complex<complextype> dot = wavefunctionfinal->Dot(temp);

	//return dot product square
	return real(conj(dot)*dot);

}


template<typename complextype>
template<typename xType>
complextype PhotonScattering<complextype>::MeanTransRate(Temperature T, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,
		complextype energyfinallevel, complextype energyinitiallevel, complextype finaldensity,complextype initialdensity, complextype power, shared_ptr<Material> material, Band band){

	complextype factor = 2*Constant::pi*pow(Constant::hbar,3)*Constant::e*Constant::e/(pow(0.067,2)*Constant::emass*Constant::emass*Constant::eps0*Constant::c);

	cout << "coupling=" << couplingSquare(wavefunctionfinal,wavefunctioninitial) << endl;
	cout << "factor=" << factor << endl;
	cout << "O=" << couplingSquare(wavefunctionfinal,wavefunctioninitial)*2/(0.067*(energyfinallevel-energyinitiallevel)) << endl;


	return factor*power*(initialdensity-finaldensity)*couplingSquare(wavefunctionfinal,wavefunctioninitial)/(pow(energyfinallevel-energyinitiallevel,2)*initialdensity*sqrt(material->getDielectric_optical()));

}



} /* namespace epital */

#endif
