/*
 * PhononScattering.cpp
 *
 *  Created on: Dec 3, 2012
 *      Author: marcel
 */

#include <thread>
#include <iostream>
#include <limits>
#include <string>
#include <cstring>
#include <cstdio>
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
#include <WaveTimeEvolution.hpp>
#include <Solution.hpp>
#include <gsl/gsl_integration.h>

#include <cmath>
#include <o2scl/funct.h>
#include <o2scl/inte_qagi_gsl.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/exception.h>

#include "PhononScattering.hpp"

#ifndef PHONONSCATTERING_CPP_
#define PHONONSCATTERING_CPP_

using namespace std;
namespace epital{

template<typename complextype>
PhononScattering<complextype>::PhononScattering() {
}

template<typename complextype>
PhononScattering<complextype>::~PhononScattering() {
}

template<typename complextype>
template<typename xType>
complextype PhononScattering<complextype>::FormFactorSquare(shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial, complextype phononwavenumber){
	//Assert invalid wavefunctions
	if(wavefunctionfinal->getGrid()!=(wavefunctioninitial->getGrid())){
			throw std::runtime_error("Incompatible grid at FormFactor on PhononScattering");
	};

	//Temporary for multipication
	DiscreteFunction<std::complex<complextype>,xType> temp(wavefunctioninitial->getGrid());

	//Multiplication
	xType step = wavefunctionfinal->getGrid().getIncrement();
	xType initialpoint = wavefunctionfinal->getGrid().getInferiorLimit();

	#pragma omp parallel for
	for(long int k = 0; k < wavefunctioninitial->getGrid().getSize(); k++){
		temp[k]=std::exp(std::complex<complextype>(0.0,-phononwavenumber*(initialpoint+k*step)))*(*wavefunctioninitial)(k);
	}

	complex<complextype> dot = wavefunctionfinal.Dot(temp);
	//return dot product square
	return real(conj(dot)*dot);

}

template<typename complextype>
template<typename xType>
complextype PhononScattering<complextype>::Integrand<xType>::thefunction(complextype phononwavenumber){

	//Assert invalid wavefunction
	if(wavefunctionfinal->getGrid()!=(wavefunctioninitial->getGrid())){
				throw std::runtime_error("Incompatible grid at integrand function on PhononScattering");
		};

	//steps for form factor...
	DiscreteFunction<std::complex<complextype>,xType> temp(wavefunctioninitial->getGrid());

	xType step = wavefunctionfinal->getGrid().getIncrement();
	xType initialpoint = wavefunctionfinal->getGrid().getInferiorLimit();

	for(long int k = 0; k < wavefunctioninitial->getGrid().getSize(); k++){
		temp[k]=std::exp(std::complex<complextype>(0.0,-phononwavenumber*(initialpoint+k*step)))*(*wavefunctioninitial)(k);
	}

    //calculate the integrand with a factor of 1e+10
	return 1.0e+15*Constant::pi*real(conj(wavefunctionfinal->Dot(temp))*wavefunctionfinal->Dot(temp))/( 1.0e+5 * (sqrt( pow(phononwavenumber,4) + 2*phononwavenumber*phononwavenumber*(2*kparalleli*kparalleli - (2*effectivemass*deltaenergy/(Constant::hbar*Constant::hbar)) ) + pow((2*effectivemass*deltaenergy/(Constant::hbar*Constant::hbar) ),2) ) ) );
}

template<typename complextype>
template<typename xType>
PhononScattering<complextype>::Integrand<xType>::Integrand(shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal,shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,complextype kparalleli,complextype deltaenergy,complextype effectivemass): kparalleli(kparalleli), deltaenergy(deltaenergy), effectivemass(effectivemass){
	if(wavefunctionfinal==nullptr||wavefunctioninitial==nullptr){
		throw std::runtime_error("Invalid wavefunction pointer (initial or final) at Integrand");
	}
	this->wavefunctionfinal=wavefunctionfinal;
	this->wavefunctioninitial=wavefunctioninitial;
}


template<typename complextype>
template<typename xType>
complextype PhononScattering<complextype>::IntegrandforMean<xType>::functionemi(complextype initialenergy){
	return 1.0e+10*(1/ ( std::exp( (initialenergy - fermiinitial) / (Constant::kb*T)  ) + 1 ) )* (1 - (1/ ( std::exp( (initialenergy - (Constant::hbar*wLO) - fermifinal) / (Constant::kb*T)  ) + 1 ) ) ) * PhononScattering<complextype>::EmissionRateLO(T,wavefunctionfinal,wavefunctioninitial,energyfinallevel,initialenergy,0.0,material,band);
}

template<typename complextype>
template<typename xType>
complextype PhononScattering<complextype>::IntegrandforMean<xType>::functionabs(complextype initialenergy){
	return 1.0e+10*(1/ ( std::exp( (initialenergy - fermiinitial) / (Constant::kb*T)  ) + 1 ) )* (1 - (1/ ( std::exp( (initialenergy + (Constant::hbar*wLO) - fermifinal) / (Constant::kb*T)  ) + 1 ) ) ) * PhononScattering<complextype>::AbsRateLO(T,wavefunctionfinal,wavefunctioninitial,energyfinallevel,initialenergy,0.0,material,band);
}

template<typename complextype>
template<typename xType>
complextype PhononScattering<complextype>::IntegrandforMean<xType>::functioncarriersub(complextype initialenergy){
	return 1.0e+10*(1/ ( std::exp( (initialenergy - fermiinitial) / (Constant::kb*T)  ) + 1 ) );
}

template<typename complextype>
template<typename xType>
PhononScattering<complextype>::IntegrandforMean<xType>::IntegrandforMean(Temperature T, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,
		complextype energyfinallevel,complextype fermifinal, complextype fermiinitial, shared_ptr<Material> material, Band band): T(T), wavefunctionfinal(wavefunctionfinal), wavefunctioninitial(wavefunctioninitial),
				energyfinallevel(energyfinallevel), fermifinal(fermifinal), fermiinitial(fermiinitial), material(material), band(band){
	wLO = material->getwLO();
}


template<typename complextype>
template<typename xType>
complextype PhononScattering<complextype>::EmissionRateLO(Temperature T, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,complextype energyfinallevel, complextype energyinitiallevel,complextype kparalleli, shared_ptr<Material> material, Band band){
	complextype integral;
	complextype integralerror;
	complextype factor;
	complextype effectivemass;
	complextype deltaenergy;
	complextype phononpop;

	bool integration_ok = false;

	if(!material)
		throw std::runtime_error("null material at EmissionRateLO");

	//Assert invalid wavefunction
	if(wavefunctionfinal->getGrid()!=(wavefunctioninitial->getGrid())){
				throw std::runtime_error("Incompatible grid at EmissionRateLO");
		};

	switch(band){
	case Band::Conduction:
		effectivemass = material->effectiveMass_gamma();
		break;
	case Band::ConductionL:
		effectivemass = material->effectiveMass_LDOS();
		break;
	case Band::ConductionX:
		effectivemass = material->effectiveMass_XDOS();
		break;
	case Band::HeavyHole:
		effectivemass = pow(material->effectiveMass_HHl()*material->effectiveMass_HHt()*material->effectiveMass_HHt(),1.0/3.0);
		break;
	case Band::LightHole:
		effectivemass = pow(material->effectiveMass_lHl()*material->effectiveMass_lHt()*material->effectiveMass_lHt(),1.0/3.0);
		break;
	default:
		throw std::runtime_error("Invalid carrier type(Band)");
		break;
	}

	//Phonon population -> bose-einstein distribuition
	phononpop = 1/(std::exp((Constant::hbar*material->getwLO()) / (Constant::kb*T)  )-1);

	deltaenergy=energyfinallevel-energyinitiallevel+Constant::hbar*material->getwLO();

	//Selection rule in energy
	if((kparalleli*kparalleli - (2*effectivemass*deltaenergy/(Constant::hbar*Constant::hbar)) ) < 0.0){
		return 0.0;
	}
	else{
		//P factor
		double Patm=(1/Constant::eps0)*(1/material->getDielectric_optical()-1/material->getDielectric_static());

		//Y'' factor
		factor=2*effectivemass*Constant::e*Constant::e*Patm*material->getwLO()/(pow(2*Constant::pi*Constant::hbar,2));

		//Integrand construction
		Integrand<xType> integ(wavefunctionfinal,wavefunctioninitial,kparalleli,deltaenergy,effectivemass);
		std::function<complextype(xType)> f = std::bind(std::mem_fn<xType(xType)> (&PhononScattering<complextype>::Integrand<xType>::thefunction),&integ,std::placeholders::_1);
		//o2scl::funct_mfptr<PhononScattering<complextype>::Integrand<xType>> f(&integ,&PhononScattering<complextype>::Integrand<xType>::thefunction);

		//Integrator definition
		o2scl::inte_qagi_gsl<> integrator;
		integrator.tol_rel = 1e-3;
		integrator.tol_abs = 0;

		while(!integration_ok){
			try{
				int ret1=integrator.integ_err(f,0.0,0.0,integral,integralerror);
					if(ret1!=0)
						throw std::runtime_error("Integration error on EmissionRateLO() error code:" + ret1);
				integration_ok = true;
			}
			catch(o2scl::exc_runtime_error& e ){
				integrator.tol_rel = integrator.tol_rel*10;
				if(integrator.tol_rel>=1)
					throw std::runtime_error("Integration error on EmissionRateLO() convergence error");
				integration_ok= false;
			}
		}


		//Correct normalization of integration
		integral/=1.0e+10;

		return (phononpop+1)*integral*factor/2;
	}
}


template<typename complextype>
template<typename xType>
complextype PhononScattering<complextype>::AbsRateLO(Temperature T, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,complextype energyfinallevel, complextype energyinitiallevel,complextype kparalleli, shared_ptr<Material> material, Band band){
	complextype integral;
	complextype integralerror;
	complextype factor;
	complextype effectivemass;
	complextype deltaenergy;
	complextype phononpop;

	bool integration_ok = false;

	//Assert invalid wavefunction
	if(wavefunctionfinal->getGrid()!=(wavefunctioninitial->getGrid())){
				throw std::runtime_error("Incompatible grid at EmissionRateLO");
		};

	switch(band){
	case Band::Conduction:
		effectivemass = material->effectiveMass_gamma();
		break;
	case Band::ConductionL:
		effectivemass = material->effectiveMass_LDOS();
		break;
	case Band::ConductionX:
		effectivemass = material->effectiveMass_XDOS();
		break;
	case Band::HeavyHole:
		effectivemass = pow(material->effectiveMass_HHl()*material->effectiveMass_HHt()*material->effectiveMass_HHt(),1.0/3.0);
		break;
	case Band::LightHole:
		effectivemass = pow(material->effectiveMass_lHl()*material->effectiveMass_lHt()*material->effectiveMass_lHt(),1.0/3.0);
		break;
	default:
		throw std::runtime_error("Invalid carrier type(Band)");
		break;
	}

	//Phonon population -> bose-einstein distribuition
	phononpop = 1/(std::exp((Constant::hbar*material->getwLO()) / (Constant::kb*T)  )-1);

	deltaenergy=energyfinallevel-energyinitiallevel-Constant::hbar*material->getwLO();

	//Selection rule in energy
	if((kparalleli*kparalleli - (2*effectivemass*deltaenergy/(Constant::hbar*Constant::hbar)) ) < 0.0){
		return 0.0;
	}
	else{
		//P factor
		double Patm=(1/Constant::eps0)*(1/material->getDielectric_optical()-1/material->getDielectric_static());

		//Y'' factor
		factor=2*effectivemass*Constant::e*Constant::e*Patm*material->getwLO()/(pow(2*Constant::pi*Constant::hbar,2));

		//Integrand construction
		Integrand<xType> integ(wavefunctionfinal,wavefunctioninitial,kparalleli,deltaenergy,effectivemass);
		std::function<complextype(xType)> f = std::bind(std::mem_fn<xType(xType)> (&PhononScattering<complextype>::Integrand<xType>::thefunction),&integ,std::placeholders::_1);
//		o2scl::funct_mfptr<PhononScattering<complextype>::Integrand<xType>> f(&integ,&PhononScattering<complextype>::Integrand<xType>::thefunction);

		//Integrator definition
		o2scl::inte_qagi_gsl<> integrator;
		integrator.tol_rel = 1e-3;
		integrator.tol_abs = 0;


		while(!integration_ok){
			try{
				int ret1=integrator.integ_err(f,0.0,0.0,integral,integralerror);
				if(ret1!=0)
					throw std::runtime_error("Integration error on AbsRateLO() error code:" + ret1);
				integration_ok = true;
			}
			catch(o2scl::exc_runtime_error& e ){
				integrator.tol_rel = integrator.tol_rel*10;
				if(integrator.tol_rel>=1)
					throw std::runtime_error("Integration error on AbsRateLO() convergence error");
				integration_ok= false;
			}
		}



		//Correct normalization of integration
		integral/=1.0e+10;

		return (phononpop)*integral*factor/2;
	}
}

template<typename complextype>
template<typename xType>
complextype PhononScattering<complextype>::MeanEmissionRateLO(Temperature T, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,complextype energyfinallevel, complextype energyinitiallevel,complextype fermifinal,complextype fermiinitial , shared_ptr<Material> material, Band band, long int limitpar){
	complextype integral;
	complextype integralcarr;
	complextype integralerror;
	complextype integralcarrerror;

	bool integration_ok = false;
	bool integrationb_ok = false;

	//Assert invalid wavefunction
	if(wavefunctionfinal->getGrid()!=(wavefunctioninitial->getGrid())){
				throw std::runtime_error("Incompatible grid at MeanEmissionRateLO");
		};

	PhononScattering<complextype>::IntegrandforMean<xType> integrand(T,wavefunctionfinal,wavefunctioninitial,energyfinallevel,fermifinal,fermiinitial,material,band);
	std::function<complextype(xType)> f = std::bind(std::mem_fn<xType(xType)> (&PhononScattering<complextype>::IntegrandforMean<xType>::functionemi),&integrand,std::placeholders::_1);
	std::function<complextype(xType)> fcarr = std::bind(std::mem_fn<xType(xType)> (&PhononScattering<complextype>::IntegrandforMean<xType>::functioncarriersub),&integrand,std::placeholders::_1);

	//o2scl::funct_mfptr<PhononScattering<complextype>::IntegrandforMean<xType>> f(&integrand,&PhononScattering<complextype>::IntegrandforMean<xType>::functionemi);
	//o2scl::funct_mfptr<PhononScattering<complextype>::IntegrandforMean<xType>> fcarr(&integrand,&PhononScattering<complextype>::IntegrandforMean<xType>::functioncarriersub);


	o2scl::inte_qag_gsl<> integrator;
	integrator.tol_rel = 1e-3;
	integrator.tol_abs = 0;


	while(!integration_ok){
		try{
			int ret1=integrator.integ_err(f,energyinitiallevel,(energyinitiallevel+limitpar*Constant::kb*T),integral,integralerror);
				if(ret1!=0)
					throw std::runtime_error("Integration error on MeanEmissionRateLO() error code:" + ret1);
			integration_ok = true;

		}
		catch(o2scl::exc_runtime_error& e ){
			integrator.tol_rel = integrator.tol_rel*10;
			if(integrator.tol_rel>=1)
				throw std::runtime_error("Integration error on MeanEmissionRateLO() convergence error");
			integration_ok= false;
		}
		catch(std::runtime_error& f){
			if(limitpar>10)
				limitpar=10;
			else
				limitpar-=2;

			if(limitpar<4)
				throw std::runtime_error(f.what());
			integrationb_ok = false;
		}
	}

	integrator.tol_rel = 1e-3;
	integrator.tol_abs = 0;

	while(!integrationb_ok){
		try{
			int ret2=integrator.integ_err(fcarr,energyinitiallevel,(energyinitiallevel+limitpar*Constant::kb*T),integralcarr,integralcarrerror);
				if(ret2!=0)
					throw std::runtime_error("Integration error on carrier of MeanEmissionRateLO() error code:" + ret2);
			integrationb_ok = true;

		}
		catch(o2scl::exc_runtime_error& e){
			integrator.tol_rel = integrator.tol_rel*10;
			if(integrator.tol_rel>=1)
				throw std::runtime_error("Integration error on carrier MeanEmissionRateLO() convergence error");
			integrationb_ok= false;
		}
		catch(std::runtime_error& f){
			if(limitpar>10)
				limitpar=10;
			else
				limitpar-=2;

			if(limitpar<4)
				throw std::runtime_error(f.what());
			integrationb_ok = false;
		}
	}


	return integral/integralcarr;


}

template<typename complextype>
template<typename xType>
complextype PhononScattering<complextype>::MeanAbsRateLO(Temperature T, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,complextype energyfinallevel, complextype energyinitiallevel,complextype fermifinal,complextype fermiinitial , shared_ptr<Material> material, Band band, long int limitpar){
	complextype integral;
	complextype integralcarr;
	complextype integralerror;
	complextype integralcarrerror;
	bool integration_ok = false;
	bool integrationb_ok = false;

	//Assert invalid wavefunction
	if(wavefunctionfinal->getGrid()!=(wavefunctioninitial->getGrid())){
				throw std::runtime_error("Incompatible grid at MeanAbsRateLO");
		};

	PhononScattering<complextype>::IntegrandforMean<xType> integrand(T,wavefunctionfinal,wavefunctioninitial,energyfinallevel,fermifinal,fermiinitial,material,band);

	std::function<complextype(xType)> f = std::bind(std::mem_fn<xType(xType)> (&PhononScattering<complextype>::IntegrandforMean<xType>::functionabs),&integrand,std::placeholders::_1);
	std::function<complextype(xType)> fcarr = std::bind(std::mem_fn<xType(xType)> (&PhononScattering<complextype>::IntegrandforMean<xType>::functioncarriersub),&integrand,std::placeholders::_1);

	//o2scl::funct_mfptr<PhononScattering<complextype>::IntegrandforMean<xType>> f(&integrand,&PhononScattering<complextype>::IntegrandforMean<xType>::functionabs);
	//o2scl::funct_mfptr<PhononScattering<complextype>::IntegrandforMean<xType>> fcarr(&integrand,&PhononScattering<complextype>::IntegrandforMean<xType>::functioncarriersub);


	o2scl::inte_qag_gsl<> integrator;
	integrator.tol_rel = 1e-3;
	integrator.tol_abs = 0;

	while(!integration_ok){
		try{
			int ret1=integrator.integ_err(f,energyinitiallevel,(energyinitiallevel+limitpar*Constant::kb*T),integral,integralerror);
			if(ret1!=0)
				throw std::runtime_error("Integration error on MeanAbsRateLO() error code:" + ret1);
			integration_ok = true;

		}
		catch(o2scl::exc_runtime_error& e ){
			integrator.tol_rel = integrator.tol_rel*10;
			if(integrator.tol_rel>=1)
				throw std::runtime_error("Integration error on MeanAbsRateLO() convergence error");
			integration_ok= false;
		}
		catch(std::runtime_error& f){
			if(limitpar>10)
				limitpar=10;
			else
				limitpar-=2;

			if(limitpar<4)
				throw std::runtime_error(f.what());
			integrationb_ok = false;
		}
	}


	integrator.tol_rel = 1e-3;
	integrator.tol_abs = 0;

	while(!integrationb_ok){
		try{
			int ret2=integrator.integ_err(fcarr,energyinitiallevel,(energyinitiallevel+limitpar*Constant::kb*T),integralcarr,integralcarrerror);
				if(ret2!=0)
					throw std::runtime_error("Integration error on MeanAbsRateLO() error code:" + ret2);
			integrationb_ok = true;

		}
		catch(o2scl::exc_runtime_error& e ){
			integrator.tol_rel = integrator.tol_rel*10;
			if(integrator.tol_rel>=1)
				throw std::runtime_error("Integration error on MeanAbsRateLO() convergence error");
			integrationb_ok= false;
		}
		catch(std::runtime_error& f){
			if(limitpar>10)
				limitpar=10;
			else
				limitpar-=2;

			if(limitpar<4)
				throw std::runtime_error(f.what());
			integrationb_ok = false;
		}
	}

	return integral/integralcarr;

}



template<typename complextype>
template<typename xType>
complextype PhononScattering<complextype>::IntersubbandMeanAbsRateLO(Temperature T, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,
		complextype energyfinallevel, complextype energyinitiallevel,complextype quasifermifinal,complextype quasifermiinitial , shared_ptr<Material> material, Band band, long int limitpar){
	complextype integral;
	complextype integralcarr;
	complextype integralerror;
	complextype integralcarrerror;
	bool integration_ok = false;
	bool integrationb_ok = false;

	//Assert invalid wavefunction
	if(wavefunctionfinal->getGrid()!=(wavefunctioninitial->getGrid())){
				throw std::runtime_error("Incompatible grid at IntersubbandMeanAbsRateLO");
		};

	PhononScattering<complextype>::IntegrandforMean<xType> integrand(T,wavefunctionfinal,wavefunctioninitial,energyfinallevel,quasifermifinal,quasifermiinitial,material,band);

	std::function<complextype(xType)> f = std::bind(std::mem_fn<xType(xType)> (&PhononScattering<complextype>::IntegrandforMean<xType>::functionabs),&integrand,std::placeholders::_1);
	std::function<complextype(xType)> fcarr = std::bind(std::mem_fn<xType(xType)> (&PhononScattering<complextype>::IntegrandforMean<xType>::functioncarriersub),&integrand,std::placeholders::_1);
	//o2scl::funct_mfptr<PhononScattering<complextype>::IntegrandforMean<xType>> f(&integrand,&PhononScattering<complextype>::IntegrandforMean<xType>::functionabs);
	//o2scl::funct_mfptr<PhononScattering<complextype>::IntegrandforMean<xType>> fcarr(&integrand,&PhononScattering<complextype>::IntegrandforMean<xType>::functioncarriersub);


	o2scl::inte_qag_gsl<> integrator;
	integrator.tol_rel = 1e-3;
	integrator.tol_abs = 0;

	while(!integration_ok){
		try{
			int ret1=integrator.integ_err(f,energyinitiallevel,(energyinitiallevel+limitpar*Constant::kb*T),integral,integralerror);
			if(ret1!=0)
				throw std::runtime_error("Integration error on IntersubbandMeanAbsRateLO() error code:" + ret1);
			integration_ok = true;

		}
		catch(o2scl::exc_runtime_error& e ){
			integrator.tol_rel = integrator.tol_rel*10;
			if(integrator.tol_rel>=1)
				throw std::runtime_error("Integration error on IntersubbandMeanAbsRateLO() convergence error");
			integration_ok= false;
		}
		catch(std::runtime_error& f){
			if(limitpar>10)
				limitpar=10;
			else
				limitpar-=2;

			if(limitpar<4)
				throw std::runtime_error(f.what());
			integrationb_ok = false;
		}
	}


	integrator.tol_rel = 1e-3;
	integrator.tol_abs = 0;

	while(!integrationb_ok){
		try{
			int ret2=integrator.integ_err(fcarr,energyinitiallevel,(energyinitiallevel+limitpar*Constant::kb*T),integralcarr,integralcarrerror);
				if(ret2!=0)
					throw std::runtime_error("Integration error on IntersubbandMeanAbsRateLO() error code:" + ret2);
			integrationb_ok = true;

		}
		catch(o2scl::exc_runtime_error& e ){
			integrator.tol_rel = integrator.tol_rel*10;
			if(integrator.tol_rel>=1)
				throw std::runtime_error("Integration error on IntersubbandMeanAbsRateLO() convergence error");
			integrationb_ok= false;
		}
		catch(std::runtime_error& f){
			if(limitpar>10)
				limitpar=10;
			else
				limitpar-=2;

			if(limitpar<4)
				throw std::runtime_error(f.what());
			integrationb_ok = false;
		}
	}

	return integral/integralcarr;

}

template<typename complextype>
template<typename xType>
complextype PhononScattering<complextype>::IntersubbandMeanEmissionRateLO(Temperature T, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,
		complextype energyfinallevel, complextype energyinitiallevel,complextype quasifermifinal,complextype quasifermiinitial , shared_ptr<Material> material, Band band, long int limitpar){
	complextype integral;
	complextype integralcarr;
	complextype integralerror;
	complextype integralcarrerror;

	bool integration_ok = false;
	bool integrationb_ok = false;

	//Assert invalid wavefunction
	if(wavefunctionfinal->getGrid()!=(wavefunctioninitial->getGrid())){
				throw std::runtime_error("Incompatible grid at MeanEmissionRateLO");
		};

	PhononScattering<complextype>::IntegrandforMean<xType> integrand(T,wavefunctionfinal,wavefunctioninitial,energyfinallevel,quasifermifinal,quasifermiinitial,material,band);

	std::function<complextype(xType)> f = std::bind(std::mem_fn<xType(xType)> (&PhononScattering<complextype>::IntegrandforMean<xType>::functionemi),&integrand,std::placeholders::_1);
	std::function<complextype(xType)> fcarr = std::bind(std::mem_fn<xType(xType)> (&PhononScattering<complextype>::IntegrandforMean<xType>::functioncarriersub),&integrand,std::placeholders::_1);

	//o2scl::funct_mfptr<PhononScattering<complextype>::IntegrandforMean<xType>> f(&integrand,&PhononScattering<complextype>::IntegrandforMean<xType>::functionemi);
	//o2scl::funct_mfptr<PhononScattering<complextype>::IntegrandforMean<xType>> fcarr(&integrand,&PhononScattering<complextype>::IntegrandforMean<xType>::functioncarriersub);


	o2scl::inte_qag_gsl<> integrator;
	integrator.tol_rel = 1e-3;
	integrator.tol_abs = 0;


	while(!integration_ok){
		try{
			int ret1=integrator.integ_err(f,energyinitiallevel,(energyinitiallevel+limitpar*Constant::kb*T),integral,integralerror);
				if(ret1!=0)
					throw std::runtime_error("Integration error on IntersubbandMeanEmissionRateLO() error code:" + ret1);
			integration_ok = true;

		}
		catch(o2scl::exc_runtime_error& e ){
			integrator.tol_rel = integrator.tol_rel*10;
			if(integrator.tol_rel>=1)
				throw std::runtime_error("Integration error on IntersubbandMeanEmissionRateLO() convergence error");
			integration_ok= false;
		}
		catch(std::runtime_error& f){
			if(limitpar>10)
				limitpar=10;
			else
				limitpar-=2;

			if(limitpar<4)
				throw std::runtime_error(f.what());
			integrationb_ok = false;
		}
	}

	integrator.tol_rel = 1e-3;
	integrator.tol_abs = 0;

	while(!integrationb_ok){
		try{
			int ret2=integrator.integ_err(fcarr,energyinitiallevel,(energyinitiallevel+limitpar*Constant::kb*T),integralcarr,integralcarrerror);
				if(ret2!=0)
					throw std::runtime_error("Integration error on carrier of MeanEmissionRateLO() error code:" + ret2);
			integrationb_ok = true;

		}
		catch(o2scl::exc_runtime_error& e){
			integrator.tol_rel = integrator.tol_rel*10;
			if(integrator.tol_rel>=1)
				throw std::runtime_error("Integration error on carrier MeanEmissionRateLO() convergence error");
			integrationb_ok= false;
		}
		catch(std::runtime_error& f){
			if(limitpar>10)
				limitpar=10;
			else
				limitpar-=2;

			if(limitpar<4)
				throw std::runtime_error(f.what());
			integrationb_ok = false;
		}
	}


	return integral/integralcarr;


}





}

#endif
