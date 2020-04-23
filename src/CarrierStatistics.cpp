/*
 * CarrierStatistics.cpp
 *
 *  Created on: Dec 5, 2012
 *      Author: marcel
 */


#include "CMaterial.hpp"
#include <vector>
#include "CarrierStatistics.hpp"
#include "Solution.hpp"

//#include <o2scl/funct.h>
//#include <o2scl/root_brent_gsl.h>
//#include <o2scl/root_cern.h>

#include "BrentsRootFinder.hpp"

#ifndef CARRIERSTATISTICS_CPP_
#define CARRIERSTATISTICS_CPP_

namespace epital{

template<typename complextype>
CarrierStatistics<complextype>::CarrierStatistics() {

}

template<typename complextype>
CarrierStatistics<complextype>::~CarrierStatistics() {

}

template<typename complextype>
template<typename xType>
complextype CarrierStatistics<complextype>::FermiEnergyWells(Solution<complextype,xType> solution, complextype carriersperperiod, shared_ptr<Material> material, Temperature T, Band banda){

	complextype effectivemass;
	complextype fermilevel;
	complextype bracket;
	complextype carriersignal;

	vector<complextype> levels;

	for(int k=0; k < solution.getLevels();k++){
		levels.push_back(solution.getEnergy(k).real());
	}


	switch(banda){
	case Band::Conduction:
		effectivemass = material->effectiveMass_gamma();
		fermilevel=material->ValenceBandEnergy() + material->gap_gamma()+material->gap_gamma()*0.05;
		bracket = levels.back();
		carriersignal=1;
		break;
	case Band::ConductionL:
		effectivemass = material->effectiveMass_LDOS();
		fermilevel=material->ValenceBandEnergy() + material->gap_L();
		bracket = levels.back();
		carriersignal=1;
		break;
	case Band::ConductionX:
		effectivemass = material->effectiveMass_XDOS();
		fermilevel=material->ValenceBandEnergy() + material->gap_X();
		bracket = levels.back();
		carriersignal=1;
		break;
	case Band::HeavyHole:
		effectivemass = pow(material->effectiveMass_HHl()*material->effectiveMass_HHt()*material->effectiveMass_HHt(),1.0/3.0);
		fermilevel=material->ValenceBandEnergy();
		bracket = levels.back();
		carriersignal=-1;
		break;
	case Band::LightHole:
		effectivemass = pow(material->effectiveMass_lHl()*material->effectiveMass_lHt()*material->effectiveMass_lHt(),1.0/3.0);
		fermilevel=material->ValenceBandEnergy();
		bracket = levels.back();
		carriersignal=-1;
		break;
	default:
		throw std::runtime_error("Invalid carrier type(Band)");
		break;
	}


	CarrierNumber function(carriersperperiod, T, effectivemass, levels, solution.getEnergy(solution.getLevels()-1).real()+carriersignal*50*Constant::kb*T,carriersignal);

//	o2scl::root_brent_gsl<> solver;
//	o2scl::funct_mfptr<CarrierStatistics<complextype>::CarrierNumber> tosolve(&function,&CarrierStatistics<complextype>::CarrierNumber::tosolvefermilevel);
//	solver.solve_bkt(fermilevel,bracket,tosolve);
	auto func = std::mem_fn(&CarrierStatistics<complextype>::CarrierNumber::tosolvefermilevel);
	std::function<double(double)> funcpointer = std::bind(func, function, std::placeholders::_1);
	fermilevel = BrentsRootFinder::zbrent(funcpointer,fermilevel,bracket,1e-25l);


	return fermilevel;
}

template<typename complextype>
template<typename xType>
complextype CarrierStatistics<complextype>::FermiEnergyWells(SolutionPW<complextype,xType> solution, complextype carriersperperiod, shared_ptr<Material> material, Temperature T, Band banda){

	complextype effectivemass;
	complextype fermilevel;
	complextype bracket;
	complextype carriersignal;

	vector<complextype> levels;

	for(int k=0; k < solution.getLevels();k++){
		levels.push_back(solution.getEnergy(k).real());
	}


	switch(banda){
	case Band::Conduction:
		effectivemass = material->effectiveMass_gamma();
		fermilevel=material->ValenceBandEnergy() + material->gap_gamma()+energy_from_eV(0.001);
		bracket = levels.back();
		carriersignal=1;
		break;
	case Band::ConductionL:
		effectivemass = material->effectiveMass_LDOS();
		fermilevel=material->ValenceBandEnergy() + material->gap_L();
		bracket = levels.back();
		carriersignal=1;
		break;
	case Band::ConductionX:
		effectivemass = material->effectiveMass_XDOS();
		fermilevel=material->ValenceBandEnergy() + material->gap_X();
		bracket = levels.back();
		carriersignal=1;
		break;
	case Band::HeavyHole:
		effectivemass = pow(material->effectiveMass_HHl()*material->effectiveMass_HHt()*material->effectiveMass_HHt(),1.0/3.0);
		fermilevel=material->ValenceBandEnergy();
		bracket = levels.back();
		carriersignal=-1;
		break;
	case Band::LightHole:
		effectivemass = pow(material->effectiveMass_lHl()*material->effectiveMass_lHt()*material->effectiveMass_lHt(),1.0/3.0);
		fermilevel=material->ValenceBandEnergy();
		bracket = levels.back();
		carriersignal=-1;
		break;
	default:
		throw std::runtime_error("Invalid carrier type(Band)");
		break;
	}


	CarrierNumber function(carriersperperiod, T, effectivemass, levels, solution.getEnergy(solution.getLevels()-1).real()+carriersignal*50*Constant::kb*T,carriersignal);

	//o2scl::root_brent_gsl<> solver;
	//o2scl::funct_mfptr<CarrierStatistics<complextype>::CarrierNumber> tosolve(&function,&CarrierStatistics<complextype>::CarrierNumber::tosolvefermilevel);
	//solver.solve_bkt(fermilevel,bracket,tosolve);
	auto func = std::mem_fn(&CarrierStatistics<complextype>::CarrierNumber::tosolvefermilevel);
	std::function<double(double)> funcpointer = std::bind(func, function, std::placeholders::_1);
/*
	///TODO paralelize
	vector<pair<double,double>> guesses;
	bool negsignal;
	if(funcpointer(fermilevel)<0)
		negsignal=true;
	else
		negsignal=false;
	if(fermilevel<bracket){
		for(double i=fermilevel+energy_from_eV(0.001); i<=bracket; i+=energy_from_eV(0.001)){
			if((funcpointer(i)>=0&&negsignal) || (funcpointer(i)<0&& (!negsignal) )  ){
				negsignal=!negsignal;
				guesses.push_back(make_pair(i-energy_from_eV(0.001),i));
			}
			if(guesses.size()==static_cast<unsigned long int>(levels.size()))
				break;
		}
	}
	else{
		for(double i=fermilevel-energy_from_eV(0.001); i>=bracket; i-=energy_from_eV(0.001)){
			if((funcpointer(i)>=0&&negsignal) || (funcpointer(i)<0&& (!negsignal) )  ){
				negsignal=!negsignal;
				guesses.push_back(make_pair(i+energy_from_eV(0.001),i));
			}
			if(guesses.size()==static_cast<unsigned long int>(levels.size()))
				break;
		}
	}

*/

	//cout << funcpointer(fermilevel) << "/" << funcpointer(levels[0]) << "/" << funcpointer(levels[1]) << "\n";
	fermilevel = BrentsRootFinder::zbrent(funcpointer,fermilevel,bracket,1e-25l);


	return fermilevel;
}



template<typename complextype>
template<typename xType>
complextype CarrierStatistics<complextype>::QuasiFermiEnergyWells(std::complex<complextype> subbandenergy, complextype subbandcarriersperperiod, shared_ptr<Material> material, Temperature T, Band banda){

	complextype effectivemass;
	complextype fermilevel;
	complextype bracket;
	complextype carriersignal;

	vector<complextype> levels;


	levels.push_back(subbandenergy.real());


	switch(banda){
	case Band::Conduction:
		effectivemass = material->effectiveMass_gamma();
		fermilevel=material->ValenceBandEnergy() + material->gap_gamma()+material->gap_gamma()*0.05;
		bracket = levels.back();
		carriersignal=1;
		break;
	case Band::ConductionL:
		effectivemass = material->effectiveMass_LDOS();
		fermilevel=material->ValenceBandEnergy() + material->gap_L();
		bracket = levels.back();
		carriersignal=1;
		break;
	case Band::ConductionX:
		effectivemass = material->effectiveMass_XDOS();
		fermilevel=material->ValenceBandEnergy() + material->gap_X();
		bracket = levels.back();
		carriersignal=1;
		break;
	case Band::HeavyHole:
		effectivemass = pow(material->effectiveMass_HHl()*material->effectiveMass_HHt()*material->effectiveMass_HHt(),1.0/3.0);
		fermilevel=material->ValenceBandEnergy();
		bracket = levels.back();
		carriersignal=-1;
		break;
	case Band::LightHole:
		effectivemass = pow(material->effectiveMass_lHl()*material->effectiveMass_lHt()*material->effectiveMass_lHt(),1.0/3.0);
		fermilevel=material->ValenceBandEnergy();
		bracket = levels.back();
		carriersignal=-1;
		break;
	default:
		throw std::runtime_error("Invalid carrier type(Band)");
		break;
	}


	CarrierNumber function(subbandcarriersperperiod, T, effectivemass, levels, subbandenergy.real()+carriersignal*30*Constant::kb*T,carriersignal);

	auto func = std::mem_fn(&CarrierStatistics<complextype>::CarrierNumber::tosolvefermilevel);
	std::function<double(double)> funcpointer = std::bind(func, function, std::placeholders::_1);
	fermilevel = BrentsRootFinder::zbrent(funcpointer,fermilevel,bracket,1e-25l);



	return fermilevel;
}


template<typename complextype>
complextype CarrierStatistics<complextype>::CarriersinSubband(Temperature T, complextype subbandminima, complextype fermienergy,shared_ptr<Material> material,Band banda, complextype energytop){
	complextype effectivemass;
	complextype carriersignal=1.0;

	switch(banda){
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
		carriersignal=-1.0;
		break;
	case Band::LightHole:
		effectivemass = pow(material->effectiveMass_lHl()*material->effectiveMass_lHt()*material->effectiveMass_lHt(),1.0/3.0);
		carriersignal=-1.0;
		break;
	default:
		throw std::runtime_error("Invalid carrier type(Band)in CarriersinSubband");
		break;
	}


	complextype factor = effectivemass*Constant::kb*T/(Constant::pi*Constant::hbar*Constant::hbar);

	return factor*( ( carriersignal*(fermienergy+energytop-subbandminima) / (Constant::kb*T)  ) - std::log( 1 + std::exp( carriersignal*(energytop)  / (Constant::kb*T) ) ) + std::log( 1 + std::exp( carriersignal*(subbandminima-fermienergy) / (Constant::kb*T)  ) ) );
}

template<typename complextype>
complextype CarrierStatistics<complextype>::CarrierNumber::tosolvefermilevel(complextype energy){

	complextype sum=0.0;
	complextype factor = effectivemass*Constant::kb*temp/(Constant::pi*Constant::hbar*Constant::hbar);
	for(complextype level : energies){
		sum+=factor*( ( carriersignal*(supenergy-level) / (Constant::kb*temp)  ) - std::log( 1 + std::exp( carriersignal*(supenergy-energy)  / (Constant::kb*temp) ) ) + std::log( 1 + std::exp( carriersignal*(level-energy) / (Constant::kb*temp)  ) ) );
	}

	return carriersperperiod - sum;

}

template<typename complextype>
CarrierStatistics<complextype>::CarrierNumber::CarrierNumber(complextype carriersperperiod,Temperature temp,complextype effectivemass, vector<complextype> energies, complextype supenergy,int carriersignal)
	:  effectivemass(effectivemass), supenergy(supenergy), temp(temp), energies(energies), carriersperperiod(carriersperperiod), carriersignal(carriersignal){

}



}
#endif
