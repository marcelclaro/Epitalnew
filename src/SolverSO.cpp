/*
 * SolverSO.cpp
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

#include "SolverSO.hpp"

#ifndef SOLVERSO_CPP_
#define SOLVERSO_CPP_

using namespace epital;

template<typename complextype,typename xType>
SolverSO<complextype,xType>::SolverSO() {

}

template<typename complextype,typename xType>
SolverSO<complextype,xType>::~SolverSO() {

}

template<typename complextype,typename xType>
template<typename yType>
Solution<complextype,xType> SolverSO<complextype,xType>::SolvePeriodic(Heterostructure<yType>& sample, Band selectedband, long int levels, complextype timedifusion, complextype timestep, long int gridpoints, std::function<std::complex<complextype>(xType)> inifuncpointer){
	if(inifuncpointer==nullptr){
		inifuncpointer = std::bind(SolverSO<complextype,xType>::oddevengaussian,std::placeholders::_1, 10*sample.getTotalWidth(), sample.Begin()+0.5*sample.getTotalWidth(), 1);
	}


	//TODO Assert invalid parameters

	Grid1D<xType> basegrid(sample.Begin(),sample.End(),gridpoints);

	vector<shared_ptr<DiscreteFunction<std::complex<complextype>,xType>>> initialfunctions(levels);
	for(long int i = 0; i< levels; i++)
		initialfunctions[i] = make_shared<DiscreteFunction<std::complex<complextype>,xType>>(basegrid,inifuncpointer);

	for(long int k = 0; k< levels; k++)
		initialfunctions[k]->normalize(1.0);

	Solution<complextype,xType> solution;
	vector<shared_ptr<DiscreteFunction<std::complex<complextype>,xType>>> templevels;

	switch(selectedband){
	case Band::Conduction:
	{
		DiscreteFunction<yType,xType> potential = sample.ConductionBand_gamma(basegrid);
		DiscreteFunction<yType,xType> effmass = sample.effectiveMass_gamma(basegrid);

		WaveTimeEvolution<std::complex<complextype>,xType> evoluter(WaveTimeEvolution<std::complex<complextype>,xType>::Carrier::ELECTRON);

		evoluter.groundlevel_periodic(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[0]),potential,effmass);
		templevels.push_back(initialfunctions[0]);
		solution.addLevel(initialfunctions[0],evoluter.Energy((*initialfunctions[0]),potential,effmass).real());

		for(long int k = 1; k < levels; k++){
			evoluter.newlevel_periodic(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[k]),potential,templevels,effmass);
			templevels.push_back(initialfunctions[k]);
			solution.addLevel(initialfunctions[k],evoluter.Energy((*initialfunctions[k]),potential,effmass).real());
		}
		break;
	}
	case Band::ConductionX:
	{
		DiscreteFunction<yType,xType> potential = sample.ConductionBand_X(basegrid);
		DiscreteFunction<yType,xType> effmass = sample.effectiveMass_Xt(basegrid);

		WaveTimeEvolution<std::complex<complextype>,xType> evoluter(WaveTimeEvolution<std::complex<complextype>,xType>::Carrier::ELECTRON);

		evoluter.groundlevel_periodic(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[0]),potential,effmass);
		templevels.push_back(initialfunctions[0]);
		solution.addLevel(initialfunctions[0],evoluter.Energy((*initialfunctions[0]),potential,effmass).real());

		for(long int k = 1; k < levels; k++){
			evoluter.newlevel_periodic(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[k]),potential,templevels,effmass);
			templevels.push_back(initialfunctions[k]);
			solution.addLevel(initialfunctions[k],evoluter.Energy((*initialfunctions[k]),potential,effmass).real());
		}
		break;
	}
	case Band::ConductionL:
	{
		DiscreteFunction<yType,xType> potential = sample.ConductionBand_L(basegrid);
		DiscreteFunction<yType,xType> effmass = sample.effectiveMass_Lt(basegrid);

		WaveTimeEvolution<std::complex<complextype>,xType> evoluter(WaveTimeEvolution<std::complex<complextype>,xType>::Carrier::ELECTRON);

		evoluter.groundlevel_periodic(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[0]),potential,effmass);
		templevels.push_back(initialfunctions[0]);
		solution.addLevel(initialfunctions[0],evoluter.Energy((*initialfunctions[0]),potential,effmass).real());

		for(long int k = 1; k < levels; k++){
			evoluter.newlevel_periodic(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[k]),potential,templevels,effmass);
			templevels.push_back(initialfunctions[k]);
			solution.addLevel(initialfunctions[k],evoluter.Energy((*initialfunctions[k]),potential,effmass).real());
		}
		break;
	}
	case Band::HeavyHole:
	{
		DiscreteFunction<yType,xType> potential = sample.ValenceBand_HH(basegrid);
		DiscreteFunction<yType,xType> effmass = sample.effectiveMass_HHt(basegrid);

		WaveTimeEvolution<std::complex<complextype>,xType> evoluter(WaveTimeEvolution<std::complex<complextype>,xType>::Carrier::HOLE);

		evoluter.groundlevel_periodic(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[0]),potential,effmass);
		templevels.push_back(initialfunctions[0]);
		solution.addLevel(initialfunctions[0],evoluter.Energy((*initialfunctions[0]),potential,effmass).real());

		for(long int k = 1; k < levels; k++){
			evoluter.newlevel_periodic(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[k]),potential,templevels,effmass);
			templevels.push_back(initialfunctions[k]);
			solution.addLevel(initialfunctions[k],evoluter.Energy((*initialfunctions[k]),potential,effmass).real());
		}

		break;
	}
	case Band::LightHole:
	{
		DiscreteFunction<yType,xType> potential = sample.ValenceBand_lH(basegrid);
		DiscreteFunction<yType,xType> effmass = sample.effectiveMass_lHt(basegrid);

		WaveTimeEvolution<std::complex<complextype>,xType> evoluter(WaveTimeEvolution<std::complex<complextype>,xType>::Carrier::HOLE);

		evoluter.groundlevel_periodic(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[0]),potential,effmass);
		templevels.push_back(initialfunctions[0]);
		solution.addLevel(initialfunctions[0],evoluter.Energy((*initialfunctions[0]),potential,effmass).real());

		for(long int k = 1; k < levels; k++){
			evoluter.newlevel_periodic(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[k]),potential,templevels,effmass);
			templevels.push_back(initialfunctions[k]);
			solution.addLevel(initialfunctions[k],evoluter.Energy((*initialfunctions[k]),potential,effmass).real());
		}

		break;
	}
	default:
		throw std::invalid_argument("invalid band");
	}


	return solution;
}


template<typename complextype,typename xType>
template<typename yType>
Solution<complextype,xType> SolverSO<complextype,xType>::SolvePeriodicFFT(Heterostructure<yType>& sample, Band selectedband, long int levels, complextype timedifusion, complextype timestep, long int gridpoints, std::function<std::complex<complextype>(xType)> inifuncpointer){
	if(inifuncpointer==nullptr){
		inifuncpointer = std::bind(SolverSO<complextype,xType>::oddevengaussian,std::placeholders::_1, 10*sample.getTotalWidth(), sample.Begin()+0.5*sample.getTotalWidth(), 1);
	}


	//TODO Assert invalid parameters

	Grid1D<xType> basegrid(sample.Begin(),sample.End(),gridpoints);

	vector<shared_ptr<DiscreteFunction<std::complex<complextype>,xType>>> initialfunctions(levels);
	for(long int i = 0; i< levels; i++)
		initialfunctions[i] = make_shared<DiscreteFunction<std::complex<complextype>,xType>>(basegrid,inifuncpointer);

	for(long int k = 0; k< levels; k++)
		initialfunctions[k]->normalize(1.0);

	Solution<complextype,xType> solution;
	vector<shared_ptr<DiscreteFunction<std::complex<complextype>,xType>>> templevels;

	switch(selectedband){
	case Band::Conduction:
	{
		DiscreteFunction<yType,xType> potential = sample.ConductionBand_gamma(basegrid);
		DiscreteFunction<yType,xType> effmass = sample.effectiveMass_gamma(basegrid);

		WaveTimeEvolution<std::complex<complextype>,xType> evoluter(WaveTimeEvolution<std::complex<complextype>,xType>::Carrier::ELECTRON);

		evoluter.groundlevel_periodicFFT(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[0]),potential,effmass);
		templevels.push_back(initialfunctions[0]);
		solution.addLevel(initialfunctions[0],evoluter.Energy((*initialfunctions[0]),potential,effmass).real());

		for(long int k = 1; k < levels; k++){
			evoluter.newlevel_periodic(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[k]),potential,templevels,effmass);
			templevels.push_back(initialfunctions[k]);
			solution.addLevel(initialfunctions[k],evoluter.Energy((*initialfunctions[k]),potential,effmass).real());
		}
		break;
	}
	case Band::ConductionX:
	{
		DiscreteFunction<yType,xType> potential = sample.ConductionBand_X(basegrid);
		DiscreteFunction<yType,xType> effmass = sample.effectiveMass_Xt(basegrid);

		WaveTimeEvolution<std::complex<complextype>,xType> evoluter(WaveTimeEvolution<std::complex<complextype>,xType>::Carrier::ELECTRON);

		evoluter.groundlevel_periodicFFT(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[0]),potential,effmass);
		templevels.push_back(initialfunctions[0]);
		solution.addLevel(initialfunctions[0],evoluter.Energy((*initialfunctions[0]),potential,effmass).real());

		for(long int k = 1; k < levels; k++){
			evoluter.newlevel_periodic(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[k]),potential,templevels,effmass);
			templevels.push_back(initialfunctions[k]);
			solution.addLevel(initialfunctions[k],evoluter.Energy((*initialfunctions[k]),potential,effmass).real());
		}
		break;
	}
	case Band::ConductionL:
	{
		DiscreteFunction<yType,xType> potential = sample.ConductionBand_L(basegrid);
		DiscreteFunction<yType,xType> effmass = sample.effectiveMass_Lt(basegrid);

		WaveTimeEvolution<std::complex<complextype>,xType> evoluter(WaveTimeEvolution<std::complex<complextype>,xType>::Carrier::ELECTRON);

		evoluter.groundlevel_periodicFFT(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[0]),potential,effmass);
		templevels.push_back(initialfunctions[0]);
		solution.addLevel(initialfunctions[0],evoluter.Energy((*initialfunctions[0]),potential,effmass).real());

		for(long int k = 1; k < levels; k++){
			evoluter.newlevel_periodic(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[k]),potential,templevels,effmass);
			templevels.push_back(initialfunctions[k]);
			solution.addLevel(initialfunctions[k],evoluter.Energy((*initialfunctions[k]),potential,effmass).real());
		}
		break;
	}
	case Band::HeavyHole:
	{
		DiscreteFunction<yType,xType> potential = sample.ValenceBand_HH(basegrid);
		DiscreteFunction<yType,xType> effmass = sample.effectiveMass_HHt(basegrid);

		WaveTimeEvolution<std::complex<complextype>,xType> evoluter(WaveTimeEvolution<std::complex<complextype>,xType>::Carrier::HOLE);

		evoluter.groundlevel_periodicFFT(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[0]),potential,effmass);
		templevels.push_back(initialfunctions[0]);
		solution.addLevel(initialfunctions[0],evoluter.Energy((*initialfunctions[0]),potential,effmass).real());

		for(long int k = 1; k < levels; k++){
			evoluter.newlevel_periodic(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[k]),potential,templevels,effmass);
			templevels.push_back(initialfunctions[k]);
			solution.addLevel(initialfunctions[k],evoluter.Energy((*initialfunctions[k]),potential,effmass).real());
		}

		break;
	}
	case Band::LightHole:
	{
		DiscreteFunction<yType,xType> potential = sample.ValenceBand_lH(basegrid);
		DiscreteFunction<yType,xType> effmass = sample.effectiveMass_lHt(basegrid);

		WaveTimeEvolution<std::complex<complextype>,xType> evoluter(WaveTimeEvolution<std::complex<complextype>,xType>::Carrier::HOLE);

		evoluter.groundlevel_periodicFFT(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[0]),potential,effmass);
		templevels.push_back(initialfunctions[0]);
		solution.addLevel(initialfunctions[0],evoluter.Energy((*initialfunctions[0]),potential,effmass).real());

		for(long int k = 1; k < levels; k++){
			evoluter.newlevel_periodic(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[k]),potential,templevels,effmass);
			templevels.push_back(initialfunctions[k]);
			solution.addLevel(initialfunctions[k],evoluter.Energy((*initialfunctions[k]),potential,effmass).real());
		}

		break;
	}
	default:
		throw std::invalid_argument("invalid band");
	}


	return solution;
}


template<typename complextype,typename xType>
template<typename yType>
Solution<complextype,xType> SolverSO<complextype,xType>::Solve(Heterostructure<yType>& sample, Band selectedband, long int levels, complextype timedifusion, complextype timestep, long int gridpoints, std::function<std::complex<complextype>(xType)> inifuncpointer){
	if(inifuncpointer==nullptr){
		inifuncpointer = std::bind(SolverSO<complextype,xType>::oddevengaussian,std::placeholders::_1, 10*sample.getTotalWidth(), sample.Begin()+500, 1);
	}


	//TODO Assert invalid parameters

	Grid1D<xType> basegrid(sample.Begin(),sample.End(),gridpoints);

	vector<shared_ptr<DiscreteFunction<std::complex<complextype>,xType>>> initialfunctions(levels);
	//initialfunctions.push_back();
	for(long int i = 0; i< levels; i++)
		initialfunctions[i] = make_shared<DiscreteFunction<std::complex<complextype>,xType>>(basegrid,inifuncpointer);

	for(long int k = 0; k< levels; k++)
		initialfunctions[k]->normalize(1.0);

	Solution<complextype,xType> solution;
	vector<shared_ptr<DiscreteFunction<std::complex<complextype>,xType>>> templevels;



	switch(selectedband){
	case Band::Conduction:
	{
		DiscreteFunction<yType,xType> potential = sample.ConductionBand_gamma(basegrid);
		DiscreteFunction<yType,xType> effmass = sample.effectiveMass_gamma(basegrid);

		WaveTimeEvolution<std::complex<complextype>,xType> evoluter(WaveTimeEvolution<std::complex<complextype>,xType>::Carrier::ELECTRON);

		evoluter.groundlevel_finite(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[0]),potential,effmass);
		templevels.push_back(initialfunctions[0]);
		solution.addLevel(initialfunctions[0],evoluter.Energy((*initialfunctions[0]),potential,effmass).real());

		for(long int k = 1; k < levels; k++){
			evoluter.newlevel_finite(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),*(initialfunctions[k]),potential,templevels,effmass);
			templevels.push_back(initialfunctions[k]);
			solution.addLevel(initialfunctions[k],evoluter.Energy((*initialfunctions[k]),potential,effmass).real());
		}
		break;
	}
	case Band::ConductionX:
	{
		DiscreteFunction<yType,xType> potential = sample.ConductionBand_X(basegrid);
		DiscreteFunction<yType,xType> effmass = sample.effectiveMass_Xt(basegrid);

		WaveTimeEvolution<std::complex<complextype>,xType> evoluter(WaveTimeEvolution<std::complex<complextype>,xType>::Carrier::ELECTRON);

		evoluter.groundlevel_finite(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[0]),potential,effmass);
		templevels.push_back(initialfunctions[0]);
		solution.addLevel(initialfunctions[0],evoluter.Energy((*initialfunctions[0]),potential,effmass).real());

		for(long int k = 1; k < levels; k++){
			evoluter.newlevel_finite(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[k]),potential,templevels,effmass);
			templevels.push_back(initialfunctions[k]);
			solution.addLevel(initialfunctions[k],evoluter.Energy((*initialfunctions[k]),potential,effmass).real());
		}
		break;
	}
	case Band::ConductionL:
	{
		DiscreteFunction<yType,xType> potential = sample.ConductionBand_L(basegrid);
		DiscreteFunction<yType,xType> effmass = sample.effectiveMass_Lt(basegrid);

		WaveTimeEvolution<std::complex<complextype>,xType> evoluter(WaveTimeEvolution<std::complex<complextype>,xType>::Carrier::ELECTRON);

		evoluter.groundlevel_finite(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[0]),potential,effmass);
		templevels.push_back(initialfunctions[0]);
		solution.addLevel(initialfunctions[0],evoluter.Energy((*initialfunctions[0]),potential,effmass).real());

		for(long int k = 1; k < levels; k++){
			evoluter.newlevel_finite(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[k]),potential,templevels,effmass);
			templevels.push_back(initialfunctions[k]);
			solution.addLevel(initialfunctions[k],evoluter.Energy((*initialfunctions[k]),potential,effmass).real());
		}
		break;
	}
	case Band::HeavyHole:
	{
		DiscreteFunction<yType,xType> potential = sample.ValenceBand_HH(basegrid);
		DiscreteFunction<yType,xType> effmass = sample.effectiveMass_HHt(basegrid);

		WaveTimeEvolution<std::complex<complextype>,xType> evoluter(WaveTimeEvolution<std::complex<complextype>,xType>::Carrier::HOLE);

		evoluter.groundlevel_finite(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[0]),potential,effmass);
		templevels.push_back(initialfunctions[0]);
		solution.addLevel(initialfunctions[0],evoluter.Energy((*initialfunctions[0]),potential,effmass).real());

		for(long int k = 1; k < levels; k++){
			evoluter.newlevel_finite(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[k]),potential,templevels,effmass);
			templevels.push_back(initialfunctions[k]);
			solution.addLevel(initialfunctions[k],evoluter.Energy((*initialfunctions[k]),potential,effmass).real());
		}

		break;
	}
	case Band::LightHole:
	{
		DiscreteFunction<yType,xType> potential = sample.ValenceBand_lH(basegrid);
		DiscreteFunction<yType,xType> effmass = sample.effectiveMass_lHt(basegrid);

		WaveTimeEvolution<std::complex<complextype>,xType> evoluter(WaveTimeEvolution<std::complex<complextype>,xType>::Carrier::HOLE);

		evoluter.groundlevel_finite(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[0]),potential,effmass);
		templevels.push_back(initialfunctions[0]);
		solution.addLevel(initialfunctions[0],evoluter.Energy((*initialfunctions[0]),potential,effmass).real());

		for(long int k = 1; k < levels; k++){
			evoluter.newlevel_finite(std::complex<complextype>(0,-timedifusion),static_cast<long int>((timedifusion/timestep)),(*initialfunctions[k]),potential,templevels,effmass);
			templevels.push_back(initialfunctions[k]);
			solution.addLevel(initialfunctions[k],evoluter.Energy((*initialfunctions[k]),potential,effmass).real());
		}

		break;
	}
	default:
		throw std::invalid_argument("invalid band");
	}


	return solution;
}

template<typename complextype,typename xType>
std::complex<double> SolverSO<complextype,xType>::commongaussian(double x,double sigma, double mean, double norm){
	std::complex<double> arg(-(x-mean)*(x-mean)/(2*sigma));
	std::complex<double> value = (exp(arg)*norm/2.0)+(exp(arg)*std::complex<double>(0.0,1.0)*norm/2.0);
	return value;
}

template<typename complextype,typename xType>
std::complex<double> SolverSO<complextype,xType>::oddevengaussian(double x,double sigma, double mean, double norm){
	return (mean-x)*commongaussian(x,sigma,mean,norm)+commongaussian(x,sigma,mean,norm*10);
}

#endif
