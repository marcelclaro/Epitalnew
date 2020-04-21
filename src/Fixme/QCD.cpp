/*
 * QCD.cpp
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
#include <vector>
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
#include <SolverSO.hpp>
#include <PhononScattering.hpp>
#include <CarrierStatistics.hpp>
#include <QCD.hpp>

#ifndef QCD_CPP_
#define QCD_CPP_

namespace epital {

template<typename complextype,typename xType>
QCD<complextype,xType>::QCD(vector<typename Heterostructure<complextype>::Epilayer>& activeregion_, vector<typename Heterostructure<complextype>::Epilayer>& contacts_, shared_ptr<Material> substrate_):  solved(false), dsolved(false),solveinitializated(false) {

	sample = new Heterostructure<double>(activeregion_,contacts_,substrate_);
	doubledsample = new Heterostructure<double>(activeregion_,contacts_,substrate_);
	doubledsample->duplicate();

}

template<typename complextype,typename xType>
QCD<complextype,xType>::~QCD() {
	delete sample;
	delete doubledsample;
}

template<typename complextype,typename xType>
complextype QCD<complextype,xType>::ResistanceArea(Temperature T){

	if(!solveinitializated)
		throw std::runtime_error("Solve not initialized for QCD<>::ResistanceArea");

	if(sample->hasContact())
		throw std::runtime_error("Invalid sample(must be without contacts) on QCD<>::ResistanceArea()");

	vector<long int> leftcascadelevels;
	vector<long int> rightcascadelevels;

	doubledsample->setBIASladder(0.75_Volts/20,15);

	if(!dsolved){
		doubledsolution = SolverSO<complextype,complextype>::Solve(doubledsample,selectedband,2*levels,timedifusion,timestep,2*gridpoints,inifuncpointer);
		dsolved=true;
	}

	if(!solved){
		simplesolution = SolverSO<complextype,complextype>::SolvePeriodic(sample,selectedband,levels,timedifusion,timestep,gridpoints,inifuncpointer);
		solved=true;
	}

	cout << "Solved System" << endl;

	std::function<std::complex<complextype>(xType)> firstfunc = std::bind(&QCD<complextype,xType>::localizefunction, std::placeholders::_1, sample->Begin(),sample->End());
	DiscreteFunction<std::complex<complextype>,xType> firstcascade(doubledsolution.getWavefunction(0)->getGrid(),firstfunc);

	cout << "Begin level separation" << endl;

	for(long int k = 0; k< (2*levels); k++){
		if(doubledsolution.getWavefunction(k)->meanLocalization() < sample->End()){
			leftcascadelevels.push_back(k);
		}
		else{
			rightcascadelevels.push_back(k);
		}
	}

	cout << "Plot" << endl;

	Grid1D<double> basegrid(doubledsample->Begin(),doubledsample->End(),2*gridpoints);
	DiscreteFunction<double,double> potential = doubledsample->ConductionBand_gamma(basegrid);

	//Create graphics
	Graphics pt(doubledsolution.getWavefunction(leftcascadelevels[0])->probabilitydensity().Real()+doubledsolution.getEnergy(leftcascadelevels[0]).real());
	Graphics pt2(doubledsolution.getWavefunction(leftcascadelevels[1])->probabilitydensity().Real()+doubledsolution.getEnergy(leftcascadelevels[1]).real());
	Graphics pt3(doubledsolution.getWavefunction(leftcascadelevels[2])->probabilitydensity().Real()+doubledsolution.getEnergy(leftcascadelevels[2]).real());
	Graphics pt4(doubledsolution.getWavefunction(leftcascadelevels[3])->probabilitydensity().Real()+doubledsolution.getEnergy(leftcascadelevels[3]).real());
	Graphics pt5(doubledsolution.getWavefunction(leftcascadelevels[4])->probabilitydensity().Real()+doubledsolution.getEnergy(leftcascadelevels[4]).real());
	Graphics pt6(doubledsolution.getWavefunction(leftcascadelevels[5])->probabilitydensity().Real()+doubledsolution.getEnergy(leftcascadelevels[5]).real());
	Graphics pt7(doubledsolution.getWavefunction(leftcascadelevels[6])->probabilitydensity().Real()+doubledsolution.getEnergy(leftcascadelevels[6]).real());
	Graphics pt8(doubledsolution.getWavefunction(leftcascadelevels[7])->probabilitydensity().Real()+doubledsolution.getEnergy(leftcascadelevels[7]).real());
	//Create graphics
	Graphics ptb(doubledsolution.getWavefunction(rightcascadelevels[0])->probabilitydensity().Real()+doubledsolution.getEnergy(rightcascadelevels[0]).real());
	Graphics pt2b(doubledsolution.getWavefunction(rightcascadelevels[1])->probabilitydensity().Real()+doubledsolution.getEnergy(rightcascadelevels[1]).real());
	Graphics pt3b(doubledsolution.getWavefunction(rightcascadelevels[2])->probabilitydensity().Real()+doubledsolution.getEnergy(rightcascadelevels[2]).real());
	Graphics pt4b(doubledsolution.getWavefunction(rightcascadelevels[3])->probabilitydensity().Real()+doubledsolution.getEnergy(rightcascadelevels[3]).real());
	Graphics pt5b(doubledsolution.getWavefunction(rightcascadelevels[4])->probabilitydensity().Real()+doubledsolution.getEnergy(rightcascadelevels[4]).real());
	Graphics pt6b(doubledsolution.getWavefunction(rightcascadelevels[5])->probabilitydensity().Real()+doubledsolution.getEnergy(rightcascadelevels[5]).real());
	Graphics pt7b(doubledsolution.getWavefunction(rightcascadelevels[6])->probabilitydensity().Real()+doubledsolution.getEnergy(rightcascadelevels[6]).real());
	Graphics pt8b(doubledsolution.getWavefunction(rightcascadelevels[7])->probabilitydensity().Real()+doubledsolution.getEnergy(rightcascadelevels[7]).real());



	Graphics gama(potential);
	Graphics all({&gama,&pt,&pt2,&pt3,&pt4,&pt5,&pt6,&pt7,&pt8});
	all.setXrange(doubledsample->Begin(),doubledsample->End());
	all.DataConvertionY(energy_to_eV);
	all.Plot();

	Graphics allb({&gama,&ptb,&pt2b,&pt3b,&pt4b,&pt5b,&pt6b,&pt7b,&pt8b});
	allb.setXrange(doubledsample->Begin(),doubledsample->End());
	allb.DataConvertionY(energy_to_eV);
	allb.Plot();


	complextype fermilevel = CarrierStatistics<complextype>::FermiEnergyWells(simplesolution,sample->getTotalDoping(),sample->getSubstrate(),T,selectedband);

	cout << "Fermilevel= " << energy_to_eV(fermilevel) <<  endl;

	complextype sum=0;
	long int leftindex = 0;

	for(long int i: leftcascadelevels){
		complextype carriers = CarrierStatistics<complextype>::CarriersinSubband(T,doubledsolution.getEnergy(i).real(),fermilevel,sample->getSubstrate(),selectedband,fermilevel+30*Constant::kb*T);
		long int rightindex = 0;
		for(long int j: rightcascadelevels){
			if(rightindex>=levels)
				break;
			if(rightindex!=leftindex&&(i>1)&&(j>(i+1))){
				cout << "Scatering between levels" << leftindex << "," << rightindex;
				complextype scateringrate = carriers*( PhononScattering<double>::MeanAbsRateLO(T,doubledsolution.getWavefunction(j),doubledsolution.getWavefunction(i),doubledsolution.getEnergy(j).real(),doubledsolution.getEnergy(i).real(),fermilevel-0.75_Volts/20,fermilevel,sample->getSubstrate(),selectedband)
						+ PhononScattering<double>::MeanEmissionRateLO(T,doubledsolution.getWavefunction(j),doubledsolution.getWavefunction(i),doubledsolution.getEnergy(j).real(),doubledsolution.getEnergy(i).real(),fermilevel-0.75_Volts/20,fermilevel,sample->getSubstrate(),selectedband) );

				sum+=scateringrate;
				cout << " Rate: " << frequency_to_SI(scateringrate*3.57116e+20) << endl;
			}
			++rightindex;
		}
		++leftindex;
	}

	return ( Constant::kb*T / (Constant::e*Constant::e*sum) );
}

template<typename complextype,typename xType>
std::complex<complextype> QCD<complextype,xType>::localizefunction(xType x, xType begin, xType end){
	if(x>begin && x<end){
		return std::complex<complextype>(1.0,0.0);
	}
	else
		return std::complex<complextype>(0.0,0.0);
}


template<typename complextype,typename xType>
Solution<complextype,xType> QCD<complextype,xType>::Solve(){

	if(!solveinitializated)
		throw std::runtime_error("Solve not initialized for QCD<>:Solve");

	if(solved)
		return simplesolution;
	else{
		simplesolution = SolverSO<complextype,xType>::SolvePeriodic(sample,selectedband,levels,timedifusion,timestep,gridpoints,inifuncpointer);
		solved=true;
		return simplesolution;
	}
}

template<typename complextype,typename xType>
Solution<complextype,xType> QCD<complextype,xType>::getSolution(){
	if(solved){
		return simplesolution;
	}
	else{
		throw std::runtime_error("QCD<>::getSolution() requires the use of QCD<>::Solve() first");
		return Solution<complextype,xType>();
	}
}

template<typename complextype,typename xType>
void QCD<complextype,xType>::setSolveParams(Band selectedband, long int levels, complextype timedifusion , complextype timestep , long int gridpoints, std::function<std::complex<complextype>(xType)> inifuncpointer ){
	this->selectedband = selectedband;
	this->levels = levels;
	this->timedifusion = timedifusion;
	this->timestep = timestep;
	this->inifuncpointer = inifuncpointer;
	this->gridpoints = gridpoints;

	solveinitializated = true;
	solved = false;
	dsolved = false;
}



} /* namespace epital */

#endif
