/*
 * SolverTM.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: marcel
 */

#include "SolverTM.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include "SolutionPW.hpp"
#include "WannierFunction.hpp"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "BrentsRootFinder.hpp"


using namespace Eigen;


#ifndef SOLVERTM_CPP_
#define SOLVERTM_CPP_

//#define PREC 50

namespace epital {

template<typename complextype, typename xType>
SolverTM<complextype,xType>::SolverTM(Heterostructure<xType>& sample, Band selectedband) {
	switch(selectedband){
	case Band::Conduction:
	{
		potential=sample.ConductionBand_gamma();
		effectivemasses=sample.effectiveMass_gamma();
		carriersignal=1.0;
		break;
	}
	case Band::ConductionX:
	{
		potential=sample.ConductionBand_L();
		effectivemasses=sample.effectiveMass_Lt();
		carriersignal=1.0;
		break;
	}
	case Band::ConductionL:
	{
		potential=sample.ConductionBand_X();
		effectivemasses=sample.effectiveMass_Xt();
		carriersignal=1.0;
		break;
	}
	case Band::HeavyHole:
	{
		potential=sample.ValenceBand_HH();
		effectivemasses=sample.effectiveMass_HHt();
		carriersignal=-1.0;
		break;
	}
	case Band::LightHole:
	{
		potential=sample.ValenceBand_lH();
		effectivemasses=sample.effectiveMass_lHt();
		carriersignal=-1.0;
		break;
	}
	default:
		throw std::invalid_argument("invalid band");
	}

	widths=sample.getWidths();
	limitslst=sample.getLimitsLists();

}

template<typename complextype, typename xType>
SolverTM<complextype,xType>::~SolverTM() {

}

template<typename complextype, typename xType>
template<typename yType>
SolutionPW<complextype,xType> SolverTM<complextype,xType>::SolvePeriodic(Heterostructure<yType>& sample, Band selectedband, long int levels, complextype levelseparation, complextype initialenergy, complextype finalenergy){


	vector<shared_ptr<PlaneWavesFunction<std::complex<complextype>,xType>>> functions;

	vector<pair<xType,xType>> limitslst_local=sample.getLimitsLists();
	///TODO Error if out of bounds

	SolverTM<complextype,complextype> solver(sample,selectedband);
	auto func = std::mem_fn(&SolverTM<complextype,xType>::Determinantz);
	std::function<double(double)> funcpointer = std::bind(func, solver, std::placeholders::_1);

	///TODO paralelize
	vector<pair<double,double>> guesses;
	bool negsignal;
	if(funcpointer(initialenergy)<0)
		negsignal=true;
	else
		negsignal=false;
	if(initialenergy<finalenergy){
		for(double i=initialenergy+levelseparation; i<=finalenergy; i+=levelseparation){
			if((funcpointer(i)>=0&&negsignal) || (funcpointer(i)<0&& (!negsignal) )  ){
				negsignal=!negsignal;
				guesses.push_back(make_pair(i-levelseparation,i));
			}
			if(guesses.size()==static_cast<unsigned long int>(levels))
				break;
		}
	}
	else{
		for(double i=initialenergy-levelseparation; i>=finalenergy; i-=levelseparation){
			if((funcpointer(i)>=0&&negsignal) || (funcpointer(i)<0&& (!negsignal) )  ){
				negsignal=!negsignal;
				guesses.push_back(make_pair(i+levelseparation,i));
			}
			if(guesses.size()==static_cast<unsigned long int>(levels))
				break;
		}
	}

	//o2scl::gsl_root_brent<> solverr;
	//o2scl::funct_mfptr<SolverTM<complextype,complextype>> tosolve(&solver,&SolverTM<complextype,complextype>::Determinantz);

	SolutionPW<complextype, xType> sol;

	for(pair<double,double> k : guesses){
		//complextype root = (BrentsRootFinder::brent_find_minima(funcpointer,k.first,k.second,58,20000L)).first;
		double root = BrentsRootFinder::zbrent(funcpointer,k.first,k.second,1e-25l);
		//solverr.solve(root,tosolve);
		sol.addLevel(solver.getWavefunction(root,Constant::pi/(2*sample.getTotalWidth())),std::complex<complextype>(root));
		//cout << "k =" << k.first << " Energy(eV)= " << energy_to_eV(root) << " (root0)= " << funcpointer(root)<< endl;
	}

	return sol;
}

template<typename complextype, typename xType>
template<typename yType>
vector<complextype> SolverTM<complextype,xType>::getminibandwidth(Heterostructure<yType>& sample, Band selectedband,SolutionPW<complextype,xType> solution, long int boundlevels){


	long int levels = solution.getLevels();

	vector<complextype> widths;


	SolverTM<complextype,complextype> solver(sample,selectedband);
	auto func = std::mem_fn(&SolverTM<complextype,xType>::Determinantnonzero);

	complextype a=1.0;
	complextype b=-1.0;

	std::function<complextype(xType)> funcpointera = std::bind(func, solver, std::placeholders::_1,a);
	std::function<complextype(xType)> funcpointerb = std::bind(func, solver, std::placeholders::_1,b);


	complextype lastresult=1.0e-6;
	complextype boundn=energy_from_eV(1.0e-6);
	complextype boundp=energy_from_eV(1.0e-6);

	for(long int i = 0; i < levels && i < boundlevels ; i++){
		long double rootp=solution.getEnergy(i).real();
		while((funcpointera(rootp-boundn) > 0.0l && funcpointera(rootp+boundp) > 0.0l) || (funcpointera(rootp-boundn) < 0.0l && funcpointera(rootp+boundp) < 0.0l)){
			boundn*=2.0;
			if((funcpointera(rootp-boundn) > 0.0l && funcpointera(rootp+boundp) > 0.0l) || (funcpointera(rootp-boundn) < 0.0l && funcpointera(rootp+boundp) < 0.0l))
				boundp*=2.0;
		}
		rootp = BrentsRootFinder::zbrent(funcpointera,(rootp-boundn),(rootp+boundn),(1.0e-19));
		long double rootm=solution.getEnergy(i).real();

		boundn=energy_from_eV(1.0e-6);
		boundp=energy_from_eV(1.0e-6);
		while((funcpointerb(rootm-boundn) > 0.0l && funcpointerb(rootm+boundp) > 0.0l) || (funcpointerb(rootm-boundn) < 0.0l && funcpointerb(rootm+boundp) < 0.0l)){
			boundn*=2.0;
			if((funcpointerb(rootm-boundn) > 0.0l && funcpointerb(rootm+boundp) > 0.0l) || (funcpointerb(rootm-boundn) < 0.0l && funcpointerb(rootm+boundp) < 0.0l))
				boundp*=2.0;
		}
		rootm = BrentsRootFinder::zbrent(funcpointerb,(rootm-boundn),(rootm+boundp),(1e-19));
		widths.push_back(rootp-rootm);
		lastresult=abs(rootp-rootm);
	}

	return widths;
}


template<typename complextype, typename xType>
template<typename yType>
SolutionWannier<complextype,xType> SolverTM<complextype,xType>::getWannierFunctions(Heterostructure<yType>& sample, Band selectedband,SolutionPW<complextype,xType> solution,long int brillouinpoints){

	long int levels = solution.getLevels();
	complextype period = sample.getLimitsLists().back().second-sample.getLimitsLists().front().first;
	complextype halfunitcell = Constant::pi/period;
	complextype blochstep = 2*halfunitcell/brillouinpoints;
	pair<long double,long double> guesse;

	SolverTM<complextype,complextype> solver(sample,selectedband);

	auto func = std::mem_fn(&SolverTM<complextype,xType>::Determinantnonzero);

	SolutionWannier<complextype, xType> sol;

//TODO Maybe is better to define a new integration function
	for(long int i = 0; i < levels ; i++){
		xType meanLoc = solution.getWavefunction(i)->meanLocalization();
		shared_ptr<WannierFunction<complextype, xType>> wannier = make_shared<WannierFunction<complextype, xType>>(period,blochstep);
		cout << "meanLoc=" << meanLoc << endl;
		for(long int j=-brillouinpoints/2; j <= brillouinpoints/2; j++){
			/*solver for each q*/
			long double rootp=solution.getEnergy(i).real();
			complextype a=j*blochstep*period;
			long double b=std::cos(a);
			std::function<double(double)> funcpointerb = std::bind(func, solver, std::placeholders::_1,b);

			bool notbrackted = true;
			bool guessneg = (funcpointerb(rootp)<0.0);
			long double expandone=rootp;
			double expandtwo=rootp;

			while(notbrackted){
				expandone+=(0.00001_eV);
				if((funcpointerb(expandone)<0.0)!=guessneg){
					notbrackted=false;
					continue;
				}
				expandtwo-=(0.00001_eV);
				if((funcpointerb(expandtwo)<0.0)!=guessneg)
					notbrackted=false;
			}

			rootp=BrentsRootFinder::zbrent(funcpointerb,expandtwo,expandone,(1e-50));
			shared_ptr<PlaneWavesFunction<complextype,xType>> function = solver.getWavefunction(rootp,a/period);
			/*phase correction*/
			//std::complex<complextype> funcvalue = (*function)(meanLoc);
			cout << "value[" << i << "," << j << "]=" << (*function)(0.0) << "/" << (*function)(sample.getTotalWidth()) <<  "-%=" << 100.0*std::real((*function)(0.0)-(*function)(sample.getTotalWidth()))/(*function)(0.0) << endl;
			//complextype teta;
			//teta = std::atan2(-imag(funcvalue),real(funcvalue));
			//function->addPhase(teta);

			wannier->addFunction(function,std::complex<complextype>(rootp));
		}
		sol.addLevel(wannier,wannier->getEnergy(brillouinpoints/2));
	}


	return sol;
}


template<typename complextype, typename xType>
double SolverTM<complextype,xType>::Determinantz(double energy){

	vector<Matrix<std::complex<complextype>, 2, 2>> matrixs(potential.size()+2);
	vector<std::complex<complextype>> ki(potential.size());

	Matrix<std::complex<complextype>, 2, 2> resultmatrix;

	for(long unsigned int i=0; i< potential.size() ; i++){
		ki[i]=std::sqrt(std::complex<complextype>(2.0)*std::complex<complextype>(effectivemasses[i])*std::complex<complextype>(carriersignal)*(std::complex<complextype>(potential[i])-std::complex<complextype>(energy)))/std::complex<complextype>(Constant::hbar);
		//cout << "k[" << i << "]="<< ki[i] << std::endl;
	}

	for(long unsigned int i=1; i< potential.size()+1 ; i++){
		(matrixs[i])(0,0)=std::cosh(-ki[i-1]*std::complex<complextype>(widths[i-1]));
		(matrixs[i])(0,1)=std::complex<complextype>(effectivemasses[i-1])*std::sinh(-ki[i-1]*std::complex<complextype>(widths[i-1]))/ki[i-1];
		(matrixs[i])(1,0)=ki[i-1]*std::sinh(-ki[i-1]*std::complex<complextype>(widths[i-1]))/std::complex<complextype>(effectivemasses[i-1]);
		(matrixs[i])(1,1)=std::cosh(-ki[i-1]*std::complex<complextype>(widths[i-1]));
	//cout << matrixs[i] << std::endl;
	}

	resultmatrix=matrixs[potential.size()];
	for(long int l = potential.size()-1 ;l>=1; l--){
		resultmatrix=matrixs[l]*resultmatrix;
	}

	return ((resultmatrix(0,0)+resultmatrix(1,1)).real()+(resultmatrix(0,0)+resultmatrix(1,1)).imag())/2.0;


	/*	vector<Matrix<std::complex<complextype>, 2, 2>> matrixs(potential.size()+2);
	vector<std::complex<complextype>> ki(potential.size());

	Matrix<std::complex<complextype>, 2, 2> resultmatrix;

	for(long unsigned int i=0; i< potential.size() ; i++){
		ki[i]=std::sqrt(std::complex<complextype>(2.0,0.0)*effectivemasses[i]*carriersignal*(potential[i]-energy))/Constant::hbar;
		//cout << "k[" << i << "]="<< ki[i] << std::endl;
	}

	for(long unsigned int i=1; i< potential.size()+1 ; i++){
		(matrixs[i])(0,0)=std::cosh(-ki[i-1]*widths[i-1]);
		(matrixs[i])(0,1)=effectivemasses[i-1]*std::sinh(-ki[i-1]*widths[i-1])/ki[i-1];
		(matrixs[i])(1,0)=ki[i-1]*std::sinh(-ki[i-1]*widths[i-1])/effectivemasses[i-1];
		(matrixs[i])(1,1)=std::cosh(-ki[i-1]*widths[i-1]);
	//cout << matrixs[i] << std::endl;
	}

	resultmatrix=matrixs[potential.size()];
	for(long int l = potential.size()-1 ;l>=1; l--){
		resultmatrix=matrixs[l]*resultmatrix;
	}

	//return std::real((resultmatrix(0,0)-1.0)*(resultmatrix(1,1)-1.0)-resultmatrix(1,0)*resultmatrix(0,1))+std::imag((resultmatrix(0,0)-1.0)*(resultmatrix(1,1)-1.0)-resultmatrix(1,0)*resultmatrix(0,1));
	return std::real((resultmatrix(0,0)+resultmatrix(1,1))/2.0)+std::imag((resultmatrix(0,0)+resultmatrix(1,1))/2.0);*/
}

template<typename complextype, typename xType>
double SolverTM<complextype,xType>::Determinantnonzero(double energy,double& value){

	vector<Matrix<std::complex<complextype>, 2, 2>> matrixs(potential.size()+2);
	vector<std::complex<complextype>> ki(potential.size());

	Matrix<std::complex<complextype>, 2, 2> resultmatrix;

	for(long unsigned int i=0; i< potential.size() ; i++){
		ki[i]=std::sqrt(std::complex<complextype>(2.0)*std::complex<complextype>(effectivemasses[i])*std::complex<complextype>(carriersignal)*(std::complex<complextype>(potential[i])-std::complex<complextype>(energy)))/std::complex<complextype>(Constant::hbar);
		//cout << "k[" << i << "]="<< ki[i] << std::endl;
	}

	for(long unsigned int i=1; i< potential.size()+1 ; i++){
		(matrixs[i])(0,0)=std::cosh(-ki[i-1]*std::complex<complextype>(widths[i-1]));
		(matrixs[i])(0,1)=std::complex<complextype>(effectivemasses[i-1])*std::sinh(-ki[i-1]*std::complex<complextype>(widths[i-1]))/ki[i-1];
		(matrixs[i])(1,0)=ki[i-1]*std::sinh(-ki[i-1]*std::complex<complextype>(widths[i-1]))/std::complex<complextype>(effectivemasses[i-1]);
		(matrixs[i])(1,1)=std::cosh(-ki[i-1]*std::complex<complextype>(widths[i-1]));
	//cout << matrixs[i] << std::endl;
	}

	resultmatrix=matrixs[potential.size()];
	for(long int l = potential.size()-1 ;l>=1; l--){
		resultmatrix=matrixs[l]*resultmatrix;
	}

	return (((resultmatrix(0,0)+resultmatrix(1,1))/2.0-value).real()+((resultmatrix(0,0)+resultmatrix(1,1))/2.0-value).imag());



	/*vector<Matrix<std::complex<complextype>, 2, 2>> matrixs(potential.size()+2);
	vector<std::complex<complextype>> ki(potential.size());

	Matrix<std::complex<complextype>, 2, 2> resultmatrix;

	for(long unsigned int i=0; i< potential.size() ; i++){
		ki[i]=std::sqrt(std::complex<complextype>(2.0,0.0)*effectivemasses[i]*carriersignal*(potential[i]-energy))/Constant::hbar;
	}

	for(long unsigned int i=1; i< potential.size()+1 ; i++){
		(matrixs[i])(0,0)=std::cosh(-ki[i-1]*widths[i-1]);
		(matrixs[i])(0,1)=effectivemasses[i-1]*std::sinh(-ki[i-1]*widths[i-1])/ki[i-1];
		(matrixs[i])(1,0)=ki[i-1]*std::sinh(-ki[i-1]*widths[i-1])/effectivemasses[i-1];
		(matrixs[i])(1,1)=std::cosh(-ki[i-1]*widths[i-1]);
	}

	resultmatrix=matrixs[potential.size()];
	for(long int l = potential.size()-1 ;l>=1; l--){
		resultmatrix=matrixs[l]*resultmatrix;
	}

	//return std::real((resultmatrix(0,0)-std::exp(std::complex<complextype>(0.0,1.0)*value))*(resultmatrix(1,1)-std::exp(std::complex<complextype>(0.0,1.0)*value))-resultmatrix(1,0)*resultmatrix(0,1))+std::imag((resultmatrix(0,0)-std::exp(std::complex<complextype>(0.0,1.0)*value))*(resultmatrix(1,1)-std::exp(std::complex<complextype>(0.0,1.0)*value))-resultmatrix(1,0)*resultmatrix(0,1));
	return std::real((resultmatrix(0,0)+resultmatrix(1,1))/2.0-value)+std::imag((resultmatrix(0,0)+resultmatrix(1,1))/2.0-value);*/
}


template<typename complextype, typename xType>
shared_ptr<PlaneWavesFunction<complextype,xType>> SolverTM<complextype,xType>::getWavefunction(double energy,complextype bloch){

	vector<Matrix<std::complex<complextype>, 2, 2>> matrixs(potential.size()+1);
		vector<Matrix<std::complex<complextype>, 2, 2>> complementmatrix(potential.size());
		vector<Matrix<std::complex<complextype>, 2, 2>> resultmatrix(potential.size());
		vector<std::complex<complextype>> ki(potential.size());
		vector<Matrix<std::complex<complextype>, 2, 1> > coefs(potential.size());
		Matrix<std::complex<complextype>, 2, 2> result;


		complextype periodsize = limitslst.back().second-limitslst.front().first;

		for(long unsigned int i=0; i< potential.size() ; i++){
			ki[i]=std::sqrt(std::complex<complextype>(2.0)*std::complex<complextype>(effectivemasses[i])*std::complex<complextype>(carriersignal)*(std::complex<complextype>(potential[i])-std::complex<complextype>(energy)))/std::complex<complextype>(Constant::hbar);
			//cout << "k[" << i << "]="<< ki[i] << std::endl;
		}

		for(long unsigned int i=1; i< potential.size() ; i++){
			(matrixs[i])(0,0)=std::cosh(-ki[i-1]*std::complex<complextype>(widths[i-1]));
			(matrixs[i])(0,1)=std::complex<complextype>(effectivemasses[i-1])*std::sinh(-ki[i-1]*std::complex<complextype>(widths[i-1]))/ki[i-1];
			(matrixs[i])(1,0)=ki[i-1]*std::sinh(-ki[i-1]*std::complex<complextype>(widths[i-1]))/std::complex<complextype>(effectivemasses[i-1]);
			(matrixs[i])(1,1)=std::cosh(-ki[i-1]*std::complex<complextype>(widths[i-1]));
		//cout << matrixs[i] << std::endl;
		}

		for(long unsigned int i=0; i< potential.size() ; i++){
			(complementmatrix[i])(0,0)=std::complex<complextype>(0.5)*std::exp(-ki[i]*std::complex<complextype>(limitslst[i].second));
			(complementmatrix[i])(0,1)=std::complex<complextype>(0.5)*std::complex<complextype>(effectivemasses[i])*std::exp(-ki[i]*std::complex<complextype>(limitslst[i].second))/ki[i];
			(complementmatrix[i])(1,0)=std::complex<complextype>(0.5)*std::exp(ki[i]*std::complex<complextype>(limitslst[i].second));
			(complementmatrix[i])(1,1)=std::complex<complextype>(-0.5)*std::complex<complextype>(effectivemasses[i])*std::exp(ki[i]*std::complex<complextype>(limitslst[i].second))/ki[i];
		//cout << matrixs[i] << std::endl;
		}

		xType z1=limitslst.back().first;
		xType z2=limitslst.back().second;

		(matrixs[potential.size()])(0,0)=std::exp(ki[potential.size()-1]*std::complex<complextype>(z1));
		(matrixs[potential.size()])(0,1)=std::exp(-ki[potential.size()-1]*std::complex<complextype>(z1));
		(matrixs[potential.size()])(1,0)=ki[potential.size()-1]*std::exp(ki[potential.size()-1]*std::complex<complextype>(z1))/std::complex<complextype>(effectivemasses[potential.size()-1]);
		(matrixs[potential.size()])(1,1)=-ki[potential.size()-1]*std::exp(-ki[potential.size()-1]*std::complex<complextype>(z1))/std::complex<complextype>(effectivemasses[potential.size()-1]);

		(matrixs[0])(0,0)=std::complex<complextype>(0.5)*std::exp(-ki[potential.size()-1]*std::complex<complextype>(z2));
		(matrixs[0])(0,1)=std::complex<complextype>(0.5)*std::complex<complextype>(effectivemasses[potential.size()-1])*std::exp(-ki[potential.size()-1]*std::complex<complextype>(z2))/ki[potential.size()-1];
		(matrixs[0])(1,0)=std::complex<complextype>(0.5)*std::exp(ki[potential.size()-1]*std::complex<complextype>(z2));
		(matrixs[0])(1,1)=std::complex<complextype>(-0.5)*std::complex<complextype>(effectivemasses[potential.size()-1])*std::exp(ki[potential.size()-1]*std::complex<complextype>(z2))/ki[potential.size()-1];

		result=matrixs[potential.size()];
		for(long int l = potential.size()-1;l>=0; l--){
			result=matrixs[l]*result;
		}


		for(long int l = potential.size()-2;l>=0; l--){
			resultmatrix[l]=matrixs[potential.size()];
			for(long int g = potential.size()-1; g>(l+1); g-- ){
				resultmatrix[l]=matrixs[g]*resultmatrix[l];
			}
			resultmatrix[l]=complementmatrix[l]*resultmatrix[l];
		}

		coefs[potential.size()-1]=Matrix<std::complex<complextype>, 2, 1>(std::complex<complextype>(1.0),std::complex<complextype>(-1.0)*(result(0,0)-std::exp(std::complex<complextype>(0.0,-1.0)*std::complex<complextype>(bloch)*std::complex<complextype>(periodsize)))/result(0,1));
		for(long int h=potential.size()-2;h>=0;h--){
			coefs[h]=resultmatrix[h]*coefs[potential.size()-1];
		}

		shared_ptr<PlaneWavesFunction<complextype,xType>> function = make_shared<PlaneWavesFunction<complextype,xType>>(potential.size(),limitslst);


		for(long unsigned int k=0; k< potential.size() ; k++){
			(function->data[k]).push_back(make_pair((coefs[k])(0),ki[k]-std::complex<complextype>(0.0,bloch)));
			(function->data[k]).push_back(make_pair((coefs[k])(1),-ki[k]-std::complex<complextype>(0.0,bloch)));
		}

		function->normalize(1.0);
		function->setPeriodic();

		return function;

	/*	vector<Matrix<std::complex<complextype>, 2, 2>> matrixs(potential.size()+1);
	vector<Matrix<std::complex<complextype>, 2, 2>> complementmatrix(potential.size());
	vector<Matrix<std::complex<complextype>, 2, 2>> resultmatrix(potential.size());
	vector<std::complex<complextype>> ki(potential.size());
	vector<Matrix<std::complex<complextype>, 2, 1> > coefs(potential.size());
	Matrix<std::complex<complextype>, 2, 2> result;


	complextype periodsize = limitslst.back().second-limitslst.front().first;

	for(long unsigned int i=0; i< potential.size() ; i++){
		ki[i]=std::sqrt(std::complex<complextype>(2.0,0.0)*effectivemasses[i]*carriersignal*(potential[i]-energy))/Constant::hbar;
		//cout << "k[" << i << "]="<< ki[i] << std::endl;
	}

	for(long unsigned int i=1; i< potential.size() ; i++){
		(matrixs[i])(0,0)=std::cosh(-ki[i-1]*widths[i-1]);
		(matrixs[i])(0,1)=effectivemasses[i-1]*std::sinh(-ki[i-1]*widths[i-1])/ki[i-1];
		(matrixs[i])(1,0)=ki[i-1]*std::sinh(-ki[i-1]*widths[i-1])/effectivemasses[i-1];
		(matrixs[i])(1,1)=std::cosh(-ki[i-1]*widths[i-1]);
	//cout << matrixs[i] << std::endl;
	}

	for(long unsigned int i=0; i< potential.size() ; i++){
		(complementmatrix[i])(0,0)=0.5*std::exp(-ki[i]*limitslst[i].second);
		(complementmatrix[i])(0,1)=0.5*effectivemasses[i]*std::exp(-ki[i]*limitslst[i].second)/ki[i];
		(complementmatrix[i])(1,0)=0.5*std::exp(ki[i]*limitslst[i].second);
		(complementmatrix[i])(1,1)=-0.5*effectivemasses[i]*std::exp(ki[i]*limitslst[i].second)/ki[i];
	//cout << matrixs[i] << std::endl;
	}

	xType z1=limitslst.back().first;
	xType z2=limitslst.back().second;

	(matrixs[potential.size()])(0,0)=std::exp(ki[potential.size()-1]*z1);
	(matrixs[potential.size()])(0,1)=std::exp(-ki[potential.size()-1]*z1);
	(matrixs[potential.size()])(1,0)=ki[potential.size()-1]*std::exp(ki[potential.size()-1]*z1)/effectivemasses[potential.size()-1];
	(matrixs[potential.size()])(1,1)=-ki[potential.size()-1]*std::exp(-ki[potential.size()-1]*z1)/effectivemasses[potential.size()-1];

	(matrixs[0])(0,0)=0.5*std::exp(-ki[potential.size()-1]*z2);
	(matrixs[0])(0,1)=0.5*effectivemasses[potential.size()-1]*std::exp(-ki[potential.size()-1]*z2)/ki[potential.size()-1];
	(matrixs[0])(1,0)=0.5*std::exp(ki[potential.size()-1]*z2);
	(matrixs[0])(1,1)=-0.5*effectivemasses[potential.size()-1]*std::exp(ki[potential.size()-1]*z2)/ki[potential.size()-1];

	result=matrixs[potential.size()];
	for(long int l = potential.size()-1;l>=0; l--){
		result=matrixs[l]*result;
	}


	for(long int l = potential.size()-2;l>=0; l--){
		resultmatrix[l]=matrixs[potential.size()];
		for(long int g = potential.size()-1; g>(l+1); g-- ){
			resultmatrix[l]=matrixs[g]*resultmatrix[l];
		}
		resultmatrix[l]=complementmatrix[l]*resultmatrix[l];
	}

	coefs[potential.size()-1]=Matrix<std::complex<complextype>, 2, 1>(std::complex<complextype>(1.0,0.0),std::complex<complextype>(-1.0,0.0)*(result(0,0)-std::exp(std::complex<complextype>(0.0,-1.0)*bloch*periodsize))/result(0,1));
	for(long int h=potential.size()-2;h>=0;h--){
		coefs[h]=resultmatrix[h]*coefs[potential.size()-1];
	}

	shared_ptr<PlaneWavesFunction<complextype,xType>> function = make_shared<PlaneWavesFunction<complextype,xType>>(potential.size(),limitslst);


	for(long unsigned int k=0; k< potential.size() ; k++){
		(function->data[k]).push_back(make_pair((coefs[k])(0),ki[k]-std::complex<complextype>(0.0,bloch)));
		(function->data[k]).push_back(make_pair((coefs[k])(1),-ki[k]-std::complex<complextype>(0.0,bloch)));
	}

	function->normalize(1.0);
	function->setPeriodic();

	return function;*/

}


/*
template<typename complextype, typename xType>
complextype SolverTM<complextype,xType>::Jp(complextype energy){

	vector<std::complex<complextype>> ki(potential.size());

	vector<std::complex<complextype>> fi(potential.size());
	vector<std::complex<complextype>> gi(potential.size());
	vector<std::complex<complextype>> hi(potential.size());

	for(long unsigned int i=0; i< potential.size() ; i++){
		ki[i]=std::sqrt(std::complex<complextype>(2.0,0.0)*effectivemasses[i]*(potential[i]-energy))/Constant::hbar;
	}

	for(long unsigned int i=0; i< potential.size() ; i++){
		fi[i]=1/cos(ki[i]*widths[i]);
		gi[i]=ki[i]*tan(ki[i]*widths[i])/effectivemasses[i];
		hi[i]=effectivemasses[i]*tan(ki[i]*widths[i])/ki[i];
	}


	std::complex<complextype> S1n(0.0,0.0);
	vector<std::complex<complextype>> S1p(potential.size());
	vector<std::complex<complextype>> Spn(potential.size());

	for(long unsigned int i=0; i< potential.size() ; i++){
		fi[i]=1/cos(ki[i]*widths[i]);
		gi[i]=ki[i]*tan(ki[i]*widths[i])/effectivemasses[i];
		hi[i]=effectivemasses[i]*tan(ki[i]*widths[i])/ki[i];
	}
*/
/*	for(long int s=0;s<=potential.size()-1;s++){
		vector<long int> i2s(s);
		vector<long int> tops(s);
		i2s[0]=p+s;
		tops[0]=i2s[0]-1
		for(long int k=1; k< s;k++){
			i2s[k]=p+s-k;
			tops[k]=i2s[k-1]-1;
		}
		long int currinc=s-1;
		bool incrementing=true;
		while(incrementing){
			std::complex<complextype> product(1.0,0.0);
			for(long int u=1;u<=s;u++){
				product*=-Lpq(i2s[u],i2s[u+1],fi,hi,gi);
			}
			sum+=product;
			if(i2s[currinc]!=tops[currinc]){
				i2s[currinc]++;
			}
			else{
				if(currinc==0)
					incrementing=false;
				currinc--;
			}
		}
	}

*//*
}

template<typename complextype, typename xType>
inline complextype SolverTM<complextype,xType>::Lpq(long int p, long int q, vector<std::complex<complextype>>& fi, vector<std::complex<complextype>>& hi , vector<std::complex<complextype>>& gi){
	auto product=1;
	for(long unsigned int i= p; i <=q ; i++)
		product*=fi[i]*fi[i];
	return hi[p]*gi[q]*product;
}

template<typename complextype, typename xType>
complextype SolverTM<complextype,xType>::Jm(complextype energy){

	return 0.0;
}
*/
} /* namespace epital */



#endif
