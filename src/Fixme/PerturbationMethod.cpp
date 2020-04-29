/*
 * PerturbationMethod.cpp
 *
 *  Created on: Apr 15, 2013
 *      Author: marcel
 */


#include "PerturbationMethod.hpp"


#ifndef PERTURBATIONMETHOD_CPP_
#define PERTURBATIONMETHOD_CPP_

#include <vector>

namespace epital {

template <typename timetype, typename complextype>
PerturbationMethod<timetype,complextype>::PerturbationMethod(PerturbationMethod::Carrier carrier) {
	carriertype=carrier;

}

template <typename timetype, typename complextype>
PerturbationMethod<timetype,complextype>::~PerturbationMethod() {

}


template <typename timetype, typename complextype>
template<typename yType, typename xType>
std::complex<complextype> PerturbationMethod<timetype,complextype>::Vab(Solution<complextype,xType>& solution, DiscreteFunction<yType,xType>& effectivemass,yType blochvector ,long int a, long int b){


	DiscreteFunction<std::complex<complextype>,xType> ket_a(solution.getWavefunction(b));
	DiscreteFunction<std::complex<complextype>,xType> ket_b(solution.getWavefunction(b));
	DiscreteFunction<std::complex<complextype>,xType> ket_c(solution.getWavefunction(b));

	auto datab = effectivemass.getdata();

	ket_a.Differential();
	auto datac = ket_a.getdata();
	#pragma omp parallel for
	for(long int i = 0; i < ket_a.getGrid().getSize(); i++){
		*(datac+i)/=(*(datab+i));
	}
	ket_a*=(-Constant::hbar*Constant::hbar*blochvector*std::complex<complextype>(0.0,1.0)/2.0);

	auto datad = ket_b.getdata();
	#pragma omp parallel for
	for(long int i = 0; i < ket_b.getGrid().getSize(); i++){
		*(datad+i)/=(*(datab+i));
	}
	ket_b*=(Constant::hbar*Constant::hbar*blochvector*blochvector/2.0);


	auto datae = ket_c.getdata();
	#pragma omp parallel for
	for(long int i = 0; i < ket_c.getGrid().getSize(); i++){
		*(datae+i)/=(*(datab+i));
	}
	ket_c.Differential();
	ket_c*=(-Constant::hbar*Constant::hbar*blochvector*std::complex<complextype>(0.0,1.0)/2.0);

	return (solution.getWavefunction(a)->Dot(ket_a)+solution.getWavefunction(a)->Dot(ket_b)+solution.getWavefunction(a)->Dot(ket_c));

}

template <typename timetype, typename complextype>
template<typename xType>
std::complex<complextype> PerturbationMethod<timetype,complextype>::Eab(Solution<complextype,xType>& solution,long int a, long int b){
	return (solution.getEnergy(a)-solution.getEnergy(b));
}


template <typename timetype, typename complextype>
template<typename yType, typename xType>
std::complex<complextype> PerturbationMethod<timetype,complextype>::EnergyBrilluoin(long int level, Solution<complextype,xType> solution, DiscreteFunction<yType,xType>& effectivemass,yType blochvector){
	std::complex<complextype> energysign;

	if(carriertype==Carrier::ELECTRON)
		energysign=std::complex<complextype>(1.0,0);
	else
		energysign=std::complex<complextype>(-1.0,0);

	std::complex<complextype> elements[200][200];
	std::complex<complextype> deltas[200][200];

	for(long int i=0;i<solution.getLevels();i++){
		for(long int j=0;j<solution.getLevels();j++)
			elements[i][j]=Vab(solution,effectivemass,blochvector,i,j);
	}

	for(long int i=0;i<solution.getLevels();i++){
		for(long int j=0;j<solution.getLevels();j++)
			deltas[i][j]=Eab(solution,i,j);
	}



	auto firstorder = elements[level][level];

	std::complex<complextype> secondorder = 0.0;
	for(long int k2 = 0; k2<solution.getLevels(); k2++  ){
		if(k2!=level){
			secondorder += (elements[k2][level]*elements[level][k2])/deltas[level][k2];
		}
	}

	std::complex<complextype> thirdorder = 0.0;
	for(long int k3 = 0; k3<solution.getLevels(); k3++){
		if(k3!=level){
			thirdorder -= firstorder*elements[k3][level]*elements[level][k3]/(deltas[level][k3]*deltas[level][k3]);
			for(long int k2 = 0; k2 < solution.getLevels(); k2++  ){
				if(k2!=level){
					thirdorder += elements[level][k3]*elements[k3][k2]*elements[k2][level]/
							(deltas[level][k2]*deltas[level][k3]);
				}
			}
		}
	}

	std::complex<complextype> forthorder = 0.0;

	for(long int k4 = 0; k4<solution.getLevels(); k4++ ){
		if(k4!=level){
			forthorder -= secondorder*elements[level][k4]*elements[k4][level]/(deltas[level][k4]*deltas[level][k4]);
			forthorder += firstorder*firstorder*elements[level][k4]*elements[k4][level]/(deltas[level][k4]*deltas[level][k4]*deltas[level][k4]);
			for(long int k3 = 0; k3<solution.getLevels(); k3++ ){
				if(k3!=level){
					forthorder -= 2.0*firstorder*elements[level][k4]*elements[k4][k3]*elements[k3][level]/
							(deltas[level][k3]*deltas[level][k3]*deltas[level][k4]);
					for(long int k2 = 0; k2<solution.getLevels(); k2++ ){
						if(k2!=level){
							forthorder += elements[level][k4]*elements[k4][k3]*elements[k3][k2]*elements[k2][level]/
									(deltas[level][k2]*deltas[level][k3]*deltas[level][k4]);
						}
					}
				}
			}
		}
	}


	std::complex<complextype> fifthorder = 0.0;

	for(long int k5 = 0; k5<solution.getLevels(); k5++ ){
		if(k5!=level){
			fifthorder -= 2.0*firstorder*2.0*secondorder*elements[level][k5]*elements[k5][level]/
					(deltas[level][k5]*deltas[level][k5]*deltas[level][k5]);
			fifthorder -= firstorder*firstorder*firstorder*elements[level][k5]*elements[k5][level]/
					(deltas[level][k5]*deltas[level][k5]*deltas[level][k5]*deltas[level][k5]);
			for(long int k4 = 0; k4<solution.getLevels(); k4++ ){
				if(k4!=level){
					fifthorder += firstorder*firstorder*2.0*elements[level][k5]*elements[k5][k4]*elements[k4][level]/
							(deltas[level][k4]*deltas[level][k4]*deltas[level][k4]*deltas[level][k5]);
					fifthorder -= 2.0*secondorder*elements[level][k5]*elements[k5][k4]*elements[k4][level]/
							(deltas[level][k4]*deltas[level][k4]*deltas[level][k5]);
					for(long int k3 = 0; k3<solution.getLevels(); k3++ ){
						if(k3!=level){
							fifthorder += firstorder*firstorder*elements[level][k5]*elements[k5][k3]*elements[k3][level]/
									(deltas[level][k5]*deltas[level][k5]*deltas[level][k3]*deltas[level][k3]);
							fifthorder -= 2.0*firstorder*elements[level][k5]*elements[k5][level]*elements[level][k3]*elements[k3][level]/
									(deltas[level][k3]*deltas[level][k3]*deltas[level][k5]*deltas[level][k5]);
							fifthorder -= 2.0*firstorder*elements[level][k5]*elements[k5][k4]*elements[k4][k3]*elements[k3][level]/
									(deltas[level][k3]*deltas[level][k3]*deltas[level][k4]*deltas[level][k5]);
							for(long int k2 = 0; k2<solution.getLevels(); k2++ ){
								if(k2!=level){
									fifthorder += elements[level][k5]*elements[k5][k4]*elements[k3][k4]*elements[k3][k2]*elements[k2][level]/
											(deltas[level][k2]*deltas[level][k3]*deltas[level][k4]*deltas[level][k5]);
									fifthorder -= 2.0*firstorder*elements[level][k5]*elements[k5][k4]*elements[k4][k2]*elements[k2][level]/
													(deltas[level][k2]*deltas[level][k4]*deltas[level][k4]*deltas[level][k5]);
									fifthorder -= elements[level][k5]*elements[k5][level]*elements[level][k3]*elements[k3][k2]*elements[k2][level]/
													(deltas[level][k5]*deltas[level][k5]*deltas[level][k2]*deltas[level][k3]);
								}
							}
						}
					}
				}
			}
		}
	}


	return (firstorder+secondorder+thirdorder+forthorder/*+fifthorder*/);
}


} /* namespace epital */

#endif
