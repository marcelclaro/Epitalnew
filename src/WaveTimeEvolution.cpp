/*
 * WaveTimeEvolution.cpp
 *
 *  Created on: Aug 21, 2012
 *      Author: marcel
 */


#include <vector>
#include <fftw3.h>
#include "Graphics.hpp"
#include "WaveTimeEvolution.hpp"

#ifndef WAVETIMEEVOLUTION_CPP_
#define WAVETIMEEVOLUTION_CPP_

namespace epital{

template <typename timetype, typename complextype>
WaveTimeEvolution<timetype,complextype>::WaveTimeEvolution(WaveTimeEvolution::Carrier carrier) {
	carriertype=carrier;
}

template <typename timetype, typename complextype>
WaveTimeEvolution<timetype,complextype>::~WaveTimeEvolution() {
	// TODO Auto-generated destructor stub
}

template <typename timetype, typename complextype>
template<typename yType, typename xType>
DiscreteFunction<std::complex<complextype>,xType>& WaveTimeEvolution<timetype,complextype>::evolute_finite(timetype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential, DiscreteFunction<yType,xType>& effectivemass){

	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	//get time propagator potential factor
	DiscreteFunction<std::complex<complextype>,xType> potpropag = potentialpropagator(methodtimestep,potential);

	//for each time step apply time propagator
	for(long int i = 0 ; i< timesteps; i++){
		initial*=potpropag;
		kinect_propag_evolution_finite<yType,xType>(methodtimestep,initial,effectivemass);
		initial*=potpropag;
	}

	return initial;
}

template <typename timetype, typename complextype>
template<typename yType, typename xType>
DiscreteFunction<std::complex<complextype>,xType>& WaveTimeEvolution<timetype,complextype>::evolute_periodic(timetype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential, DiscreteFunction<yType,xType>& effectivemass){

	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	//get time propagator potential factor
	DiscreteFunction<std::complex<complextype>,xType> potpropag = potentialpropagator(methodtimestep,potential);

	//for each time step apply time propagator
	for(long int i = 0 ; i< timesteps; i++){
		initial*=potpropag;
		kinect_propag_evolution_periodic<yType,xType>(methodtimestep,initial,effectivemass);
		initial*=potpropag;
	}

	return initial;
}


template <typename timetype, typename complextype>
template<typename yType, typename xType>
DiscreteFunction<std::complex<complextype>,xType>& WaveTimeEvolution<timetype,complextype>::groundlevel_finite(timetype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential, DiscreteFunction<yType,xType>& effectivemass){
	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	std::complex<complextype> modulus = initial.Modulus();

	DiscreteFunction<std::complex<complextype>,xType> potpropag = potentialpropagator(methodtimestep,potential);

	for(long int i = 0 ; i< timesteps; i++){
		initial*=potpropag;
		kinect_propag_evolution_finite<yType,xType>(methodtimestep,initial,effectivemass);
		initial*=potpropag;
		initial*=modulus/initial.Modulus();
	}

	return initial;
}

template <typename timetype, typename complextype>
template<typename yType, typename xType>
DiscreteFunction<std::complex<complextype>,xType>& WaveTimeEvolution<timetype,complextype>::groundlevel_periodic(timetype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential, DiscreteFunction<yType,xType>& effectivemass){
	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	std::complex<complextype> modulus = initial.Modulus();

	DiscreteFunction<std::complex<complextype>,xType> potpropag = potentialpropagator(methodtimestep,potential);
	//DiscreteFunction<std::complex<complextype>,xType> bloch = fftmulti(blochvector,methodtimestep,initial.getGrid(),effectivemass(0));

	for(long int i = 0 ; i< timesteps; i++){
		initial*=potpropag;
		kinect_propag_evolution_periodic<yType,xType>(methodtimestep,initial,effectivemass);
		initial*=potpropag;
		initial*=modulus/initial.Modulus();
	}

	return initial;
}

template <typename timetype, typename complextype>
template<typename yType, typename xType>
DiscreteFunction<std::complex<complextype>,xType>& WaveTimeEvolution<timetype,complextype>::groundlevel_periodicFFT(timetype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential, DiscreteFunction<yType,xType>& effectivemass){
	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	std::complex<complextype> modulus = initial.Modulus();

	DiscreteFunction<std::complex<complextype>,xType> potpropag = potentialpropagator(methodtimestep,potential);
	DiscreteFunction<std::complex<complextype>,xType> bloch = fftmulti(0.0,methodtimestep,initial.getGrid(),effectivemass(0));

	for(long int i = 0 ; i< timesteps; i++){
		initial*=potpropag;
		initial.FFT();
		//kinect_propag_evolution_periodic<yType,xType>(methodtimestep,initial,effectivemass);
		initial*=bloch;
		initial.FFT();
		initial*=potpropag;
		initial*=modulus/initial.Modulus();
	}

	return initial;
}


template <typename timetype, typename complextype>
template<typename yType, typename xType>
DiscreteFunction<std::complex<complextype>,xType>& WaveTimeEvolution<timetype,complextype>::newlevel_finite(timetype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential, vector<shared_ptr<DiscreteFunction<std::complex<complextype>,xType>>> otherslevels, DiscreteFunction<yType,xType>& effectivemass){
	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	std::complex<complextype> modulus = initial.Modulus();

	DiscreteFunction<std::complex<complextype>,xType> potpropag = potentialpropagator(methodtimestep,potential);

	for(long int i = 0 ; i< timesteps; i++){

		for(auto k : otherslevels){
			std::complex<complextype> dotprod = (*k).Dot(initial);
			long int pos = 0;
			for(auto& i : initial){
				i-=dotprod*((*k)(pos));
				++pos;
			}
		}


		initial*=potpropag;
		kinect_propag_evolution_finite<yType,xType>(methodtimestep,initial,effectivemass);
		initial*=potpropag;
		initial*=modulus/initial.Modulus();
	}

	return initial;
}

template <typename timetype, typename complextype>
template<typename yType, typename xType>
DiscreteFunction<std::complex<complextype>,xType>& WaveTimeEvolution<timetype,complextype>::newlevel_periodicFFT(timetype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential, vector<shared_ptr<DiscreteFunction<std::complex<complextype>,xType>>> otherslevels, DiscreteFunction<yType,xType>& effectivemass){
	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	std::complex<complextype> modulus = initial.Modulus();


	DiscreteFunction<std::complex<complextype>,xType> potpropag = potentialpropagator(methodtimestep,potential);
	DiscreteFunction<std::complex<complextype>,xType> bloch = fftmulti(0.0,methodtimestep,initial.getGrid(),effectivemass(0));

	for(long int i = 0 ; i< timesteps; i++){

		for(auto k : otherslevels){
			std::complex<complextype> dotprod = (*k).Dot(initial);
			long int pos = 0;
			for(auto& i : initial){
				i-=dotprod*((*k)(pos));
				++pos;
			}
		}

		initial*=potpropag;
		initial.FFT();
		//kinect_propag_evolution_periodic<yType,xType>(methodtimestep,initial,effectivemass);
		initial*=bloch;
		initial.FFT();
		initial*=potpropag;
		initial*=modulus/initial.Modulus();
	}


	return initial;
}


template <typename timetype, typename complextype>
template<typename yType, typename xType>
DiscreteFunction<std::complex<complextype>,xType>& WaveTimeEvolution<timetype,complextype>::newlevel_periodic(timetype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential, vector<shared_ptr<DiscreteFunction<std::complex<complextype>,xType>>> otherslevels, DiscreteFunction<yType,xType>& effectivemass){
	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	std::complex<complextype> modulus = initial.Modulus();


	DiscreteFunction<std::complex<complextype>,xType> potpropag = potentialpropagator(methodtimestep,potential);

	for(long int i = 0 ; i< timesteps; i++){

		for(auto k : otherslevels){
			std::complex<complextype> dotprod = (*k).Dot(initial);
			long int pos = 0;
			for(auto& i : initial){
				i-=dotprod*((*k)(pos));
				++pos;
			}
		}

		initial*=potpropag;
		kinect_propag_evolution_periodic<yType,xType>(methodtimestep,initial,effectivemass);
		initial*=potpropag;
		initial*=modulus/initial.Modulus();
	}


	return initial;
}


template <typename timetype, typename complextype>
template<typename yType, typename xType>
DiscreteFunction<std::complex<complextype>,complextype> WaveTimeEvolution<timetype,complextype>::correlationfunction_periodic(complextype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential, DiscreteFunction<yType,xType>& effectivemass){


	DiscreteFunction<std::complex<complextype>,xType> toevolute = initial;
	std::shared_ptr<Grid1D<complextype>> grid(new Grid1D<complextype>(-3.141592654/(time/timesteps),3.141592654/(time/timesteps),timesteps));
	DiscreteFunction<std::complex<complextype>,complextype> corr_func(*grid);

	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	DiscreteFunction<std::complex<complextype>,xType> potpropag = potentialpropagator(methodtimestep,potential);

	for(long int i = 0 ; i< timesteps; i++){
		toevolute*=potpropag;
		kinect_propag_evolution_periodic<yType,xType>(methodtimestep,toevolute,effectivemass);
		toevolute*=potpropag;
		corr_func[i]=initial.Dot(toevolute);
	}

	return corr_func;
}


template <typename timetype, typename complextype>
template<typename yType, typename xType>
DiscreteFunction<std::complex<complextype>,complextype> WaveTimeEvolution<timetype,complextype>::correlationfunction_finite(complextype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential, DiscreteFunction<yType,xType>& effectivemass){

	DiscreteFunction<std::complex<complextype>,xType> toevolute = initial;
	std::shared_ptr<Grid1D<complextype>> grid(new Grid1D<complextype>(-3.141592654/(time/timesteps),3.141592654/(time/timesteps),timesteps));
	DiscreteFunction<std::complex<complextype>,complextype> corr_func(*grid);

	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	DiscreteFunction<std::complex<complextype>,xType> potpropag = potentialpropagator(methodtimestep,potential);

	for(long int i = 0 ; i< timesteps; i++){
		toevolute*=potpropag;
		kinect_propag_evolution_finite<yType,xType>(methodtimestep,toevolute,effectivemass);
		toevolute*=potpropag;
		corr_func[i]=initial.Dot(toevolute);
	}

	return corr_func;
}


template <typename timetype, typename complextype>
template<typename yType, typename xType>
DiscreteFunction<std::complex<complextype>,xType>& WaveTimeEvolution<timetype,complextype>::kinect_propag_evolution_finite(timetype methodtimestep, DiscreteFunction<std::complex<complextype>,xType>& wavefunction, DiscreteFunction<yType,xType>& effectivemass){


	/*
	 * R
	 * E
	 * D
	 * O
	 * */


	//defining beta...

	xType gridstep = wavefunction.getGrid().getIncrement();
	DiscreteFunction<std::complex<complextype>,xType> beta(wavefunction.getGrid());

	for(long int i = 0; i< beta.getGrid().getSize();i++){
		beta[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*methodtimestep/(4.0*effectivemass(i)*gridstep*gridstep);
	}

	DiscreteFunction<std::complex<complextype>,xType> copy(wavefunction);

	Grid::Position last = wavefunction.getGrid().getSize()-1;


	//Multiply vector
	wavefunction[0]=beta(0)*copy(1)+(1.0-2.0*beta(0))*copy(0);
	wavefunction[last]=beta(last)*copy(last-1)+(1.0-2.0*beta(last))*copy(last);

	#pragma omp parallel for
	for(long int k = 1; k < last; k++){
		wavefunction[k]=beta(k)*copy(k-1)+(1.0-2.0*beta(k))*copy(k)+beta(k)*copy(k+1);
	}

	DiscreteFunction<std::complex<complextype>,xType> as(wavefunction.getGrid());
	DiscreteFunction<std::complex<complextype>,xType> bs(wavefunction.getGrid());
	DiscreteFunction<std::complex<complextype>,xType> cs(wavefunction.getGrid());


	as[0]=0;
	bs[0]=2.0*beta(0)+1.0;
	cs[0]=-beta(0);

    #pragma omp parallel for
	for(long int i = 1; i < last ; i++){
		as[i]=-beta(i);
		bs[i]=2.0*beta(i)+1.0;
		cs[i]=-beta(i);
	}

	as[last]=-beta(last);
	bs[last]=2.0*beta(last)+1.0;
	cs[last]=0;

	//Thomas Algorithm

    for (long int i = 1; i <= last; i++)
    {
            std::complex<complextype> m = as(i)/bs(i-1);
            bs[i]-=m*cs(i-1);
            wavefunction[i]-=m*wavefunction(i-1);
    }

    wavefunction[last]/=bs[last];
    for (long int i = (last - 1); i >= 0; --i)
            wavefunction[i]=(wavefunction(i)-cs(i)*wavefunction(i+1))/bs(i);

    return wavefunction;

}

template <typename timetype, typename complextype>
template<typename yType, typename xType>
DiscreteFunction<std::complex<complextype>,xType>& WaveTimeEvolution<timetype,complextype>::kinect_propag_evolution_periodic(timetype methodtimestep, DiscreteFunction<std::complex<complextype>,xType>& wavefunction, DiscreteFunction<yType,xType>& effectivemass){


	/*
	 * R
	 * E
	 * D
	 * O
	 * */



	//defining beta...

	xType gridstep = wavefunction.getGrid().getIncrement();
	DiscreteFunction<std::complex<complextype>,xType> beta(wavefunction.getGrid());

	for(long int i = 0; i< beta.getGrid().getSize();i++){
		beta[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*methodtimestep/(4.0*effectivemass(i)*gridstep*gridstep);
	}

	Grid::Position last = wavefunction.getGrid().getSize()-1;

	DiscreteFunction<std::complex<complextype>,xType> copy(wavefunction);
	DiscreteFunction<std::complex<complextype>,xType> uvec(wavefunction.getGrid());


	std::complex<complextype> gammaconst=-beta(0);


	//Multiply vector
	wavefunction[0]=beta(0)*(1.0+effectivemass(0)*(1.0/effectivemass(1) - 1.0/effectivemass(last) )/4.0)*copy(1)+(1.0-(2.0)*beta(0))*copy(0)+beta(0)*(1.0-effectivemass(0)*(1.0/effectivemass(1) - 1.0/effectivemass(last))/4.0)*copy(last);
	wavefunction[last]=beta(last)*(1.0-effectivemass(last)*(1.0/effectivemass(0) - 1.0/effectivemass(last-1))/4.0)*copy(last-1)+(1.0-(2.0)*beta(last))*copy(last)+beta(last)*(1.0+effectivemass(last)*(1.0/effectivemass(0) - 1.0/effectivemass(last-1))/4.0)*copy(0);

	#pragma omp parallel for
	for(long int k = 1; k < last; k++){
		wavefunction[k]=beta(k)*(1.0-effectivemass(k)*(1.0/effectivemass(k+1) - 1.0/effectivemass(k-1))/4.0)*copy(k-1)+(1.0-(2.0)*beta(k))*copy(k)+beta(k)*(1.0+effectivemass(k)*(1.0/effectivemass(k+1) - 1.0/effectivemass(k-1))/4.0)*copy(k+1);
	}

	DiscreteFunction<std::complex<complextype>,xType> as(wavefunction.getGrid());
	DiscreteFunction<std::complex<complextype>,xType> bs(wavefunction.getGrid());
	DiscreteFunction<std::complex<complextype>,xType> cs(wavefunction.getGrid());


	as[0]=0;
	bs[0]=(2.0)*beta(0)+1.0-gammaconst;
	cs[0]=-beta(0)*(1.0+effectivemass(0)*(1.0/effectivemass(1) - 1.0/effectivemass(last))/4.0);

    #pragma omp parallel for
	for(long int i = 1; i < last ; i++){
		as[i]=beta(i)*(-1.0+effectivemass(i)*(1.0/effectivemass(i+1) - 1.0/effectivemass(i-1))/4.0);
		bs[i]=(2.0)*beta(i)+1.0;
		cs[i]=-beta(i)*(1.0+effectivemass(i)*(1.0/effectivemass(i+1) - 1.0/effectivemass(i-1))/4.0);
	}

	as[last]=beta(last)*(-1.0+effectivemass(last)*(1.0/effectivemass(0) - 1.0/effectivemass(last-1))/4.0);
	bs[last]=(2.0)*beta(last)+1.0
			-(beta(0)*beta(last)*(1.0 - (effectivemass(last)*(1.0/effectivemass(0) - 1.0/effectivemass(last-1))/4.0)*(effectivemass(0)*(1.0/effectivemass(1) - 1.0/effectivemass(last))/4.0) ))/gammaconst;
	cs[last]=0;

	//Thomas Algorithm 1

    for (long int i = 1; i <= last; i++)
    {
            std::complex<complextype> m = as(i)/bs(i-1);
            bs[i]-=m*cs(i-1);
            wavefunction[i]-=m*wavefunction(i-1);
    }

    wavefunction[last]/=bs[last];
    for (long int i = (last - 1); i >= 0; i--)
            wavefunction[i]=(wavefunction(i)-cs(i)*wavefunction(i+1))/bs(i);

    // ************

	uvec*=0;
	uvec[0]=gammaconst;
	uvec[last]=-beta(last)*(1.0+effectivemass(last)*(1.0/effectivemass(0) - 1.0/effectivemass(last-1))/4.0);

	as[0]=0;
	bs[0]=(2.0)*beta(0)+1.0-gammaconst;
	cs[0]=-beta(0)*(1.0+effectivemass(0)*(1.0/effectivemass(1) - 1.0/effectivemass(last))/4.0);

    #pragma omp parallel for
	for(long int i = 1; i < last ; i++){
		as[i]=beta(i)*(-1.0+effectivemass(i)*(1.0/effectivemass(i+1) - 1.0/effectivemass(i-1))/4.0);
		bs[i]=(2.0)*beta(i)+1.0;
		cs[i]=-beta(i)*(1.0+effectivemass(i)*(1.0/effectivemass(i+1) - 1.0/effectivemass(i-1))/4.0);
	}

	as[last]=beta(last)*(-1.0+effectivemass(last)*(1.0/effectivemass(0) - 1.0/effectivemass(last-1))/4.0);
	bs[last]=(2.0)*beta(last)+1.0
			-(beta(0)*beta(last)*(1.0 - (effectivemass(last)*(1.0/effectivemass(0) - 1.0/effectivemass(last-1))/4.0)*(effectivemass(0)*(1.0/effectivemass(1) - 1.0/effectivemass(last))/4.0) ))/gammaconst;
	cs[last]=0;

	//Thomas Algorithm 2

    for (long int i = 1; i <= last; i++)
    {
            std::complex<complextype> m2 = as(i)/bs(i-1);
            bs[i]-=m2*cs(i-1);
            uvec[i]-=m2*uvec(i-1);
    }

    uvec[last]/=bs[last];
    for (long int i = (last - 1); i >= 0; i--)
            uvec[i]=(uvec(i)-cs(i)*uvec(i+1))/bs(i);


    std::complex<complextype> afactor=(wavefunction(0)+wavefunction(last)*/*std::conj*/(beta(0)*(-1.0+effectivemass(0)*(1.0/effectivemass(1) - 1.0/effectivemass(last))/4.0)/gammaconst))/(1.0+uvec(0)+uvec(last)*/*std::conj*/(beta(0)*(-1.0+effectivemass(0)*(1.0/effectivemass(1) - 1.0/effectivemass(last))/4.0)/gammaconst));
    for(long int k = 0; k <= last; k++){
    	wavefunction[k]=wavefunction(k)-uvec(k)*afactor;
    }

    return wavefunction;
}

template <typename timetype, typename complextype>
template<typename yType, typename xType>
std::complex<complextype> WaveTimeEvolution<timetype,complextype>::Energy(DiscreteFunction<std::complex<complextype>,xType>& wavefunction, DiscreteFunction<yType,xType>& potential, DiscreteFunction<yType,xType>& effectivemass){
	std::complex<complextype> energysign;

		if(carriertype==Carrier::ELECTRON)
			energysign=std::complex<complextype>(1.0,0);
		else
			energysign=std::complex<complextype>(-1.0,0);

	DiscreteFunction<std::complex<complextype>,xType> normalized(wavefunction);
	normalized.normalize(1.0);

	DiscreteFunction<std::complex<complextype>,xType> copy(normalized);
	DiscreteFunction<std::complex<complextype>,xType> copy2(normalized);

	auto data1 = copy2.getdata();
	auto data2 = potential.getdata();
	#pragma omp parallel for
	for(long int i = 0; i < copy2.getGrid().getSize(); i++){
		*(data1+i)=*(data1+i)*(*(data2+i))*energysign;
	}

	copy.Differential();
	auto dataa = copy.getdata();
	auto datab = effectivemass.getdata();
	#pragma omp parallel for
	for(long int i = 0; i < copy.getGrid().getSize(); i++){
		*(dataa+i)/=(*(datab+i));
	}

	copy.Differential();
	copy*=(-Constant::hbar*Constant::hbar/2.0);

	return energysign*(normalized.Dot(copy)+normalized.Dot(copy2));

}



template <typename timetype, typename complextype>
template<typename yType, typename xType>
DiscreteFunction<std::complex<complextype>,xType> WaveTimeEvolution<timetype,complextype>::potentialpropagator(timetype methodtimestep, DiscreteFunction<yType,xType> potential){
	std::complex<complextype> energysign;

	if(carriertype==Carrier::ELECTRON)
		energysign=std::complex<complextype>(1.0,0);
	else
		energysign=std::complex<complextype>(-1.0,0);

	DiscreteFunction<std::complex<complextype>,xType> result(potential.getGrid());
	#pragma omp parallel for
	for(long int k = 0; k < potential.getGrid().getSize(); k++){
		result[k]=exp(std::complex<complextype>(0.0,-1.0)*energysign*potential(k)*methodtimestep/(2.0*Constant::hbar));
	}

	return result;

}

template <typename timetype, typename complextype>
template<typename yType, typename xType>
DiscreteFunction<std::complex<complextype>,xType> WaveTimeEvolution<timetype,complextype>::fftmulti(yType blochvector,timetype methodtimestep,Grid1D<xType> grid,yType effmass){


	DiscreteFunction<std::complex<complextype>,xType> result(grid);
	long int sizehalf = grid.getSize()/2;
	#pragma omp parallel for
	for(long int k = 0; k < sizehalf; k++){
		result[k]=std::exp( std::complex<complextype>(0.0,-1.0)*(blochvector+2*Constant::pi*k/(grid.getSuperiorLimit()-grid.getInferiorLimit()))*(blochvector+2*Constant::pi*k/(grid.getSuperiorLimit()-grid.getInferiorLimit()))*methodtimestep/(2*Constant::hbar*effmass) );

	}
	for(long int k = sizehalf; k < grid.getSize(); k++){
		result[k]=std::exp( std::complex<complextype>(0.0,-1.0)*(blochvector+2*Constant::pi*(k-2*sizehalf)/(grid.getSuperiorLimit()-grid.getInferiorLimit()))*(blochvector+2*Constant::pi*(k-2*sizehalf)/(grid.getSuperiorLimit()-grid.getInferiorLimit()))*methodtimestep/(2*Constant::hbar*effmass) );
	}

	return result;

}

/*
template <typename timetype, typename complextype>
template<typename yType, typename xType>
DiscreteFunction<std::complex<complextype>,xType> WaveTimeEvolution<timetype,complextype>::blochwave(yType blochvector,Grid1D<xType> grid){

	DiscreteFunction<std::complex<complextype>,xType> result(grid);

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		result[k]=std::exp( std::complex<complextype>(0.0,1.0)*blochvector*(k*grid.getIncrement()+grid.getInferiorLimit()) );
	}

	return result;
}
*/

}

#endif /* WAVETIMEEVOLUTION_CPP_ */
