/*
 * WaveTimeEvolution2D.hpp
 *
 *  Created on: Nov 7, 2014
 *      Author: marcel
 */

#include <thread>
#include <iostream>
#include <limits>
#include <string>
#include <cstdio>
#include <vector>
#include <complex>
#include <armadillo>

#include "CustomTypes.hpp"
#include "Units.hpp"
#include "Grid.hpp"
#include "DFunction2D.hpp"
//#include "Heterostructure2D.hpp"


#ifndef WAVETIMEEVOLUTION2D_HPP_
#define WAVETIMEEVOLUTION2D_HPP_


namespace epital {

/**
 *\brief  This class implements the evolution in time (real or complex) of wavefunction on a grid
 */
template <typename timetype, typename complextype>
class WaveTimeEvolution2D {
public:

	/**
	 * Type of carrier of wavefunction
	 */
	enum Carrier{
		HOLE=-1,  //!< HOLE
		ELECTRON=1//!< ELECTRON
	} carriertype;


	WaveTimeEvolution2D(WaveTimeEvolution2D::Carrier carrier = Carrier::ELECTRON);
	virtual ~WaveTimeEvolution2D();

	template<typename T, typename gridtype>
	DiscreteFunction2D<std::complex<complextype>,gridtype>& groundlevel_periodicFFT(timetype time, long int timesteps, DiscreteFunction2D<std::complex<complextype>,gridtype>& initial, DiscreteFunction2D<T,gridtype>& potential, T effectivemass);

	template<typename T, typename gridtype>
	DiscreteFunction2D<std::complex<complextype>,gridtype>& newlevel_periodicFFT(timetype time, long int timesteps, DiscreteFunction2D<std::complex<complextype>,gridtype>& initial, DiscreteFunction2D<T,gridtype>& potential,std::vector<std::shared_ptr<DiscreteFunction2D<std::complex<complextype>,gridtype>>> otherslevels, T effectivemass);

	template<typename T, typename gridtype>
	DiscreteFunction2D<std::complex<complextype>,gridtype>& groundlevel_finite(timetype time, long int timesteps, DiscreteFunction2D<std::complex<complextype>,gridtype>& initial, DiscreteFunction2D<T,gridtype>& potential, DiscreteFunction2D<T,gridtype>& effectivemass);

	template<typename T, typename gridtype>
	DiscreteFunction2D<std::complex<complextype>,gridtype>& groundlevel_finite_valence(timetype time, long int timesteps, DiscreteFunction2D<std::complex<complextype>,gridtype>& initial, DiscreteFunction2D<T,gridtype>& potential, DiscreteFunction2D<T,gridtype>& effectivemasst,DiscreteFunction2D<T,gridtype>& effectivemassl);

	template<typename T, typename gridtype>
	DiscreteFunction2D<std::complex<complextype>,gridtype>& newlevel_finite(timetype time, long int timesteps, DiscreteFunction2D<std::complex<complextype>,gridtype>& initial, DiscreteFunction2D<T,gridtype>& potential,std::vector<std::shared_ptr<DiscreteFunction2D<std::complex<complextype>,gridtype>>> otherslevels, DiscreteFunction2D<T,gridtype>& effectivemass);

/*
	template<typename T, typename gridtype>
	DiscreteFunction2D<std::complex<complextype>,gridtype>& kinect_propag_evolution_finite(timetype methodtimestep, DiscreteFunction2D<std::complex<complextype>,gridtype>& wavefunction, DiscreteFunction2D<T,gridtype>& effectivemass);
*/

	template<typename T, typename gridtype>
	std::complex<complextype> Energy(DiscreteFunction2D<std::complex<complextype>,gridtype>& wavefunction, DiscreteFunction2D<T,gridtype>& potential, DiscreteFunction2D<T,gridtype>& effectivemass);

	template<typename T, typename gridtype>
	std::complex<complextype> Energy(DiscreteFunction2D<std::complex<complextype>,gridtype>& wavefunction, DiscreteFunction2D<T,gridtype>& potential, DiscreteFunction2D<T,gridtype>& effectivemasst,DiscreteFunction2D<T,gridtype>& effectivemassl);

private:
	/**
	 * Potential factor of time propagator
	 * @param methodtimestep time step of the method
	 * @param potential potental on the grid
	 * @return the potential factor
	 */
	template<typename T, typename gridtype>
	DiscreteFunction2D<std::complex<complextype>,gridtype> potentialpropagator(timetype methodtimestep, DiscreteFunction2D<T,gridtype> potential);

	template<typename T, typename gridtype>
	DiscreteFunction2D<std::complex<complextype>,gridtype> fftmulti(timetype methodtimestep,Grid2D<gridtype> grid,T effmass);

};

} /* namespace epital */

#include "WaveTimeEvolution2D.cpp"

#endif /* WAVETIMEEVOLUTION2D_HPP_ */
