/*
 * WaveTimeEvolution.hpp
 *
 *  Created on: Aug 21, 2012
 *      Author: marcel
 */

#ifndef WAVETIMEEVOLUTION_HPP_
#define WAVETIMEEVOLUTION_HPP_

#include <thread>
#include <iostream>
#include <limits>
#include <string>
#include <cstdio>
#include <vector>
#include <complex>

#include "CustomTypes.hpp"
#include "Units.hpp"
#include "Grid.hpp"
#include "DFunction.hpp"
#include "Heterostructure.hpp"

namespace epital{

/**
 *\brief  This class implements the evolution in time (real or complex) of wavefunction on a grid
 */
template <typename timetype, typename complextype>
class WaveTimeEvolution {
public:
	/**
	 * Type of carrier of wavefunction
	 */
	enum Carrier{
		HOLE=-1,  //!< HOLE
		ELECTRON=1//!< ELECTRON
	} carriertype;

	/**
	 * contructor
	 * @param timesteps number of time steps
	 */
	WaveTimeEvolution(WaveTimeEvolution::Carrier carrier = Carrier::ELECTRON);
	virtual ~WaveTimeEvolution();


	/**
	 * Evolute in time the wavefunction on a finite grid.
	 * @param time total time to evolute.
	 * @param timesteps time step of the method.
	 * @param initial initial function.
	 * @param potential potential on the grid.
	 * @param effectivemass effective mass on the grid.
	 * @return reference to evolved wavefunction.
	 */
	template<typename yType, typename xType>
	DiscreteFunction<std::complex<complextype>,xType>& evolute_finite(timetype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential, DiscreteFunction<yType,xType>& effectivemass);


	/**
	 * Find ground state from a initial wavefunction.
	 * @param time time to function stabilization.
	 * @param timestep time step of the method.
	 * @param initial initial wavefunction (normally a gaussian wavefunction).
	 * @param potential potential on the grid.
	 * @param effectivemass effective mass on the grid.
	 * @return reference to ground state wavefunction (the same that initial wavefunction)
	 */
	template<typename yType, typename xType>
	DiscreteFunction<std::complex<complextype>,xType>& groundlevel_finite(timetype time, long int timestep, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential, DiscreteFunction<yType,xType>& effectivemass);


	/**
	 * Find upperstate levels from a initial wavefunction
	 * @param time time to function stabilization.
	 * @param timestep time step of the method.
	 * @param initial initial wavefunction (normally a gaussian wavefunction).
	 * @param potential potential on the grid.
	 * @param otherslevels vector of pointer to ortogonal wavefunctions already taken.
	 * @param effectivemass effective mass on the grid
	 * @return reference to founded wavefunction (the same that initial wavefunction)
	 */
	template<typename yType, typename xType>
	DiscreteFunction<std::complex<complextype>,xType>& newlevel_finite(timetype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential,vector<shared_ptr<DiscreteFunction<std::complex<complextype>,xType>>> otherslevels, DiscreteFunction<yType,xType>& effectivemass);


	/**
	 * Correlation function of the wavefunction in evolution and initial function
	 * @param time time to function stabilization.
	 * @param timestep time step of the method.
	 * @param initial initial wavefunction (normally a gaussian wavefunction).
	 * @param potential potential on the grid.
	 * @param otherslevels vector of pointer to ortogonal wavefunctions already taken.
	 * @param effectivemass effective mass on the grid
	 * @return the correlation funtion
	 */
	template<typename yType, typename xType>
	DiscreteFunction<std::complex<complextype>,complextype> correlationfunction_finite(complextype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential, DiscreteFunction<yType,xType>& effectivemass);


	/**
	 * Evolute in time the wavefunction on a periodic potential (width of the grid).
	 * @param time total time to evolute.
	 * @param timesteps time step of the method.
	 * @param initial initial function.
	 * @param potential potential on the grid.
	 * @param effectivemass effective mass on the grid.
	 * @return reference to evolved wavefunction.
	 */
	template<typename yType, typename xType>
	DiscreteFunction<std::complex<complextype>,xType>& evolute_periodic(timetype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential, DiscreteFunction<yType,xType>& effectivemass);


	/**
	 * Find ground state from a initial wavefunction on a periodic potential (width of the grid).
	 * @param time time to function stabilization.
	 * @param timestep time step of the method.
	 * @param initial initial wavefunction (normally a gaussian wavefunction).
	 * @param potential potential on the grid.
	 * @param otherslevels vector of pointer to ortogonal wavefunctions already taken.
	 * @param effectivemass effective mass on the grid
	 * @return reference to ground state wavefunction (the same that initial wavefunction)
	 */
	template<typename yType, typename xType>
	DiscreteFunction<std::complex<complextype>,xType>& groundlevel_periodic(timetype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential, DiscreteFunction<yType,xType>& effectivemass);


	/**
	 * Find ground state from a initial wavefunction on a periodic potential (width of the grid).
	 * @param time time to function stabilization.
	 * @param timestep time step of the method.
	 * @param initial initial wavefunction (normally a gaussian wavefunction).
	 * @param potential potential on the grid.
	 * @param otherslevels vector of pointer to ortogonal wavefunctions already taken.
	 * @param effectivemass effective mass on the grid
	 * @return reference to ground state wavefunction (the same that initial wavefunction)
	 */
	template<typename yType, typename xType>
	DiscreteFunction<std::complex<complextype>,xType>& groundlevel_periodicFFT(timetype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential, DiscreteFunction<yType,xType>& effectivemass);


	/**
	 * Find upperstate levels from a initial wavefunction on a periodic potential (width of the grid).
	 * @param time time to function stabilization.
	 * @param timestep time step of the method.
	 * @param initial initial wavefunction (normally a gaussian wavefunction).
	 * @param potential potential on the grid.
	 * @param otherslevels vector of pointer to ortogonal wavefunctions already taken.
	 * @param effectivemass effective mass on the grid
	 * @return reference to founded wavefunction (the same that initial wavefunction)
	 */
	template<typename yType, typename xType>
	DiscreteFunction<std::complex<complextype>,xType>& newlevel_periodic(timetype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential,vector<shared_ptr<DiscreteFunction<std::complex<complextype>,xType>>> otherslevels, DiscreteFunction<yType,xType>& effectivemass);

	/**
	 * Find upperstate levels from a initial wavefunction on a periodic potential (width of the grid).
	 * @param time time to function stabilization.
	 * @param timestep time step of the method.
	 * @param initial initial wavefunction (normally a gaussian wavefunction).
	 * @param potential potential on the grid.
	 * @param otherslevels vector of pointer to ortogonal wavefunctions already taken.
	 * @param effectivemass effective mass on the grid
	 * @return reference to founded wavefunction (the same that initial wavefunction)
	 */
	template<typename yType, typename xType>
	DiscreteFunction<std::complex<complextype>,xType>& newlevel_periodicFFT(timetype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential,vector<shared_ptr<DiscreteFunction<std::complex<complextype>,xType>>> otherslevels, DiscreteFunction<yType,xType>& effectivemass);




	/**
	 * Correlation function of the wavefunction in evolution and initial function  on a periodic potential (width of the grid).
	 * @param time time to function stabilization.
	 * @param timestep time step of the method.
	 * @param initial initial wavefunction (normally a gaussian wavefunction).
	 * @param potential potential on the grid.
	 * @param otherslevels vector of pointer to ortogonal wavefunctions already taken.
	 * @param effectivemass effective mass on the grid
	 * @return the correlation funtion
	 */
	template<typename yType, typename xType>
	DiscreteFunction<std::complex<complextype>,complextype> correlationfunction_periodic(complextype time, long int timesteps, DiscreteFunction<std::complex<complextype>,xType>& initial, DiscreteFunction<yType,xType>& potential, DiscreteFunction<yType,xType>& effectivemass);


	/**
	 * Total energy of the wave function.
	 * @param wavefunction the wavefunction reference
	 * @param potential potential on the grid
	 * @param effectivemass effective mass on the grid
	 * @return Energy value ( valence band top is 0).
	 */
	template<typename yType, typename xType>
	std::complex<complextype> Energy(DiscreteFunction<std::complex<complextype>,xType>& wavefunction, DiscreteFunction<yType,xType>& potential, DiscreteFunction<yType,xType>& effectivemass);


private:

	/**
	 * Potential factor of time propagator
	 * @param methodtimestep time step of the method
	 * @param potential potental on the grid
	 * @return the potential factor
	 */
	template<typename yType, typename xType>
	DiscreteFunction<std::complex<complextype>,xType> potentialpropagator(timetype methodtimestep, DiscreteFunction<yType,xType> potential);


	template<typename yType, typename xType>
	DiscreteFunction<std::complex<complextype>,xType> fftmulti(yType blochvector,timetype methodtimestep,Grid1D<xType> grid,yType effmass);


	/*Blochwave used in tests*/
	//template<typename yType, typename xType>
	//DiscreteFunction<std::complex<complextype>,xType> blochwave(yType blochvector,Grid1D<xType> grid);


	/**
	 * Apply kinect factor of time propagator application ( for finite potentials).
	 * @param methodtimestep time step of the method
	 * @param wavefunction wavefunction to apply kin. factor.
	 * @param effectivemass effective mass on the grid
	 * @return reference to wavefunction after kinect factor of time propagator application
	 */
	template<typename yType, typename xType>
	DiscreteFunction<std::complex<complextype>,xType>& kinect_propag_evolution_finite(timetype methodtimestep, DiscreteFunction<std::complex<complextype>,xType>& wavefunction, DiscreteFunction<yType,xType>& effectivemass);

	/**
	 * Wavefunction after kinect factor of time propagator application ( for periodic potentials).
	 * @param methodtimestep time step of the method
	 * @param wavefunction wavefunction to apply kin. factor.
	 * @param effectivemass effective mass on the grid
	 * @return reference to wavefunction after kinect factor of time propagator application
	 */
	template<typename yType, typename xType>
	DiscreteFunction<std::complex<complextype>,xType>& kinect_propag_evolution_periodic(timetype methodtimestep, DiscreteFunction<std::complex<complextype>,xType>& wavefunction, DiscreteFunction<yType,xType>& effectivemass);

};

}

#include "WaveTimeEvolution.cpp"


#endif /* WAVETIMEEVOLUTION_HPP_ */
