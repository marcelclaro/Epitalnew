/*
 * WannierFunction.hpp
 *
 *  Created on: Aug 2, 2013
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
#include <fstream>
#include <utility>

#include "Graphics.hpp"
#include "CustomTypes.hpp"
#include "CMaterial.hpp"
#include "Materials.hpp"
#include "Units.hpp"
#include "Alloy.hpp"
#include "Grid.hpp"
#include "DFunction.hpp"
#include "Heterostructure.hpp"
#include "WaveTimeEvolution.hpp"



#ifndef WANNIERFUNCTION_HPP_
#define WANNIERFUNCTION_HPP_


#include "PWFunction.hpp"

namespace epital{


/**
 * \brief This class are used as result of a solve method, it has the wavefunctions and energies.
 */
template<typename complextype, typename xType>
class WannierFunction{
public:
	WannierFunction(xType periodsize, xType vectorstep);
	virtual ~WannierFunction();

	//WannierFunction(string filename);

	/**
	 * Add a level to solution, append in order of energy
	 * @param wavefunction wavefunction to add
	 * @param energy corresponding energy
	 */
	void addFunction(shared_ptr<PlaneWavesFunction<complextype,xType>> wavefunction,std::complex<complextype> energy);

	/**
	 * Get wavefunction
	 * @param number number of level 0=Ground State
	 * @return Function pointer
	 */
	shared_ptr<PlaneWavesFunction<complextype,xType>> getWavefunction(long int number);

	/**
	 * Get energy of level
	 * @param number number of the level 0=Ground state
	 * @return
	 */
	std::complex<complextype> getEnergy(long int number);


	/**
	 * Get number of levels.
	 * @return return number of levels.
	 */
	long int getLevels();

	std::complex<complextype> getValue(long int period, xType position);


	xType getPeriod();

	xType getBlochVectorStep();

	DiscreteFunction<std::complex<complextype>,xType> getDiscrete(long int period, long int periodsback, long int periodsfar, long int points, bool normalize = true);



	//bool saveToFile(string filename);

protected:
	long int levels; ///< Number of levels.
	xType periodsize;
	xType vectorstep;
	string file;

	vector<shared_ptr<PlaneWavesFunction<complextype,xType>>> waves; ///<Pointers to wavefunctions
	vector<std::complex<complextype>> energies; ///< Energies of founded wavefunctions

};

} /* namespace epital */

#include "WannierFunction.cpp"

#endif /* WANNIERFUNCTION_HPP_ */
