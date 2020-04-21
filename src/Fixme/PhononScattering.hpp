/*
 * PhononScattering.hpp
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
#include <list>
#include <complex>
#include <chrono>
#include <ratio>

#include <CustomTypes.hpp>
#include <CMaterial.hpp>
#include "Materials.hpp"
#include "Units.hpp"
#include "Alloy.hpp"
#include "Grid.hpp"
#include "DFunction.hpp"
#include "Graphics.hpp"
#include "Heterostructure.hpp"
#include "WaveTimeEvolution.hpp"
#include "Solution.hpp"

#ifndef PHONONSCATTERING_HPP_
#define PHONONSCATTERING_HPP_

namespace epital{

/**\brief This class do the calculations for phonon scattering in quantum wells
 *
 */
template<typename complextype>
class PhononScattering {
public:


	/**Empty constructor
	 *
	 */
	PhononScattering();

	/**
	 * Empty destructor
	 */
	virtual ~PhononScattering();

	/**
	 * The Form factor of phonon interaction
	 * @param wavefunctionfinal Final wavefunction
	 * @param wavefunctioninitial Initial wavefunction
	 * @param phononwavenumber z component of phonon wavenumber
	 * @return
	 */
	template<typename xType>
	static complextype FormFactorSquare(shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial, complextype phononwavenumber);

	/**
	 * Transition rate of phonon emission process in quantum wells
	 * @param T Temperature of lattice
	 * @param wavefunctionfinal Final wave function
	 * @param wavefunctioninitial Initial wavefunction
	 * @param energyfinallevel Energy of final state (subband minima)
	 * @param energyinitiallevel Energy of initial state (subband minima or total energy[ + k//])
	 * @param kparalleli Parallel modulos of wavenumber
	 * @param material Material of lattice
	 * @param band Carrier type
	 * @return Transistion rate
	 */
	template<typename xType>
	static complextype EmissionRateLO(Temperature T, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,
			complextype energyfinallevel, complextype energyinitiallevel,complextype kparalleli, shared_ptr<Material> material, Band band);

	/**
	 * Transition rate of phonon absorption process in quantum wells
	 * @param T Temperature of lattice
	 * @param wavefunctionfinal Final wave function
	 * @param wavefunctioninitial Initial wavefunction
	 * @param energyfinallevel Energy of final state (subband minima)
	 * @param energyinitiallevel Energy of initial state (subband minima or total energy[ + k//])
	 * @param kparalleli Parallel modulos of carrier wavenumber
	 * @param material Material of lattice
	 * @param band Carrier type
	 * @return Transistion rate
	 */
	template<typename xType>
	static complextype AbsRateLO(Temperature T, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,
			complextype energyfinallevel, complextype energyinitiallevel,complextype kparalleli, shared_ptr<Material> material, Band band);


	template<typename xType>
	static complextype MeanEmissionRateLO(Temperature T, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,
			complextype energyfinallevel, complextype energyinitiallevel,complextype fermifinal,complextype fermiinitial , shared_ptr<Material> material, Band band, long int limitpar =  10);


	template<typename xType>
	static complextype MeanAbsRateLO(Temperature T, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,
			complextype energyfinallevel, complextype energyinitiallevel,complextype fermifinal,complextype fermiinitial , shared_ptr<Material> material, Band band, long int limitpar =  10);



	template<typename xType>
	static complextype IntersubbandMeanEmissionRateLO(Temperature T, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,
			complextype energyfinallevel, complextype energyinitiallevel,complextype quasifermifinal,complextype quasifermiinitial , shared_ptr<Material> material, Band band, long int limitpar =  10);


	template<typename xType>
	static complextype IntersubbandMeanAbsRateLO(Temperature T, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,
			complextype energyfinallevel, complextype energyinitiallevel,complextype quasifermifinal,complextype quasifermiinitial , shared_ptr<Material> material, Band band, long int limitpar =  10);



private:

	/**
	 * Integrand for transition rate in quantum wells (x1e+10)
	 */
	template<typename xType>
	class Integrand{
	public:
		shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal; ///< Final wavefunction pointer
		shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial; ///< Initial wavefunction pointer
		complextype kparalleli;   ///<parallel modulos of carrier wavenumber
		complextype deltaenergy; ///< Ef-Ei+-hwLO
		complextype effectivemass;  ///< Carrier effective mass

		/**
		 * Constructor
		 * @param wavefunctionfinal
		 * @param wavefunctioninitial
		 * @param kparalleli
		 * @param deltaenergy
		 * @param effectivemass
		 */
		Integrand(shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal,shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,complextype kparalleli,complextype deltaenergy,complextype effectivemass);

		/**
		 * Integrand function (x1e+10)
		 * @param phononwavenumber z component of phonon wavenumber
		 * @return Function value * 1e+10
		 */
		complextype thefunction(complextype phononwavenumber);
	};

	/**
	 * Integrand for mean transition rate in quantum wells
	 */
	template<typename xType>
	class IntegrandforMean{
	public:
		Temperature T;
		shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal;
		shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial;
		complextype energyfinallevel;
		complextype fermifinal;
		complextype fermiinitial;
		shared_ptr<Material> material;
		Band band;
		complextype wLO;

		IntegrandforMean(Temperature T, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctionfinal, shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunctioninitial,
				complextype energyfinallevel,complextype fermifinal, complextype fermiinitial, shared_ptr<Material> material, Band band);


		complextype functionemi(complextype initialenergy);
		complextype functionabs(complextype initialenergy);
		complextype functioncarriersub(complextype initialenergy);
	};

};

}

#include "PhononScattering.cpp"


#endif /* PHONONSCATTERING_HPP_ */
