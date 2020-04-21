/*
 * CarrierStatistics.hpp
 *
 *  Created on: Dec 5, 2012
 *      Author: marcel
 */

#include <Solution.hpp>
#include <SolutionPW.hpp>
#include <CMaterial.hpp>
#include <vector>

#ifndef CARRIERSTATISTICS_HPP_
#define CARRIERSTATISTICS_HPP_



namespace epital{

/**
 * \brief This class do the calculation of Fermi levels and carrier statistics
 */
template<typename complextype>
class CarrierStatistics {
public:
	/**
	 * Empty constructor.
	 */
	CarrierStatistics();
	/**
	 * Empty destructor.
	 */
	virtual ~CarrierStatistics();

	/**
	 * Approximated calculation of fermi level in quantum well
	 * @param solution Solution with energy levels
	 * @param carriersperperiod Total 2d density of carrier in periods
	 * @param material Material of structure
	 * @param T Temperature of electrons (lattice)
	 * @param banda Type of carriers
	 * @return
	 */
	template<typename xType>
	static complextype FermiEnergyWells(Solution<complextype,xType> solution, complextype carriersperperiod, shared_ptr<Material> material, Temperature T, Band banda);

	template<typename xType>
	static complextype FermiEnergyWells(SolutionPW<complextype,xType> solution, complextype carriersperperiod, shared_ptr<Material> material, Temperature T, Band banda);



	template<typename xType>
	static complextype QuasiFermiEnergyWells(std::complex<complextype> subbandenergy, complextype subbandcarriersperperiod, shared_ptr<Material> material, Temperature T, Band banda);


	/**
	 * Calculate the 2d density of carriers in a subband
	 * @param T Temperature
	 * @param subbandminima Subband minima
	 * @param fermienergy Fermi energy
	 * @param material Material
	 * @param banda Band of subband
	 * @param energytop Max energy in integration
	 * @return 2d density in subband
	 */
	static complextype CarriersinSubband(Temperature T, complextype subbandminima, complextype fermienergy,shared_ptr<Material> material,Band banda, complextype energytop);


private:

	/**
	 * \brief Object used for root finding in FermiEnergyWells
	 */
	class CarrierNumber{
	public:
		complextype effectivemass;
		complextype supenergy;
		Temperature temp;
		vector<complextype> energies;
		complextype carriersperperiod;
		complextype carriersignal;

		/**
		 * Initialization constructor
		 * @param carriersperperiod
		 * @param temp
		 * @param effectivemass
		 * @param energies
		 * @param supenergy
		 */
		CarrierNumber(complextype carriersperperiod,Temperature temp,complextype effectivemass, vector<complextype> energies, complextype supenergy,int carriersignal);

		/**
		 * Function to be nulled
		 * @param energy Fermi energy
		 * @return Difference between calculated carriers based in Fermi level and real carrier
		 */
		complextype tosolvefermilevel(complextype energy);
	};
};

}
#include "CarrierStatistics.cpp"

#endif /* CARRIERSTATISTICS_HPP_ */
