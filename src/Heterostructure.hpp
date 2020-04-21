/*
 * Heterostructure<TYPE>.hpp
 *
 *  Created on: Jul 31, 2012
 *      Author: marcel
 */

#ifndef Heterostructure_HPP_
#define Heterostructure_HPP_

#include <memory>
#include <vector>
#include <stdexcept>
#include <complex>
#include "CustomTypes.hpp"
#include "CMaterial.hpp"
#include "Materials.hpp"
#include "Units.hpp"
#include "Alloy.hpp"
#include "Grid.hpp"
#include "DFunction.hpp"

using namespace std;

namespace epital{

/**
 * @class Heterostructure
 * @brief Class for heterostructure
 */
template <typename TYPE>
class Heterostructure {
public:
	/**
	 * @class Epilayer
	 * @brief Define a layer object that forming heterstructure
	 */
	class Epilayer{
	public:
		/**
		 * Constructor
		 * @param material material that compose epilayer
		 * @param layerwidth layer width
		 * @param doping doping of the layer
		 */
		Epilayer(shared_ptr<Material> material, TYPE layerwidth, TYPE doping = 1.0e+14_under_cubic_cm);

		/**
		 * Destructor
		 */
		virtual ~Epilayer();

		/**
		 * Set temperature of epilayer
		 * @param temp Temperature
		 */
		void setTemperature(Temperature temp);

		/**
		 * Get material of layer
		 * @return shared pointer to material
		 */
		shared_ptr<Material> getMaterial();


		/**
		 * Get width of layer
		 * @return width
		 */
		TYPE getWidth() const;

		/**
		 * Get doping from layer
		 * @return doping density
		 */
		TYPE getDoping() const;

		friend class Heterostructure;

	protected:
		shared_ptr<Material> mat;  ///< Material
		TYPE width;  ///< Layer width
		TYPE doping; ///< Doping density

	};

	/**
	 * Constructor
	 * @param activeregion_ Vector with list of epilayers forming active region
	 * @param contacts_ Vector with list of epilayers forming contacts region
	 * @param substrate_ Shared pointer to Material that form substrate.
	 */
	Heterostructure<TYPE>(vector<Epilayer>& activeregion_, vector<Epilayer>& contacts_, shared_ptr<Material> substrate_=make_shared<GaAs>(77),int roughness=0);

	virtual ~Heterostructure<TYPE>();

	/**
	 * Duplicate the heterostructure
	 */
	void duplicate();



	/**
	 * Get total doping ( a 2D density)
	 * @return width
	 */
	TYPE getTotalDoping();


	/**
	 * Say if the structure has or not a contact
	 * @return true if has, false otherwise
	 */
	bool hasContact();


	shared_ptr<Material> getSubstrate();


	/**
	 * Set temperature in the all heterostructure
	 * @param temp Temperature
	 */
	void setTemperature(Temperature temp);

	/**
	* Set bias applied in contacts or heterostructure border
	* @param BIAS applied bias
	*/
	void setBIAS(TYPE BIAS);

	void setBIASladder(TYPE BIAS, long int layersperiod);


	/**
	 * Get total width of the sample
	 * @return width
	 */
	TYPE getTotalWidth();


	/**
	 * Get begin of heterostructure position
	 * @return begin position
	 */
	inline TYPE Begin() const;


	/**
	 * Get end of heterostructure position
	 * @return end position
	 */
	inline TYPE End() const;

	/* TODO try to use a guess based on the last position to accelerate functions */
	/* TODO Paralellize function*/

	/**
	 * Get Conduction band gamma valley energy
	 * @param position position inquired
	 * @return energy value
	 */
	inline TYPE ConductionBand_gamma(TYPE position) const;
	/**
	 * Get Conduction band gamma valley energy function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> ConductionBand_gamma(Grid1D<T>& basegrid) const;
//TODO
	vector<TYPE> ConductionBand_gamma() const;

	/**
	 * Get Conduction band L valley energy
	 * @param position position inquired
	 * @return energy value
	 */
	inline TYPE ConductionBand_L(TYPE position) const;
	/**
	 * Get Conduction band L valley energy function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> ConductionBand_L(Grid1D<T>& basegrid) const;

	vector<TYPE> ConductionBand_L() const;

	/**
	 * Get Conduction band X valley energy
	 * @param position position inquired
	 * @return energy value
	 */
	inline TYPE ConductionBand_X(TYPE position) const;
	/**
	 * Get Conduction band X valley energy function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> ConductionBand_X(Grid1D<T>& basegrid) const;

	vector<TYPE> ConductionBand_X() const;

	/**
	 * Get heavy-hole Valence band energy
	 * @param position position inquired
	 * @return energy value
	 */
	inline TYPE ValenceBand_HH(TYPE position) const;
	/**
	 * Get heavy-hole Valence band energy function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> ValenceBand_HH(Grid1D<T>& basegrid) const;

	vector<TYPE> ValenceBand_HH() const;

	/**
	 * Get ligth-hole Valence band energy
	 * @param position position inquired
	 * @return energy value
	 */
	inline TYPE ValenceBand_lH(TYPE position) const;
	/**
	 * Get ligth-hole Valence band energy function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> ValenceBand_lH(Grid1D<T>& basegrid) const;

	vector<TYPE> ValenceBand_lH() const;

	/**
	 * Get gap at gamma valey
	 * @param position position inquired
	 * @return gap value
	 */
	inline TYPE gap_gamma(TYPE position) const;
	/**
	 * Get gap at gamma valey function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> gap_gammaDiscrete(Grid1D<T>& basegrid) const;

	/**
	 * Get gap at L valey
	 * @param position position inquired
	 * @return gap value
	 */
	inline TYPE gap_L(TYPE position) const;
	/**
	 * Get gap at L valey function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> gap_LDiscrete(Grid1D<T>& basegrid) const;


	/**
	 * Get gap at X valey
	 * @param position position inquired
	 * @return gap value
	 */
	inline TYPE gap_X(TYPE position) const;
	/**
	 * Get gap at X valey function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> gap_XDiscrete(Grid1D<T>& basegrid) const;


	/**
	 * Effective mass of conduction band electron at gamma valley
	 * @param position position inquired
	 * @return mass value
	 */
	inline TYPE effectiveMass_gamma(TYPE position) const;
	/**
	 * Get Effective mass of conduction band electron at gamma valley function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> effectiveMass_gamma(Grid1D<T>& basegrid) const;

	vector<TYPE> effectiveMass_gamma() const;

	/**
	 * Effective mass of conduction band electron at L valley (001) direction
	 * @param position position inquired
	 * @return mass value
	 */
	inline TYPE effectiveMass_Lt(TYPE position) const;
	/**
	 * Get Effective mass of conduction band electron at L valley (001) direction function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> effectiveMass_Lt(Grid1D<T>& basegrid) const;

	vector<TYPE> effectiveMass_Lt() const;

	/**
	* Effective mass of conduction band electron at L valley (110) direction
	* @param position position inquired
	* @return mass value
	*/
	inline TYPE effectiveMass_Ll(TYPE position) const;
	/**
	 * Get Effective mass of conduction band electron at L valley (110) direction function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> effectiveMass_Ll(Grid1D<T>& basegrid) const;


	/**
	 * Effective mass of conduction band electron at X valley (001) direction
	 * @param position position inquired
	 * @return mass value
	 */
	inline TYPE effectiveMass_Xt(TYPE position) const;
	/**
	 * Get Effective mass of conduction band electron at X valley (001) direction function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> effectiveMass_Xt(Grid1D<T>& basegrid) const;

	vector<TYPE> effectiveMass_Xt() const;

	/**
	 * Effective mass of conduction band electron at X valley (110) direction
	 * @param position position inquired
	 * @return mass value
	 */
	inline TYPE effectiveMass_Xl(TYPE position) const;
	/**
	 * Get Effective mass of conduction band electron at X valley (110) direction function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> effectiveMass_Xl(Grid1D<T>& basegrid) const;


	/**
	 * Density of states effective mass of conduction band electron at L valley direction
	 * @param position position inquired
	 * @return mass value
	 */
	inline TYPE effectiveMass_LDOS(TYPE position) const;
	/**
	 * Density of states effective mass of conduction band electron at L valley direction function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> effectiveMass_LDOS(Grid1D<T>& basegrid) const;

	/**
	 * Density of states effective mass of conduction band electron at X valley direction
	 * @param position position inquired
	 * @return mass value
	 */
	inline TYPE effectiveMass_XDOS(TYPE position) const;
	/**
	 * Density of states effective mass of conduction band electron at X valley direction function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> effectiveMass_XDOS(Grid1D<T>& basegrid) const;

	/**
	 * Effective mass of heavy-hole valence band electron at (001) direction
	 * @param position position inquired
	 * @return mass value
	 */
	inline TYPE effectiveMass_HHt(TYPE position) const;
	/**
	 * Effective mass of heavy-hole valence band electron at (001) direction function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> effectiveMass_HHt(Grid1D<T>& basegrid) const;

	vector<TYPE> effectiveMass_HHt() const;

	/**
	 * Effective mass of heavy-hole valence band electron at (110) direction
	 * @param position position inquired
	 * @return mass value
	 */
	inline TYPE effectiveMass_HHl(TYPE position) const;
	/**
	 * Effective mass of heavy-hole valence band electron at (110) direction function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> effectiveMass_HHl(Grid1D<T>& basegrid) const;

	/**
	 * Effective mass of light-hole valence band electron at (001) direction
	 * @param position position inquired
	 * @return mass value
	 */
	inline TYPE effectiveMass_lHt(TYPE position) const;
	/**
	 * Effective mass of light-hole valence band electron at (001) direction function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> effectiveMass_lHt(Grid1D<T>& basegrid) const;

	vector<TYPE> effectiveMass_lHt() const;

	/**
	 * Effective mass of light-hole valence band electron at (110) direction
	 * @param position position inquired
	 * @return mass value
	 */
	inline TYPE effectiveMass_lHl(TYPE position) const;
	/**
	 * Effective mass of light-hole valence band electron at (110) direction function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> effectiveMass_lHl(Grid1D<T>& basegrid) const;

	/**
	 * First Luttinger parameter
	 * @param position position inquired
	 * @return parameter value
	 */
	inline TYPE luttingerParam_1(TYPE position) const;
	/**
	 * First Luttinger parameter
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> luttingerParam_1(Grid1D<T>& basegrid) const;

	/**
	 * Second Luttinger parameter
	 * @param position position inquired
	 * @return parameter value
	 */
	inline TYPE luttingerParam_2(TYPE position) const;
	/**
	 * Second Luttinger parameter function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> luttingerParam_2(Grid1D<T>& basegrid) const;

	/**
	 * Third Luttinger parameter
	 * @param position position inquired
	 * @return parameter value
	 */
	inline TYPE luttingerParam_3(TYPE position) const;
	/**
	 * Get Third Luttinger parameter function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> luttingerParam_3(Grid1D<T>& basegrid) const;

	/**
	 * Lattice parameter
	 * @param position position inquired
	 * @return parameter value
	 */
	inline TYPE latticeparam(TYPE position) const;
	/**
	 * Get Lattice parameter function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> latticeparam(Grid1D<T>& basegrid) const;

	/**
	 * Doping density
	 * @param position position inquired
	 * @return parameter value
	 */
	inline TYPE doping(TYPE position) const;
	/**
	 * Get Doping density function
	 * @param basegrid The base grid to the function
	 * @return shared pointer to function
	 */
	template <typename T>
	DiscreteFunction<TYPE,T> doping(Grid1D<T>& basegrid) const;

	vector<TYPE> getWidths() const;

	vector<pair<TYPE,TYPE>> getLimitsLists() const;

protected:
	vector<Epilayer>& activeregion;  ///<Active region layers vector
	vector<Epilayer>& contacts; ///<Contact region layers vector
	shared_ptr<Material> substrate; ///< Substrate material shared pointer

	vector<shared_ptr<Material> > matlst; ///<Ordered material list
	vector<pair<TYPE,TYPE> > limitslst; ///<Orderred layer limits
	vector<TYPE> widthlst; ///< Ordered layer width list
	vector<TYPE> dopinglst; ///< Doping density list
	vector<TYPE> dEc; ///< Strain induced conduction band shift
	vector<TYPE> dEv_lH; ///< Strain induced light-hole valence band shift
	vector<TYPE> dEv_HH; ///< Strain induced heavy-hole band shift

	TYPE bias;
	TYPE ladderbias;
	vector<TYPE> linearelectricfield;
	vector<TYPE> layerbeginenergy;
	vector<TYPE> electricfieldladder;

	bool hascontact;


};


}

#include "Heterostructure.cpp"


#endif /* Heterostructure_HPP_ */
