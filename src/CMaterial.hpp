/*
 * Material.hpp
 *
 *  Created on: Jun 20, 2012
 *      Author: marcel
 */

/**
 * @file CMaterial.hpp
 * @author Marcel S. Claro <marcelclaro@gmail.com>
 *
 * @section License
 *
 * Copyright (c) 2012, Marcel Santos Claro
 *All rights reserved.
 *
 *Redistribution and use in source and binary forms, with or without
 *modification, are permitted provided that the following conditions are met:
 *1. Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *2. Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *3. All advertising materials mentioning features or use of this software
 *   must display the following acknowledgement:
 *   This product includes software developed by the Marcel Santos Claro.
 *4. Neither the name of the Marcel Santos Claro nor the
 *   names of its contributors may be used to endorse or promote products
 *   derived from this software without specific prior written permission.
 *
 *THIS SOFTWARE IS PROVIDED BY Marcel Santos Claro ''AS IS'' AND ANY
 *EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 *WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *DISCLAIMED. IN NO EVENT SHALL Marcel Santos Claro BE LIABLE FOR ANY
 *DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 *ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 */


#include "CustomTypes.hpp"
#include <string>

#ifndef CMATERIAL_HPP_
#define CMATERIAL_HPP_

namespace epital{
/**
 * @class Material
 * @brief This class define a abstract class for materials(alloys) used in this library
 *
 * This class define a abstract class for definition of materials(alloys) used in this library. The material parameter like effective mass, lattice constant,
 * etc must be defined in constructor of derived class. Another parameter must be implemented in this class first.
 *
 */
class Material {
public:
	/**
	 * Not applicable
	 */
	Material() = default;
	/**
    * Not applicable
	 */
	virtual ~Material();


	inline std::string getmaterialname() const;

	//Material properties declarations

	/**
	 *	Effective mass of conduction band at gamma valley in electron masses
	 * @return Effective mass
	 */
	inline double effectiveMass_gamma() const;
	/**
	 *	Effective mass of conduction band at L valley in (001) direction (in electron masses)
	 * @return Effective mass
	 */
	inline double effectiveMass_Lt() const;
	/**
	 *	Effective mass of conduction band at L valley in (110) directions (in electron masses)
	 * @return Effective mass
	 */
	inline double effectiveMass_Ll() const;
	/**
	 *	Effective mass of conduction band at X valley in (001) direction (in electron masses)
	 * @return Effective mass
	 */
	inline double effectiveMass_Xt() const;
	/**
	 *	Effective mass of conduction band at X valley in (110) directions (in electron masses)
	 * @return Effective mass
	 */
	inline double effectiveMass_Xl() const;
	/**
	 *	Density of states effective mass of conduction band at L valley (in electron masses)
	 * @return DOS Effective mass
	 */
	inline double effectiveMass_LDOS() const;
	/**
	 *	Density of states effective mass of conduction band at X valley (in electron masses)
	 * @return DOS Effective mass
	 */
	inline double effectiveMass_XDOS() const;
	/**
	 *	Effective mass of heavy holes in (001) direction (in electron masses)
	 * @return Effective mass
	 */
	inline double effectiveMass_HHt() const;
	/**
	 *	Effective mass of heavy holes in (110) directions (in electron masses)
	 * @return Effective mass
	 */
	inline double effectiveMass_HHl() const;
	/**
	 *	Effective mass of light holes in (001) direction (in electron masses)
	 * @return Effective mass
	 */
	inline double effectiveMass_lHt() const;
	/**
	 *	Effective mass of light holes in (110) directions (in electron masses)
	 * @return Effective mass
	 */
	inline double effectiveMass_lHl() const;
	/**
	 *	First Luttinger parameter value
	 * @return Luttinger Parameter value
	 */
	inline double luttingerParam_1() const;
	/**
	 *	Second Luttinger parameter value
	 * @return Luttinger Parameter value
	 */
	inline double luttingerParam_2() const;
	/**
	 *	Third Luttinger parameter value
	 * @return Luttinger Parameter value
	 */
	inline double luttingerParam_3() const;
	/**
	 *  Gap in gamma valley (in atomic units)
	 * @return gap value
	 */
	inline double gap_gamma() const;
	/**
	 *  Gap in L valley (in atomic units)
	 * @return gap value
	 */
	inline double gap_L() const;
	/**
	 *  Gap in X valley (in atomic units)
	 * @return gap value
	 */
	inline double gap_X() const;
	/**
	 * Lattice parameter of material (in atomic units)
	 * @return lattice parameter value
	 */
	inline double latticeparam() const;

	/**
	 * Get the fraction of minority (intermediate) element of a Ternary (Quaternary) alloy.
	 * @return
	 */
	inline double getComposition() const;

	/**
	* Get the fraction of minority element of a Quaternary alloy.
    * @return
	*/
	inline double getCompositionQuate() const;

	/**
	 * Get valence band energy relative to GaAs valence band
	 * @return Valence band energy (GaAs=0)
	 */
	inline double ValenceBandEnergy() const;

	/**
	 * Get strain parameter ac
	 * @return ac
	 */
	inline double getStrain_ac() const;

	/**
	 * Get strain parameter av
	 * @return av
	 */
	inline double getStrain_av() const;

	/**
	 * Get strain parameter b
	 * @return b
	 */
	inline double getStrain_b() const;

	/**
	 * Get elastic constant C11
	 * @return C11
	 */
	inline double getC11() const;

	/**
	 * Get elastic constant C12
	 * @return C12
	 */
	inline double getC12() const;

	/**
	 * Get elastic constant C44
	 * @return C44
	 */
	inline double getC44() const;

	/**
	 * Change material temperature. This is a pure abstract method and needs definition in derived classes. In derived classes this method
	 * must recalculated parameters that change with temperature like gap.
	 * @param temp New temperature
	 */
	virtual void setTemperature(Temperature temp) = 0;

	/**
	 * Get LO phonon frequency
	 * @return frequency
	 */
	inline double getwLO() const;


	/**
	 * Get TO phonon frequency
	 * @return frequency
	 */
	inline double getwTO() const;


	/**
	 * Get static dielectric constant
	 * @return relative dielectric constant
	 */
	inline double getDielectric_static() const;

	/**
	 * Get optical dielectric constant
	 * @return relative dielectric constant
	 */
	inline double getDielectric_optical() const;


	/**
	 * Get current material temperature
	 * @return Current temperature
	 */
	Temperature getTemperature() const;

protected:
	/**
	 * Varshini equation for calculation of band gap in variable temperature
	 * @param T Temperature
	 * @param egzero gap at 0K
	 * @param alpha alpha parameter
	 * @param beta beta parameter
	 * @return gap value
	 */
	double varshini(Temperature T,double egzero, double alpha, double beta);

	/**
	 * Bose-Einstein equation for calculation of band gap in variable temperature
	 * @param T Temperature
	 * @param egzero gap at 0K
	 * @param alpha alpha parameter
	 * @param beta beta parameter
	 * @return gap value
	 */
	double boseeinsteingap(Temperature T, double egzero, double a, double b);

	 /**
	 * Calculate effective mass of heavy hole in (001) direction from Luttinger parameter
	 * @param luttinger1 First Luttinger parameter
	 * @param luttinger2 Second Luttinger parameter
	 * @param luttinger3 Third Luttinger parameter
	 * @return Effective mass
	 */
	double HHmass_t(double luttinger1,double luttinger2, double luttinger3);
	/**
	 * Calculate effective mass of light hole in (001) direction from Luttinger parameter
	 * @param luttinger1 First Luttinger parameter
	 * @param luttinger2 Second Luttinger parameter
	 * @param luttinger3 Third Luttinger parameter
	 * @return Effective mass
	 */
	double lHmass_t(double luttinger1,double luttinger2, double luttinger3);
	/**
	 * Calculate effective mass of heavy hole in (110) directions from Luttinger parameter
	 * @param luttinger1 First Luttinger parameter
	 * @param luttinger2 Second Luttinger parameter
	 * @param luttinger3 Third Luttinger parameter
	 * @return Effective mass
	 */
	double HHmass_l(double luttinger1,double luttinger2, double luttinger3);
	/**
	 * Calculate effective mass of light hole in (110) directions from Luttinger parameter
	 * @param luttinger1 First Luttinger parameter
	 * @param luttinger2 Second Luttinger parameter
	 * @param luttinger3 Third Luttinger parameter
	 * @return Effective mass
	 */
	double lHmass_l(double luttinger1,double luttinger2, double luttinger3);



	Temperature temp; ///< Current temperature
	float compositionA; ///< Compositional fraction of secondary element in a ternary(quarternary) alloy
	float compositionB; ///< Compositional fraction of another-secondary element in a quarternary alloy

	double gapgamma; ///<Gap at gamma valley
	double gapL; ///<Gap at L valley
	double gapX; ///<Gap at X valley

	double lparam; ///<Lattice parameter

	double emass_gama; ///<Effective mass in conduction band at gamma valley
	double emass_Lt; ///<Effective mass in conduction band at L valley in (001) direction
	double emass_Ll; ///<Effective mass in conduction band at L valley in (110) directions
	double emass_Xt; ///<Effective mass in conduction band at X valley in (001) direction
	double emass_Xl; ///<Effective mass in conduction band at L valley in (110) directions
	double emass_Ldos; ///<Density of states effective mass in conduction band at L valley
	double emass_Xdos; ///<Density of states effective mass in conduction band at X valley
	double emassHH_l; ///<Effective mass of heavy holes in (110) directions
	double emassHH_t; ///<Effective mass of heavy holes in (001) direction
	double emasslH_l; ///<Effective mass of light holes in (110) directions
	double emasslH_t; ///<Effective mass of light holes in (001) direction

	double lutt1; ///<First Luttinger parameter
	double lutt2; ///<Second Luttinger parameter
	double lutt3; ///<Third Luttinger parameter

	double alpha_gamma; ///<Varshini parameter alpha for gamma valley
	double beta_gamma; ///<Varshini parameter beta for gamma valley
	double alpha_L; ///<Varshini parameter alpha for L valley
	double beta_L; ///<Varshini parameter beta for L valley
	double alpha_X; ///<Varshini parameter alpha for X valley
	double beta_X; ///<Varshini parameter beta for X valley
	double gapgammatzero; ///<Gap at 0 K for gamma valley
	double gapLtzero; ///<Gap at 0 K for L valley
	double gapXtzero; ///<Gap at 0 K for X valley
	double gapa;  ///bose-einstein gap parameter a
	double gapb;  ///bose-einstein gap parameter teta_b

	double valencebandpos; ///< Energy of valence band relative to GaAs.

	double ac;  ///< strain ac parameter.
	double av;  ///< strain av parameter.
	double b;   ///< strain b parameter.

	double c11; ///< Elastic constant c11
	double c12; ///< Elastic constant c12
	double c44; ///< Elastic constant c44

	double wLO; ///< LO phonon frequency
	double wTO; ///< TO phonon frequency

	double dielectric0; ///< Static dielectric constant
	double dielectricinf; ///<Optical dielectric constant

	std::string materialname;

};

inline double Material::effectiveMass_gamma() const{
	return emass_gama;
};

inline double Material::effectiveMass_Lt() const {
	return emass_Lt;
};

inline double Material::effectiveMass_Ll() const {
	return emass_Ll;
};

inline double Material::effectiveMass_Xt() const {
	return emass_Xt;
};

inline double Material::effectiveMass_Xl() const {
	return emass_Xl;
};

inline double Material::effectiveMass_LDOS() const {
	return emass_Ldos;
};

inline double Material::effectiveMass_XDOS() const {
	return emass_Xdos;
};

inline double Material::gap_gamma() const {
	return gapgamma;
};

inline double Material::gap_L() const {
	return gapL;
};

inline double Material::gap_X() const {
	return gapX;
};

inline double Material::latticeparam() const {
	return lparam;
};

inline Temperature Material::getTemperature() const {
	return temp;
};

inline double Material::effectiveMass_HHt() const {
	return emassHH_t;
};

inline double Material::effectiveMass_HHl() const {
	return emassHH_l;
};

inline double Material::effectiveMass_lHt() const {
	return emasslH_t;
};

inline double Material::effectiveMass_lHl() const {
	return emasslH_l;
};

inline double Material::luttingerParam_1() const {
	return lutt1;
};

inline double Material::luttingerParam_2() const {
	return lutt2;
};

inline double Material::luttingerParam_3() const {
	return lutt3;
};

inline double Material::getComposition() const{
	return compositionA;
}

inline double Material::getCompositionQuate() const{
	return compositionB;
}

inline double Material::ValenceBandEnergy() const{
	return valencebandpos;
}

inline double Material::getStrain_ac() const{
	return ac;
}

inline double Material::getStrain_av() const{
	return av;
}

inline double Material::getStrain_b() const{
	return b;
}

inline double Material::getC11() const{
	return c11;
}

inline double Material::getC12() const{
	return c12;
}

inline double Material::getC44() const{
	return c44;
}


inline double Material::getwLO() const{
	return wLO;
}


inline double Material::getwTO() const{
	return wTO;
}



inline double Material::getDielectric_static() const{
	return dielectric0;
}


inline double Material::getDielectric_optical() const{
	return dielectricinf;
}

inline std::string Material::getmaterialname() const{
	return materialname;
}


}
#endif /* CMATERIAL_HPP_ */
