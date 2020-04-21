/*
 * Alloy.hpp
 *
 *  Created on: Jun 21, 2012
 *      Author: marcel
 */

/**
 * @file Alloy.hpp
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

#ifndef ALLOY_HPP_
#define ALLOY_HPP_



#include "Materials.hpp"

namespace epital{
/**
 * @class Alloy
 * @brief This class define a material made by 2 materials with variable composition \f$(A_{1-x}B_xC)\f$
 *
 * Used when a material is a alloy that composition varies continuously in space or with other parameter
 *
 */
class Alloy {
public:
	/**
	 * Constructor for AlGaAs alloy
	 * @param gaas A GaAs material
	 * @param alas A AlAs material
	 */
	Alloy(const GaAs& gaas,const AlAs& alas);
	/**
	 * Constructor for InGaAs alloy
	 * @param gaas A GaAs material
	 * @param inas A InAs material
	 */
	Alloy(const GaAs& gaas,const  InAs& inas);
	/**
	 * Constructor for InAlAs alloy
	 * @param alas A AlAs material
	 * @param inas A InAs material
	 */
	Alloy(const AlAs& alas,const  InAs& inas);

	/**
	 * Standard destructor
	 */
	virtual ~Alloy();


	/**
	 * Effective mass of conduction band at gamma valley in electron masses
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return  Effective mass
	 */
	inline double effectiveMass_gamma(double composition) const;

	/**
	 * Effective mass of conduction band at L valley in (001) direction (in electron masses)
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return Effective mass
	 */
	inline double effectiveMass_Lt(double composition) const;

	/**
	 * Effective mass of conduction band at L valley in (110) directions (in electron masses)
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return Effective mass
	 */
	inline double effectiveMass_Ll(double composition) const;

	/**
	 * Effective mass of conduction band at X valley in (001) direction (in electron masses)
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return Effective mass
	 */
	inline double effectiveMass_Xt(double composition) const;

	/**
	 * Effective mass of conduction band at X valley in (110) directions (in electron masses)
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return Effective mass
	 */
	inline double effectiveMass_Xl(double composition) const;

	/**
	 * Density of states effective mass of conduction band at L valley (in electron masses)
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return DOS Effective mass
	 */
	inline double effectiveMass_LDOS(double composition) const;

	/**
	 * 	Density of states effective mass of conduction band at X valley (in electron masses)
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return DOS Effective mass
	 */
	inline double effectiveMass_XDOS(double composition) const;

	/**
	 * Effective mass of heavy holes in (001) direction (in electron masses)
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return Effective mass
	 */
	inline double effectiveMass_HHt(double composition) const;

	/**
	 * Effective mass of heavy holes in (110) directions (in electron masses)
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return Effective mass
	 */
	inline double effectiveMass_HHl(double composition) const;

	/**
	 * Effective mass of light holes in (001) direction (in electron masses)
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return Effective mass
	 */
	inline double effectiveMass_lHt(double composition) const;

	/**
	 * Effective mass of light holes in (110) directions (in electron masses)
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return Effective mass
	 */
	inline double effectiveMass_lHl(double composition) const;

	/**
	 * First Luttinger parameter value
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return Luttinger Parameter value
	 */
	inline double luttingerParam_1(double composition) const;

	/**
	 * Second Luttinger parameter value
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return Luttinger Parameter value
	 */
	inline double luttingerParam_2(double composition) const;

	/**
	 * Third Luttinger parameter value
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return Luttinger Parameter value
	 */
	inline double luttingerParam_3(double composition) const;

	/**
	 * Gap in gamma valley (in atomic units)
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return gap value
	 */
	inline double gap_gamma(double composition) const;

	/**
	 * Gap in L valley (in atomic units)
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return gap value
	 */
	inline double gap_L(double composition) const;

	/**
	 * Gap in X valley (in atomic units)
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return gap value
	 */
	inline double gap_X(double composition) const;

	/**
	 * Lattice parameter of material (in atomic units)
	 * @param composition Composition of alloy in % of secondary material (x value)
	 * @return lattice parameter value
	 */
	inline double latticeparam(double composition) const;
protected:
	//Bowing parameters
	double bowgap_gamma1; ///< First bowing parameter for Gap in gamma valley.
	double bowgap_gamma2; ///< Second bowing parameter for Gap in gamma valley.
	double bowgap_X; ///<bowing parameter for Gap in X valley.
	double bowgap_L; ///<bowing parameter for Gap in L valley.
	double boweffectivemass_e; ///<bowing parameter for electron effective mass
	double boweffectivemass_lh; ///<bowing parameter for light-hole effective mass
	double boweffectivemass_hh; ///<bowing parameter for heavy-hole effective mass

private:
	const Material& principal; ///< Primary material from alloy
	const Material& secundary; ///< Secondary material from alloy
};

inline double Alloy::effectiveMass_gamma(double composition) const{
	return (1-composition)*principal.effectiveMass_gamma()+composition*secundary.effectiveMass_gamma()-composition*(1-composition)*boweffectivemass_e;
};

inline double Alloy::effectiveMass_Lt(double composition) const {
	return (1-composition)*principal.effectiveMass_Lt()+composition*secundary.effectiveMass_Lt();
};

inline double Alloy::effectiveMass_Ll(double composition) const {
	return (1-composition)*principal.effectiveMass_Ll()+composition*secundary.effectiveMass_Ll();
};

inline double Alloy::effectiveMass_Xt(double composition) const {
	return (1-composition)*principal.effectiveMass_Xt()+composition*secundary.effectiveMass_Xt();
};

inline double Alloy::effectiveMass_Xl(double composition) const {
	return (1-composition)*principal.effectiveMass_Xl()+composition*secundary.effectiveMass_Xl();
};

inline double Alloy::effectiveMass_LDOS(double composition) const {
	return (1-composition)*principal.effectiveMass_LDOS()+composition*secundary.effectiveMass_LDOS();
};

inline double Alloy::effectiveMass_XDOS(double composition) const {
	return (1-composition)*principal.effectiveMass_XDOS()+composition*secundary.effectiveMass_XDOS();
};

inline double Alloy::gap_gamma(double composition) const {
	return (1-composition)*principal.gap_gamma()+composition*secundary.gap_gamma()-composition*(1-composition)*(bowgap_gamma1+bowgap_gamma2*composition);
};

inline double Alloy::gap_L(double composition) const {
	return (1-composition)*principal.gap_L()+composition*secundary.gap_L()-composition*(1-composition)*bowgap_L;
};

inline double Alloy::gap_X(double composition) const {
	return (1-composition)*principal.gap_X()+composition*secundary.gap_X()-composition*(1-composition)*bowgap_X;
};

inline double Alloy::latticeparam(double composition) const {
	return (1-composition)*principal.latticeparam()+composition*secundary.latticeparam();
};


inline double Alloy::effectiveMass_HHt(double composition) const {
	return (1-composition)*principal.effectiveMass_HHt()+composition*secundary.effectiveMass_HHt()-composition*(1-composition)*boweffectivemass_hh;
};

inline double Alloy::effectiveMass_HHl(double composition) const {
	return (1-composition)*principal.effectiveMass_HHl()+composition*secundary.effectiveMass_HHl();
};

inline double Alloy::effectiveMass_lHt(double composition) const {
	return (1-composition)*principal.effectiveMass_lHt()+composition*secundary.effectiveMass_lHt()-composition*(1-composition)*boweffectivemass_lh;
};

inline double Alloy::effectiveMass_lHl(double composition) const {
	return (1-composition)*principal.effectiveMass_lHl()+composition*secundary.effectiveMass_lHl();
};

inline double Alloy::luttingerParam_1(double composition) const {
	return (1-composition)*principal.luttingerParam_1()+composition*secundary.luttingerParam_1();
};

inline double Alloy::luttingerParam_2(double composition) const {
	return (1-composition)*principal.luttingerParam_2()+composition*secundary.luttingerParam_2();
};

inline double Alloy::luttingerParam_3(double composition) const {
	return (1-composition)*principal.luttingerParam_3()+composition*secundary.luttingerParam_3();
};

}

#endif /* ALLOY_HPP_ */
