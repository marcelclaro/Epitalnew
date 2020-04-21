/*
 * Units.cpp
 *
 *  Created on: Jun 20, 2012
 *      Author: marcel
 */

/**
 * @file Units.cpp
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

#include "Units.hpp"
#include "Constants.hpp"

#ifndef UNITS_CPP_
#define UNITS_CPP_

namespace epital{

using namespace Constant;

inline double mass_to_SI(double au){
	return au*Constant::emass_SI;
}

inline double mass_from_SI(double si){
	return si/Constant::emass_SI;
}

inline double energy_to_SI(double au){
	return au*4.35974394e-18;
}

inline double energy_from_SI(double si){
	return si/4.35974394e-18;
}

inline double energy_to_eV(double au){
	return au*27.21138386;
}

inline double energy_from_eV(double ev){
	return ev/27.21138386;
}

inline double length_to_SI(double au){
	return au*5.2917721092e-11;
}

inline double length_from_SI(double si){
	return si/5.2917721092e-11;
}

inline double time_to_SI(double au){
	return  au*2.418884326502e-17;
}

inline double time_from_SI(double si){
	return si/2.418884326502e-17;
}

inline double potential_to_SI(double au){
	return au*27.21138505;
}

inline double potential_from_SI(double si){
	return si/27.21138505;
}

inline double electricfield_to_SI(double au){
	return au*5.14220651e+11;
}

inline double electricfield_from_SI(double si){
	return si/5.14220651e+11;
}

inline double charge_to_SI(double au){
	return au*Constant::e_SI;
}

inline double charge_from_SI(double si){
	return si/Constant::e_SI;
}

inline double velocity_to_SI(double au){
	return au*2.18769126379e+06;
}

inline double velocity_from_SI(double si){
	return si/2.18769126379e+06;
}

inline double force_to_SI(double au){
	return au*8.23872278e-08;
}

inline double force_from_SI(double si){
	return si/8.23872278e-08;
}

inline double current_to_SI(double au){
	return au*6.62361795e-03;
}

inline double current_from_SI(double si){
	return si/6.62361795e-03;
}

inline double permittivity_to_SI(double au){
	return au*1.112650056e-10;
}

inline double permittivity_from_SI(double si){
	return si/1.112650056e-10;
}

inline double frequency_to_SI(double au){
	return au/2.418884326502e-17;
}

inline double frequency_from_SI(double si){
	return si*2.418884326502e-17;
}

inline double power_to_SI(double au){
	return au*1.80237814e-1;
}

inline double power_from_SI(double si){
	return si/1.80237814e-1;
}

inline double powerperarea_to_SI(double au){
	return au*6.43640931e19;
}

inline double powerperarea_from_SI(double si){
	return si/6.43640931e19;
}


inline double uni_density_to_SI(double au){
	return au/148.184711486e-33;
}

inline double uni_density_from_SI(double si){
	return si*148.184711486e-33;
}

inline double pressure_to_SI(double au){
	return au*2.9421912e+13;
}

inline double pressure_from_SI(double si){
	return si/2.9421912e+13;
}


constexpr double operator "" _m(long double value){
	return value/5.2917721092e-11;
}

constexpr double operator "" _Angs(long double value){
	return value/5.2917721092e-01;
}

constexpr double operator "" _eV(long double value){
	return value/27.21138386;
}

constexpr double operator "" _meV(long double value){
	return (value*0.001)/27.21138386;
}

constexpr double operator "" _Volts(long double value){
	return value/27.21138505;
}

constexpr double operator "" _under_cubic_cm(long double value){
	return value*148.184711486e-27;
}

constexpr double operator "" _Pa(long double value){
	return value/2.9421912e+13;
}

constexpr double operator "" _fs(long double value){
	return value/2.418884326502e-2;
}

constexpr double operator "" _THz(long double value){
	return value*1.0e12*2.418884326502e-17;
}

}

#endif
