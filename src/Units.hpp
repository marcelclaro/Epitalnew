/*
 * Units.hpp
 *
 *  Created on: Jun 20, 2012
 *      Author: marcel
 */

/**
 * @file Units.hpp
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

#ifndef UNITS_HPP_
#define UNITS_HPP_

namespace epital{
/**
 * Conversion from atomic units of mass to Kg (SI)
 * @param au Value in atomic units
 * @return Value in Kg (SI)
 */
inline double mass_to_SI(double au);
/**
 * Conversion from Kg to atomic units of mass
 * @param si Value in Kg
 * @return Value in atomic units
 */
inline double mass_from_SI(double si);

/**
 * Conversion from atomic units of energy (Hartree) to Joule (SI)
 * @param au Value in atomic units (Hartree)
 * @return Value in Joules (SI)
 */
inline double energy_to_SI(double au);
/**
 * Conversion from Joule to atomic units of energy (Hartree)
 * @param si Value in Joules
 * @return Value in Hartree
 */
inline double energy_from_SI(double si);

/**
 * Conversion from atomic units (Hartree) of energy to eV (SI)
 * @param au Value in atomic units
 * @return Value in eV (SI)
 */
inline double energy_to_eV(double au);
/**
 * Conversion from eV (SI) to atomic units of energy (Hartree)
 * @param ev Value in eV (SI)
 * @return Value in atomic units (Hartree)
 */
inline double energy_from_eV(double ev);

/**
 * Conversion from atomic units of length to meters (SI)
 * @param au Value in atomic units
 * @return Value in m (SI)
 */
inline double length_to_SI(double au);
/**
 * Conversion from meters (SI) to atomic units of length
 * @param si Value in m (SI)
 * @return Value in atomic units
 */
inline double length_from_SI(double si);

/**
 * Conversion from atomic units of time to seconds (SI)
 * @param au Value in atomic units
 * @return Value in seconds (SI)
 */
inline double time_to_SI(double au);
/**
 * Conversion from seconds to atomic units of time
 * @param si Value in seconds
 * @return Value in atomic units of time
 */
inline double time_from_SI(double si);

/**
 * Conversion from atomic units of electric potential to Volts (SI)
 * @param au Value in atomic units
 * @return Value in V (SI)
 */
inline double potential_to_SI(double au);
/**
 * Conversion from Volts to atomic units of electric potential
 * @param si Value in volts
 * @return Value in atomic units of electric potential
 */
inline double potential_from_SI(double si);


/**
 * Conversion from atomic units of electric field to V/m (SI)
 * @param au Value in atomic units
 * @return Value in V/m (SI)
 */
inline double electricfield_to_SI(double au);
/**
 * Conversion from V/m to atomic units of electric field
 * @param si Value in V/m
 * @return Value in atomic units of electric field
 */
inline double electricfield_from_SI(double si);


/**
 * Conversion from atomic units of charge to Coulombs (SI)
 * @param au Value in atomic units of charge
 * @return Value in C (SI)
 */
inline double charge_to_SI(double au);
/**
 * Conversion from Coulombs to atomic units of charge
 * @param si Value in C (SI)
 * @return Value in atomic units of charge
 */
inline double charge_from_SI(double si);


/**
 * Conversion from atomic units of velocity to m/s (SI)
 * @param au Value in atomic units
 * @return Value in m/s (SI)
 */
inline double velocity_to_SI(double au);
/**
 * Conversion from Kg to atomic units of velocity
 * @param si Value in m/s
 * @return Value in atomic units of velocity
 */
inline double velocity_from_SI(double si);

/**
 * Conversion from atomic units of force to Newtons (SI)
 * @param au Value in atomic units of force
 * @return Value in atomic units of force
 */
inline double force_to_SI(double au);
/**
 * Conversion from newtons to atomic units of force
 * @param si Value in N (SI)
 * @return Value in atomic units of force
 */
inline double force_from_SI(double si);

/**
 * Conversion from atomic units of current to Amperes (SI)
 * @param au Value in atomic units of current
 * @return Value in A (SI)
 */
inline double current_to_SI(double au);
/**
 * Conversion from Ampere to atomic units of current
 * @param si Value in A (SI)
 * @return Value in atomic units of current
 */
inline double current_from_SI(double si);

/**
 * Conversion from atomic units of permittivity to F/m (SI)
 * @param au Value in atomic units of permittivity
 * @return Value in F/m (SI)
 */
inline double permittivity_to_SI(double au);
/**
 * Conversion from F/m (SI) to atomic units of permittivity
 * @param si Value in F/m
 * @return Value in atomic units of permittivity
 */
inline double permittivity_from_SI(double si);


/**
 * Conversion from atomic units of frequency to Hertz (SI)
 * @param au Value in atomic units
 * @return Value in Hz (SI)
 */
inline double frequency_to_SI(double au);
/**
 * Conversion from Hertz to atomic units of frequency
 * @param si Value in Hz (SI)
 * @return Value in atomic units of frequency
 */
inline double frequency_from_SI(double si);


/**
 * Conversion from atomic units of power to Watts (SI)
 * @param au Value in atomic units of power
 * @return Value in W (SI)
 */
inline double power_to_SI(double au);
/**
 * Conversion from Watts (SI)  to atomic units of power
 * @param si Value in W (SI)
 * @return Value in atomic units of power
 */
inline double power_from_SI(double si);

/**
 * Conversion from atomic units of power area density to W/m² (SI)
 * @param au Value in atomic units of power area density
 * @return Value in W/m² (SI)
 */
inline double powerperarea_to_SI(double au);
/**
 * Conversion from W/m² (SI)  to atomic units of power area density
 * @param si Value in W/m² (SI)
 * @return Value in atomic units of power
 */
inline double powerperarea_from_SI(double si);


/**
 * Conversion from atomic units of density to 1/m³ (SI)
 * @param au Value in atomic units of density
 * @return Value in 1/m³ (SI)
 */
inline double uni_density_to_SI(double au);
/**
 * Conversion from 1/m³ to atomic units of density
 * @param si Value in 1/m³ (SI)
 * @return Value in atomic units of density
 */
inline double uni_density_from_SI(double si);

/**
 * Conversion from atomic units of pressure to Pa (SI)
 * @param au Value in atomic units of pressure
 * @return Value in 1/m³ (SI)
 */
inline double pressure_to_SI(double au);
/**
 * Conversion from Pa to atomic units of Pressure
 * @param si Value in Pa (SI)
 * @return Value in atomic units of pressure
 */
inline double pressure_from_SI(double si);


/**
 * Literal for conversion from meters to atomic units
 * @param value Value in meters
 */
constexpr double operator "" _m(long double value);
/**
 * Literal for conversion from Angstrons to atomic units
 * @param value Value in Angstrons
 */
constexpr double operator "" _Angs(long double value);
/**
 * Literal for conversion from eV (electron volts) to atomic units
 * @param value Value in eV
 */
constexpr double operator "" _eV(long double value);
/**
 * Literal for conversion from meV (mili-electron volts) to atomic units
 * @param value Value in meV
 */
constexpr double operator "" _meV(long double value);
/**
 * Literal for conversion from Volts to atomic units
 * @param value Value in Volts
 */
constexpr double operator "" _Volts(long double value);

/**
 * Literal for conversion from 1/cm³ to atomic units
 * @param value Value in 1/cm³
 */
constexpr double operator "" _under_cubic_cm(long double value);

/**
 * Literal for conversion from Pa to atomic units of pressure
 * @param value Value in Pascal
 */
constexpr double operator "" _Pa(long double value);

/**
 * Literal for conversion from fs to atomic units of time
 * @param value Value in fs
 */
constexpr double operator "" _fs(long double value);

/**
 * Literal for conversion from THz to atomic units of time
 * @param value Value in THz
 */
constexpr double operator "" _THz(long double value);

}

#include "Units.cpp"

#endif /* UNITS_HPP_ */
