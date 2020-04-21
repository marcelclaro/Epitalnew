/*
 * Constants.hpp
 *
 *  Created on: Jun 20, 2012
 *      Author: marcel
 */

/**
 * @file Constants.hpp
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

#ifndef CONSTANTS_HPP_
#define CONSTANTS_HPP_

namespace epital{
/**
 * Physical constants by CODATA 2010
 * @brief Physical constants
 */
namespace Constant{

const double pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899;

//In atomic system of units
const double emass = 1e+0; ///<Electron Mass constant in Hartree Atomic Units
const double e = 1e+0; ///<Electron charge constant in Hartree Atomic Units
const double h = 6.283185307179586476; ///<Planck constant in Hartree Atomic Units
const double hbar = 1e+0; ///<Planck constant divided by 2*pi in Hartree Atomic Units
const double kb = 3.1668114e-6; ///<Boltzmann constant in Hartree Atomic Units
const double c = 137.035999074; ///<Light velocity constant in Hartree Atomic Units
const double eps0 = 0.0795775/*8.854187817620e-12/1.112650056e-10*/; ///<Vacuum Permittivity

// In SI system of unit
const double emass_SI = 9.10938291e-31; ///<Electron Mass constant in International System of Units
const double e_SI = 1.602176565e-19; ///<Electron charge constant in International System of Units
const double h_SI = 6.62606957e-34; ///<Planck constant in International System of Units
const double hbar_SI = 1.054571726e-34; ///<\f$(\hbar)\f$ : Planck constant divided by \f$(2\pi)\f$  in International System of Units
const double kb_SI = 1.3806504e-23; ///<Boltzmann constant in International System of Units
const double c_SI = 299792458; ///<Light velocity constant in International System of Units
const double eps0_SI = 8.854187817620e-12;

//In eV(SI) system of units
const double h_eV = 4.135667516e-15; ///<Planck constant in International System of Units using electron volts (eV)
const double hbar_eV = 6.58211928e-16; ///<\f$(\hbar)\f$ : Planck constant divided by \f$(2\pi)\f$  in International System of Units using electron volts (eV)
const double kb_eV= 8.617343e-5; ///<Boltzmann constant in International System of Units using electron volts (eV)
}

}

#endif /* CONSTANTS_HPP_ */
