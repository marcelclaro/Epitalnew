/*
 * Alloy.cpp
 *
 *  Created on: Jun 21, 2012
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

#include "Alloy.hpp"
#include "Units.hpp"

using namespace epital;

Alloy::Alloy(const GaAs& gaas,const AlAs& alas): principal(gaas), secundary(alas) {
	bowgap_gamma1=-0.127_eV;
	bowgap_gamma2=1.310_eV;
	bowgap_X=0.055_eV;
	bowgap_L=0;
	boweffectivemass_e=0;
	boweffectivemass_lh=0;
	boweffectivemass_hh=0;
}

Alloy::Alloy(const GaAs& gaas,const  InAs& inas): principal(gaas), secundary(inas){
	bowgap_gamma1=-0.477_eV;
	bowgap_gamma2=0;
	bowgap_X=1.4_eV;
	bowgap_L=0.33_eV;
	boweffectivemass_e=0.0091;
	boweffectivemass_lh=0.0202;
	boweffectivemass_hh=-0.145;
}
Alloy::Alloy(const AlAs& alas,const  InAs& inas): principal(alas), secundary(inas){
	bowgap_gamma1=-0.70_eV;
	bowgap_gamma2=0;
	bowgap_X=0;
	bowgap_L=0;
	boweffectivemass_e=0.049;
	boweffectivemass_lh=0;
	boweffectivemass_hh=0;
}

Alloy::~Alloy() {

}
