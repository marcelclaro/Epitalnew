/*
 * Material.cpp
 *
 *  Created on: Jun 20, 2012
 *      Author: marcel
 */

/**
 * @file CMaterial.cpp
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

#include <cmath>
#include "CMaterial.hpp"
#include "CustomTypes.hpp"

using namespace epital;

Material::~Material() {};

double Material::varshini(Temperature T,double egzero, double alpha, double beta){
	return egzero-(alpha*pow(T,2.0))/(T+beta);
};

double Material::boseeinsteingap(Temperature T,double egzero, double a, double b){
	return egzero-2*a/(exp(b/T)-1);
};

double Material::HHmass_t(double luttinger1,double luttinger2, double luttinger3){
	return 1/(luttinger1-2*luttinger2);
};

double Material::lHmass_t(double luttinger1,double luttinger2, double luttinger3){
	return 1/(luttinger1+2*luttinger2);
};

double Material::HHmass_l(double luttinger1,double luttinger2, double luttinger3){
	return 1/(0.5*(2*luttinger1-luttinger2-3*luttinger3));
}


double Material::lHmass_l(double luttinger1,double luttinger2, double luttinger3){
	return 1/(0.5*(2*luttinger1+luttinger2+3*luttinger3));
}
