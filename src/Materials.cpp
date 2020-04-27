/*
 * Materials.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: marcel
 */

/**
 * @file Materials.cpp
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
#include "Materials.hpp"

using namespace epital;



// it need holes effective masses
GaAs::GaAs(Temperature temper){
	compositionA=0;
	compositionB=0;
	emass_gama=0.067;
	emass_Lt=0.0754;
	emass_Ll=1.9;
	emass_Xt=0.23;
	emass_Xl=1.3;
	emass_Ldos=0.56;
	emass_Xdos=0.85;

	lutt1=6.98;
	lutt2=2.06;
	lutt3=2.93;

	emassHH_t=HHmass_t(lutt1,lutt2,lutt3);
	emassHH_l=HHmass_l(lutt1,lutt2,lutt3);
	emasslH_t=lHmass_t(lutt1,lutt2,lutt3);
	emasslH_l=lHmass_l(lutt1,lutt2,lutt3);
	//emassHH_t=0.395;
	//emassHH_l=0.934;
	//emasslH_t=0.089;
	//emasslH_l=0.072;

	alpha_gamma=0.5405_meV;
	beta_gamma=204;
	alpha_L=0.605_meV;
	beta_L=204;
	alpha_X=0.460_meV;
	beta_X=204;
	gapgammatzero=1.519_eV;
	gapLtzero=1.815_eV;
	gapXtzero=1.981_eV;
	gapa=0;
	gapb=0;
	valencebandpos=0;

	temp=temper;
	gapgamma=varshini(temp,gapgammatzero,alpha_gamma,beta_gamma);
	gapL=varshini(temp,gapLtzero,alpha_L,beta_L);
	gapX=varshini(temp,gapXtzero,alpha_X,beta_X);
	lparam=5.65325_Angs+(temp-300)*3.88e-05_Angs;

	ac=-7.17_eV;
	av=-1.16_eV;
	b=-2.0_eV;

	c11=122.1e+9_Pa;
	c12=56.6e+9_Pa;
	c44=60.0e+9_Pa;

	wLO=54.9934_THz;
	wTO=50.6113_THz;

	dielectric0 = 12.90;
	dielectricinf = 10.92;

	materialname = "GaAs";
}
void GaAs::setTemperature(Temperature temper){
	temp=temper;
	gapgamma=varshini(temp,gapgammatzero,alpha_gamma,beta_gamma);
	gapL=varshini(temp,gapLtzero,alpha_L,beta_L);
	gapX=varshini(temp,gapXtzero,alpha_X,beta_X);
	lparam=5.65325_Angs+(temp-300)*3.88e-05_Angs;
}


AlAs::AlAs(Temperature temper){
	compositionA=0;
	compositionB=0;
	emass_gama=0.150;
	emass_Lt=0.15;
	emass_Ll=1.32;
	emass_Xt=0.22;
	emass_Xl=0.97;
	emass_Ldos=0.78;
	emass_Xdos=0.76;

	lutt1=3.76;
	lutt2=0.82;
	lutt3=1.42;

	emassHH_t=HHmass_t(lutt1,lutt2,lutt3);
	emassHH_l=HHmass_l(lutt1,lutt2,lutt3);
	emasslH_t=lHmass_t(lutt1,lutt2,lutt3);
	emasslH_l=lHmass_l(lutt1,lutt2,lutt3);
	//emassHH_t=0.51;
	//emassHH_l=1.09;
	//emasslH_t=0.18;
	//emasslH_l=0.15;

	alpha_gamma=0.885_meV;
	beta_gamma=530;
	alpha_L=0.605_meV;
	beta_L=204;
	alpha_X=0.70_meV;
	beta_X=530;
	gapgammatzero=3.099_eV;
	gapLtzero=2.46_eV;
	gapXtzero=2.24_eV;
	valencebandpos=-0.53_eV;
	gapa=0;
	gapb=0;

	ac=-5.64_eV;
	av=-2.47_eV;
	b=-2.3_eV;

	c11=125.0e+9_Pa;
	c12=53.4e+9_Pa;
	c44=64.2e+9_Pa;

	wLO=76.0749_THz;
	wTO=67.9818_THz;


	dielectric0 = 10.06;
	dielectricinf = 8.16;

	temp=temper;
	gapgamma=varshini(temp,gapgammatzero,alpha_gamma,beta_gamma);
	gapL=varshini(temp,gapLtzero,alpha_L,beta_L);
	gapX=varshini(temp,gapXtzero,alpha_X,beta_X);
	lparam=5.6611_Angs+(temp-300)*2.90e-05_Angs;

	materialname = "AlAs";
}
void AlAs::setTemperature(Temperature temper){
	temp=temper;
	gapgamma=varshini(temp,gapgammatzero,alpha_gamma,beta_gamma);
	gapL=varshini(temp,gapLtzero,alpha_L,beta_L);
	gapX=varshini(temp,gapXtzero,alpha_X,beta_X);
	lparam=5.6611_Angs+(temp-300)*2.90e-05_Angs;
}


InAs::InAs(Temperature temper){
	compositionA=0;
	compositionB=0;
	emass_gama=0.026;
	emass_Lt=0.05;
	emass_Ll=0.64;
	emass_Xt=0.16;
	emass_Xl=1.13;
	emass_Ldos=0.29;
	emass_Xdos=0.64;

	lutt1=20;
	lutt2=8.5;
	lutt3=9.2;

	emassHH_t=HHmass_t(lutt1,lutt2,lutt3);
	emassHH_l=HHmass_l(lutt1,lutt2,lutt3);
	emasslH_t=lHmass_t(lutt1,lutt2,lutt3);
	emasslH_l=lHmass_l(lutt1,lutt2,lutt3);
	//emassHH_t=0.35;
	//emassHH_l=0.85; //Needs review
	//emasslH_t=0.025;
	//emasslH_l=0.025; //copied emasslH_t


	alpha_gamma=0.276_meV;
	beta_gamma=93;
	alpha_L=0.276_meV;
	beta_L=93;
	alpha_X=0.276_meV;
	beta_X=93;
	gapgammatzero=0.417_eV;
	gapLtzero=1.133_eV;
	gapXtzero=1.433_eV;
	valencebandpos=0.21_eV;
	gapa=0;
	gapb=0;


	ac=-5.08_eV;
	av=-1.00_eV;
	b=-1.8_eV;

	c11=83.29e+9_Pa;
	c12=45.26e+9_Pa;
	c44=39.59e+9_Pa;

	wLO=45.2028_THz;
	wTO=41.0576_THz;

	dielectric0 = 15.15;
	dielectricinf = 15.25;

	temp=temper;
	gapgamma=varshini(temp,gapgammatzero,alpha_gamma,beta_gamma);
	gapL=varshini(temp,gapLtzero,alpha_L,beta_L);
	gapX=varshini(temp,gapXtzero,alpha_X,beta_X);
	lparam=6.0583_Angs+(temp-300)*2.74e-05_Angs;

	materialname = "InAs";
}
void InAs::setTemperature(Temperature temper){
	temp=temper;
	gapgamma=varshini(temp,gapgammatzero,alpha_gamma,beta_gamma);
	gapL=varshini(temp,gapLtzero,alpha_L,beta_L);
	gapX=varshini(temp,gapXtzero,alpha_X,beta_X);
	lparam=6.0583_Angs+(temp-300)*2.74e-05_Angs;
}


AlGaAs::AlGaAs(Temperature temper, double composition){
	compositionA=composition;
	compositionB=0;
	emass_gama=(1-compositionA)*0.067+compositionA*0.150;
	emass_Lt=(1-compositionA)*0.0754+compositionA*0.150;
	emass_Ll=(1-compositionA)*1.9+compositionA*1.32;
	emass_Xt=(1-compositionA)*0.23+compositionA*0.22;
	emass_Xl=(1-compositionA)*1.3+compositionA*0.97;
	emass_Ldos=(1-compositionA)*0.56+compositionA*0.78;
	emass_Xdos=(1-compositionA)*0.85+compositionA*0.76;

	double lutt1GaAs=6.98;
	double lutt2GaAs=2.06;
	double lutt3GaAs=2.93;

	double lutt1AlAs=3.76;
	double lutt2AlAs=0.82;
	double lutt3AlAs=1.42;

	lutt1=(1-compositionA)*lutt1GaAs+compositionA*lutt1AlAs;
	lutt2=(1-compositionA)*lutt2GaAs+compositionA*lutt2AlAs;
	lutt3=(1-compositionA)*lutt3GaAs+compositionA*lutt3AlAs;

	emassHH_t=HHmass_t(lutt1,lutt2,lutt3);
	emassHH_l=HHmass_l(lutt1,lutt2,lutt3);
	emasslH_t=lHmass_t(lutt1,lutt2,lutt3);
	emasslH_l=lHmass_l(lutt1,lutt2,lutt3);

	temp=temper;
	double gapgammaGaAs=varshini(temp,1.519_eV,0.5405_meV,204);
	double gapgammaAlAs=varshini(temp,3.099_eV,0.885_meV,530);
	double	gapLGaAs=varshini(temp,1.815_eV,0.605_meV,204);
	double	gapLAlAs=varshini(temp,2.46_eV,0.605_meV,204);
	double gapXGaAs=varshini(temp,1.981_eV,0.460_meV,204);
	double gapXAlAs=varshini(temp,2.24_eV,0.70_meV,530);
	valencebandpos=compositionA*(-0.53_eV);
	gapa=0;
	gapb=0;

	ac=(1-compositionA)*(-7.17_eV)+compositionA*(-5.64_eV);
	av=(1-compositionA)*(-1.16_eV)+compositionA*(-2.47_eV);
	b=(1-compositionA)*(-2.0_eV)+compositionA*(-2.3_eV);

	c11=(1-compositionA)*(122.1e+9_Pa)+compositionA*(125.0e+9_Pa);
	c12=(1-compositionA)*(56.6e+9_Pa)+compositionA*(53.4e+9_Pa);
	c44=(1-compositionA)*(60.0e+9_Pa)+compositionA*(64.2e+9_Pa);

	wLO=(1-compositionA)*(54.9934_THz)+compositionA*(76.0749_THz);
	wTO=(1-compositionA)*(50.6113_THz)+compositionA*(67.9818_THz);

	dielectric0 = (1-compositionA)*(12.90)+compositionA*(10.06);
	dielectricinf = (1-compositionA)*(10.92)+compositionA*(8.16);

	gapgamma=(1-compositionA)*gapgammaGaAs+compositionA*gapgammaAlAs-compositionA*(1-compositionA)*(-0.127_eV+1.310_eV*compositionA);
	gapL=(1-compositionA)*gapLGaAs+compositionA*gapLAlAs;
	gapX=(1-compositionA)*gapXGaAs+compositionA*gapXAlAs-compositionA*(1-compositionA)*0.055_eV;


	lparam=5.65325_Angs+0.0078_Angs*compositionA+(temp-300)*(3.88e-05_Angs*(1-compositionA)+2.90e-05_Angs*compositionA);

	materialname = "AlGaAs";


}
void AlGaAs::setTemperature(Temperature temper){
	temp=temper;
	double gapgammaGaAs=varshini(temp,1.519_eV,0.5405_meV,204);
	double gapgammaAlAs=varshini(temp,3.099_eV,0.885_meV,530);
	double	gapLGaAs=varshini(temp,1.815_eV,0.605_meV,204);
	double	gapLAlAs=varshini(temp,2.46_eV,0.605_meV,204);
	double gapXGaAs=varshini(temp,1.981_eV,0.460_meV,204);
	double gapXAlAs=varshini(temp,2.24_eV,0.70_meV,530);

	gapgamma=(1-compositionA)*gapgammaGaAs+compositionA*gapgammaAlAs-compositionA*(1-compositionA)*(-0.127_eV+1.310_eV*compositionA);
	gapL=(1-compositionA)*gapLGaAs+compositionA*gapLAlAs;
	gapX=(1-compositionA)*gapXGaAs+compositionA*gapXAlAs-compositionA*(1-compositionA)*0.055_eV;


	lparam=5.65325_Angs+0.0078_Angs*compositionA+(temp-300)*(3.88e-05_Angs*(1-compositionA)+2.90e-05_Angs*compositionA);
}


InGaAs::InGaAs(Temperature temper, double composition){
	compositionA=composition;
	compositionB=0;
	emass_gama=(1-compositionA)*0.067+compositionA*0.026-compositionA*(1-compositionA)*0.0091;
	emass_Lt=(1-compositionA)*0.0754+compositionA*0.05;
	emass_Ll=(1-compositionA)*1.9+compositionA*0.64;
	emass_Xt=(1-compositionA)*0.23+compositionA*0.16;
	emass_Xl=(1-compositionA)*1.3+compositionA*1.13;
	emass_Ldos=(1-compositionA)*0.56+compositionA*0.29;
	emass_Xdos=(1-compositionA)*0.85+compositionA*0.64;

	double lutt1GaAs=6.98;
	double lutt2GaAs=2.06;
	double lutt3GaAs=2.93;

	double lutt1InAs=20.0;
	double lutt2InAs=8.5;
	double lutt3InAs=9.2;

	lutt1=(1-compositionA)*lutt1GaAs+compositionA*lutt1InAs;
	lutt2=(1-compositionA)*lutt2GaAs+compositionA*lutt2InAs;
	lutt3=(1-compositionA)*lutt3GaAs+compositionA*lutt3InAs;

	double emassHH_tGaAs=HHmass_t(lutt1GaAs,lutt2GaAs,lutt3GaAs);
	double emassHH_lGaAs=HHmass_l(lutt1GaAs,lutt2GaAs,lutt3GaAs);
	double emasslH_tGaAs=lHmass_t(lutt1GaAs,lutt2GaAs,lutt3GaAs);
	double emasslH_lGaAs=lHmass_l(lutt1GaAs,lutt2GaAs,lutt3GaAs);

	double emassHH_tInAs=HHmass_t(lutt1InAs,lutt2InAs,lutt3InAs);
	double emassHH_lInAs=HHmass_l(lutt1InAs,lutt2InAs,lutt3InAs);
	double emasslH_tInAs=lHmass_t(lutt1InAs,lutt2InAs,lutt3InAs);
	double emasslH_lInAs=lHmass_l(lutt1InAs,lutt2InAs,lutt3InAs);

	emassHH_t=(1-compositionA)*emassHH_tGaAs+compositionA*emassHH_tInAs-compositionA*(1-compositionA)*(-0.145);
	emassHH_l=(1-compositionA)*emassHH_lGaAs+compositionA*emassHH_lInAs;
	emasslH_t=(1-compositionA)*emasslH_tGaAs+compositionA*emasslH_tInAs-compositionA*(1-compositionA)*0.0202;
	emasslH_l=(1-compositionA)*emasslH_lGaAs+compositionA*emasslH_lInAs;

	ac=(1-compositionA)*(-7.17_eV)+compositionA*(-5.08_eV);
	av=(1-compositionA)*(-1.16_eV)+compositionA*(-1.00_eV);
	b=(1-compositionA)*(-2.0_eV)+compositionA*(-1.8_eV);

	c11=(1-compositionA)*(122.1e+9_Pa)+compositionA*(83.29e+9_Pa);
	c12=(1-compositionA)*(56.6e+9_Pa)+compositionA*(45.26e+9_Pa);
	c44=(1-compositionA)*(60.0e+9_Pa)+compositionA*(39.59e+9_Pa);


	wLO=(1-compositionA)*(54.9934_THz)+compositionA*(45.2028_THz);
	wTO=(1-compositionA)*(50.6113_THz)+compositionA*(41.0576_THz);

	dielectric0 = (1-compositionA)*(12.90)+compositionA*(15.15);
	dielectricinf = (1-compositionA)*(10.92)+compositionA*(12.25);

	temp=temper;
	double gapgammaGaAs=varshini(temp,1.519_eV,0.5405_meV,204);
	double gapgammaInAs=varshini(temp,0.417_eV,0.276_meV,93);
	double	gapLGaAs=varshini(temp,1.815_eV,0.605_meV,204);
	double	gapLInAs=varshini(temp,1.133_eV,0.276_meV,93);
	double gapXGaAs=varshini(temp,1.981_eV,0.460_meV,204);
	double gapXInAs=varshini(temp,1.433_eV,0.276_meV,93);
	valencebandpos=compositionA*(0.21_eV)-compositionA*(1-compositionA)*(-0.38_eV);
	gapa=0;
	gapb=0;

	gapgamma=(1-compositionA)*gapgammaGaAs+compositionA*gapgammaInAs-compositionA*(1-compositionA)*0.477_eV;
	gapL=(1-compositionA)*gapLGaAs+compositionA*gapLInAs-compositionA*(1-compositionA)*1.4_eV;
	gapX=(1-compositionA)*gapXGaAs+compositionA*gapXInAs-compositionA*(1-compositionA)*0.33_eV;

	lparam=5.65325_Angs+0.40505_Angs*compositionA+(temp-300)*(3.88e-05_Angs*(1-compositionA)+2.74e-05_Angs*compositionA);

	materialname = "InGaAs";

}
void InGaAs::setTemperature(Temperature temper){
	temp=temper;
	double gapgammaGaAs=varshini(temp,1.519_eV,0.5405_meV,204);
	double gapgammaInAs=varshini(temp,0.417_eV,0.276_meV,93);
	double	gapLGaAs=varshini(temp,1.815_eV,0.605_meV,204);
	double	gapLInAs=varshini(temp,1.133_eV,0.276_meV,93);
	double gapXGaAs=varshini(temp,1.981_eV,0.460_meV,204);
	double gapXInAs=varshini(temp,1.433_eV,0.276_meV,93);

	gapgamma=(1-compositionA)*gapgammaGaAs+compositionA*gapgammaInAs-compositionA*(1-compositionA)*0.477_eV;
	gapL=(1-compositionA)*gapLGaAs+compositionA*gapLInAs-compositionA*(1-compositionA)*1.4_eV;
	gapX=(1-compositionA)*gapXGaAs+compositionA*gapXInAs-compositionA*(1-compositionA)*0.33_eV;

	lparam=5.65325_Angs+0.40505_Angs*compositionA+(temp-300)*(3.88e-05_Angs*(1-compositionA)+2.74e-05_Angs*compositionA);
}

InAlAs::InAlAs(Temperature temper, double composition){
	compositionA=composition;
	compositionB=0;
	emass_gama=(1-compositionA)*0.150+compositionA*0.026-compositionA*(1-compositionA)*0.049;
	emass_Lt=(1-compositionA)*.150+compositionA*0.05;
	emass_Ll=(1-compositionA)*1.32+compositionA*0.64;
	emass_Xt=(1-compositionA)*.22+compositionA*0.16;
	emass_Xl=(1-compositionA)*0.97+compositionA*1.13;
	emass_Ldos=(1-compositionA)*0.56+compositionA*0.29;
	emass_Xdos=(1-compositionA)*0.85+compositionA*0.64;

	double lutt1InAs=20.0;
	double lutt2InAs=8.5;
	double lutt3InAs=9.2;

	double lutt1AlAs=3.76;
	double lutt2AlAs=0.82;
	double lutt3AlAs=1.42;

	lutt1=(1-compositionA)*lutt1AlAs+compositionA*lutt1InAs;
	lutt2=(1-compositionA)*lutt2AlAs+compositionA*lutt2InAs;
	lutt3=(1-compositionA)*lutt3AlAs+compositionA*lutt3InAs;

	double emassHH_tAlAs=HHmass_t(lutt1AlAs,lutt2AlAs,lutt3AlAs);
	double emassHH_lAlAs=HHmass_l(lutt1AlAs,lutt2AlAs,lutt3AlAs);
	double emasslH_tAlAs=lHmass_t(lutt1AlAs,lutt2AlAs,lutt3AlAs);
	double emasslH_lAlAs=lHmass_l(lutt1AlAs,lutt2AlAs,lutt3AlAs);

	double emassHH_tInAs=HHmass_t(lutt1InAs,lutt2InAs,lutt3InAs);
	double emassHH_lInAs=HHmass_l(lutt1InAs,lutt2InAs,lutt3InAs);
	double emasslH_tInAs=lHmass_t(lutt1InAs,lutt2InAs,lutt3InAs);
	double emasslH_lInAs=lHmass_l(lutt1InAs,lutt2InAs,lutt3InAs);

	emassHH_t=(1-compositionA)*emassHH_tAlAs+compositionA*emassHH_tInAs;
	emassHH_l=(1-compositionA)*emassHH_lAlAs+compositionA*emassHH_lInAs;
	emasslH_t=(1-compositionA)*emasslH_tAlAs+compositionA*emasslH_tInAs;
	emasslH_l=(1-compositionA)*emasslH_lAlAs+compositionA*emasslH_lInAs;

	ac=(1-compositionA)*(-5.64_eV)+compositionA*(-5.08_eV)-(1-compositionA)*compositionA*(-1.4_eV);
	av=(1-compositionA)*(-2.47_eV)+compositionA*(-1.00_eV);
	b=(1-compositionA)*(-2.3_eV)+compositionA*(-1.8_eV);

	c11=(1-compositionA)*(125.0e+9_Pa)+compositionA*(83.29e+9_Pa);
	c12=(1-compositionA)*(53.4e+9_Pa)+compositionA*(45.26e+9_Pa);
	c44=(1-compositionA)*(54.2e+9_Pa)+compositionA*(39.59e+9_Pa);


	wLO=(1-compositionA)*(76.0749_THz)+compositionA*(45.2028_THz);
	wTO=(1-compositionA)*(67.9818_THz)+compositionA*(41.0576_THz);


	dielectric0 = (1-compositionA)*(10.06)+compositionA*(15.15);
	dielectricinf = (1-compositionA)*(8.16)+compositionA*(12.25);

	temp=temper;
	double gapgammaInAs=varshini(temp,0.417_eV,0.276_meV,93);
	double	gapLInAs=varshini(temp,1.133_eV,0.276_meV,93);
	double gapXInAs=varshini(temp,1.433_eV,0.276_meV,93);
	double gapgammaAlAs=varshini(temp,3.099_eV,0.885_meV,530);
	double	gapLAlAs=varshini(temp,2.46_eV,0.605_meV,204);
	double gapXAlAs=varshini(temp,2.24_eV,0.70_meV,530);
	valencebandpos=(1-compositionA)*(-0.53_eV)+compositionA*(0.21_eV)-compositionA*(1-compositionA)*(-0.64_eV);
	gapa=0;
	gapb=0;

	gapgamma=(1-compositionA)*gapgammaAlAs+compositionA*gapgammaInAs-compositionA*(1-compositionA)*0.70_eV;
	gapL=(1-compositionA)*gapLAlAs+compositionA*gapLInAs-compositionA*(1-compositionA)*0.0_eV;
	gapX=(1-compositionA)*gapXAlAs+compositionA*gapXInAs-compositionA*(1-compositionA)*0.0_eV;

	lparam=5.6611_Angs+0.3972_Angs*compositionA+(temp-300)*(2.90e-05_Angs*(1-compositionA)+2.74e-05_Angs*compositionA);

	materialname = "InAlAs";

}
void InAlAs::setTemperature(Temperature temper){
	temp=temper;
	double gapgammaInAs=varshini(temp,0.417_eV,0.276_meV,93);
	double	gapLInAs=varshini(temp,1.133_eV,0.276_meV,93);
	double gapXInAs=varshini(temp,1.433_eV,0.276_meV,93);
	double gapgammaAlAs=varshini(temp,3.099_eV,0.885_meV,530);
	double	gapLAlAs=varshini(temp,2.46_eV,0.605_meV,204);
	double gapXAlAs=varshini(temp,2.24_eV,0.70_meV,530);
	valencebandpos=(1-compositionA)*(-0.53_eV)+compositionA*(0.21_eV)-compositionA*(1-compositionA)*(-0.64_eV);

	gapgamma=(1-compositionA)*gapgammaAlAs+compositionA*gapgammaInAs-compositionA*(1-compositionA)*0.70_eV;
	gapL=(1-compositionA)*gapLAlAs+compositionA*gapLInAs-compositionA*(1-compositionA)*0.0_eV;
	gapX=(1-compositionA)*gapXAlAs+compositionA*gapXInAs-compositionA*(1-compositionA)*0.0_eV;

	lparam=5.6611_Angs+0.3972_Angs*compositionA+(temp-300)*(2.90e-05_Angs*(1-compositionA)+2.74e-05_Angs*compositionA);
}


InAlGaAs::InAlGaAs(Temperature temper, double compositionin, double compositional){
	compositionA=compositionin/(1-compositional);
	compositionB=compositional;
	emass_gama=(1-compositionA)*0.067+compositionA*0.026-compositionA*(1-compositionA)*0.0091;
	emass_Lt=(1-compositionA)*0.0754+compositionA*0.05;
	emass_Ll=(1-compositionA)*1.9+compositionA*0.64;
	emass_Xt=(1-compositionA)*0.23+compositionA*0.16;
	emass_Xl=(1-compositionA)*1.3+compositionA*1.13;
	emass_Ldos=(1-compositionA)*0.56+compositionA*0.29;
	emass_Xdos=(1-compositionA)*0.85+compositionA*0.64;

	emass_gama=(1-compositionB)*emass_gama+compositionB*0.150;
	emass_Lt=(1-compositionB)*emass_Lt+compositionB*0.150;
	emass_Ll=(1-compositionB)*emass_Ll+compositionB*1.32;
	emass_Xt=(1-compositionB)*emass_Xt+compositionB*0.22;
	emass_Xl=(1-compositionB)*emass_Xl+compositionB*0.97;
	emass_Ldos=(1-compositionB)*emass_Ldos+compositionB*0.78;
	emass_Xdos=(1-compositionB)*emass_Xdos+compositionB*0.76;


	double lutt1GaAs=6.98;
	double lutt2GaAs=2.06;
	double lutt3GaAs=2.93;

	double lutt1InAs=20.0;
	double lutt2InAs=8.5;
	double lutt3InAs=9.2;

	double lutt1AlAs=3.76;
	double lutt2AlAs=0.82;
	double lutt3AlAs=1.42;

	lutt1=(1-compositionA)*lutt1GaAs+compositionA*lutt1InAs;
	lutt2=(1-compositionA)*lutt2GaAs+compositionA*lutt2InAs;
	lutt3=(1-compositionA)*lutt3GaAs+compositionA*lutt3InAs;

	lutt1=(1-compositionB)*lutt1+compositionB*lutt1AlAs;
	lutt2=(1-compositionB)*lutt2+compositionB*lutt2AlAs;
	lutt3=(1-compositionB)*lutt3+compositionB*lutt3AlAs;

	double emassHH_tGaAs=HHmass_t(lutt1GaAs,lutt2GaAs,lutt3GaAs);
	double emassHH_lGaAs=HHmass_l(lutt1GaAs,lutt2GaAs,lutt3GaAs);
	double emasslH_tGaAs=lHmass_t(lutt1GaAs,lutt2GaAs,lutt3GaAs);
	double emasslH_lGaAs=lHmass_l(lutt1GaAs,lutt2GaAs,lutt3GaAs);

	double emassHH_tInAs=HHmass_t(lutt1InAs,lutt2InAs,lutt3InAs);
	double emassHH_lInAs=HHmass_l(lutt1InAs,lutt2InAs,lutt3InAs);
	double emasslH_tInAs=lHmass_t(lutt1InAs,lutt2InAs,lutt3InAs);
	double emasslH_lInAs=lHmass_l(lutt1InAs,lutt2InAs,lutt3InAs);

	double emassHH_tAlAs=HHmass_t(lutt1AlAs,lutt2AlAs,lutt3AlAs);
	double emassHH_lAlAs=HHmass_l(lutt1AlAs,lutt2AlAs,lutt3AlAs);
	double emasslH_tAlAs=lHmass_t(lutt1AlAs,lutt2AlAs,lutt3AlAs);
	double emasslH_lAlAs=lHmass_l(lutt1AlAs,lutt2AlAs,lutt3AlAs);

	emassHH_t=(1-compositionA)*emassHH_tGaAs+compositionA*emassHH_tInAs-compositionA*(1-compositionA)*(-0.145);
	emassHH_l=(1-compositionA)*emassHH_lGaAs+compositionA*emassHH_lInAs;
	emasslH_t=(1-compositionA)*emasslH_tGaAs+compositionA*emasslH_tInAs-compositionA*(1-compositionA)*0.0202;
	emasslH_l=(1-compositionA)*emasslH_lGaAs+compositionA*emasslH_lInAs;

	emassHH_t=(1-compositionB)*emassHH_t+compositionB*emassHH_tAlAs;
	emassHH_l=(1-compositionB)*emassHH_l+compositionB*emassHH_lAlAs;
	emasslH_t=(1-compositionB)*emasslH_t+compositionB*emasslH_tAlAs;
	emasslH_l=(1-compositionB)*emasslH_l+compositionB*emasslH_lAlAs;


	ac=(1-compositionA)*(-7.17_eV)+compositionA*(-5.08_eV);
	av=(1-compositionA)*(-1.16_eV)+compositionA*(-1.00_eV);
	b=(1-compositionA)*(-2.0_eV)+compositionA*(-1.8_eV);

	ac=(1-compositionB)*ac+compositionB*(-5.64_eV);
	av=(1-compositionB)*av+compositionB*(-2.47_eV);
	b=(1-compositionB)*b+compositionB*(-2.3_eV);


	c11=(1-compositionA)*(122.1e+9_Pa)+compositionA*(83.29e+9_Pa);
	c12=(1-compositionA)*(56.6e+9_Pa)+compositionA*(45.26e+9_Pa);
	c44=(1-compositionA)*(60.0e+9_Pa)+compositionA*(39.59e+9_Pa);

	c11=(1-compositionB)*c11+compositionB*(125.0e+9_Pa);
	c12=(1-compositionB)*c12+compositionB*(53.4e+9_Pa);
	c44=(1-compositionB)*c44+compositionB*(64.2e+9_Pa);


	wLO=(1-compositionA)*(54.9934_THz)+compositionA*(45.2028_THz);
	wTO=(1-compositionA)*(50.6113_THz)+compositionA*(41.0576_THz);

	wLO=(1-compositionB)*wLO+compositionB*(76.0749_THz);
	wTO=(1-compositionB)*wTO+compositionB*(67.9818_THz);

	dielectric0 = (1-compositionA)*(12.90)+compositionA*(15.15);
	dielectricinf = (1-compositionA)*(10.92)+compositionA*(12.25);

	dielectric0 = (1-compositionB)*dielectric0+compositionB*(10.06);
	dielectricinf = (1-compositionB)*dielectricinf+compositionB*(8.16);

	temp=temper;
	double gapgammaGaAs=varshini(temp,1.519_eV,0.5405_meV,204);
	double gapgammaInAs=varshini(temp,0.417_eV,0.276_meV,93);
	double gapgammaAlAs=varshini(temp,3.099_eV,0.885_meV,530);
	double	gapLGaAs=varshini(temp,1.815_eV,0.605_meV,204);
	double	gapLInAs=varshini(temp,1.133_eV,0.276_meV,93);
	double	gapLAlAs=varshini(temp,2.46_eV,0.605_meV,204);
	double gapXGaAs=varshini(temp,1.981_eV,0.460_meV,204);
	double gapXInAs=varshini(temp,1.433_eV,0.276_meV,93);
	double gapXAlAs=varshini(temp,2.24_eV,0.70_meV,530);
	valencebandpos=compositionA*(0.21_eV)-compositionA*(1-compositionA)*(-0.38_eV);
	valencebandpos=(1-compositional)*valencebandpos+compositional*(-0.53_eV);
	gapa=0;
	gapb=0;

	gapgamma=(1-compositionin-compositional)*gapgammaGaAs+compositionin*gapgammaInAs-compositionin*(1-compositionin)*0.477_eV
			+compositional*gapgammaAlAs-compositional*(1-compositional)*(-0.127_eV+1.310_eV*compositional)+compositionin*compositional*0.55_eV;
	gapL=(1-compositionin-compositionB)*gapLGaAs+compositionin*gapLInAs-compositionin*(1-compositionin)*1.4_eV+compositionB*gapLAlAs;
	gapX=(1-compositionin-compositionB)*gapXGaAs+compositionin*gapXInAs-compositionin*(1-compositionin)*0.33_eV+compositionB*gapXAlAs-compositionB*(1-compositionB)*0.055_eV;

	lparam=5.65325_Angs+0.40505_Angs*compositionin+0.0078_Angs*compositionB+(temp-300)*(3.88e-05_Angs*(1-compositionin-compositionB)+2.74e-05_Angs*compositionin+2.90e-05_Angs*compositionB);

	materialname = "InAlGaAs";


}
void InAlGaAs::setTemperature(Temperature temper){
	temp=temper;
	double gapgammaGaAs=varshini(temp,1.519_eV,0.5405_meV,204);
	double gapgammaInAs=varshini(temp,0.417_eV,0.276_meV,93);
	double gapgammaAlAs=varshini(temp,3.099_eV,0.885_meV,530);
	double	gapLGaAs=varshini(temp,1.815_eV,0.605_meV,204);
	double	gapLInAs=varshini(temp,1.133_eV,0.276_meV,93);
	double	gapLAlAs=varshini(temp,2.46_eV,0.605_meV,204);
	double gapXGaAs=varshini(temp,1.981_eV,0.460_meV,204);
	double gapXInAs=varshini(temp,1.433_eV,0.276_meV,93);
	double gapXAlAs=varshini(temp,2.24_eV,0.70_meV,530);
	valencebandpos=compositionA*(0.21_eV)-compositionA*(1-compositionA)*(-0.38_eV);
	valencebandpos=(1-compositionB)*valencebandpos+compositionB*(-0.53_eV);

	gapgamma=(1-compositionA-compositionB)*gapgammaGaAs+compositionA*gapgammaInAs-compositionA*(1-compositionA)*0.477_eV
			+compositionB*gapgammaAlAs-compositionB*(1-compositionB)*(-0.127_eV+1.310_eV*compositionB)+compositionA*compositionB*0.55_eV;
	gapL=(1-compositionA-compositionB)*gapLGaAs+compositionA*gapLInAs-compositionA*(1-compositionA)*1.4_eV+compositionB*gapLAlAs;
	gapX=(1-compositionA-compositionB)*gapXGaAs+compositionA*gapXInAs-compositionA*(1-compositionA)*0.33_eV+compositionB*gapXAlAs-compositionB*(1-compositionB)*0.055_eV;

	lparam=5.65325_Angs+0.40505_Angs*compositionA+0.0078_Angs*compositionB+(temp-300)*(3.88e-05_Angs*(1-compositionA-compositionB)+2.74e-05_Angs*compositionA+2.90e-05_Angs*compositionB);
}

ZnSe::ZnSe(Temperature temper){
	compositionA=0;
	compositionB=0;
	emass_gama=0.150;  //tamargo
	emass_Lt=0.160;  //not real
	emass_Ll=0.160;  //not real
	emass_Xt=0.160;  //not real
	emass_Xl=0.160;  //not real
	emass_Ldos=0.160;  //not real
	emass_Xdos=0.160;  //not real

	lutt1=3.94;  //adachi
	lutt2=1.0;  //
	lutt3=1.52;  //

	emassHH_t=0.66;  //tamargo
	emassHH_l=1.11;  //(adachi)
	emasslH_t=0.145;
	emasslH_l=0.143; //adachi

	alpha_gamma=0.73_meV;
	beta_gamma=295;
	alpha_L=0.73_meV;
	beta_L=295;
	alpha_X=0.73_meV;
	beta_X=295;
	gapgammatzero=2.800_eV;  //tamargo
	gapLtzero=3.8_eV;  //not real
	gapXtzero=3.4_eV;  //not real
	valencebandpos=-1.06_eV; //it need review
	gapa=73.0_meV;  //tamargo
	gapb=260;  //tamargo

	temp=temper;
	gapgamma=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapL=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapX=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	lparam=5.667_Angs+(temp-300)*7.4e-06_Angs;  //it need review

	ac=-7.51_eV;  //origin
	av=-8.6_eV;
	b=-1.24_eV;

	c11=90.3e+9_Pa; //landolt
	c12=53.6e+9_Pa;
	c44=39.4e+9_Pa;

	wLO=7.59_THz;  //landolt
	wTO=6.39_THz;

	dielectric0 = 8.6;   //landolt
	dielectricinf = 5.7;

	materialname = "ZnSe";

}

void ZnSe::setTemperature(Temperature temper){
	temp=temper;Material::
	gapgamma=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapL=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapX=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	lparam=5.667_Angs+(temp-300)*7.4e-06_Angs;  //it need review
}

ZnCdSe::ZnCdSe(Temperature temper,double compositionzn){
	compositionA=compositionzn;
	compositionB=0;
	emass_gama=(1-compositionA)*0.12+compositionA*0.150;  //landolt
	emass_Lt=(1-compositionA)*0.12+compositionA*0.150;  //not real
	emass_Ll=(1-compositionA)*0.12+compositionA*0.150;  //not real
	emass_Xt=(1-compositionA)*0.12+compositionA*0.150;  //not real
	emass_Xl=(1-compositionA)*0.12+compositionA*0.150;  //not real
	emass_Ldos=(1-compositionA)*0.12+compositionA*0.150;  //not real
	emass_Xdos=(1-compositionA)*0.12+compositionA*0.150;  //not real

	lutt1=(1-compositionA)*5.51+compositionA*3.94;  //not real adachi
	lutt2=(1-compositionA)*1.24+compositionA*1.0;  //not real
	lutt3=(1-compositionA)*2.14+compositionA*1.52;  //not real


	emassHH_t=(1-compositionA)*0.9+compositionA*0.66;
	emassHH_l=(1-compositionA)*1.7+compositionA*1.11;
	emasslH_t=(1-compositionA)*0.18+compositionA*0.145;
	emasslH_l=(1-compositionA)*0.16+compositionA*0.143;

	alpha_gamma=(1-compositionA)*0.424_meV+compositionA*0.73_meV;
	beta_gamma=(1-compositionA)*118+compositionA*295;
	alpha_L=(1-compositionA)*0.424_meV+compositionA*0.73_meV;
	beta_L=(1-compositionA)*118+compositionA*295;
	alpha_X=(1-compositionA)*0.424_meV+compositionA*0.73_meV;
	beta_X=(1-compositionA)*118+compositionA*295;
	gapgammatzero=(1-compositionA)*1.830_eV+compositionA*2.800_eV-compositionA*(1-compositionA)*0.387_eV;  //tamargo
	gapLtzero=(1-compositionA)*3.47_eV+compositionA*3.8_eV;  //not real
	gapXtzero=(1-compositionA)*3.47_eV+compositionA*3.4_eV;  //not real
	valencebandpos=(1-compositionA)*-0.934_eV+compositionA*-1.06_eV; //it need review 0.80?
	gapa=(1-compositionA)*36.0_meV+compositionA*73.0_meV;  //tamargo
	gapb=(1-compositionA)*179+compositionA*260;

	temp=temper;
	gapgamma=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapL=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapX=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	lparam=((1-compositionA)*6.077_Angs+compositionA*5.667_Angs)+(temp-300)*((1-compositionA)*7.4e-06_Angs+compositionA*7.4e-06_Angs);  //it need review

	ac=(1-compositionA)*(-11.0_eV)+compositionA*(-4.17_eV);  //origin
	av=(1-compositionA)*(-8.9_eV)+compositionA*(-8.6_eV);
	b=(1-compositionA)*(-0.80_eV)+compositionA*(-1.24_eV);

    c11=(1-compositionA)*74.6e+9_Pa+compositionA*90.3e+9_Pa; //landolt
	c12=(1-compositionA)*46.1e+9_Pa+compositionA*53.6e+9_Pa;
	c44=(1-compositionA)*13.0e+9_Pa+compositionA*39.4e+9_Pa;

	wLO=(1-compositionA)*6.35_THz+compositionA*7.59_THz;  //landolt
	wTO=(1-compositionA)*5.06_THz+compositionA*6.39_THz;

	dielectric0 = (1-compositionA)*9.29+compositionA*8.6;  //landolt
	dielectricinf = (1-compositionA)*6.20+compositionA*5.7;

	materialname = "ZnCdSe";

}

void ZnCdSe::setTemperature(Temperature temper){
	temp=temper;
	gapgamma=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapL=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapX=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	lparam=((1-compositionA)*6.077_Angs+compositionA*5.667_Angs)+(temp-300)*((1-compositionA)*7.4e-06_Angs+compositionA*7.4e-06_Angs);  //it need review
}

MgZnCdSe::MgZnCdSe(Temperature temper,double compositionzn,double compositionmg){
	compositionA=compositionzn/(1-compositionmg);
	compositionB=compositionmg;
	emass_gama=(1-compositionA)*0.12+compositionA*0.150;  //landolt
	emass_Lt=(1-compositionA)*0.12+compositionA*0.150;  //not real
	emass_Ll=(1-compositionA)*0.12+compositionA*0.150;  //not real
	emass_Xt=(1-compositionA)*0.12+compositionA*0.150;  //not real
	emass_Xl=(1-compositionA)*0.12+compositionA*0.150;  //not real
	emass_Ldos=(1-compositionA)*0.12+compositionA*0.150;  //not real
	emass_Xdos=(1-compositionA)*0.12+compositionA*0.150;  //not real

	emass_gama=(1-compositionB)*emass_gama+compositionB*0.20;  //landolt
	emass_Lt=(1-compositionB)*emass_gama+compositionB*0.20;  //not real
	emass_Ll=(1-compositionB)*emass_gama+compositionB*0.20;  //not real
	emass_Xt=(1-compositionB)*emass_gama+compositionB*0.20;  //not real
	emass_Xl=(1-compositionB)*emass_gama+compositionB*0.20;  //not real
	emass_Ldos=(1-compositionB)*emass_gama+compositionB*0.20;  //not real
	emass_Xdos=(1-compositionB)*emass_gama+compositionB*0.20;  //not real


	lutt1=(1-compositionA)*5.51+compositionA*3.94;  //not real adachi
	lutt2=(1-compositionA)*1.24+compositionA*1.0;  //not real
	lutt3=(1-compositionA)*2.14+compositionA*1.52;  //not real

	lutt1=(1-compositionB)*lutt1+compositionB*2.84;  //not real adachi
	lutt2=(1-compositionB)*lutt2+compositionB*0.43;  //not real
	lutt3=(1-compositionB)*lutt3+compositionB*1.0;  //not real

	emassHH_t=(1-compositionA)*0.9+compositionA*0.66;
	emassHH_l=(1-compositionA)*1.7+compositionA*1.11;
	emasslH_t=(1-compositionA)*0.18+compositionA*0.145;
	emasslH_l=(1-compositionA)*0.16+compositionA*0.143;

	emassHH_t=(1-compositionB)*0.9+compositionB*0.78;
	emassHH_l=(1-compositionB)*1.7+compositionB*HHmass_l(2.84,0.43,1.0);
	emasslH_t=(1-compositionB)*0.18+compositionB*0.33;
	emasslH_l=(1-compositionB)*0.16+compositionB*lHmass_l(2.84,0.43,1.0);

	alpha_gamma=(1-compositionA)*0.424_meV+compositionA*0.73_meV;
	beta_gamma=(1-compositionA)*118+compositionA*295;
	alpha_L=(1-compositionA)*0.424_meV+compositionA*0.73_meV;
	beta_L=(1-compositionA)*118+compositionA*295;
	alpha_X=(1-compositionA)*0.424_meV+compositionA*0.73_meV;
	beta_X=(1-compositionA)*118+compositionA*295;
	gapgammatzero=(1-compositionA)*1.830_eV+compositionA*2.800_eV-compositionA*(1-compositionA)*0.387_eV;  //tamargo
	gapLtzero=(1-compositionA)*3.47_eV+compositionA*3.8_eV;  //not real
	gapXtzero=(1-compositionA)*3.47_eV+compositionA*3.4_eV;  //not real
	valencebandpos=(1-compositionA)*-0.934_eV+compositionA*-1.06_eV; //it need review 0.80?
	gapa=(1-compositionA)*36.0_meV+compositionA*73.0_meV;  //tamargo
	gapb=(1-compositionA)*179+compositionA*260;

	gapgammatzero=(1-compositionB)*gapgammatzero+compositionB*4.05_eV-compositionB*(1-compositionB)*0.400_eV;  //tamargo
	gapLtzero=(1-compositionB)*gapLtzero+compositionB*3.4_eV;  //not real
	gapXtzero=(1-compositionB)*gapXtzero+compositionB*3.4_eV;  //not real
	valencebandpos=(1-compositionB)*valencebandpos+compositionB*-1.351_eV;

	temp=temper;
	gapgamma=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapL=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapX=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	lparam=(1-compositionA)*6.077_Angs+compositionA*5.667_Angs;  //it need review
	lparam=((1-compositionB)*lparam+compositionB*5.904_Angs)+(temp-300)*7.4e-06_Angs;  //it need review

	c11=75.8e+9_Pa; //landolt
	c12=48.6e+9_Pa;
	c44=31.7e+9_Pa;

	wLO=10.19_THz;  //landolt
	wTO=7.10_THz;

	dielectric0 = 3.8;  //landolt
	dielectricinf = 3.8;

	ac=(1-compositionA)*-11.0_eV+compositionA*-7.51_eV;  //origin
	av=(1-compositionA)*-8.9_eV+compositionA*-8.6_eV;
	b=(1-compositionA)*-0.80_eV+compositionA*-1.24_eV;

	ac=(1-compositionB)*ac+compositionB*-7.5_eV;  //origin
	av=(1-compositionB)*av+compositionB*-8.6_eV;
	b=(1-compositionB)*b+compositionB*-1.27_eV;

	c11=(1-compositionA)*74.6e+9_Pa+compositionA*90.3e+9_Pa; //landolt
	c12=(1-compositionA)*46.1e+9_Pa+compositionA*53.6e+9_Pa;
	c44=(1-compositionA)*13.0e+9_Pa+compositionA*39.4e+9_Pa;

	c11=(1-compositionB)*c11+compositionB*75.8e+9_Pa; //landolt
	c12=(1-compositionB)*c12+compositionB*48.6e+9_Pa;
	c44=(1-compositionB)*c44+compositionB*31.7e+9_Pa;

	wLO=(1-compositionA)*6.35_THz+compositionA*7.59_THz;  //landolt
	wTO=(1-compositionA)*5.06_THz+compositionA*6.39_THz;

	wLO=(1-compositionB)*wLO+compositionB*10.19_THz;  //landolt
	wTO=(1-compositionB)*wTO+compositionB*7.10_THz;

	dielectric0 = (1-compositionA)*9.29+compositionA*8.6;  //landolt
	dielectricinf = (1-compositionA)*6.20+compositionA*5.7;

	dielectric0 = (1-compositionB)*dielectric0+compositionB*3.8;  //landolt
	dielectricinf = (1-compositionB)*dielectricinf+compositionB*3.8;

	materialname = "MgZnCdSe";
}

void MgZnCdSe::setTemperature(Temperature temper){
	temp=temper;
	gapgamma=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapL=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapX=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	lparam=((1-compositionA)*6.077_Angs+compositionA*5.667_Angs)+(temp-300)*((1-compositionA)*7.4e-06_Angs+compositionA*7.4e-06_Angs);  //it need review
}


CdSe::CdSe(Temperature temper){
	compositionA=0;
	compositionB=0;
	emass_gama=0.12;  //landolt
	emass_Lt=0.12;  //not real
	emass_Ll=0.12;  //not real
	emass_Xt=0.12;  //not real
	emass_Xl=0.12;  //not real
	emass_Ldos=0.12;  //not real
	emass_Xdos=0.12;  //not real

	lutt1=5.51;  //not real adachi
	lutt2=1.24;  //not real
	lutt3=2.14;  //not real

	emassHH_t=0.9;  //landolt
	emassHH_l=1.7;
	emasslH_t=0.18;
	emasslH_l=0.16;

	alpha_gamma=0.424_meV;
	beta_gamma=118;
	alpha_L=0.424_meV;
	beta_L=118;
	alpha_X=0.424_meV;
	beta_X=118;
	gapgammatzero=1.830_eV;  //tamargo
	gapLtzero=3.47_eV;  //not real
	gapXtzero=3.47_eV;  //not real
	valencebandpos=-0.934_eV; //it need review 0.80?
	gapa=36.0_meV;  //tamargo
	gapb=179;

	temp=temper;
	gapgamma=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapL=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapX=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	lparam=6.077_Angs+(temp-300)*7.4e-06_Angs;  //it need review

	ac=-11.0_eV;  //origin
	av=-8.9_eV;
	b=-0.80_eV;

	c11=74.6e+9_Pa; //landolt
	c12=46.1e+9_Pa;
	c44=13.0e+9_Pa;

	wLO=6.35_THz;  //landolt
	wTO=5.06_THz;

	dielectric0 = 9.29;  //landolt
	dielectricinf = 6.20;

	materialname = "CdSe";

}

void CdSe::setTemperature(Temperature temper){
	temp=temper;
	gapgamma=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapL=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapX=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	lparam=6.077_Angs+(temp-300)*7.4e-06_Angs;  //it need review
}


MgSe::MgSe(Temperature temper){
	compositionA=0;
	compositionB=0;
	emass_gama=0.2;  //adachi
	emass_Lt=0.2;  //not real
	emass_Ll=0.2;  //not real
	emass_Xt=0.2;  //not real
	emass_Xl=0.2;  //not real
	emass_Ldos=0.2;  //not real
	emass_Xdos=0.2;  //not real

	lutt1=2.84;  //adachi
	lutt2=0.43;  //
	lutt3=1.0;  //

	//emassHH_t=HHmass_t(lutt1,lutt2,lutt3);
	emassHH_l=HHmass_l(lutt1,lutt2,lutt3);
	//emasslH_t=lHmass_t(lutt1,lutt2,lutt3);
	emasslH_l=lHmass_l(lutt1,lutt2,lutt3);

	emassHH_t=0.78;  //tamargo
	//emassHH_l=1.7;
	emasslH_t=0.33;
	//emasslH_l=0.16;

	alpha_gamma=0.424_meV;
	beta_gamma=118;
	alpha_L=0.424_meV;
	beta_L=118;
	alpha_X=0.424_meV;
	beta_X=118;
	gapgammatzero=4.05_eV;  //landolt
	gapLtzero=3.4_eV;  //not real
	gapXtzero=3.4_eV;  //not real
	valencebandpos=-1.351_eV; //it need review
	gapa=36.0_meV;  //tamargo
	gapb=179;
	temp=temper;
	gapgamma=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapL=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapX=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	lparam=5.904_Angs+(temp-300)*7.4e-06_Angs;  //it need review

	ac=-7.5_eV;  //origin
	av=-8.6_eV;
	b=-1.27_eV;

	c11=75.8e+9_Pa; //landolt
	c12=48.6e+9_Pa;
	c44=31.7e+9_Pa;

	wLO=10.19_THz;  //landolt
	wTO=7.10_THz;

	dielectric0 = 3.8;  //landolt
	dielectricinf = 3.8;

	materialname = "MgSe";

}

void MgSe::setTemperature(Temperature temper){
	temp=temper;
	gapgamma=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapL=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapX=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	lparam=5.904_Angs+(temp-300)*7.4e-06_Angs;  //it need review
}

CdTe::CdTe(Temperature temper){
	compositionA=0;
	compositionB=0;
	emass_gama=0.094;  //adachi
	emass_Lt=0.094;  //not real
	emass_Ll=0.094;  //not real
	emass_Xt=0.094;  //not real
	emass_Xl=0.094;  //not real
	emass_Ldos=0.094;  //not real
	emass_Xdos=0.094;  //not real

	lutt1=4.14;  //adachi
	lutt2=1.09;  //
	lutt3=1.62;  //

	emassHH_t=HHmass_t(lutt1,lutt2,lutt3);
	emassHH_l=HHmass_l(lutt1,lutt2,lutt3);
	emasslH_t=lHmass_t(lutt1,lutt2,lutt3);
	emasslH_l=lHmass_l(lutt1,lutt2,lutt3);

	//emassHH_l=0.81; //landolt
	//emasslH_l=0.12;

	alpha_gamma=0.46_meV;
	beta_gamma=160;
	alpha_L=0.46_meV;
	beta_L=160;
	alpha_X=0.46_meV;
	beta_X=160;
	gapgammatzero=1.597_eV;  //origin
	gapLtzero=3.48_eV;  //not real
	gapXtzero=3.48_eV;
	valencebandpos=0.12_eV; //it need review
	gapa=36.0_meV;  //tamargo
	gapb=179;

	temp=temper;
	gapgamma=varshini(temp,gapgammatzero,alpha_gamma,beta_gamma);
	gapL=varshini(temp,gapLtzero,alpha_L,beta_L);
	gapX=varshini(temp,gapXtzero,alpha_X,beta_X);
	lparam=6.481_Angs+(temp-300)*4.67e-06_Angs;  //origin/landolt

	ac=-4.8_eV;  //origin
	av=-3.96_eV;
	b=-1.18_eV;

	c11=75.8e+9_Pa; //landolt
	c12=48.6e+9_Pa;
	c44=31.7e+9_Pa;

	wLO=5.08_THz;  //landolt
	wTO=4.20_THz;

	dielectric0 = 10.4;  //landolt
	dielectricinf = 7.1;

	materialname = "CdTe";

}

void CdTe::setTemperature(Temperature temper){
	temp=temper;
	gapgamma=varshini(temp,gapgammatzero,alpha_gamma,beta_gamma);
	gapL=varshini(temp,gapLtzero,alpha_L,beta_L);
	gapX=varshini(temp,gapXtzero,alpha_X,beta_X);
	lparam=6.481_Angs+(temp-300)*4.67e-06_Angs;  //origin/landolt
}


ZnTe::ZnTe(Temperature temper){
	compositionA=0;
	compositionB=0;
	emass_gama=0.123;  //landolt
	emass_Lt=0.123;  //not real
	emass_Ll=0.123;  //not real
	emass_Xt=0.123;  //not real
	emass_Xl=0.123;  //not real
	emass_Ldos=0.123;  //not real
	emass_Xdos=0.123;  //not real

	lutt1=3.96;  //adachi
	lutt2=0.86;  //
	lutt3=1.39;  //

	emassHH_t=HHmass_t(lutt1,lutt2,lutt3);
	emassHH_l=HHmass_l(lutt1,lutt2,lutt3);
	emasslH_t=lHmass_t(lutt1,lutt2,lutt3);
	emasslH_l=lHmass_l(lutt1,lutt2,lutt3);

	//emassHH_t=0.45;  //adachi
	//emassHH_l=0.85;
	//emasslH_t=0.176;
	//emasslH_l=0.148;

	alpha_gamma=0.424_meV;
	beta_gamma=118;
	alpha_L=0.424_meV;
	beta_L=118;
	alpha_X=0.424_meV;
	beta_X=118;
	gapgammatzero=2.3945_eV;  //landolt
	gapLtzero=3.05_eV;  //not real
	gapXtzero=3.05_eV;
	valencebandpos=-0.2_eV; //it need review
	gapa=37.5_meV;  //tamargo
	gapb=163;

	temp=temper;
	gapgamma=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapL=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapX=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	lparam=6.0882_Angs+(temp-300)*8.33e-06_Angs;  //origin/landolt

	ac=-5.19_eV;  //origin
	av=-0.07_eV;
	b=-1.26_eV;

	c11=72.2e+9_Pa; //landolt
	c12=40.9e+9_Pa;
	c44=30.8e+9_Pa;

	wLO=6.20_THz;  //landolt
	wTO=5.30_THz;

	dielectric0 = 10.3;  //landolt
	dielectricinf = 7.28;

	materialname = "ZnTe";

}

void ZnTe::setTemperature(Temperature temper){
	temp=temper;
	gapgamma=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapL=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapX=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	lparam=6.0882_Angs+(temp-300)*8.33e-06_Angs;  //origin/landolt
}



/*
Bi2Se3::Bi2Se3(Temperature temper){
	compositionA=0;
	compositionB=0;
	emass_gama=0.02;  //landolt
	emass_Lt=0.123;  //not real
	emass_Ll=0.123;  //not real
	emass_Xt=0.123;  //not real
	emass_Xl=0.123;  //not real
	emass_Ldos=0.123;  //not real
	emass_Xdos=0.123;  //not real

	lutt1=3.96;  //adachi
	lutt2=0.86;  //
	lutt3=1.39;  //

	emassHH_t=HHmass_t(lutt1,lutt2,lutt3);
	emassHH_l=HHmass_l(lutt1,lutt2,lutt3);
	emasslH_t=lHmass_t(lutt1,lutt2,lutt3);
	emasslH_l=lHmass_l(lutt1,lutt2,lutt3);

	//emassHH_t=0.45;  //adachi
	//emassHH_l=0.85;
	//emasslH_t=0.176;
	//emasslH_l=0.148;

	alpha_gamma=0.424_meV;
	beta_gamma=118;
	alpha_L=0.424_meV;
	beta_L=118;
	alpha_X=0.424_meV;
	beta_X=118;
	gapgammatzero=0.160_eV;  //landolt
	gapLtzero=3.48_eV;  //not real
	gapXtzero=3.48_eV;
	valencebandpos=-0.2_eV; //it need review
	gapa=36_meV;  //tamargo
	gapb=179;

	temp=temper;
	gapgamma=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapL=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	gapX=boseeinsteingap(temp,gapgammatzero,gapa,gapb);
	lparam=6.0882_Angs+(temp-300)*5.15e-05_Angs;  //origin/landolt

	ac=-5.19_eV;  //origin
	av=-0.07_eV;
	b=-1.26_eV;

	c11=72.2e+9_Pa; //landolt
	c12=40.9e+9_Pa;
	c44=30.8e+9_Pa;

	wLO=6.20_THz;  //landolt
	wTO=5.30_THz;

	dielectric0 = 10.3;  //landolt
	dielectricinf = 7.28;

}

*/
