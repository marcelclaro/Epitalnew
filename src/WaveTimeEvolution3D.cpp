/*
 * WaveTimeEvolution3D.cpp
 *
 *  Created on: Nov 7, 2014
 *      Author: marcel
 */

#include "WaveTimeEvolution3D.hpp"

#ifndef WAVETIMEEVOLUTION3D_CPP_
#define WAVETIMEEVOLUTION3D_CPP_

namespace epital {

template <typename timetype, typename complextype>
WaveTimeEvolution3D<timetype,complextype>::WaveTimeEvolution3D(WaveTimeEvolution3D::Carrier carrier): carriertype(carrier) {
	// TODO Auto-generated constructor stub

}

template <typename timetype, typename complextype>
WaveTimeEvolution3D<timetype,complextype>::~WaveTimeEvolution3D() {
	// TODO Auto-generated destructor stub
}

template <typename timetype, typename complextype>
template<typename T, typename gridtype>
DiscreteFunction3D<std::complex<complextype>,gridtype> WaveTimeEvolution3D<timetype,complextype>::potentialpropagator(timetype methodtimestep, DiscreteFunction3D<T,gridtype> potential){
	std::complex<complextype> energysign;

	if(carriertype==Carrier::ELECTRON)
		energysign=std::complex<complextype>(1.0,0);
	else
		energysign=std::complex<complextype>(-1.0,0);

	DiscreteFunction3D<std::complex<complextype>,gridtype> result(potential.getGrid());


	for(long int k=0; k<potential.getGrid().getSizeX();k++)
		for(long int l=0; l<potential.getGrid().getSizeY();l++)
			for(long int m=0; m<potential.getGrid().getSizeZ();m++){
				result(k,l,m)=exp(std::complex<complextype>(0.0,-1.0)*energysign*methodtimestep*potential(k,l,m)/(2.0*Constant::hbar));
			}

	return result;
}

template <typename timetype, typename complextype>
template<typename T, typename gridtype>
DiscreteFunction3D<std::complex<complextype>,gridtype>& WaveTimeEvolution3D<timetype,complextype>::groundlevel_periodicFFT(timetype time, long int timesteps, DiscreteFunction3D<std::complex<complextype>,gridtype>& initial, DiscreteFunction3D<T,gridtype>& potential, T effectivemass){
	timetype methodtimestep=time/static_cast<timetype>(timesteps);


	initial/=initial.Modulus();

	DiscreteFunction3D<std::complex<complextype>,gridtype> potpropag = potentialpropagator(methodtimestep,potential);
	DiscreteFunction3D<std::complex<complextype>,gridtype> bloch = fftmulti(methodtimestep,initial.getGrid(),effectivemass);

	for(long int i = 0 ; i< timesteps; i++){
		//std::cout << "step: " << i  << " ..." << std::endl;
		initial*=potpropag;
		initial.FFT();
		initial*=bloch;
		initial.iFFT();
		initial*=potpropag;
		std::complex<complextype> tonorm=initial.Modulus();
		initial*=1.0/tonorm;


	}

	return initial;
}

template <typename timetype, typename complextype>
template<typename T, typename gridtype>
DiscreteFunction3D<std::complex<complextype>,gridtype>& WaveTimeEvolution3D<timetype,complextype>::groundlevel_finite_valence(timetype time, long int timesteps, DiscreteFunction3D<std::complex<complextype>,gridtype>& initial, DiscreteFunction3D<T,gridtype>& potential, DiscreteFunction3D<T,gridtype>& effectivemasst,DiscreteFunction3D<T,gridtype>& effectivemassl){
	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	initial/=initial.Modulus();

	DiscreteFunction3D<std::complex<complextype>,gridtype> potpropag = potentialpropagator(methodtimestep,potential);

	timetype adi_methodtimestep = methodtimestep/3.0;

	gridtype gridstepx = initial.getGrid().getxIncrement();
	gridtype gridstepy = initial.getGrid().getyIncrement();
	gridtype gridstepz = initial.getGrid().getzIncrement();


	long int xmax=initial.getGrid().getSizeX();
	long int ymax=initial.getGrid().getSizeY();
	long int zmax=initial.getGrid().getSizeZ();

	long int xlast=xmax-1;
	long int ylast=ymax-1;
	long int zlast=zmax-1;


	DiscreteFunction3D<std::complex<complextype>,gridtype> as(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> bs(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> cs(initial.getGrid());

	DiscreteFunction3D<std::complex<complextype>,gridtype> betax(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> betay(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> betaz(initial.getGrid());

	#pragma omp parallel for
	for(long int i = 0; i< betax.getGrid().getSize();i++){
		betax.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*adi_methodtimestep/(4.0*effectivemassl.getdata()[i]*gridstepx*gridstepx);
	}

	#pragma omp parallel for
	for(long int i = 0; i< betay.getGrid().getSize();i++){
		betay.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*adi_methodtimestep/(4.0*effectivemassl.getdata()[i]*gridstepy*gridstepy);
	}

	#pragma omp parallel for
	for(long int i = 0; i< betaz.getGrid().getSize();i++){
		betaz.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*adi_methodtimestep/(4.0*effectivemasst.getdata()[i]*gridstepz*gridstepz);
	}


	DiscreteFunction3D<std::complex<complextype>,gridtype> x(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> y(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> z(initial.getGrid());

	gridtype dx = gridstepx;
	gridtype dy = gridstepy;
	gridtype dz = gridstepz;

	for(long int l=0;l<ymax;++l){
		for(long int m=0;m<zmax;++m){
			//edge: order h and h²
			x(0,l,m) = (effectivemassl(0,l,m)*dx/2.0)*(1.0/effectivemassl(1,l,m)-1.0/effectivemassl(0,l,m))/dx;
			x(1,l,m) = (effectivemassl(1,l,m)*dx/2.0)*(1.0/effectivemassl(2,l,m)-1.0/effectivemassl(0,l,m))/(2*dx);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int k = 2; k < (xmax-2); k++){
				x(k,l,m)=(effectivemassl(k,l,m)*dx/2.0)*((1.0/effectivemassl(k+1,l,m))*8.0-(1.0/effectivemassl(k+2,l,m))-(1.0/effectivemassl(k-1,l,m))*8.0+(1.0/effectivemassl(k-2,l,m)))/(12.0*dx);
			}

			//edge: order h and h²
			x(xmax-2,l,m)= (effectivemassl(xmax-2,l,m)*dx/2.0)*(1.0/effectivemassl(xmax-1,l,m)-1.0/effectivemassl(xmax-3,l,m))/(2.0*dx);
			x(xmax-1,l,m)= (effectivemassl(xmax-1,l,m)*dx/2.0)*(1.0/effectivemassl(xmax-1,l,m)-1.0/effectivemassl(xmax-2,l,m))/dx;
		}
	}

	for(long int k=0;k<xmax;++k){
		for(long int m=0;m<zmax;++m){
			//edge: order h and h²
			y(k,0,m) = (effectivemassl(k,0,m)*dy/2.0)*(1.0/effectivemassl(k,1,m)-1.0/effectivemassl(k,0,m))/dy;
			y(k,1,m) = (effectivemassl(k,1,m)*dy/2.0)*(1.0/effectivemassl(k,2,m)-1.0/effectivemassl(k,0,m))/(2*dy);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int l = 2; l < (ymax-2); l++){
				y(k,l,m)=(effectivemassl(k,l,m)*dy/2.0)*((1.0/effectivemassl(k,l+1,m))*8.0-(1.0/effectivemassl(k,l+2,m))-(1.0/effectivemassl(k,l-1,m))*8.0+(1.0/effectivemassl(k,l-2,m)))/(12.0*dy);
			}

			//edge: order h and h²
			y(k,ymax-2,m)= (effectivemassl(k,ymax-2,m)*dy/2.0)*(1.0/effectivemassl(k,ymax-1,m)-1.0/effectivemassl(k,ymax-3,m))/(2.0*dy);
			y(k,ymax-1,m)= (effectivemassl(k,ymax-1,m)*dy/2.0)*(1.0/effectivemassl(k,ymax-1,m)-1.0/effectivemassl(k,ymax-2,m))/dy;
		}
	}

	for(long int k=0;k<xmax;++k){
		for(long int l=0;l<ymax;++l){
			//edge: order h and h²
			z(k,l,0) = (effectivemasst(k,l,0)*dz/2.0)*(1.0/effectivemasst(k,l,1)-1.0/effectivemasst(k,l,0))/dz;
			z(k,l,1) = (effectivemasst(k,l,1)*dz/2.0)*(1.0/effectivemasst(k,l,2)-1.0/effectivemasst(k,l,0))/(2*dz);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int m = 2; m < (zmax-2); m++){
				z(k,l,m)=(effectivemasst(k,l,m)*dz/2.0)*((1.0/effectivemasst(k,l,m+1))*8.0-(1.0/effectivemasst(k,l,m+2))-(1.0/effectivemasst(k,l,m-1))*8.0+(1.0/effectivemasst(k,l,m-2)))/(12.0*dz);
			}

			//edge: order h and h²
			z(k,l,zmax-2)= (effectivemasst(k,l,zmax-2)*dz/2.0)*(1.0/effectivemasst(k,l,zmax-1)-1.0/effectivemasst(k,l,zmax-3))/(2.0*dz);
			z(k,l,zmax-1)= (effectivemasst(k,l,zmax-1)*dz/2.0)*(1.0/effectivemasst(k,l,zmax-1)-1.0/effectivemasst(k,l,zmax-2))/dz;
		}
	}



	for(long int i = 0 ; i< timesteps; i++){
		std::cout << "step: " << i  << " ..." << std::endl;
		initial*=potpropag;

		{
			/* x direction*/

			DiscreteFunction3D<std::complex<complextype>,gridtype> copy(initial);

			#pragma omp parallel for
			for(long int l=0;l<ymax;++l)
				for(long int m=0;m<zmax;++m){
					//Multiply vector
					initial(0,l,m)=betax(0,l,m)*copy(1,l,m)*(1.0+x(0,l,m))+(1.0-2.0*betax(0,l,m))*copy(0,l,m);
					initial(xlast,l,m)=betax(xlast,l,m)*copy(xlast-1,l,m)*(1.0-x(xlast,l,m))+(1.0-2.0*betax(xlast,l,m))*copy(xlast,l,m);

					for(long int k = 1; k < xlast; k++){
							initial(k,l,m)=betax(k,l,m)*copy(k-1,l,m)*(1.0-x(k,l,m))+(1.0-2.0*betax(k,l,m))*copy(k,l,m)+betax(k,l,m)*copy(k+1,l,m)*(1.0+x(k,l,m));
					}
			}


			#pragma omp parallel for
			for(long int l=0;l<ymax;++l)
				for(long int m=0;m<zmax;++m){
					as(0,l,m)=0;
					bs(0,l,m)=2.0*betax(0,l,m)+1.0;
					cs(0,l,m)=-betax(0,l,m)*(1.0+x(0,l,m));

					for(long int i = 1; i < xlast ; i++){
						as(i,l,m)=-betax(i,l,m)*(1.0-x(i,l,m));
						bs(i,l,m)=2.0*betax(i,l,m)+1.0;
						cs(i,l,m)=-betax(i,l,m)*(1.0+x(i,l,m));
					}

					as(xlast,l,m)=-betax(xlast,l,m)*(1.0-x(xlast,l,m));
					bs(xlast,l,m)=2.0*betax(xlast,l,m)+1.0;
					cs(xlast,l,m)=0;

					//Thomas Algorithm

					for (long int i = 1; i <= xlast; i++){
						std::complex<complextype> n = as(i,l,m)/bs(i-1,l,m);
						bs(i,l,m)-=n*cs(i-1,l,m);
						initial(i,l,m)-=n*initial(i-1,l,m);
					}

					initial(xlast,l,m)/=bs(xlast,l,m);
					for (long int i = (xlast - 1); i >= 0; --i)
						initial(i,l,m)=(initial(i,l,m)-cs(i,l,m)*initial(i+1,l,m))/bs(i,l,m);

				}

			/* y direction*/

			copy=initial;

			#pragma omp parallel for
			for(long int k=0;k<xmax;++k)
				for(long int m=0;m<zmax;++m){
					//Multiply vector
					initial(k,0,m)=betay(k,0,m)*copy(k,1,m)*(1.0+y(k,0,m))+(1.0-2.0*betay(k,0,m))*copy(k,0,m);
					initial(k,ylast,m)=betay(k,ylast,m)*copy(k,ylast-1,m)*(1.0-y(k,ylast,m))+(1.0-2.0*betay(k,ylast,m))*copy(k,ylast,m);

					for(long int l = 1; l < ylast; l++){
							initial(k,l,m)=betay(k,l,m)*copy(k,l-1,m)*(1.0-y(k,l,m))+(1.0-2.0*betay(k,l,m))*copy(k,l,m)+betay(k,l,m)*copy(k,l+1,m)*(1.0+y(k,l,m));
					}
				}


			#pragma omp parallel for
			for(long int k=0;k<xmax;++k)
				for(long int m=0;m<zmax;++m){
					as(k,0,m)=0;
					bs(k,0,m)=2.0*betay(k,0,m)+1.0;
					cs(k,0,m)=-betay(k,0,m)*(1.0+y(k,0,m));

					for(long int i = 1; i < ylast ; i++){
						as(k,i,m)=-betay(k,i,m)*(1.0-y(k,i,m));
						bs(k,i,m)=2.0*betay(k,i,m)+1.0;
						cs(k,i,m)=-betay(k,i,m)*(1.0+y(k,i,m));
					}

					as(k,ylast,m)=-betay(k,ylast,m)*(1.0-y(k,ylast,m));
					bs(k,ylast,m)=2.0*betay(k,ylast,m)+1.0;
					cs(k,ylast,m)=0;

					//Thomas Algorithm

					for (long int i = 1; i <= ylast; i++){
						std::complex<complextype> n = as(k,i,m)/bs(k,i-1,m);
						bs(k,i,m)-=n*cs(k,i-1,m);
						initial(k,i,m)-=n*initial(k,i-1,m);
					}

					initial(k,ylast,m)/=bs(k,ylast,m);
					for (long int i = (ylast - 1); i >= 0; --i)
						initial(k,i,m)=(initial(k,i,m)-cs(k,i,m)*initial(k,i+1,m))/bs(k,i,m);

				}


			/* z direction*/

			copy=initial;

			#pragma omp parallel for
			for(long int k=0;k<xmax;++k)
				for(long int l=0;l<ymax;++l){
					//Multiply vector
					initial(k,l,0)=betaz(k,l,0)*copy(k,l,1)*(1.0+z(k,l,0))+(1.0-2.0*betaz(k,l,0))*copy(k,l,0);
					initial(k,l,zlast)=betaz(k,l,zlast)*copy(k,l,zlast-1)*(1.0-z(k,l,zlast))+(1.0-2.0*betaz(k,l,zlast))*copy(k,l,zlast);

					for(long int m = 1; m < zlast; m++){
							initial(k,l,m)=betaz(k,l,m)*copy(k,l,m-1)*(1.0-z(k,l,m))+(1.0-2.0*betaz(k,l,m))*copy(k,l,m)+betaz(k,l,m)*copy(k,l,m+1)*(1.0+z(k,l,m));
					}
				}


			#pragma omp parallel for
			for(long int k=0;k<xmax;++k)
				for(long int l=0;l<ymax;++l){
					as(k,l,0)=0;
					bs(k,l,0)=2.0*betaz(k,l,0)+1.0;
					cs(k,l,0)=-betaz(k,l,0)*(1.0+z(k,l,0));

					for(long int i = 1; i < zlast ; i++){
						as(k,l,i)=-betaz(k,l,i)*(1.0-z(k,l,i));
						bs(k,l,i)=2.0*betaz(k,l,i)+1.0;
						cs(k,l,i)=-betaz(k,l,i)*(1.0+z(k,l,i));
					}

					as(k,l,zlast)=-betaz(k,l,zlast)*(1.0-z(k,l,zlast));
					bs(k,l,zlast)=2.0*betaz(k,l,zlast)+1.0;
					cs(k,l,zlast)=0;

					//Thomas Algorithm

					for (long int i = 1; i <= zlast; i++){
						std::complex<complextype> n = as(k,l,i)/bs(k,l,i-1);
						bs(k,l,i)-=n*cs(k,l,i-1);
						initial(k,l,i)-=n*initial(k,l,i-1);
					}

					initial(k,l,zlast)/=bs(k,l,zlast);
					for (long int i = (zlast - 1); i >= 0; --i)
						initial(k,l,i)=(initial(k,l,i)-cs(k,l,i)*initial(k,l,i+1))/bs(k,l,i);

				}
		}



		initial*=potpropag;
		std::complex<complextype> tonorm=initial.Modulus();
		initial*=1.0/tonorm;

	}

	return initial;
}


template <typename timetype, typename complextype>
template<typename T, typename gridtype>
DiscreteFunction3D<std::complex<complextype>,gridtype>& WaveTimeEvolution3D<timetype,complextype>::groundlevel_finite(timetype time, long int timesteps, DiscreteFunction3D<std::complex<complextype>,gridtype>& initial, DiscreteFunction3D<T,gridtype>& potential, DiscreteFunction3D<T,gridtype>& effectivemass){
	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	initial/=initial.Modulus();

	DiscreteFunction3D<std::complex<complextype>,gridtype> potpropag = potentialpropagator(methodtimestep,potential);

	timetype adi_methodtimestep = methodtimestep/3.0;

	gridtype gridstepx = initial.getGrid().getxIncrement();
	gridtype gridstepy = initial.getGrid().getyIncrement();
	gridtype gridstepz = initial.getGrid().getzIncrement();


	long int xmax=initial.getGrid().getSizeX();
	long int ymax=initial.getGrid().getSizeY();
	long int zmax=initial.getGrid().getSizeZ();

	long int xlast=xmax-1;
	long int ylast=ymax-1;
	long int zlast=zmax-1;


	DiscreteFunction3D<std::complex<complextype>,gridtype> as(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> bs(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> cs(initial.getGrid());

	DiscreteFunction3D<std::complex<complextype>,gridtype> betax(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> betay(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> betaz(initial.getGrid());

	#pragma omp parallel for
	for(long int i = 0; i< betax.getGrid().getSize();i++){
		betax.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*adi_methodtimestep/(4.0*effectivemass.getdata()[i]*gridstepx*gridstepx);
	}

	#pragma omp parallel for
	for(long int i = 0; i< betay.getGrid().getSize();i++){
		betay.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*adi_methodtimestep/(4.0*effectivemass.getdata()[i]*gridstepy*gridstepy);
	}

	#pragma omp parallel for
	for(long int i = 0; i< betaz.getGrid().getSize();i++){
		betaz.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*adi_methodtimestep/(4.0*effectivemass.getdata()[i]*gridstepz*gridstepz);
	}


	DiscreteFunction3D<std::complex<complextype>,gridtype> x(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> y(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> z(initial.getGrid());

	gridtype dx = gridstepx;
	gridtype dy = gridstepy;
	gridtype dz = gridstepz;

	for(long int l=0;l<ymax;++l){
		for(long int m=0;m<zmax;++m){
			//edge: order h and h²
			x(0,l,m) = (effectivemass(0,l,m)*dx/2.0)*(1.0/effectivemass(1,l,m)-1.0/effectivemass(0,l,m))/dx;
			x(1,l,m) = (effectivemass(1,l,m)*dx/2.0)*(1.0/effectivemass(2,l,m)-1.0/effectivemass(0,l,m))/(2*dx);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int k = 2; k < (xmax-2); k++){
				x(k,l,m)=(effectivemass(k,l,m)*dx/2.0)*((1.0/effectivemass(k+1,l,m))*8.0-(1.0/effectivemass(k+2,l,m))-(1.0/effectivemass(k-1,l,m))*8.0+(1.0/effectivemass(k-2,l,m)))/(12.0*dx);
			}

			//edge: order h and h²
			x(xmax-2,l,m)= (effectivemass(xmax-2,l,m)*dx/2.0)*(1.0/effectivemass(xmax-1,l,m)-1.0/effectivemass(xmax-3,l,m))/(2.0*dx);
			x(xmax-1,l,m)= (effectivemass(xmax-1,l,m)*dx/2.0)*(1.0/effectivemass(xmax-1,l,m)-1.0/effectivemass(xmax-2,l,m))/dx;
		}
	}

	for(long int k=0;k<xmax;++k){
		for(long int m=0;m<zmax;++m){
			//edge: order h and h²
			y(k,0,m) = (effectivemass(k,0,m)*dy/2.0)*(1.0/effectivemass(k,1,m)-1.0/effectivemass(k,0,m))/dy;
			y(k,1,m) = (effectivemass(k,1,m)*dy/2.0)*(1.0/effectivemass(k,2,m)-1.0/effectivemass(k,0,m))/(2*dy);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int l = 2; l < (ymax-2); l++){
				y(k,l,m)=(effectivemass(k,l,m)*dy/2.0)*((1.0/effectivemass(k,l+1,m))*8.0-(1.0/effectivemass(k,l+2,m))-(1.0/effectivemass(k,l-1,m))*8.0+(1.0/effectivemass(k,l-2,m)))/(12.0*dy);
			}

			//edge: order h and h²
			y(k,ymax-2,m)= (effectivemass(k,ymax-2,m)*dy/2.0)*(1.0/effectivemass(k,ymax-1,m)-1.0/effectivemass(k,ymax-3,m))/(2.0*dy);
			y(k,ymax-1,m)= (effectivemass(k,ymax-1,m)*dy/2.0)*(1.0/effectivemass(k,ymax-1,m)-1.0/effectivemass(k,ymax-2,m))/dy;
		}
	}

	for(long int k=0;k<xmax;++k){
		for(long int l=0;l<ymax;++l){
			//edge: order h and h²
			z(k,l,0) = (effectivemass(k,l,0)*dz/2.0)*(1.0/effectivemass(k,l,1)-1.0/effectivemass(k,l,0))/dz;
			z(k,l,1) = (effectivemass(k,l,1)*dz/2.0)*(1.0/effectivemass(k,l,2)-1.0/effectivemass(k,l,0))/(2*dz);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int m = 2; m < (zmax-2); m++){
				z(k,l,m)=(effectivemass(k,l,m)*dz/2.0)*((1.0/effectivemass(k,l,m+1))*8.0-(1.0/effectivemass(k,l,m+2))-(1.0/effectivemass(k,l,m-1))*8.0+(1.0/effectivemass(k,l,m-2)))/(12.0*dz);
			}

			//edge: order h and h²
			z(k,l,zmax-2)= (effectivemass(k,l,zmax-2)*dz/2.0)*(1.0/effectivemass(k,l,zmax-1)-1.0/effectivemass(k,l,zmax-3))/(2.0*dz);
			z(k,l,zmax-1)= (effectivemass(k,l,zmax-1)*dz/2.0)*(1.0/effectivemass(k,l,zmax-1)-1.0/effectivemass(k,l,zmax-2))/dz;
		}
	}



	for(long int i = 0 ; i< timesteps; i++){
		//std::cout << "step: " << i  << " ..." << std::endl;
		initial*=potpropag;

		{
			/* x direction*/

			DiscreteFunction3D<std::complex<complextype>,gridtype> copy(initial);

			#pragma omp parallel for
			for(long int l=0;l<ymax;++l)
				for(long int m=0;m<zmax;++m){
					//Multiply vector
					initial(0,l,m)=betax(0,l,m)*copy(1,l,m)*(1.0+x(0,l,m))+(1.0-2.0*betax(0,l,m))*copy(0,l,m);
					initial(xlast,l,m)=betax(xlast,l,m)*copy(xlast-1,l,m)*(1.0-x(xlast,l,m))+(1.0-2.0*betax(xlast,l,m))*copy(xlast,l,m);

					for(long int k = 1; k < xlast; k++){
							initial(k,l,m)=betax(k,l,m)*copy(k-1,l,m)*(1.0-x(k,l,m))+(1.0-2.0*betax(k,l,m))*copy(k,l,m)+betax(k,l,m)*copy(k+1,l,m)*(1.0+x(k,l,m));
					}
			}


			#pragma omp parallel for
			for(long int l=0;l<ymax;++l)
				for(long int m=0;m<zmax;++m){
					as(0,l,m)=0;
					bs(0,l,m)=2.0*betax(0,l,m)+1.0;
					cs(0,l,m)=-betax(0,l,m)*(1.0+x(0,l,m));

					for(long int i = 1; i < xlast ; i++){
						as(i,l,m)=-betax(i,l,m)*(1.0-x(i,l,m));
						bs(i,l,m)=2.0*betax(i,l,m)+1.0;
						cs(i,l,m)=-betax(i,l,m)*(1.0+x(i,l,m));
					}

					as(xlast,l,m)=-betax(xlast,l,m)*(1.0-x(xlast,l,m));
					bs(xlast,l,m)=2.0*betax(xlast,l,m)+1.0;
					cs(xlast,l,m)=0;

					//Thomas Algorithm

					for (long int i = 1; i <= xlast; i++){
						std::complex<complextype> n = as(i,l,m)/bs(i-1,l,m);
						bs(i,l,m)-=n*cs(i-1,l,m);
						initial(i,l,m)-=n*initial(i-1,l,m);
					}

					initial(xlast,l,m)/=bs(xlast,l,m);
					for (long int i = (xlast - 1); i >= 0; --i)
						initial(i,l,m)=(initial(i,l,m)-cs(i,l,m)*initial(i+1,l,m))/bs(i,l,m);

				}

			/* y direction*/

			copy=initial;

			#pragma omp parallel for
			for(long int k=0;k<xmax;++k)
				for(long int m=0;m<zmax;++m){
					//Multiply vector
					initial(k,0,m)=betay(k,0,m)*copy(k,1,m)*(1.0+y(k,0,m))+(1.0-2.0*betay(k,0,m))*copy(k,0,m);
					initial(k,ylast,m)=betay(k,ylast,m)*copy(k,ylast-1,m)*(1.0-y(k,ylast,m))+(1.0-2.0*betay(k,ylast,m))*copy(k,ylast,m);

					for(long int l = 1; l < ylast; l++){
							initial(k,l,m)=betay(k,l,m)*copy(k,l-1,m)*(1.0-y(k,l,m))+(1.0-2.0*betay(k,l,m))*copy(k,l,m)+betay(k,l,m)*copy(k,l+1,m)*(1.0+y(k,l,m));
					}
				}


			#pragma omp parallel for
			for(long int k=0;k<xmax;++k)
				for(long int m=0;m<zmax;++m){
					as(k,0,m)=0;
					bs(k,0,m)=2.0*betay(k,0,m)+1.0;
					cs(k,0,m)=-betay(k,0,m)*(1.0+y(k,0,m));

					for(long int i = 1; i < ylast ; i++){
						as(k,i,m)=-betay(k,i,m)*(1.0-y(k,i,m));
						bs(k,i,m)=2.0*betay(k,i,m)+1.0;
						cs(k,i,m)=-betay(k,i,m)*(1.0+y(k,i,m));
					}

					as(k,ylast,m)=-betay(k,ylast,m)*(1.0-y(k,ylast,m));
					bs(k,ylast,m)=2.0*betay(k,ylast,m)+1.0;
					cs(k,ylast,m)=0;

					//Thomas Algorithm

					for (long int i = 1; i <= ylast; i++){
						std::complex<complextype> n = as(k,i,m)/bs(k,i-1,m);
						bs(k,i,m)-=n*cs(k,i-1,m);
						initial(k,i,m)-=n*initial(k,i-1,m);
					}

					initial(k,ylast,m)/=bs(k,ylast,m);
					for (long int i = (ylast - 1); i >= 0; --i)
						initial(k,i,m)=(initial(k,i,m)-cs(k,i,m)*initial(k,i+1,m))/bs(k,i,m);

				}


			/* z direction*/

			copy=initial;

			#pragma omp parallel for
			for(long int k=0;k<xmax;++k)
				for(long int l=0;l<ymax;++l){
					//Multiply vector
					initial(k,l,0)=betaz(k,l,0)*copy(k,l,1)*(1.0+z(k,l,0))+(1.0-2.0*betaz(k,l,0))*copy(k,l,0);
					initial(k,l,zlast)=betaz(k,l,zlast)*copy(k,l,zlast-1)*(1.0-z(k,l,zlast))+(1.0-2.0*betaz(k,l,zlast))*copy(k,l,zlast);

					for(long int m = 1; m < zlast; m++){
							initial(k,l,m)=betaz(k,l,m)*copy(k,l,m-1)*(1.0-z(k,l,m))+(1.0-2.0*betaz(k,l,m))*copy(k,l,m)+betaz(k,l,m)*copy(k,l,m+1)*(1.0+z(k,l,m));
					}
				}


			#pragma omp parallel for
			for(long int k=0;k<xmax;++k)
				for(long int l=0;l<ymax;++l){
					as(k,l,0)=0;
					bs(k,l,0)=2.0*betaz(k,l,0)+1.0;
					cs(k,l,0)=-betaz(k,l,0)*(1.0+z(k,l,0));

					for(long int i = 1; i < zlast ; i++){
						as(k,l,i)=-betaz(k,l,i)*(1.0-z(k,l,i));
						bs(k,l,i)=2.0*betaz(k,l,i)+1.0;
						cs(k,l,i)=-betaz(k,l,i)*(1.0+z(k,l,i));
					}

					as(k,l,zlast)=-betaz(k,l,zlast)*(1.0-z(k,l,zlast));
					bs(k,l,zlast)=2.0*betaz(k,l,zlast)+1.0;
					cs(k,l,zlast)=0;

					//Thomas Algorithm

					for (long int i = 1; i <= zlast; i++){
						std::complex<complextype> n = as(k,l,i)/bs(k,l,i-1);
						bs(k,l,i)-=n*cs(k,l,i-1);
						initial(k,l,i)-=n*initial(k,l,i-1);
					}

					initial(k,l,zlast)/=bs(k,l,zlast);
					for (long int i = (zlast - 1); i >= 0; --i)
						initial(k,l,i)=(initial(k,l,i)-cs(k,l,i)*initial(k,l,i+1))/bs(k,l,i);

				}
		}



		initial*=potpropag;
		std::complex<complextype> tonorm=initial.Modulus();
		initial*=1.0/tonorm;

	}

	return initial;
}




template <typename timetype, typename complextype>
template<typename T, typename gridtype>
DiscreteFunction3D<std::complex<complextype>,gridtype>& WaveTimeEvolution3D<timetype,complextype>::newlevel_finite(timetype time, long int timesteps, DiscreteFunction3D<std::complex<complextype>,gridtype>& initial, DiscreteFunction3D<T,gridtype>& potential,std::vector<std::shared_ptr<DiscreteFunction3D<std::complex<complextype>,gridtype>>> otherslevels, DiscreteFunction3D<T,gridtype>& effectivemass){
	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	initial/=initial.Modulus();

	DiscreteFunction3D<std::complex<complextype>,gridtype> potpropag = potentialpropagator(methodtimestep,potential);

	timetype adi_methodtimestep = methodtimestep/3.0;

	gridtype gridstepx = initial.getGrid().getxIncrement();
	gridtype gridstepy = initial.getGrid().getyIncrement();
	gridtype gridstepz = initial.getGrid().getzIncrement();


	long int xmax=initial.getGrid().getSizeX();
	long int ymax=initial.getGrid().getSizeY();
	long int zmax=initial.getGrid().getSizeZ();

	long int xlast=xmax-1;
	long int ylast=ymax-1;
	long int zlast=zmax-1;


	DiscreteFunction3D<std::complex<complextype>,gridtype> as(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> bs(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> cs(initial.getGrid());

	DiscreteFunction3D<std::complex<complextype>,gridtype> betax(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> betay(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> betaz(initial.getGrid());

	#pragma omp parallel for
	for(long int i = 0; i< betax.getGrid().getSize();i++){
		betax.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*adi_methodtimestep/(4.0*effectivemass.getdata()[i]*gridstepx*gridstepx);
	}

	#pragma omp parallel for
	for(long int i = 0; i< betay.getGrid().getSize();i++){
		betay.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*adi_methodtimestep/(4.0*effectivemass.getdata()[i]*gridstepy*gridstepy);
	}

	#pragma omp parallel for
	for(long int i = 0; i< betaz.getGrid().getSize();i++){
		betaz.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*adi_methodtimestep/(4.0*effectivemass.getdata()[i]*gridstepz*gridstepz);
	}


	DiscreteFunction3D<std::complex<complextype>,gridtype> x(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> y(initial.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> z(initial.getGrid());

	gridtype dx = gridstepx;
	gridtype dy = gridstepy;
	gridtype dz = gridstepz;

	for(long int l=0;l<ymax;++l){
		for(long int m=0;m<zmax;++m){
			//edge: order h and h²
			x(0,l,m) = (effectivemass(0,l,m)*dx/2.0)*(1.0/effectivemass(1,l,m)-1.0/effectivemass(0,l,m))/dx;
			x(1,l,m) = (effectivemass(1,l,m)*dx/2.0)*(1.0/effectivemass(2,l,m)-1.0/effectivemass(0,l,m))/(2*dx);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int k = 2; k < (xmax-2); k++){
				x(k,l,m)=(effectivemass(k,l,m)*dx/2.0)*((1.0/effectivemass(k+1,l,m))*8.0-(1.0/effectivemass(k+2,l,m))-(1.0/effectivemass(k-1,l,m))*8.0+(1.0/effectivemass(k-2,l,m)))/(12.0*dx);
			}

			//edge: order h and h²
			x(xmax-2,l,m)= (effectivemass(xmax-2,l,m)*dx/2.0)*(1.0/effectivemass(xmax-1,l,m)-1.0/effectivemass(xmax-3,l,m))/(2.0*dx);
			x(xmax-1,l,m)= (effectivemass(xmax-1,l,m)*dx/2.0)*(1.0/effectivemass(xmax-1,l,m)-1.0/effectivemass(xmax-2,l,m))/dx;
		}
	}

	for(long int k=0;k<xmax;++k){
		for(long int m=0;m<zmax;++m){
			//edge: order h and h²
			y(k,0,m) = (effectivemass(k,0,m)*dy/2.0)*(1.0/effectivemass(k,1,m)-1.0/effectivemass(k,0,m))/dy;
			y(k,1,m) = (effectivemass(k,1,m)*dy/2.0)*(1.0/effectivemass(k,2,m)-1.0/effectivemass(k,0,m))/(2*dy);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int l = 2; l < (ymax-2); l++){
				y(k,l,m)=(effectivemass(k,l,m)*dy/2.0)*((1.0/effectivemass(k,l+1,m))*8.0-(1.0/effectivemass(k,l+2,m))-(1.0/effectivemass(k,l-1,m))*8.0+(1.0/effectivemass(k,l-2,m)))/(12.0*dy);
			}

			//edge: order h and h²
			y(k,ymax-2,m)= (effectivemass(k,ymax-2,m)*dy/2.0)*(1.0/effectivemass(k,ymax-1,m)-1.0/effectivemass(k,ymax-3,m))/(2.0*dy);
			y(k,ymax-1,m)= (effectivemass(k,ymax-1,m)*dy/2.0)*(1.0/effectivemass(k,ymax-1,m)-1.0/effectivemass(k,ymax-2,m))/dy;
		}
	}

	for(long int k=0;k<xmax;++k){
		for(long int l=0;l<ymax;++l){
			//edge: order h and h²
			z(k,l,0) = (effectivemass(k,l,0)*dz/2.0)*(1.0/effectivemass(k,l,1)-1.0/effectivemass(k,l,0))/dz;
			z(k,l,1) = (effectivemass(k,l,1)*dz/2.0)*(1.0/effectivemass(k,l,2)-1.0/effectivemass(k,l,0))/(2*dz);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int m = 2; m < (zmax-2); m++){
				z(k,l,m)=(effectivemass(k,l,m)*dz/2.0)*((1.0/effectivemass(k,l,m+1))*8.0-(1.0/effectivemass(k,l,m+2))-(1.0/effectivemass(k,l,m-1))*8.0+(1.0/effectivemass(k,l,m-2)))/(12.0*dz);
			}

			//edge: order h and h²
			z(k,l,zmax-2)= (effectivemass(k,l,zmax-2)*dz/2.0)*(1.0/effectivemass(k,l,zmax-1)-1.0/effectivemass(k,l,zmax-3))/(2.0*dz);
			z(k,l,zmax-1)= (effectivemass(k,l,zmax-1)*dz/2.0)*(1.0/effectivemass(k,l,zmax-1)-1.0/effectivemass(k,l,zmax-2))/dz;
		}
	}



	for(long int i = 0 ; i< timesteps; i++){
		//std::cout << "step: " << i  << " ..." << std::endl;

		for(auto k : otherslevels){
			std::complex<complextype> dotprod = (*k).Dot(initial);
			initial-=((*k)*dotprod);
		}

		initial*=potpropag;

		{
			/* x direction*/

			DiscreteFunction3D<std::complex<complextype>,gridtype> copy(initial);

			#pragma omp parallel for
			for(long int l=0;l<ymax;++l)
				for(long int m=0;m<zmax;++m){
					//Multiply vector
					initial(0,l,m)=betax(0,l,m)*copy(1,l,m)*(1.0+x(0,l,m))+(1.0-2.0*betax(0,l,m))*copy(0,l,m);
					initial(xlast,l,m)=betax(xlast,l,m)*copy(xlast-1,l,m)*(1.0-x(xlast,l,m))+(1.0-2.0*betax(xlast,l,m))*copy(xlast,l,m);

					for(long int k = 1; k < xlast; k++){
							initial(k,l,m)=betax(k,l,m)*copy(k-1,l,m)*(1.0-x(k,l,m))+(1.0-2.0*betax(k,l,m))*copy(k,l,m)+betax(k,l,m)*copy(k+1,l,m)*(1.0+x(k,l,m));
					}
			}


			#pragma omp parallel for
			for(long int l=0;l<ymax;++l)
				for(long int m=0;m<zmax;++m){
					as(0,l,m)=0;
					bs(0,l,m)=2.0*betax(0,l,m)+1.0;
					cs(0,l,m)=-betax(0,l,m)*(1.0+x(0,l,m));

					for(long int i = 1; i < xlast ; i++){
						as(i,l,m)=-betax(i,l,m)*(1.0-x(i,l,m));
						bs(i,l,m)=2.0*betax(i,l,m)+1.0;
						cs(i,l,m)=-betax(i,l,m)*(1.0+x(i,l,m));
					}

					as(xlast,l,m)=-betax(xlast,l,m)*(1.0-x(xlast,l,m));
					bs(xlast,l,m)=2.0*betax(xlast,l,m)+1.0;
					cs(xlast,l,m)=0;

					//Thomas Algorithm

					for (long int i = 1; i <= xlast; i++){
						std::complex<complextype> n = as(i,l,m)/bs(i-1,l,m);
						bs(i,l,m)-=n*cs(i-1,l,m);
						initial(i,l,m)-=n*initial(i-1,l,m);
					}

					initial(xlast,l,m)/=bs(xlast,l,m);
					for (long int i = (xlast - 1); i >= 0; --i)
						initial(i,l,m)=(initial(i,l,m)-cs(i,l,m)*initial(i+1,l,m))/bs(i,l,m);

				}

			/* y direction*/

			copy=initial;

			#pragma omp parallel for
			for(long int k=0;k<xmax;++k)
				for(long int m=0;m<zmax;++m){
					//Multiply vector
					initial(k,0,m)=betay(k,0,m)*copy(k,1,m)*(1.0+y(k,0,m))+(1.0-2.0*betay(k,0,m))*copy(k,0,m);
					initial(k,ylast,m)=betay(k,ylast,m)*copy(k,ylast-1,m)*(1.0-y(k,ylast,m))+(1.0-2.0*betay(k,ylast,m))*copy(k,ylast,m);

					for(long int l = 1; l < ylast; l++){
							initial(k,l,m)=betay(k,l,m)*copy(k,l-1,m)*(1.0-y(k,l,m))+(1.0-2.0*betay(k,l,m))*copy(k,l,m)+betay(k,l,m)*copy(k,l+1,m)*(1.0+y(k,l,m));
					}
				}


			#pragma omp parallel for
			for(long int k=0;k<xmax;++k)
				for(long int m=0;m<zmax;++m){
					as(k,0,m)=0;
					bs(k,0,m)=2.0*betay(k,0,m)+1.0;
					cs(k,0,m)=-betay(k,0,m)*(1.0+y(k,0,m));

					for(long int i = 1; i < ylast ; i++){
						as(k,i,m)=-betay(k,i,m)*(1.0-y(k,i,m));
						bs(k,i,m)=2.0*betay(k,i,m)+1.0;
						cs(k,i,m)=-betay(k,i,m)*(1.0+y(k,i,m));
					}

					as(k,ylast,m)=-betay(k,ylast,m)*(1.0-y(k,ylast,m));
					bs(k,ylast,m)=2.0*betay(k,ylast,m)+1.0;
					cs(k,ylast,m)=0;

					//Thomas Algorithm

					for (long int i = 1; i <= ylast; i++){
						std::complex<complextype> n = as(k,i,m)/bs(k,i-1,m);
						bs(k,i,m)-=n*cs(k,i-1,m);
						initial(k,i,m)-=n*initial(k,i-1,m);
					}

					initial(k,ylast,m)/=bs(k,ylast,m);
					for (long int i = (ylast - 1); i >= 0; --i)
						initial(k,i,m)=(initial(k,i,m)-cs(k,i,m)*initial(k,i+1,m))/bs(k,i,m);

				}


			/* z direction*/

			copy=initial;

			#pragma omp parallel for
			for(long int k=0;k<xmax;++k)
				for(long int l=0;l<ymax;++l){
					//Multiply vector
					initial(k,l,0)=betaz(k,l,0)*copy(k,l,1)*(1.0+z(k,l,0))+(1.0-2.0*betaz(k,l,0))*copy(k,l,0);
					initial(k,l,zlast)=betaz(k,l,zlast)*copy(k,l,zlast-1)*(1.0-z(k,l,zlast))+(1.0-2.0*betaz(k,l,zlast))*copy(k,l,zlast);

					for(long int m = 1; m < zlast; m++){
							initial(k,l,m)=betaz(k,l,m)*copy(k,l,m-1)*(1.0-z(k,l,m))+(1.0-2.0*betaz(k,l,m))*copy(k,l,m)+betaz(k,l,m)*copy(k,l,m+1)*(1.0+z(k,l,m));
					}
				}


			#pragma omp parallel for
			for(long int k=0;k<xmax;++k)
				for(long int l=0;l<ymax;++l){
					as(k,l,0)=0;
					bs(k,l,0)=2.0*betaz(k,l,0)+1.0;
					cs(k,l,0)=-betaz(k,l,0)*(1.0+z(k,l,0));

					for(long int i = 1; i < zlast ; i++){
						as(k,l,i)=-betaz(k,l,i)*(1.0-z(k,l,i));
						bs(k,l,i)=2.0*betaz(k,l,i)+1.0;
						cs(k,l,i)=-betaz(k,l,i)*(1.0+z(k,l,i));
					}

					as(k,l,zlast)=-betaz(k,l,zlast)*(1.0-z(k,l,zlast));
					bs(k,l,zlast)=2.0*betaz(k,l,zlast)+1.0;
					cs(k,l,zlast)=0;

					//Thomas Algorithm

					for (long int i = 1; i <= zlast; i++){
						std::complex<complextype> n = as(k,l,i)/bs(k,l,i-1);
						bs(k,l,i)-=n*cs(k,l,i-1);
						initial(k,l,i)-=n*initial(k,l,i-1);
					}

					initial(k,l,zlast)/=bs(k,l,zlast);
					for (long int i = (zlast - 1); i >= 0; --i)
						initial(k,l,i)=(initial(k,l,i)-cs(k,l,i)*initial(k,l,i+1))/bs(k,l,i);

				}
		}



		initial*=potpropag;
		std::complex<complextype> tonorm=initial.Modulus();
		initial*=1.0/tonorm;

	}

	return initial;
}

template <typename timetype, typename complextype>
template<typename T, typename gridtype>
DiscreteFunction3D<std::complex<complextype>,gridtype>& WaveTimeEvolution3D<timetype,complextype>::newlevel_periodicFFT(timetype time, long int timesteps, DiscreteFunction3D<std::complex<complextype>,gridtype>& initial, DiscreteFunction3D<T,gridtype>& potential,std::vector<std::shared_ptr<DiscreteFunction3D<std::complex<complextype>,gridtype>>> otherslevels, T effectivemass){
	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	initial/=initial.Modulus();


	DiscreteFunction3D<std::complex<complextype>,gridtype> potpropag = potentialpropagator(methodtimestep,potential);
	DiscreteFunction3D<std::complex<complextype>,gridtype> bloch = fftmulti(methodtimestep,initial.getGrid(),effectivemass);

	for(long int i = 0 ; i< timesteps; i++){
		//std::cout << "step: " << i  << " ..." << std::endl;
		for(auto k : otherslevels){
			std::complex<complextype> dotprod = (*k).Dot(initial);
			initial-=((*k)*dotprod);
		}
		initial*=potpropag;
		initial.FFT();
		initial*=bloch;
		initial.iFFT();
		initial*=potpropag;
		std::complex<complextype> tonorm=initial.Modulus();
		initial*=1.0/tonorm;
	}

	return initial;
}

template <typename timetype, typename complextype>
template<typename T, typename gridtype>
DiscreteFunction3D<std::complex<complextype>,gridtype> WaveTimeEvolution3D<timetype,complextype>::fftmulti(timetype methodtimestep,Grid3D<gridtype> grid,T effmass){
	DiscreteFunction3D<std::complex<complextype>,gridtype> result(grid);
	long int sizehalfx = grid.getSizeX()/2;
	long int sizehalfy = grid.getSizeY()/2;
	long int sizehalfz = grid.getSizeZ()/2;

	for(long int k = 0; k < sizehalfx; k++){
		for(long int l = 0; l < sizehalfy; l++){
			for(long int m = 0; m < sizehalfz; m++){
				result(k,l,m)=std::exp( std::complex<complextype>(0.0,-1.0)*((pow(2.0*Constant::pi*k/(grid.getxSuperiorLimit()-grid.getxInferiorLimit()),2)+pow(2.0*Constant::pi*l/(grid.getySuperiorLimit()-grid.getyInferiorLimit()),2)+pow(2.0*Constant::pi*m/(grid.getzSuperiorLimit()-grid.getzInferiorLimit()),2)))*methodtimestep/(2*Constant::hbar*effmass) );
			}
			for(long int m = sizehalfz; m < grid.getSizeZ(); m++){
				result(k,l,m)=std::exp( std::complex<complextype>(0.0,-1.0)*((pow(2.0*Constant::pi*k/(grid.getxSuperiorLimit()-grid.getxInferiorLimit()),2)+pow(2.0*Constant::pi*l/(grid.getySuperiorLimit()-grid.getyInferiorLimit()),2)+pow(2.0*Constant::pi*(m-2*sizehalfz)/(grid.getzSuperiorLimit()-grid.getzInferiorLimit()),2)))*methodtimestep/(2*Constant::hbar*effmass) );
			}
		}
		for(long int l = sizehalfy; l < grid.getSizeY(); l++){
			for(long int m = 0; m < sizehalfz; m++){
				result(k,l,m)=std::exp( std::complex<complextype>(0.0,-1.0)*((pow(2.0*Constant::pi*k/(grid.getxSuperiorLimit()-grid.getxInferiorLimit()),2)+pow(2.0*Constant::pi*(l-2*sizehalfy)/(grid.getySuperiorLimit()-grid.getyInferiorLimit()),2)+pow(2.0*Constant::pi*m/(grid.getzSuperiorLimit()-grid.getzInferiorLimit()),2)))*methodtimestep/(2*Constant::hbar*effmass) );
			}
			for(long int m = sizehalfz; m < grid.getSizeZ(); m++){
				result(k,l,m)=std::exp( std::complex<complextype>(0.0,-1.0)*((pow(2.0*Constant::pi*k/(grid.getxSuperiorLimit()-grid.getxInferiorLimit()),2)+pow(2.0*Constant::pi*(l-2*sizehalfy)/(grid.getySuperiorLimit()-grid.getyInferiorLimit()),2)+pow(2.0*Constant::pi*(m-2*sizehalfz)/(grid.getzSuperiorLimit()-grid.getzInferiorLimit()),2)))*methodtimestep/(2*Constant::hbar*effmass) );
			}
		}
	}
	for(long int k = sizehalfx; k < grid.getSizeX(); k++){
		for(long int l = 0; l < sizehalfy; l++){
			for(long int m = 0; m < sizehalfz; m++){
				result(k,l,m)=std::exp( std::complex<complextype>(0.0,-1.0)*((pow(2.0*Constant::pi*(k-2*sizehalfx)/(grid.getxSuperiorLimit()-grid.getxInferiorLimit()),2)+pow(2.0*Constant::pi*l/(grid.getySuperiorLimit()-grid.getyInferiorLimit()),2)+pow(2.0*Constant::pi*m/(grid.getzSuperiorLimit()-grid.getzInferiorLimit()),2)))*methodtimestep/(2*Constant::hbar*effmass) );
			}
			for(long int m = sizehalfz; m < grid.getSizeZ(); m++){
				result(k,l,m)=std::exp( std::complex<complextype>(0.0,-1.0)*((pow(2.0*Constant::pi*(m-2*sizehalfz)/(grid.getzSuperiorLimit()-grid.getzInferiorLimit()),2)+pow(2.0*Constant::pi*(k-2*sizehalfx)/(grid.getxSuperiorLimit()-grid.getxInferiorLimit()),2)+pow(2.0*Constant::pi*l/(grid.getySuperiorLimit()-grid.getyInferiorLimit()),2)))*methodtimestep/(2*Constant::hbar*effmass) );
			}
		}
		for(long int l = sizehalfy; l < grid.getSizeY(); l++){
			for(long int m = 0; m < sizehalfz; m++){
				result(k,l,m)=std::exp( std::complex<complextype>(0.0,-1.0)*((pow(2.0*Constant::pi*(l-2*sizehalfy)/(grid.getySuperiorLimit()-grid.getyInferiorLimit()),2)+pow(2.0*Constant::pi*(k-2*sizehalfx)/(grid.getxSuperiorLimit()-grid.getxInferiorLimit()),2)+pow(2.0*Constant::pi*m/(grid.getzSuperiorLimit()-grid.getzInferiorLimit()),2)))*methodtimestep/(2*Constant::hbar*effmass) );
			}
			for(long int m = sizehalfz; m < grid.getSizeZ(); m++){
				result(k,l,m)=std::exp( std::complex<complextype>(0.0,-1.0)*((pow(2.0*Constant::pi*(l-2*sizehalfy)/(grid.getySuperiorLimit()-grid.getyInferiorLimit()),2)+pow(2.0*Constant::pi*(m-2*sizehalfz)/(grid.getzSuperiorLimit()-grid.getzInferiorLimit()),2)+pow(2.0*Constant::pi*(k-2*sizehalfx)/(grid.getxSuperiorLimit()-grid.getxInferiorLimit()),2)))*methodtimestep/(2*Constant::hbar*effmass) );
			}
		}
	}
	return result;
}
/*
template <typename timetype, typename complextype>
template<typename T, typename gridtype>
DiscreteFunction3D<std::complex<complextype>,gridtype>& WaveTimeEvolution3D<timetype,complextype>::kinect_propag_evolution_finite(timetype methodtimestep, DiscreteFunction3D<std::complex<complextype>,gridtype>& wavefunction, DiscreteFunction3D<T,gridtype>& effectivemass){
	//defining beta...

	methodtimestep/=3;

	gridtype gridstepx = wavefunction.getGrid().getxIncrement();
	gridtype gridstepy = wavefunction.getGrid().getyIncrement();
	gridtype gridstepz = wavefunction.getGrid().getzIncrement();


	long int xmax=wavefunction.getGrid().getSizeX();
	long int ymax=wavefunction.getGrid().getSizeY();
	long int zmax=wavefunction.getGrid().getSizeZ();

	long int xlast=xmax-1;
	long int ylast=ymax-1;
	long int zlast=zmax-1;

	DiscreteFunction3D<std::complex<complextype>,gridtype> betax(wavefunction.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> betay(wavefunction.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> betaz(wavefunction.getGrid());

	DiscreteFunction3D<std::complex<complextype>,gridtype> as(wavefunction.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> bs(wavefunction.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> cs(wavefunction.getGrid());

	#pragma omp parallel for
	for(long int i = 0; i< betax.getGrid().getSize();i++){
		betax.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*methodtimestep/(4.0*gridstepx*gridstepx);
	}

	#pragma omp parallel for
	for(long int i = 0; i< betay.getGrid().getSize();i++){
		betay.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*methodtimestep/(4.0*gridstepy*gridstepy);
	}

	#pragma omp parallel for
	for(long int i = 0; i< betaz.getGrid().getSize();i++){
		betaz.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*methodtimestep/(4.0*gridstepz*gridstepz);
	}


	/* x direction*//*

	DiscreteFunction3D<std::complex<complextype>,gridtype> copy(wavefunction);

	#pragma omp parallel for
	for(long int l=0;l<ymax;++l)
		for(long int m=0;m<zmax;++m){
			//Multiply vector
			wavefunction(0,l,m)=betax(0,l,m)*copy(1,l,m)/effectivemass(1,l,m)+(1.0-2.0*betax(0,l,m)/effectivemass(0,l,m))*copy(0,l,m);
			wavefunction(xlast,l,m)=betax(xlast,l,m)*copy(xlast-1,l,m)/effectivemass+(1.0-2.0*betax(xlast,l,m)/effectivemass)*copy(xlast,l,m);

			for(long int k = 1; k < xlast; k++){
					wavefunction(k,l,m)=betax(k,l,m)*copy(k-1,l,m)+(1.0-2.0*betax(k,l,m))*copy(k,l,m)+betax(k,l,m)*copy(k+1,l,m);
			}
	}


	#pragma omp parallel for
	for(long int l=0;l<ymax;++l)
		for(long int m=0;m<zmax;++m){
			as(0,l,m)=0;
			bs(0,l,m)=2.0*betax(0,l,m)+1.0;
			cs(0,l,m)=-betax(0,l,m);

			for(long int i = 1; i < xlast ; i++){
				as(i,l,m)=-betax(i,l,m);
				bs(i,l,m)=2.0*betax(i,l,m)+1.0;
				cs(i,l,m)=-betax(i,l,m);
			}

			as(xlast,l,m)=-betax(xlast,l,m);
			bs(xlast,l,m)=2.0*betax(xlast,l,m)+1.0;
			cs(xlast,l,m)=0;

			//Thomas Algorithm

			for (long int i = 1; i <= xlast; i++){
				std::complex<complextype> n = as(i,l,m)/bs(i-1,l,m);
				bs(i,l,m)-=n*cs(i-1,l,m);
				wavefunction(i,l,m)-=n*wavefunction(i-1,l,m);
			}

			wavefunction(xlast,l,m)/=bs(xlast,l,m);
			for (long int i = (xlast - 1); i >= 0; --i)
				wavefunction(i,l,m)=(wavefunction(i,l,m)-cs(i,l,m)*wavefunction(i+1,l,m))/bs(i,l,m);

		}

	/* y direction*//*

	copy=wavefunction;

	#pragma omp parallel for
	for(long int k=0;k<xmax;++k)
		for(long int m=0;m<zmax;++m){
			//Multiply vector
			wavefunction(k,0,m)=betay(k,0,m)*copy(k,1,m)+(1.0-2.0*betay(k,0,m))*copy(k,0,m);
			wavefunction(k,ylast,m)=betay(k,ylast,m)*copy(k,ylast-1,m)+(1.0-2.0*betay(k,ylast,m))*copy(k,ylast,m);

			for(long int l = 1; l < ylast; l++){
					wavefunction(k,l,m)=betay(k,l,m)*copy(k,l-1,m)+(1.0-2.0*betay(k,l,m))*copy(k,l,m)+betay(k,l,m)*copy(k,l+1,m);
			}
		}


	#pragma omp parallel for
	for(long int k=0;k<xmax;++k)
		for(long int m=0;m<zmax;++m){
			as(k,0,m)=0;
			bs(k,0,m)=2.0*betay(k,0,m)+1.0;
			cs(k,0,m)=-betay(k,0,m);

			for(long int i = 1; i < ylast ; i++){
				as(k,i,m)=-betay(k,i,m);
				bs(k,i,m)=2.0*betay(k,i,m)+1.0;
				cs(k,i,m)=-betay(k,i,m);
			}

			as(k,ylast,m)=-betay(k,ylast,m);
			bs(k,ylast,m)=2.0*betay(k,ylast,m)+1.0;
			cs(k,ylast,m)=0;

			//Thomas Algorithm

			for (long int i = 1; i <= ylast; i++){
				std::complex<complextype> n = as(k,i,m)/bs(k,i-1,m);
				bs(k,i,m)-=n*cs(k,i-1,m);
				wavefunction(k,i,m)-=n*wavefunction(k,i-1,m);
			}

			wavefunction(k,ylast,m)/=bs(k,ylast,m);
			for (long int i = (ylast - 1); i >= 0; --i)
				wavefunction(k,i,m)=(wavefunction(k,i,m)-cs(k,i,m)*wavefunction(k,i+1,m))/bs(k,i,m);

		}


	/* z direction*//*

	copy=wavefunction;

	#pragma omp parallel for
	for(long int k=0;k<xmax;++k)
		for(long int l=0;l<ymax;++l){
			//Multiply vector
			wavefunction(k,l,0)=betaz(k,l,0)*copy(k,l,1)+(1.0-2.0*betaz(k,l,0))*copy(k,l,0);
			wavefunction(k,l,zlast)=betaz(k,l,zlast)*copy(k,l,zlast-1)+(1.0-2.0*betaz(k,l,zlast))*copy(k,l,zlast);

			for(long int m = 1; m < zlast; m++){
					wavefunction(k,l,m)=betaz(k,l,m)*copy(k,l,m-1)+(1.0-2.0*betaz(k,l,m))*copy(k,l,m)+betaz(k,l,m)*copy(k,l,m+1);
			}
		}


	#pragma omp parallel for
	for(long int k=0;k<xmax;++k)
		for(long int l=0;l<ymax;++l){
			as(k,l,0)=0;
			bs(k,l,0)=2.0*betaz(k,l,0)+1.0;
			cs(k,l,0)=-betaz(k,l,0);

			for(long int i = 1; i < zlast ; i++){
				as(k,l,i)=-betaz(k,l,i);
				bs(k,l,i)=2.0*betaz(k,l,i)+1.0;
				cs(k,l,i)=-betaz(k,l,i);
			}

			as(k,l,zlast)=-betaz(k,l,zlast);
			bs(k,l,zlast)=2.0*betaz(k,l,zlast)+1.0;
			cs(k,l,zlast)=0;

			//Thomas Algorithm

			for (long int i = 1; i <= zlast; i++){
				std::complex<complextype> n = as(k,l,i)/bs(k,l,i-1);
				bs(k,l,i)-=n*cs(k,l,i-1);
				wavefunction(k,l,i)-=n*wavefunction(k,l,i-1);
			}

			wavefunction(k,l,zlast)/=bs(k,l,zlast);
			for (long int i = (zlast - 1); i >= 0; --i)
				wavefunction(k,l,i)=(wavefunction(k,l,i)-cs(k,l,i)*wavefunction(k,l,i+1))/bs(k,l,i);

		}

	return wavefunction;

}
*/


template <typename timetype, typename complextype>
template<typename T, typename gridtype>
std::complex<complextype> WaveTimeEvolution3D<timetype,complextype>::Energy(DiscreteFunction3D<std::complex<complextype>,gridtype>& wavefunction, DiscreteFunction3D<T,gridtype>& potential,DiscreteFunction3D<T,gridtype>& effectivemass){
	std::complex<complextype> energysign;

		if(carriertype==Carrier::ELECTRON)
			energysign=std::complex<complextype>(1.0,0);
		else
			energysign=std::complex<complextype>(-1.0,0);

	wavefunction.normalize(1.0);
	DiscreteFunction3D<std::complex<complextype>,gridtype> normalized(wavefunction);


	DiscreteFunction3D<std::complex<complextype>,gridtype> x(wavefunction.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> y(wavefunction.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> z(wavefunction.getGrid());

	DiscreteFunction3D<std::complex<complextype>,gridtype> xdois(wavefunction.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> ydois(wavefunction.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> zdois(wavefunction.getGrid());

	long int xmax=wavefunction.getGrid().getSizeX();
	long int ymax=wavefunction.getGrid().getSizeY();
	long int zmax=wavefunction.getGrid().getSizeZ();

	DiscreteFunction3D<std::complex<complextype>,gridtype> copy2(wavefunction.getGrid());

	gridtype dx = wavefunction.getGrid().getxIncrement();
	gridtype dy = wavefunction.getGrid().getyIncrement();
	gridtype dz = wavefunction.getGrid().getzIncrement();

	for(long int l=0;l<ymax;++l){
		for(long int m=0;m<zmax;++m){
			//edge: order h and h²
			x(0,l,m) = (wavefunction(1,l,m)-wavefunction(0,l,m))/dx;
			x(1,l,m) = (wavefunction(2,l,m)-wavefunction(0,l,m))/(2*dx);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int k = 2; k < (xmax-2); k++){
				x(k,l,m)=(wavefunction(k+1,l,m)*8.0-wavefunction(k+2,l,m)-wavefunction(k-1,l,m)*8.0+wavefunction(k-2,l,m))/(12.0*dx);
			}

			//edge: order h and h²
			x(xmax-2,l,m)= (wavefunction(xmax-1,l,m)-wavefunction(xmax-3,l,m))/(2.0*dx);
			x(xmax-1,l,m)= (wavefunction(xmax-1,l,m)-wavefunction(xmax-2,l,m))/dx;
		}
	}

	for(long int k=0;k<xmax;++k){
		for(long int m=0;m<zmax;++m){
			//edge: order h and h²
			y(k,0,m) = (wavefunction(k,1,m)-wavefunction(k,0,m))/dy;
			y(k,1,m) = (wavefunction(k,2,m)-wavefunction(k,0,m))/(2*dy);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int l = 2; l < (ymax-2); l++){
				y(k,l,m)=(wavefunction(k,l+1,m)*8.0-wavefunction(k,l+2,m)-wavefunction(k,l-1,m)*8.0+wavefunction(k,l-2,m))/(12.0*dy);
			}

			//edge: order h and h²
			y(k,ymax-2,m)= (wavefunction(k,ymax-1,m)-wavefunction(k,ymax-3,m))/(2.0*dy);
			y(k,ymax-1,m)= (wavefunction(k,ymax-1,m)-wavefunction(k,ymax-2,m))/dy;
		}
	}

	for(long int k=0;k<xmax;++k){
		for(long int l=0;l<ymax;++l){
			//edge: order h and h²
			z(k,l,0) = (wavefunction(k,l,1)-wavefunction(k,l,0))/dz;
			z(k,l,1) = (wavefunction(k,l,2)-wavefunction(k,l,0))/(2*dz);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int m = 2; m < (zmax-2); m++){
				z(k,l,m)=(wavefunction(k,l,m+1)*8.0-wavefunction(k,l,m+2)-wavefunction(k,l,m-1)*8.0+wavefunction(k,l,m-2))/(12.0*dz);
			}

			//edge: order h and h²
			z(k,l,zmax-2)= (wavefunction(k,l,zmax-1)-wavefunction(k,l,zmax-3))/(2.0*dz);
			z(k,l,zmax-1)= (wavefunction(k,l,zmax-1)-wavefunction(k,l,zmax-2))/dz;
		}
	}

	#pragma omp parallel for
	for(long int n=0; n < wavefunction.getGrid().getSize(); ++n){
		x.getdata()[n]/=effectivemass.getdata()[n];
		y.getdata()[n]/=effectivemass.getdata()[n];
		z.getdata()[n]/=effectivemass.getdata()[n];
	}

	for(long int l=0;l<ymax;++l){
		for(long int m=0;m<zmax;++m){
			//edge: order h and h²
			xdois(0,l,m) = (x(1,l,m)-x(0,l,m))/dx;
			xdois(1,l,m) = (x(2,l,m)-x(0,l,m))/(2*dx);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int k = 2; k < (xmax-2); k++){
				xdois(k,l,m)=(x(k+1,l,m)*8.0-x(k+2,l,m)-x(k-1,l,m)*8.0+x(k-2,l,m))/(12.0*dx);
			}

			//edge: order h and h²
			xdois(xmax-2,l,m)= (x(xmax-1,l,m)-x(xmax-3,l,m))/(2.0*dx);
			xdois(xmax-1,l,m)= (x(xmax-1,l,m)-x(xmax-2,l,m))/dx;
		}
	}

	for(long int k=0;k<xmax;++k){
		for(long int m=0;m<zmax;++m){
			//edge: order h and h²
			ydois(k,0,m) = (y(k,1,m)-y(k,0,m))/dy;
			ydois(k,1,m) = (y(k,2,m)-y(k,0,m))/(2*dy);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int l = 2; l < (ymax-2); l++){
				ydois(k,l,m)=(y(k,l+1,m)*8.0-y(k,l+2,m)-y(k,l-1,m)*8.0+y(k,l-2,m))/(12.0*dy);
			}

			//edge: order h and h²
			ydois(k,ymax-2,m)= (y(k,ymax-1,m)-y(k,ymax-3,m))/(2.0*dy);
			ydois(k,ymax-1,m)= (y(k,ymax-1,m)-y(k,ymax-2,m))/dy;
		}
	}

	for(long int k=0;k<xmax;++k){
		for(long int l=0;l<ymax;++l){
			//edge: order h and h²
			zdois(k,l,0) = (z(k,l,1)-z(k,l,0))/dz;
			zdois(k,l,1) = (z(k,l,2)-z(k,l,0))/(2*dz);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int m = 2; m < (zmax-2); m++){
				zdois(k,l,m)=(z(k,l,m+1)*8.0-z(k,l,m+2)-z(k,l,m-1)*8.0+z(k,l,m-2))/(12.0*dz);
			}

			//edge: order h and h²
			zdois(k,l,zmax-2)= (z(k,l,zmax-1)-z(k,l,zmax-3))/(2.0*dz);
			zdois(k,l,zmax-1)= (z(k,l,zmax-1)-z(k,l,zmax-2))/dz;
		}
	}

	#pragma omp parallel for
	for(long int n=0; n < wavefunction.getGrid().getSize(); ++n){
		copy2.getdata()[n]=xdois.getdata()[n]+ydois.getdata()[n]+zdois.getdata()[n];
	}
	copy2*=(-Constant::hbar*Constant::hbar/(2.0));

	DiscreteFunction3D<std::complex<complextype>,gridtype> copy(wavefunction);
	#pragma omp parallel for
	for(long int n=0; n < wavefunction.getGrid().getSize(); ++n){
		copy.getdata()[n]*=potential.getdata()[n]*energysign;
	}

	/*auto kinect = normalized.Dot(copy2);
	auto potener = normalized.Dot(copy);*/

	//cout << "potential: " << potener << " kinect: " << kinect << std::endl;

	return energysign*(normalized.Dot(copy)+normalized.Dot(copy2));
}


template <typename timetype, typename complextype>
template<typename T, typename gridtype>
std::complex<complextype> WaveTimeEvolution3D<timetype,complextype>::Energy(DiscreteFunction3D<std::complex<complextype>,gridtype>& wavefunction, DiscreteFunction3D<T,gridtype>& potential,DiscreteFunction3D<T,gridtype>& effectivemasst,DiscreteFunction3D<T,gridtype>& effectivemassl){
	std::complex<complextype> energysign;

		if(carriertype==Carrier::ELECTRON)
			energysign=std::complex<complextype>(1.0,0);
		else
			energysign=std::complex<complextype>(-1.0,0);

	wavefunction.normalize(1.0);
	DiscreteFunction3D<std::complex<complextype>,gridtype> normalized(wavefunction);


	DiscreteFunction3D<std::complex<complextype>,gridtype> x(wavefunction.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> y(wavefunction.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> z(wavefunction.getGrid());

	DiscreteFunction3D<std::complex<complextype>,gridtype> xdois(wavefunction.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> ydois(wavefunction.getGrid());
	DiscreteFunction3D<std::complex<complextype>,gridtype> zdois(wavefunction.getGrid());

	long int xmax=wavefunction.getGrid().getSizeX();
	long int ymax=wavefunction.getGrid().getSizeY();
	long int zmax=wavefunction.getGrid().getSizeZ();

	DiscreteFunction3D<std::complex<complextype>,gridtype> copy2(wavefunction.getGrid());

	gridtype dx = wavefunction.getGrid().getxIncrement();
	gridtype dy = wavefunction.getGrid().getyIncrement();
	gridtype dz = wavefunction.getGrid().getzIncrement();

	for(long int l=0;l<ymax;++l){
		for(long int m=0;m<zmax;++m){
			//edge: order h and h²
			x(0,l,m) = (wavefunction(1,l,m)-wavefunction(0,l,m))/dx;
			x(1,l,m) = (wavefunction(2,l,m)-wavefunction(0,l,m))/(2*dx);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int k = 2; k < (xmax-2); k++){
				x(k,l,m)=(wavefunction(k+1,l,m)*8.0-wavefunction(k+2,l,m)-wavefunction(k-1,l,m)*8.0+wavefunction(k-2,l,m))/(12.0*dx);
			}

			//edge: order h and h²
			x(xmax-2,l,m)= (wavefunction(xmax-1,l,m)-wavefunction(xmax-3,l,m))/(2.0*dx);
			x(xmax-1,l,m)= (wavefunction(xmax-1,l,m)-wavefunction(xmax-2,l,m))/dx;
		}
	}

	for(long int k=0;k<xmax;++k){
		for(long int m=0;m<zmax;++m){
			//edge: order h and h²
			y(k,0,m) = (wavefunction(k,1,m)-wavefunction(k,0,m))/dy;
			y(k,1,m) = (wavefunction(k,2,m)-wavefunction(k,0,m))/(2*dy);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int l = 2; l < (ymax-2); l++){
				y(k,l,m)=(wavefunction(k,l+1,m)*8.0-wavefunction(k,l+2,m)-wavefunction(k,l-1,m)*8.0+wavefunction(k,l-2,m))/(12.0*dy);
			}

			//edge: order h and h²
			y(k,ymax-2,m)= (wavefunction(k,ymax-1,m)-wavefunction(k,ymax-3,m))/(2.0*dy);
			y(k,ymax-1,m)= (wavefunction(k,ymax-1,m)-wavefunction(k,ymax-2,m))/dy;
		}
	}

	for(long int k=0;k<xmax;++k){
		for(long int l=0;l<ymax;++l){
			//edge: order h and h²
			z(k,l,0) = (wavefunction(k,l,1)-wavefunction(k,l,0))/dz;
			z(k,l,1) = (wavefunction(k,l,2)-wavefunction(k,l,0))/(2*dz);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int m = 2; m < (zmax-2); m++){
				z(k,l,m)=(wavefunction(k,l,m+1)*8.0-wavefunction(k,l,m+2)-wavefunction(k,l,m-1)*8.0+wavefunction(k,l,m-2))/(12.0*dz);
			}

			//edge: order h and h²
			z(k,l,zmax-2)= (wavefunction(k,l,zmax-1)-wavefunction(k,l,zmax-3))/(2.0*dz);
			z(k,l,zmax-1)= (wavefunction(k,l,zmax-1)-wavefunction(k,l,zmax-2))/dz;
		}
	}

	#pragma omp parallel for
	for(long int n=0; n < wavefunction.getGrid().getSize(); ++n){
		x.getdata()[n]/=effectivemassl.getdata()[n];
		y.getdata()[n]/=effectivemassl.getdata()[n];
		z.getdata()[n]/=effectivemasst.getdata()[n];
	}

	for(long int l=0;l<ymax;++l){
		for(long int m=0;m<zmax;++m){
			//edge: order h and h²
			xdois(0,l,m) = (x(1,l,m)-x(0,l,m))/dx;
			xdois(1,l,m) = (x(2,l,m)-x(0,l,m))/(2*dx);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int k = 2; k < (xmax-2); k++){
				xdois(k,l,m)=(x(k+1,l,m)*8.0-x(k+2,l,m)-x(k-1,l,m)*8.0+x(k-2,l,m))/(12.0*dx);
			}

			//edge: order h and h²
			xdois(xmax-2,l,m)= (x(xmax-1,l,m)-x(xmax-3,l,m))/(2.0*dx);
			xdois(xmax-1,l,m)= (x(xmax-1,l,m)-x(xmax-2,l,m))/dx;
		}
	}

	for(long int k=0;k<xmax;++k){
		for(long int m=0;m<zmax;++m){
			//edge: order h and h²
			ydois(k,0,m) = (y(k,1,m)-y(k,0,m))/dy;
			ydois(k,1,m) = (y(k,2,m)-y(k,0,m))/(2*dy);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int l = 2; l < (ymax-2); l++){
				ydois(k,l,m)=(y(k,l+1,m)*8.0-y(k,l+2,m)-y(k,l-1,m)*8.0+y(k,l-2,m))/(12.0*dy);
			}

			//edge: order h and h²
			ydois(k,ymax-2,m)= (y(k,ymax-1,m)-y(k,ymax-3,m))/(2.0*dy);
			ydois(k,ymax-1,m)= (y(k,ymax-1,m)-y(k,ymax-2,m))/dy;
		}
	}

	for(long int k=0;k<xmax;++k){
		for(long int l=0;l<ymax;++l){
			//edge: order h and h²
			zdois(k,l,0) = (z(k,l,1)-z(k,l,0))/dz;
			zdois(k,l,1) = (z(k,l,2)-z(k,l,0))/(2*dz);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int m = 2; m < (zmax-2); m++){
				zdois(k,l,m)=(z(k,l,m+1)*8.0-z(k,l,m+2)-z(k,l,m-1)*8.0+z(k,l,m-2))/(12.0*dz);
			}

			//edge: order h and h²
			zdois(k,l,zmax-2)= (z(k,l,zmax-1)-z(k,l,zmax-3))/(2.0*dz);
			zdois(k,l,zmax-1)= (z(k,l,zmax-1)-z(k,l,zmax-2))/dz;
		}
	}

	#pragma omp parallel for
	for(long int n=0; n < wavefunction.getGrid().getSize(); ++n){
		copy2.getdata()[n]=xdois.getdata()[n]+ydois.getdata()[n]+zdois.getdata()[n];
	}
	copy2*=(-Constant::hbar*Constant::hbar/(2.0));

	DiscreteFunction3D<std::complex<complextype>,gridtype> copy(wavefunction);
	#pragma omp parallel for
	for(long int n=0; n < wavefunction.getGrid().getSize(); ++n){
		copy.getdata()[n]*=potential.getdata()[n]*energysign;
	}

	/*auto kinect = normalized.Dot(copy2);
	auto potener = normalized.Dot(copy);*/

	//cout << "potential: " << potener << " kinect: " << kinect << std::endl;

	return energysign*(normalized.Dot(copy)+normalized.Dot(copy2));
}

} /* namespace epital */

#endif
