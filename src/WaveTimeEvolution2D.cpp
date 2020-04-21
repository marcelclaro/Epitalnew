/*
 * WaveTimeEvolution2D.cpp
 *
 *  Created on: Nov 7, 2014
 *      Author: marcel
 */

#include "WaveTimeEvolution2D.hpp"

#ifndef WAVETIMEEVOLUTION2D_CPP_
#define WAVETIMEEVOLUTION2D_CPP_

namespace epital {

template <typename timetype, typename complextype>
WaveTimeEvolution2D<timetype,complextype>::WaveTimeEvolution2D(WaveTimeEvolution2D::Carrier carrier): carriertype(carrier) {
	// TODO Auto-generated constructor stub

}

template <typename timetype, typename complextype>
WaveTimeEvolution2D<timetype,complextype>::~WaveTimeEvolution2D() {
	// TODO Auto-generated destructor stub
}

template <typename timetype, typename complextype>
template<typename T, typename gridtype>
DiscreteFunction2D<std::complex<complextype>,gridtype> WaveTimeEvolution2D<timetype,complextype>::potentialpropagator(timetype methodtimestep, DiscreteFunction2D<T,gridtype> potential){
	std::complex<complextype> energysign;

	if(carriertype==Carrier::ELECTRON)
		energysign=std::complex<complextype>(1.0,0);
	else
		energysign=std::complex<complextype>(-1.0,0);

	DiscreteFunction2D<std::complex<complextype>,gridtype> result(potential.getGrid());


	for(long int k=0; k<potential.getGrid().getSizeX();k++)
		for(long int l=0; l<potential.getGrid().getSizeY();l++){
				result(k,l)=exp(std::complex<complextype>(0.0,-1.0)*energysign*methodtimestep*potential(k,l)/(2.0*Constant::hbar));
		}

	return result;
}

template <typename timetype, typename complextype>
template<typename T, typename gridtype>
DiscreteFunction2D<std::complex<complextype>,gridtype>& WaveTimeEvolution2D<timetype,complextype>::groundlevel_periodicFFT(timetype time, long int timesteps, DiscreteFunction2D<std::complex<complextype>,gridtype>& initial, DiscreteFunction2D<T,gridtype>& potential, T effectivemass){
	timetype methodtimestep=time/static_cast<timetype>(timesteps);


	initial/=initial.Modulus();

	DiscreteFunction2D<std::complex<complextype>,gridtype> potpropag = potentialpropagator(methodtimestep,potential);
	DiscreteFunction2D<std::complex<complextype>,gridtype> bloch = fftmulti(methodtimestep,initial.getGrid(),effectivemass);

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
DiscreteFunction2D<std::complex<complextype>,gridtype>& WaveTimeEvolution2D<timetype,complextype>::groundlevel_finite_valence(timetype time, long int timesteps, DiscreteFunction2D<std::complex<complextype>,gridtype>& initial, DiscreteFunction2D<T,gridtype>& potential, DiscreteFunction2D<T,gridtype>& effectivemasst,DiscreteFunction2D<T,gridtype>& effectivemassl){
	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	initial/=initial.Modulus();

	DiscreteFunction2D<std::complex<complextype>,gridtype> potpropag = potentialpropagator(methodtimestep,potential);

	timetype adi_methodtimestep = methodtimestep/2.0;

	gridtype gridstepx = initial.getGrid().getxIncrement();
	gridtype gridstepy = initial.getGrid().getyIncrement();



	long int xmax=initial.getGrid().getSizeX();
	long int ymax=initial.getGrid().getSizeY();


	long int xlast=xmax-1;
	long int ylast=ymax-1;



	DiscreteFunction2D<std::complex<complextype>,gridtype> as(initial.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> bs(initial.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> cs(initial.getGrid());

	DiscreteFunction2D<std::complex<complextype>,gridtype> betax(initial.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> betay(initial.getGrid());


	#pragma omp parallel for
	for(long int i = 0; i< betax.getGrid().getSize();i++){
		betax.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*adi_methodtimestep/(4.0*effectivemasst.getdata()[i]*gridstepx*gridstepx);
	}

	#pragma omp parallel for
	for(long int i = 0; i< betay.getGrid().getSize();i++){
		betay.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*adi_methodtimestep/(4.0*effectivemassl.getdata()[i]*gridstepy*gridstepy);
	}



	DiscreteFunction2D<std::complex<complextype>,gridtype> x(initial.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> y(initial.getGrid());


	gridtype dx = gridstepx;
	gridtype dy = gridstepy;


	for(long int l=0;l<ymax;++l){
			//edge: order h and h²
			x(0,l) = (effectivemasst(0,l)*dx/2.0)*(1.0/effectivemasst(1,l)-1.0/effectivemasst(0,l))/dx;
			x(1,l) = (effectivemasst(1,l)*dx/2.0)*(1.0/effectivemasst(2,l)-1.0/effectivemasst(0,l))/(2*dx);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int k = 2; k < (xmax-2); k++){
				x(k,l)=(effectivemasst(k,l)*dx/2.0)*((1.0/effectivemasst(k+1,l))*8.0-(1.0/effectivemasst(k+2,l))-(1.0/effectivemasst(k-1,l))*8.0+(1.0/effectivemasst(k-2,l)))/(12.0*dx);
			}

			//edge: order h and h²
			x(xmax-2,l)= (effectivemasst(xmax-2,l)*dx/2.0)*(1.0/effectivemasst(xmax-1,l)-1.0/effectivemasst(xmax-3,l))/(2.0*dx);
			x(xmax-1,l)= (effectivemasst(xmax-1,l)*dx/2.0)*(1.0/effectivemasst(xmax-1,l)-1.0/effectivemasst(xmax-2,l))/dx;
	}

	for(long int k=0;k<xmax;++k){

			//edge: order h and h²
			y(k,0) = (effectivemassl(k,0)*dy/2.0)*(1.0/effectivemassl(k,1)-1.0/effectivemassl(k,0))/dy;
			y(k,1) = (effectivemassl(k,1)*dy/2.0)*(1.0/effectivemassl(k,2)-1.0/effectivemassl(k,0))/(2*dy);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int l = 2; l < (ymax-2); l++){
				y(k,l)=(effectivemassl(k,l)*dy/2.0)*((1.0/effectivemassl(k,l+1))*8.0-(1.0/effectivemassl(k,l+2))-(1.0/effectivemassl(k,l-1))*8.0+(1.0/effectivemassl(k,l-2)))/(12.0*dy);
			}

			//edge: order h and h²
			y(k,ymax-2)= (effectivemassl(k,ymax-2)*dy/2.0)*(1.0/effectivemassl(k,ymax-1)-1.0/effectivemassl(k,ymax-3))/(2.0*dy);
			y(k,ymax-1)= (effectivemassl(k,ymax-1)*dy/2.0)*(1.0/effectivemassl(k,ymax-1)-1.0/effectivemassl(k,ymax-2))/dy;

	}


	for(long int i = 0 ; i< timesteps; i++){
		std::cout << "step: " << i  << " ..." << std::endl;
		initial*=potpropag;


			/* x direction*/

			DiscreteFunction2D<std::complex<complextype>,gridtype> copy(initial);

			#pragma omp parallel for
			for(long int l=0;l<ymax;++l){
					//Multiply vector
					initial(0,l)=betax(0,l)*copy(1,l)*(1.0+x(0,l))+(1.0-2.0*betax(0,l))*copy(0,l);
					initial(xlast,l)=betax(xlast,l)*copy(xlast-1,l)*(1.0-x(xlast,l))+(1.0-2.0*betax(xlast,l))*copy(xlast,l);

					for(long int k = 1; k < xlast; k++){
							initial(k,l)=betax(k,l)*copy(k-1,l)*(1.0-x(k,l))+(1.0-2.0*betax(k,l))*copy(k,l)+betax(k,l)*copy(k+1,l)*(1.0+x(k,l));
					}
			}


			#pragma omp parallel for
			for(long int l=0;l<ymax;++l){
					as(0,l)=0;
					bs(0,l)=2.0*betax(0,l)+1.0;
					cs(0,l)=-betax(0,l)*(1.0+x(0,l));

					for(long int i = 1; i < xlast ; i++){
						as(i,l)=-betax(i,l)*(1.0-x(i,l));
						bs(i,l)=2.0*betax(i,l)+1.0;
						cs(i,l)=-betax(i,l)*(1.0+x(i,l));
					}

					as(xlast,l)=-betax(xlast,l)*(1.0-x(xlast,l));
					bs(xlast,l)=2.0*betax(xlast,l)+1.0;
					cs(xlast,l)=0;

					//Thomas Algorithm

					for (long int i = 1; i <= xlast; i++){
						std::complex<complextype> n = as(i,l)/bs(i-1,l);
						bs(i,l)-=n*cs(i-1,l);
						initial(i,l)-=n*initial(i-1,l);
					}

					initial(xlast,l)/=bs(xlast,l);
					for (long int i = (xlast - 1); i >= 0; --i)
						initial(i,l)=(initial(i,l)-cs(i,l)*initial(i+1,l))/bs(i,l);

				}

			/* y direction*/

			copy=initial;

			#pragma omp parallel for
			for(long int k=0;k<xmax;++k){
					//Multiply vector
					initial(k,0)=betay(k,0)*copy(k,1)*(1.0+y(k,0))+(1.0-2.0*betay(k,0))*copy(k,0);
					initial(k,ylast)=betay(k,ylast)*copy(k,ylast-1)*(1.0-y(k,ylast))+(1.0-2.0*betay(k,ylast))*copy(k,ylast);

					for(long int l = 1; l < ylast; l++){
							initial(k,l)=betay(k,l)*copy(k,l-1)*(1.0-y(k,l))+(1.0-2.0*betay(k,l))*copy(k,l)+betay(k,l)*copy(k,l+1)*(1.0+y(k,l));
					}
				}


			#pragma omp parallel for
			for(long int k=0;k<xmax;++k){
					as(k,0)=0;
					bs(k,0)=2.0*betay(k,0)+1.0;
					cs(k,0)=-betay(k,0)*(1.0+y(k,0));

					for(long int i = 1; i < ylast ; i++){
						as(k,i)=-betay(k,i)*(1.0-y(k,i));
						bs(k,i)=2.0*betay(k,i)+1.0;
						cs(k,i)=-betay(k,i)*(1.0+y(k,i));
					}

					as(k,ylast)=-betay(k,ylast)*(1.0-y(k,ylast));
					bs(k,ylast)=2.0*betay(k,ylast)+1.0;
					cs(k,ylast)=0;

					//Thomas Algorithm

					for (long int i = 1; i <= ylast; i++){
						std::complex<complextype> n = as(k,i)/bs(k,i-1);
						bs(k,i)-=n*cs(k,i-1);
						initial(k,i)-=n*initial(k,i-1);
					}

					initial(k,ylast)/=bs(k,ylast);
					for (long int i = (ylast - 1); i >= 0; --i)
						initial(k,i)=(initial(k,i)-cs(k,i)*initial(k,i+1))/bs(k,i);

				}


		initial*=potpropag;
		std::complex<complextype> tonorm=initial.Modulus();
		initial*=1.0/tonorm;

	}

	return initial;
}


template <typename timetype, typename complextype>
template<typename T, typename gridtype>
DiscreteFunction2D<std::complex<complextype>,gridtype>& WaveTimeEvolution2D<timetype,complextype>::groundlevel_finite(timetype time, long int timesteps, DiscreteFunction2D<std::complex<complextype>,gridtype>& initial, DiscreteFunction2D<T,gridtype>& potential, DiscreteFunction2D<T,gridtype>& effectivemass){
	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	initial/=initial.Modulus();

	DiscreteFunction2D<std::complex<complextype>,gridtype> potpropag = potentialpropagator(methodtimestep,potential);

	timetype adi_methodtimestep = methodtimestep/2.0;

	gridtype gridstepx = initial.getGrid().getxIncrement();
	gridtype gridstepy = initial.getGrid().getyIncrement();


	long int xmax=initial.getGrid().getSizeX();
	long int ymax=initial.getGrid().getSizeY();


	long int xlast=xmax-1;
	long int ylast=ymax-1;



	DiscreteFunction2D<std::complex<complextype>,gridtype> as(initial.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> bs(initial.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> cs(initial.getGrid());

	DiscreteFunction2D<std::complex<complextype>,gridtype> betax(initial.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> betay(initial.getGrid());


	#pragma omp parallel for
	for(long int i = 0; i< betax.getGrid().getSize();i++){
		betax.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*adi_methodtimestep/(4.0*effectivemass.getdata()[i]*gridstepx*gridstepx);
	}

	#pragma omp parallel for
	for(long int i = 0; i< betay.getGrid().getSize();i++){
		betay.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*adi_methodtimestep/(4.0*effectivemass.getdata()[i]*gridstepy*gridstepy);
	}



	DiscreteFunction2D<std::complex<complextype>,gridtype> x(initial.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> y(initial.getGrid());


	gridtype dx = gridstepx;
	gridtype dy = gridstepy;

	for(long int l=0;l<ymax;++l){

			//edge: order h and h²
			x(0,l) = (effectivemass(0,l)*dx/2.0)*(1.0/effectivemass(1,l)-1.0/effectivemass(0,l))/dx;
			x(1,l) = (effectivemass(1,l)*dx/2.0)*(1.0/effectivemass(2,l)-1.0/effectivemass(0,l))/(2*dx);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int k = 2; k < (xmax-2); k++){
				x(k,l)=(effectivemass(k,l)*dx/2.0)*((1.0/effectivemass(k+1,l))*8.0-(1.0/effectivemass(k+2,l))-(1.0/effectivemass(k-1,l))*8.0+(1.0/effectivemass(k-2,l)))/(12.0*dx);
			}

			//edge: order h and h²
			x(xmax-2,l)= (effectivemass(xmax-2,l)*dx/2.0)*(1.0/effectivemass(xmax-1,l)-1.0/effectivemass(xmax-3,l))/(2.0*dx);
			x(xmax-1,l)= (effectivemass(xmax-1,l)*dx/2.0)*(1.0/effectivemass(xmax-1,l)-1.0/effectivemass(xmax-2,l))/dx;

	}

	for(long int k=0;k<xmax;++k){

			//edge: order h and h²
			y(k,0) = (effectivemass(k,0)*dy/2.0)*(1.0/effectivemass(k,1)-1.0/effectivemass(k,0))/dy;
			y(k,1) = (effectivemass(k,1)*dy/2.0)*(1.0/effectivemass(k,2)-1.0/effectivemass(k,0))/(2*dy);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int l = 2; l < (ymax-2); l++){
				y(k,l)=(effectivemass(k,l)*dy/2.0)*((1.0/effectivemass(k,l+1))*8.0-(1.0/effectivemass(k,l+2))-(1.0/effectivemass(k,l-1))*8.0+(1.0/effectivemass(k,l-2)))/(12.0*dy);
			}

			//edge: order h and h²
			y(k,ymax-2)= (effectivemass(k,ymax-2)*dy/2.0)*(1.0/effectivemass(k,ymax-1)-1.0/effectivemass(k,ymax-3))/(2.0*dy);
			y(k,ymax-1)= (effectivemass(k,ymax-1)*dy/2.0)*(1.0/effectivemass(k,ymax-1)-1.0/effectivemass(k,ymax-2))/dy;

	}


	for(long int i = 0 ; i< timesteps; i++){
		//std::cout << "step: " << i  << " ..." << std::endl;
		initial*=potpropag;

		{
			/* x direction*/

			DiscreteFunction2D<std::complex<complextype>,gridtype> copy(initial);

			#pragma omp parallel for
			for(long int l=0;l<ymax;++l){

					//Multiply vector
					initial(0,l)=betax(0,l)*copy(1,l)*(1.0+x(0,l))+(1.0-2.0*betax(0,l))*copy(0,l);
					initial(xlast,l)=betax(xlast,l)*copy(xlast-1,l)*(1.0-x(xlast,l))+(1.0-2.0*betax(xlast,l))*copy(xlast,l);

					for(long int k = 1; k < xlast; k++){
							initial(k,l)=betax(k,l)*copy(k-1,l)*(1.0-x(k,l))+(1.0-2.0*betax(k,l))*copy(k,l)+betax(k,l)*copy(k+1,l)*(1.0+x(k,l));

			}


			#pragma omp parallel for
			for(long int l=0;l<ymax;++l){
					as(0,l)=0;
					bs(0,l)=2.0*betax(0,l)+1.0;
					cs(0,l)=-betax(0,l)*(1.0+x(0,l));

					for(long int i = 1; i < xlast ; i++){
						as(i,l)=-betax(i,l)*(1.0-x(i,l));
						bs(i,l)=2.0*betax(i,l)+1.0;
						cs(i,l)=-betax(i,l)*(1.0+x(i,l));
					}

					as(xlast,l)=-betax(xlast,l)*(1.0-x(xlast,l));
					bs(xlast,l)=2.0*betax(xlast,l)+1.0;
					cs(xlast,l)=0;

					//Thomas Algorithm

					for (long int i = 1; i <= xlast; i++){
						std::complex<complextype> n = as(i,l)/bs(i-1,l);
						bs(i,l)-=n*cs(i-1,l);
						initial(i,l)-=n*initial(i-1,l);
					}

					initial(xlast,l)/=bs(xlast,l);
					for (long int i = (xlast - 1); i >= 0; --i)
						initial(i,l)=(initial(i,l)-cs(i,l)*initial(i+1,l))/bs(i,l);

				}

			/* y direction*/

			copy=initial;

			#pragma omp parallel for
			for(long int k=0;k<xmax;++k){

					//Multiply vector
					initial(k,0)=betay(k,0)*copy(k,1)*(1.0+y(k,0))+(1.0-2.0*betay(k,0))*copy(k,0);
					initial(k,ylast)=betay(k,ylast)*copy(k,ylast-1)*(1.0-y(k,ylast))+(1.0-2.0*betay(k,ylast))*copy(k,ylast);

					for(long int l = 1; l < ylast; l++){
							initial(k,l)=betay(k,l)*copy(k,l-1)*(1.0-y(k,l))+(1.0-2.0*betay(k,l))*copy(k,l)+betay(k,l)*copy(k,l+1)*(1.0+y(k,l));
						}
				}


			#pragma omp parallel for
			for(long int k=0;k<xmax;++k){
					as(k,0)=0;
					bs(k,0)=2.0*betay(k,0)+1.0;
					cs(k,0)=-betay(k,0)*(1.0+y(k,0));

					for(long int i = 1; i < ylast ; i++){
						as(k,i)=-betay(k,i)*(1.0-y(k,i));
						bs(k,i)=2.0*betay(k,i)+1.0;
						cs(k,i)=-betay(k,i)*(1.0+y(k,i));
					}

					as(k,ylast)=-betay(k,ylast)*(1.0-y(k,ylast));
					bs(k,ylast)=2.0*betay(k,ylast)+1.0;
					cs(k,ylast)=0;

					//Thomas Algorithm

					for (long int i = 1; i <= ylast; i++){
						std::complex<complextype> n = as(k,i)/bs(k,i-1);
						bs(k,i)-=n*cs(k,i-1);
						initial(k,i)-=n*initial(k,i-1);
					}

					initial(k,ylast)/=bs(k,ylast);
					for (long int i = (ylast - 1); i >= 0; --i)
						initial(k,i)=(initial(k,i)-cs(k,i)*initial(k,i+1))/bs(k,i);

				}


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
DiscreteFunction2D<std::complex<complextype>,gridtype>& WaveTimeEvolution2D<timetype,complextype>::newlevel_finite(timetype time, long int timesteps, DiscreteFunction2D<std::complex<complextype>,gridtype>& initial, DiscreteFunction2D<T,gridtype>& potential,std::vector<std::shared_ptr<DiscreteFunction2D<std::complex<complextype>,gridtype>>> otherslevels, DiscreteFunction2D<T,gridtype>& effectivemass){
	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	initial/=initial.Modulus();

	DiscreteFunction2D<std::complex<complextype>,gridtype> potpropag = potentialpropagator(methodtimestep,potential);

	timetype adi_methodtimestep = methodtimestep/2.0;

	gridtype gridstepx = initial.getGrid().getxIncrement();
	gridtype gridstepy = initial.getGrid().getyIncrement();



	long int xmax=initial.getGrid().getSizeX();
	long int ymax=initial.getGrid().getSizeY();


	long int xlast=xmax-1;
	long int ylast=ymax-1;



	DiscreteFunction2D<std::complex<complextype>,gridtype> as(initial.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> bs(initial.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> cs(initial.getGrid());

	DiscreteFunction2D<std::complex<complextype>,gridtype> betax(initial.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> betay(initial.getGrid());


	#pragma omp parallel for
	for(long int i = 0; i< betax.getGrid().getSize();i++){
		betax.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*adi_methodtimestep/(4.0*effectivemass.getdata()[i]*gridstepx*gridstepx);
	}

	#pragma omp parallel for
	for(long int i = 0; i< betay.getGrid().getSize();i++){
		betay.getdata()[i]=std::complex<complextype>(0.0,1.0)*Constant::hbar*adi_methodtimestep/(4.0*effectivemass.getdata()[i]*gridstepy*gridstepy);
	}



	DiscreteFunction2D<std::complex<complextype>,gridtype> x(initial.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> y(initial.getGrid());


	gridtype dx = gridstepx;
	gridtype dy = gridstepy;


	for(long int l=0;l<ymax;++l){

			//edge: order h and h²
			x(0,l) = (effectivemass(0,l)*dx/2.0)*(1.0/effectivemass(1,l)-1.0/effectivemass(0,l))/dx;
			x(1,l) = (effectivemass(1,l)*dx/2.0)*(1.0/effectivemass(2,l)-1.0/effectivemass(0,l))/(2*dx);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int k = 2; k < (xmax-2); k++){
				x(k,l)=(effectivemass(k,l)*dx/2.0)*((1.0/effectivemass(k+1,l))*8.0-(1.0/effectivemass(k+2,l))-(1.0/effectivemass(k-1,l))*8.0+(1.0/effectivemass(k-2,l)))/(12.0*dx);
			}

			//edge: order h and h²
			x(xmax-2,l)= (effectivemass(xmax-2,l)*dx/2.0)*(1.0/effectivemass(xmax-1,l)-1.0/effectivemass(xmax-3,l))/(2.0*dx);
			x(xmax-1,l)= (effectivemass(xmax-1,l)*dx/2.0)*(1.0/effectivemass(xmax-1,l)-1.0/effectivemass(xmax-2,l))/dx;

	}

	for(long int k=0;k<xmax;++k){

			//edge: order h and h²
			y(k,0) = (effectivemass(k,0)*dy/2.0)*(1.0/effectivemass(k,1)-1.0/effectivemass(k,0))/dy;
			y(k,1) = (effectivemass(k,1)*dy/2.0)*(1.0/effectivemass(k,2)-1.0/effectivemass(k,0))/(2*dy);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int l = 2; l < (ymax-2); l++){
				y(k,l)=(effectivemass(k,l)*dy/2.0)*((1.0/effectivemass(k,l+1))*8.0-(1.0/effectivemass(k,l+2))-(1.0/effectivemass(k,l-1))*8.0+(1.0/effectivemass(k,l-2)))/(12.0*dy);
			}

			//edge: order h and h²
			y(k,ymax-2)= (effectivemass(k,ymax-2)*dy/2.0)*(1.0/effectivemass(k,ymax-1)-1.0/effectivemass(k,ymax-3))/(2.0*dy);
			y(k,ymax-1)= (effectivemass(k,ymax-1)*dy/2.0)*(1.0/effectivemass(k,ymax-1)-1.0/effectivemass(k,ymax-2))/dy;

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

			DiscreteFunction2D<std::complex<complextype>,gridtype> copy(initial);

			#pragma omp parallel for
			for(long int l=0;l<ymax;++l){
					//Multiply vector
					initial(0,l)=betax(0,l)*copy(1,l)*(1.0+x(0,l))+(1.0-2.0*betax(0,l))*copy(0,l);
					initial(xlast,l)=betax(xlast,l)*copy(xlast-1,l)*(1.0-x(xlast,l))+(1.0-2.0*betax(xlast,l))*copy(xlast,l);

					for(long int k = 1; k < xlast; k++){
							initial(k,l)=betax(k,l)*copy(k-1,l)*(1.0-x(k,l))+(1.0-2.0*betax(k,l))*copy(k,l)+betax(k,l)*copy(k+1,l)*(1.0+x(k,l));
					}
			}


			#pragma omp parallel for
			for(long int l=0;l<ymax;++l){
					as(0,l)=0;
					bs(0,l)=2.0*betax(0,l)+1.0;
					cs(0,l)=-betax(0,l)*(1.0+x(0,l));

					for(long int i = 1; i < xlast ; i++){
						as(i,l)=-betax(i,l)*(1.0-x(i,l));
						bs(i,l)=2.0*betax(i,l)+1.0;
						cs(i,l)=-betax(i,l)*(1.0+x(i,l));
					}

					as(xlast,l)=-betax(xlast,l)*(1.0-x(xlast,l));
					bs(xlast,l)=2.0*betax(xlast,l)+1.0;
					cs(xlast,l)=0;

					//Thomas Algorithm

					for (long int i = 1; i <= xlast; i++){
						std::complex<complextype> n = as(i,l)/bs(i-1,l);
						bs(i,l)-=n*cs(i-1,l);
						initial(i,l)-=n*initial(i-1,l);
					}

					initial(xlast,l)/=bs(xlast,l);
					for (long int i = (xlast - 1); i >= 0; --i)
						initial(i,l)=(initial(i,l)-cs(i,l)*initial(i+1,l))/bs(i,l);

				}

			/* y direction*/

			copy=initial;

			#pragma omp parallel for
			for(long int k=0;k<xmax;++k){
					//Multiply vector
					initial(k,0)=betay(k,0)*copy(k,1)*(1.0+y(k,0))+(1.0-2.0*betay(k,0))*copy(k,0);
					initial(k,ylast)=betay(k,ylast)*copy(k,ylast-1)*(1.0-y(k,ylast))+(1.0-2.0*betay(k,ylast))*copy(k,ylast);

					for(long int l = 1; l < ylast; l++){
							initial(k,l)=betay(k,l)*copy(k,l-1)*(1.0-y(k,l))+(1.0-2.0*betay(k,l))*copy(k,l)+betay(k,l)*copy(k,l+1)*(1.0+y(k,l));
					}
				}


			#pragma omp parallel for
			for(long int k=0;k<xmax;++k){
					as(k,0)=0;
					bs(k,0)=2.0*betay(k,0)+1.0;
					cs(k,0)=-betay(k,0)*(1.0+y(k,0));

					for(long int i = 1; i < ylast ; i++){
						as(k,i)=-betay(k,i)*(1.0-y(k,i));
						bs(k,i)=2.0*betay(k,i)+1.0;
						cs(k,i)=-betay(k,i)*(1.0+y(k,i));
					}

					as(k,ylast)=-betay(k,ylast)*(1.0-y(k,ylast));
					bs(k,ylast)=2.0*betay(k,ylast)+1.0;
					cs(k,ylast)=0;

					//Thomas Algorithm

					for (long int i = 1; i <= ylast; i++){
						std::complex<complextype> n = as(k,i)/bs(k,i-1);
						bs(k,i)-=n*cs(k,i-1);
						initial(k,i)-=n*initial(k,i-1);
					}

					initial(k,ylast)/=bs(k,ylast);
					for (long int i = (ylast - 1); i >= 0; --i)
						initial(k,i)=(initial(k,i)-cs(k,i)*initial(k,i+1))/bs(k,i);

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
DiscreteFunction2D<std::complex<complextype>,gridtype>& WaveTimeEvolution2D<timetype,complextype>::newlevel_periodicFFT(timetype time, long int timesteps, DiscreteFunction2D<std::complex<complextype>,gridtype>& initial, DiscreteFunction2D<T,gridtype>& potential,std::vector<std::shared_ptr<DiscreteFunction2D<std::complex<complextype>,gridtype>>> otherslevels, T effectivemass){
	timetype methodtimestep=time/static_cast<timetype>(timesteps);

	initial/=initial.Modulus();


	DiscreteFunction2D<std::complex<complextype>,gridtype> potpropag = potentialpropagator(methodtimestep,potential);
	DiscreteFunction2D<std::complex<complextype>,gridtype> bloch = fftmulti(methodtimestep,initial.getGrid(),effectivemass);

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
DiscreteFunction2D<std::complex<complextype>,gridtype> WaveTimeEvolution2D<timetype,complextype>::fftmulti(timetype methodtimestep,Grid2D<gridtype> grid,T effmass){
	DiscreteFunction2D<std::complex<complextype>,gridtype> result(grid);
	long int sizehalfx = grid.getSizeX()/2;
	long int sizehalfy = grid.getSizeY()/2;

	for(long int k = 0; k < sizehalfx; k++){
		for(long int l = 0; l < sizehalfy; l++){
				result(k,l)=std::exp( std::complex<complextype>(0.0,-1.0)*((pow(2.0*Constant::pi*k/(grid.getxSuperiorLimit()-grid.getxInferiorLimit()),2)+pow(2.0*Constant::pi*l/(grid.getySuperiorLimit()-grid.getyInferiorLimit()),2)))*methodtimestep/(2*Constant::hbar*effmass) );
		}
		for(long int l = sizehalfy; l < grid.getSizeY(); l++){
				result(k,l)=std::exp( std::complex<complextype>(0.0,-1.0)*((pow(2.0*Constant::pi*k/(grid.getxSuperiorLimit()-grid.getxInferiorLimit()),2)+pow(2.0*Constant::pi*(l-2*sizehalfy)/(grid.getySuperiorLimit()-grid.getyInferiorLimit()),2)))*methodtimestep/(2*Constant::hbar*effmass) );
		}
	}
	for(long int k = sizehalfx; k < grid.getSizeX(); k++){
		for(long int l = 0; l < sizehalfy; l++){
				result(k,l)=std::exp( std::complex<complextype>(0.0,-1.0)*((pow(2.0*Constant::pi*(k-2*sizehalfx)/(grid.getxSuperiorLimit()-grid.getxInferiorLimit()),2)+pow(2.0*Constant::pi*l/(grid.getySuperiorLimit()-grid.getyInferiorLimit()),2)))*methodtimestep/(2*Constant::hbar*effmass) );
		}
		for(long int l = sizehalfy; l < grid.getSizeY(); l++){
				result(k,l)=std::exp( std::complex<complextype>(0.0,-1.0)*((pow(2.0*Constant::pi*(l-2*sizehalfy)/(grid.getySuperiorLimit()-grid.getyInferiorLimit()),2)+pow(2.0*Constant::pi*(k-2*sizehalfx)/(grid.getxSuperiorLimit()-grid.getxInferiorLimit()),2)))*methodtimestep/(2*Constant::hbar*effmass) );
		}
	}
	return result;
}
/*
template <typename timetype, typename complextype>
template<typename T, typename gridtype>
DiscreteFunction2D<std::complex<complextype>,gridtype>& WaveTimeEvolution2D<timetype,complextype>::kinect_propag_evolution_finite(timetype methodtimestep, DiscreteFunction2D<std::complex<complextype>,gridtype>& wavefunction, DiscreteFunction2D<T,gridtype>& effectivemass){
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

	DiscreteFunction2D<std::complex<complextype>,gridtype> betax(wavefunction.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> betay(wavefunction.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> betaz(wavefunction.getGrid());

	DiscreteFunction2D<std::complex<complextype>,gridtype> as(wavefunction.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> bs(wavefunction.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> cs(wavefunction.getGrid());

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

	DiscreteFunction2D<std::complex<complextype>,gridtype> copy(wavefunction);

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
std::complex<complextype> WaveTimeEvolution2D<timetype,complextype>::Energy(DiscreteFunction2D<std::complex<complextype>,gridtype>& wavefunction, DiscreteFunction2D<T,gridtype>& potential,DiscreteFunction2D<T,gridtype>& effectivemass){
	std::complex<complextype> energysign;

		if(carriertype==Carrier::ELECTRON)
			energysign=std::complex<complextype>(1.0,0);
		else
			energysign=std::complex<complextype>(-1.0,0);

	wavefunction.normalize(1.0);
	DiscreteFunction2D<std::complex<complextype>,gridtype> normalized(wavefunction);


	DiscreteFunction2D<std::complex<complextype>,gridtype> x(wavefunction.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> y(wavefunction.getGrid());


	DiscreteFunction2D<std::complex<complextype>,gridtype> xdois(wavefunction.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> ydois(wavefunction.getGrid());


	long int xmax=wavefunction.getGrid().getSizeX();
	long int ymax=wavefunction.getGrid().getSizeY();


	DiscreteFunction2D<std::complex<complextype>,gridtype> copy2(wavefunction.getGrid());

	gridtype dx = wavefunction.getGrid().getxIncrement();
	gridtype dy = wavefunction.getGrid().getyIncrement();


	for(long int l=0;l<ymax;++l){

			//edge: order h and h²
			x(0,l) = (wavefunction(1,l)-wavefunction(0,l))/dx;
			x(1,l) = (wavefunction(2,l)-wavefunction(0,l))/(2*dx);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int k = 2; k < (xmax-2); k++){
				x(k,l)=(wavefunction(k+1,l)*8.0-wavefunction(k+2,l)-wavefunction(k-1,l)*8.0+wavefunction(k-2,l))/(12.0*dx);
			}

			//edge: order h and h²
			x(xmax-2,l)= (wavefunction(xmax-1,l)-wavefunction(xmax-3,l))/(2.0*dx);
			x(xmax-1,l)= (wavefunction(xmax-1,l)-wavefunction(xmax-2,l))/dx;

	}

	for(long int k=0;k<xmax;++k){

			//edge: order h and h²
			y(k,0) = (wavefunction(k,1)-wavefunction(k,0))/dy;
			y(k,1) = (wavefunction(k,2)-wavefunction(k,0))/(2*dy);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int l = 2; l < (ymax-2); l++){
				y(k,l)=(wavefunction(k,l+1)*8.0-wavefunction(k,l+2)-wavefunction(k,l-1)*8.0+wavefunction(k,l-2))/(12.0*dy);
			}

			//edge: order h and h²
			y(k,ymax-2)= (wavefunction(k,ymax-1)-wavefunction(k,ymax-3))/(2.0*dy);
			y(k,ymax-1)= (wavefunction(k,ymax-1)-wavefunction(k,ymax-2))/dy;

	}



	#pragma omp parallel for
	for(long int n=0; n < wavefunction.getGrid().getSize(); ++n){
		x.getdata()[n]/=effectivemass.getdata()[n];
		y.getdata()[n]/=effectivemass.getdata()[n];
	}

	for(long int l=0;l<ymax;++l){

			//edge: order h and h²
			xdois(0,l) = (x(1,l)-x(0,l))/dx;
			xdois(1,l) = (x(2,l)-x(0,l))/(2*dx);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int k = 2; k < (xmax-2); k++){
				xdois(k,l)=(x(k+1,l)*8.0-x(k+2,l)-x(k-1,l)*8.0+x(k-2,l))/(12.0*dx);
			}

			//edge: order h and h²
			xdois(xmax-2,l)= (x(xmax-1,l)-x(xmax-3,l))/(2.0*dx);
			xdois(xmax-1,l)= (x(xmax-1,l)-x(xmax-2,l))/dx;

	}

	for(long int k=0;k<xmax;++k){

			//edge: order h and h²
			ydois(k,0) = (y(k,1)-y(k,0))/dy;
			ydois(k,1) = (y(k,2)-y(k,0))/(2*dy);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int l = 2; l < (ymax-2); l++){
				ydois(k,l)=(y(k,l+1)*8.0-y(k,l+2)-y(k,l-1)*8.0+y(k,l-2))/(12.0*dy);
			}

			//edge: order h and h²
			ydois(k,ymax-2)= (y(k,ymax-1)-y(k,ymax-3))/(2.0*dy);
			ydois(k,ymax-1)= (y(k,ymax-1)-y(k,ymax-2))/dy;

	}


	#pragma omp parallel for
	for(long int n=0; n < wavefunction.getGrid().getSize(); ++n){
		copy2.getdata()[n]=xdois.getdata()[n]+ydois.getdata()[n];
	}
	copy2*=(-Constant::hbar*Constant::hbar/(2.0));

	DiscreteFunction2D<std::complex<complextype>,gridtype> copy(wavefunction);
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
std::complex<complextype> WaveTimeEvolution2D<timetype,complextype>::Energy(DiscreteFunction2D<std::complex<complextype>,gridtype>& wavefunction, DiscreteFunction2D<T,gridtype>& potential,DiscreteFunction2D<T,gridtype>& effectivemasst,DiscreteFunction2D<T,gridtype>& effectivemassl){
	std::complex<complextype> energysign;

		if(carriertype==Carrier::ELECTRON)
			energysign=std::complex<complextype>(1.0,0);
		else
			energysign=std::complex<complextype>(-1.0,0);

	wavefunction.normalize(1.0);
	DiscreteFunction2D<std::complex<complextype>,gridtype> normalized(wavefunction);


	DiscreteFunction2D<std::complex<complextype>,gridtype> x(wavefunction.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> y(wavefunction.getGrid());


	DiscreteFunction2D<std::complex<complextype>,gridtype> xdois(wavefunction.getGrid());
	DiscreteFunction2D<std::complex<complextype>,gridtype> ydois(wavefunction.getGrid());


	long int xmax=wavefunction.getGrid().getSizeX();
	long int ymax=wavefunction.getGrid().getSizeY();


	DiscreteFunction2D<std::complex<complextype>,gridtype> copy2(wavefunction.getGrid());

	gridtype dx = wavefunction.getGrid().getxIncrement();
	gridtype dy = wavefunction.getGrid().getyIncrement();


	for(long int l=0;l<ymax;++l){

			//edge: order h and h²
			x(0,l) = (wavefunction(1,l)-wavefunction(0,l))/dx;
			x(1,l) = (wavefunction(2,l)-wavefunction(0,l))/(2*dx);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int k = 2; k < (xmax-2); k++){
				x(k,l)=(wavefunction(k+1,l)*8.0-wavefunction(k+2,l)-wavefunction(k-1,l)*8.0+wavefunction(k-2,l))/(12.0*dx);
			}

			//edge: order h and h²
			x(xmax-2,l)= (wavefunction(xmax-1,l)-wavefunction(xmax-3,l))/(2.0*dx);
			x(xmax-1,l)= (wavefunction(xmax-1,l)-wavefunction(xmax-2,l))/dx;

	}

	for(long int k=0;k<xmax;++k){

			//edge: order h and h²
			y(k,0) = (wavefunction(k,1)-wavefunction(k,0))/dy;
			y(k,1) = (wavefunction(k,2)-wavefunction(k,0))/(2*dy);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int l = 2; l < (ymax-2); l++){
				y(k,l)=(wavefunction(k,l+1)*8.0-wavefunction(k,l+2)-wavefunction(k,l-1)*8.0+wavefunction(k,l-2))/(12.0*dy);
			}

			//edge: order h and h²
			y(k,ymax-2)= (wavefunction(k,ymax-1)-wavefunction(k,ymax-3))/(2.0*dy);
			y(k,ymax-1)= (wavefunction(k,ymax-1)-wavefunction(k,ymax-2))/dy;

	}



	#pragma omp parallel for
	for(long int n=0; n < wavefunction.getGrid().getSize(); ++n){
		x.getdata()[n]/=effectivemasst.getdata()[n];
		y.getdata()[n]/=effectivemassl.getdata()[n];
	}

	for(long int l=0;l<ymax;++l){

			//edge: order h and h²
			xdois(0,l) = (x(1,l)-x(0,l))/dx;
			xdois(1,l) = (x(2,l)-x(0,l))/(2*dx);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int k = 2; k < (xmax-2); k++){
				xdois(k,l)=(x(k+1,l)*8.0-x(k+2,l)-x(k-1,l)*8.0+x(k-2,l))/(12.0*dx);
			}

			//edge: order h and h²
			xdois(xmax-2,l)= (x(xmax-1,l)-x(xmax-3,l))/(2.0*dx);
			xdois(xmax-1,l)= (x(xmax-1,l)-x(xmax-2,l))/dx;

	}

	for(long int k=0;k<xmax;++k){

			//edge: order h and h²
			ydois(k,0) = (y(k,1)-y(k,0))/dy;
			ydois(k,1) = (y(k,2)-y(k,0))/(2*dy);

			//central difference order h⁴
			#pragma omp parallel for
			for(long int l = 2; l < (ymax-2); l++){
				ydois(k,l)=(y(k,l+1)*8.0-y(k,l+2)-y(k,l-1)*8.0+y(k,l-2))/(12.0*dy);
			}

			//edge: order h and h²
			ydois(k,ymax-2)= (y(k,ymax-1)-y(k,ymax-3))/(2.0*dy);
			ydois(k,ymax-1)= (y(k,ymax-1)-y(k,ymax-2))/dy;

	}


	#pragma omp parallel for
	for(long int n=0; n < wavefunction.getGrid().getSize(); ++n){
		copy2.getdata()[n]=xdois.getdata()[n]+ydois.getdata()[n];
	}
	copy2*=(-Constant::hbar*Constant::hbar/(2.0));

	DiscreteFunction2D<std::complex<complextype>,gridtype> copy(wavefunction);
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
