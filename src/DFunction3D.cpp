/*
 * DFunction3D.cpp
 *
 *  Created on: Nov 7, 2014
 *      Author: marcel
 */

#include "DFunction3D.hpp"
#include <fftw3.h>

#ifndef DFUNCTION3D_CPP_
#define DFUNCTION3D_CPP_

namespace epital {


template <typename T,typename gridtype>
DiscreteFunction3D<T,gridtype>::DiscreteFunction3D(Grid3D<gridtype> grid): grid(grid)  {
	data = new T[grid.getSize()];
}


template <typename T,typename gridtype>
DiscreteFunction3D<T,gridtype>::~DiscreteFunction3D() {
	delete [] data;
}

template <typename T,typename gridtype>
DiscreteFunction3D<T,gridtype>::DiscreteFunction3D(const DiscreteFunction3D<T,gridtype>& tocopy): grid(tocopy.getGrid()){
	//Create function data
	data = new T[grid.getSize()];

	for(long int size=0; size<grid.getSize(); ++size){
		data[size]= tocopy.data[size];
	}
}

template <typename T,typename gridtype>
DiscreteFunction3D<T,gridtype>::DiscreteFunction3D(DiscreteFunction3D<T,gridtype>& tocopy): grid(tocopy.getGrid()){
	//Create function data
	data = new T[grid.getSize()];
	T* tocopydata = tocopy.getdata();

	for(long int size=0; size<grid.getSize(); ++size){
		data[size]= tocopydata[size];
	}
}

template <typename T,typename gridtype>
DiscreteFunction3D<T,gridtype>::DiscreteFunction3D(std::shared_ptr<DiscreteFunction3D<T,gridtype>> tocopy): grid(tocopy->getGrid()){
	//Create function data
	data = new T[grid->getSize()];
	T* tocopydata = tocopy->getdata();

	for(long int size=0; size<grid.getSize(); ++size){
		data[size]= tocopydata[size];
	}
}

template <typename T,typename gridtype>
DiscreteFunction3D<T,gridtype>::DiscreteFunction3D(DiscreteFunction3D<T,gridtype>&& tomove): grid(tomove.getGrid()){
	data=tomove.getdata();
	tomove.data=nullptr;
}

template <typename T,typename gridtype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator=(DiscreteFunction3D<T,gridtype>&& tomove){
	if(this == &tomove)
		return *this;
	if(grid.getSize()!=tomove.getGrid().getSize() || grid.getSizeX()!=tomove.getGrid().getSizeX() || grid.getSizeY()!=tomove.getGrid().getSizeY() ){
		throw std::runtime_error("Incompatible grid at = operation in DiscreteFunction3D");
	};

	delete [] data;
	data = nullptr;
	data = tomove.getdata();
	tomove.data=nullptr;
	return *this;
}

template <typename T,typename gridtype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator=(DiscreteFunction3D<T,gridtype>& tocopy){
	if(this == &tocopy)
		return *this;

	if(grid.getSize()!=tocopy.getGrid().getSize() || grid.getSizeX()!=tocopy.getGrid().getSizeX() || grid.getSizeY()!=tocopy.getGrid().getSizeY()){
		throw std::runtime_error("Incompatible grid at = operation in DiscreteFunction");
	};

	T* tocopydata = tocopy.getdata();

	for(long int size=0; size<grid.getSize(); ++size){
		data[size]= tocopydata[size];
	}

	grid=tocopy.getGrid();
	return *this;
}


template <typename T,typename gridtype>
inline T& DiscreteFunction3D<T,gridtype>::operator()(Grid::Position3D position) {
	return data[position[0]+position[1]*grid.getSizeX()+position[2]*grid.getSizeX()*grid.getSizeY()];
}

template <typename T,typename gridtype>
inline T& DiscreteFunction3D<T,gridtype>::operator()(long int x, long int y, long int z) {
	return data[x+y*grid.getSizeX()+z*grid.getSizeX()*grid.getSizeY()];
}



template <typename T,typename gridtype>
inline T* DiscreteFunction3D<T,gridtype>::getdata(){
	return data;
}


template <typename T,typename gridtype>
const Grid3D<gridtype>& DiscreteFunction3D<T,gridtype>::getGrid() const{
    return grid;
}

template <typename T,typename gridtype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::FFT(){
	fftw_complex *in;
	fftw_plan p;
	fftw_init_threads();

	fftw_plan_with_nthreads(4);
	in = reinterpret_cast<fftw_complex*>(this->getdata());
    p = fftw_plan_dft_3d(grid.getSizeX(),grid.getSizeY(),grid.getSizeZ(), in, in, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p); /* repeat as needed */

    fftw_destroy_plan(p);
    fftw_cleanup_threads();
    return *this;
}


template <typename T,typename gridtype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::iFFT(){
	fftw_complex *in;
	fftw_plan p;
	fftw_init_threads();

	fftw_plan_with_nthreads(4);
	in = reinterpret_cast<fftw_complex*>(this->getdata());
    p = fftw_plan_dft_3d(grid.getSizeX(),grid.getSizeY(),grid.getSizeZ(), in, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(p); /* repeat as needed */

    fftw_destroy_plan(p);
    fftw_cleanup_threads();
    return *this;
}

template <typename T,typename gridtype>
T DiscreteFunction3D<T,gridtype>::Modulus(){
	long int xmax=grid.getSizeX()-1;
	long int ymax=grid.getSizeY()-1;
	long int zmax=grid.getSizeZ()-1;

	T sum=0;

	/*vertex sum*/
	sum+=(*this)(0,0,0)*conj((*this)(0,0,0))/8.0;
	sum+=(*this)(0,ymax,0)*conj((*this)(0,ymax,0))/8.0;
	sum+=(*this)(0,ymax,zmax)*conj((*this)(0,ymax,zmax))/8.0;
	sum+=(*this)(0,0,zmax)*conj((*this)(0,0,zmax))/8.0;
	sum+=(*this)(xmax,0,0)*conj((*this)(xmax,0,0))/8.0;
	sum+=(*this)(xmax,ymax,0)*conj((*this)(xmax,ymax,0))/8.0;
	sum+=(*this)(xmax,ymax,zmax)*conj((*this)(xmax,ymax,zmax))/8.0;
	sum+=(*this)(xmax,0,zmax)*conj((*this)(xmax,0,zmax))/8.0;


	/*inside + edges sum*/
	#pragma omp parallel shared(sum)
	{
	T partialsum = 0.0;

	#pragma omp for
	for(long int k=1; k<(xmax); k++){
		partialsum+=(*this)(k,0,0)*conj((*this)(k,0,0))/4.0;
		partialsum+=(*this)(k,ymax,0)*conj((*this)(k,ymax,0))/4.0;
		partialsum+=(*this)(k,ymax,zmax)*conj((*this)(k,ymax,zmax))/4.0;
		partialsum+=(*this)(k,0,zmax)*conj((*this)(k,0,zmax))/4.0;
		for(long int l=1; l<(ymax); l++){
			for(long int m=1; m<(zmax); m++){
				partialsum+=(*this)(k,l,m)*conj((*this)(k,l,m));
			}
		}
	}

	#pragma omp critical
	{
		sum+=partialsum;
	}

	}

	/*faces sum*/
	for(long int l=1; l<(ymax); l++){
		for(long int m=1; m<(zmax); m++){
			sum+=(*this)(0,l,m)*conj((*this)(0,l,m))/2.0;
			sum+=(*this)(xmax,l,m)*conj((*this)(xmax,l,m))/2.0;
		}
	}

	for(long int k=1; k<(xmax); k++){
		for(long int m=1; m<(zmax); m++){
			sum+=(*this)(k,ymax,m)*conj((*this)(k,ymax,m))/2.0;
			sum+=(*this)(k,0,m)*conj((*this)(k,0,m))/2.0;
		}
	}

	for(long int k=1; k<(xmax); k++){
		for(long int l=1; l<(ymax); l++){
			sum+=(*this)(k,l,0)*conj((*this)(k,l,0))/2.0;
			sum+=(*this)(k,l,zmax)*conj((*this)(k,l,zmax))/2.0;
		}
	}


	/*edges sum*/
	for(long int l=1; l<(ymax); l++){
		sum+=(*this)(0,l,0)*conj((*this)(0,l,0))/4.0;
		sum+=(*this)(0,l,zmax)*conj((*this)(0,l,zmax))/4.0;
		sum+=(*this)(xmax,l,0)*conj((*this)(xmax,l,0))/4.0;
		sum+=(*this)(xmax,l,zmax)*conj((*this)(xmax,l,zmax))/4.0;
	}
	for(long int m=1; m<(zmax); m++){
		sum+=(*this)(0,0,m)*conj((*this)(0,0,m))/4.0;
		sum+=(*this)(xmax,0,m)*conj((*this)(xmax,0,m))/4.0;
		sum+=(*this)(xmax,ymax,m)*conj((*this)(xmax,ymax,m))/4.0;
		sum+=(*this)(0,ymax,m)*conj((*this)(0,ymax,m))/4.0;
	}



	sum*=(xmax*ymax*zmax*grid.getxIncrement()*grid.getyIncrement()*grid.getzIncrement());
	sum=sqrt(sum);

	return sum;

}


template <typename T,typename gridtype>
void DiscreteFunction3D<T,gridtype>::saveHDF5(std::string filename){
	long int xmax=grid.getSizeX();
	long int ymax=grid.getSizeY();
	long int zmax=grid.getSizeZ();

	arma::Cube<T> cube(xmax,ymax,zmax);

	#pragma omp for
	for(long int k=0; k<xmax; k++){
		for(long int l=0; l<ymax; l++){
			for(long int m=0; m<zmax; m++){
				cube(k,l,m)=(*this)(k,l,m);
			}
		}
	}

	cube.save(filename,arma::hdf5_binary);
}

template <typename T,typename gridtype>
void DiscreteFunction3D<T,gridtype>::fillGaussian(T sigma, T norm){
	for(long int k=0; k<grid.getSizeX();k++)
		for(long int l=0; l<grid.getSizeY();l++)
			for(long int m=0; m<grid.getSizeZ();m++){
				T arg(-((k-grid.getSizeX()/2)*(k-grid.getSizeX()/2.0)+(l-grid.getSizeY()/2)*(l-grid.getSizeY()/2.0)+(m-grid.getSizeZ()/2)*(m-grid.getSizeZ()/2.0))/(2.0*sigma));
				(*this)(k,l,m)=(exp(arg)*norm/2.0)+(exp(arg)*std::complex<double>(0.0,1.0)*norm/2.0);
			}
}

template <typename T,typename gridtype>
void DiscreteFunction3D<T,gridtype>::fillGaussian2(T sigma, T norm){
	for(long int k=0; k<grid.getSizeX();k++)
		for(long int l=0; l<grid.getSizeY();l++)
			for(long int m=0; m<grid.getSizeZ();m++){
				T arg(-((k-grid.getSizeX()/2)*(k-grid.getSizeX()/2.0)+(l-grid.getSizeY()/2)*(l-grid.getSizeY()/2.0)+(m-grid.getSizeZ()/2)*(m-grid.getSizeZ()/2.0))/(2.0*sigma));
				(*this)(k,l,m)=(k-grid.getSizeX()/2.0)*(l-grid.getSizeY()/2.0)*(m-grid.getSizeZ()/2.0)*(exp(arg)*norm/2.0)+(exp(arg)*std::complex<double>(0.0,1.0)*norm/2.0);
			}
}


template <typename T,typename gridtype>
void DiscreteFunction3D<T,gridtype>::normalize(T modulus){
	T actualmod = Modulus();
	(*this)/=(actualmod/modulus);
}



template <typename T,typename gridtype>
void DiscreteFunction3D<T,gridtype>::Gradient(DiscreteFunction3D<T,gridtype>& x, DiscreteFunction3D<T,gridtype>& y,DiscreteFunction3D<T,gridtype>& z){
	long int xmax=grid.getSizeX();
	long int ymax=grid.getSizeY();
	long int zmax=grid.getSizeZ();


	for(long int l=0;l<ymax;++l){
		for(long int m=0;m<zmax;++m){
			//edge: order h and h²
			x(0,l,m) = ((*this)(1,l,m)-(*this)(0,l,m))/grid.getxIncrement();
			x(1,l,m) = ((*this)(2,l,m)-(*this)(0,l,m))/(2*grid.getxIncrement());

			//central difference order h⁴
			#pragma omp parallel for
			for(long int k = 2; k < (xmax-2); k++){
				x(k,l,m)=((*this)(k+1,l,m)*8.0-(*this)(k+2,l,m)-(*this)(k-1,l,m)*8.0+(*this)(k-2,l,m))/(12.0*grid.getxIncrement());
			}

			//edge: order h and h²
			x(xmax-2,l,m)= ((*this)(xmax-1,l,m)-(*this)(xmax-3,l,m))/(2.0*grid.getxIncrement());
			x(xmax-1,l,m)= ((*this)(xmax-1,l,m)-(*this)(xmax-2,l,m))/grid.getxIncrement();
		}
	}

	for(long int k=0;k<xmax;++k){
		for(long int m=0;m<zmax;++m){
			//edge: order h and h²
			y(k,0,m) = ((*this)(k,1,m)-(*this)(k,0,m))/grid.getyIncrement();
			y(k,1,m) = ((*this)(k,2,m)-(*this)(k,0,m))/(2*grid.getyIncrement());

			//central difference order h⁴
			#pragma omp parallel for
			for(long int l = 2; l < (ymax-2); l++){
				y(k,l,m)=((*this)(k,l+1,m)*8.0-(*this)(k,l+2,m)-(*this)(k,l-1,m)*8.0+(*this)(k,l-2,m))/(12.0*grid.getyIncrement());
			}

			//edge: order h and h²
			y(k,ymax-2,m)= ((*this)(k,ymax-1,m)-(*this)(k,ymax-3,m))/(2.0*grid.getyIncrement());
			y(k,ymax-1,m)= ((*this)(k,ymax-1,m)-(*this)(k,ymax-2,m))/grid.getyIncrement();
		}
	}

	for(long int k=0;k<xmax;++k){
		for(long int l=0;l<ymax;++l){
			//edge: order h and h²
			z(k,l,0) = ((*this)(k,l,1)-(*this)(k,l,0))/grid.getzIncrement();
			z(k,l,1) = ((*this)(k,l,2)-(*this)(k,l,0))/(2*grid.getzIncrement());

			//central difference order h⁴
			#pragma omp parallel for
			for(long int m = 2; m < (zmax-2); m++){
				z(k,l,m)=((*this)(k,l,m+1)*8.0-(*this)(k,l,m+2)-(*this)(k,l,m-1)*8.0+(*this)(k,l,m-2))/(12.0*grid.getzIncrement());
			}

			//edge: order h and h²
			z(k,l,zmax-2)= ((*this)(k,l,zmax-1)-(*this)(k,l,zmax-3))/(2.0*grid.getzIncrement());
			z(k,l,zmax-1)= ((*this)(k,l,zmax-1)-(*this)(k,l,zmax-2))/grid.getzIncrement();
		}
	}

}



template <>
template <typename TC>
DiscreteFunction3D<double,double> DiscreteFunction3D<std::complex<double>,double>::getDensity(){
	DiscreteFunction3D<double,double> density(grid);
	for(long int k=0; k<grid.getSizeX();k++)
		for(long int l=0; l<grid.getSizeY();l++)
			for(long int m=0; m<grid.getSizeZ();m++){
				density(k,l,m)=std::real((*this)(k,l,m)*conj((*this)(k,l,m)));
			}
	return density;
}

template <typename T,typename gridtype>
T DiscreteFunction3D<T,gridtype>::Dot(DiscreteFunction3D<T,gridtype>& rhs){
	long int xmax=grid.getSizeX()-1;
	long int ymax=grid.getSizeY()-1;
	long int zmax=grid.getSizeZ()-1;

	T sum=0;

	/*vertex sum*/
	sum+=conj((*this)(0,0,0))*rhs(0,0,0)/8.0;
	sum+=conj((*this)(0,ymax,0))*rhs(0,ymax,0)/8.0;
	sum+=conj((*this)(0,ymax,zmax))*rhs(0,ymax,zmax)/8.0;
	sum+=conj((*this)(0,0,zmax))*rhs(0,0,zmax)/8.0;
	sum+=conj((*this)(xmax,0,0))*rhs(xmax,0,0)/8.0;
	sum+=conj((*this)(xmax,ymax,0))*rhs(xmax,ymax,0)/8.0;
	sum+=conj((*this)(xmax,ymax,zmax))*rhs(xmax,ymax,zmax)/8.0;
	sum+=conj((*this)(xmax,0,zmax))*rhs(xmax,0,zmax)/8.0;


	/*inside + edges sum*/
	#pragma omp parallel shared(sum)
	{
	T partialsum = 0.0;

	#pragma omp for
	for(long int k=1; k<(xmax); k++){
		partialsum+=conj((*this)(k,0,0))*rhs(k,0,0)/4.0;
		partialsum+=conj((*this)(k,ymax,0))*rhs(k,ymax,0)/4.0;
		partialsum+=conj((*this)(k,ymax,zmax))*rhs(k,ymax,zmax)/4.0;
		partialsum+=conj((*this)(k,0,zmax))*rhs(k,0,zmax)/4.0;
		for(long int l=1; l<(ymax); l++){
			for(long int m=1; m<(zmax); m++){
				partialsum+=conj((*this)(k,l,m))*rhs(k,l,m);
			}
		}
	}

	#pragma omp critical
	{
		sum+=partialsum;
	}

	}

	/*faces sum*/
	for(long int l=1; l<(ymax); l++){
		for(long int m=1; m<(zmax); m++){
			sum+=conj((*this)(0,l,m))*rhs(0,l,m)/2.0;
			sum+=conj((*this)(xmax,l,m))*rhs(xmax,l,m)/2.0;
		}
	}

	for(long int k=1; k<(xmax); k++){
		for(long int m=1; m<(zmax); m++){
			sum+=conj((*this)(k,ymax,m))*rhs(k,ymax,m)/2.0;
			sum+=conj((*this)(k,0,m))*rhs(k,0,m)/2.0;
		}
	}

	for(long int k=1; k<(xmax); k++){
		for(long int l=1; l<(ymax); l++){
			sum+=conj((*this)(k,l,0))*rhs(k,l,0)/2.0;
			sum+=conj((*this)(k,l,zmax))*rhs(k,l,zmax)/2.0;
		}
	}


	/*edges sum*/
	for(long int l=1; l<(ymax); l++){
		sum+=conj((*this)(0,l,0))*rhs(0,l,0)/4.0;
		sum+=conj((*this)(0,l,zmax))*rhs(0,l,zmax)/4.0;
		sum+=conj((*this)(xmax,l,0))*rhs(xmax,l,0)/4.0;
		sum+=conj((*this)(xmax,l,zmax))*rhs(xmax,l,zmax)/4.0;
	}
	for(long int m=1; m<(zmax); m++){
		sum+=conj((*this)(0,0,m))*rhs(0,0,m)/4.0;
		sum+=conj((*this)(xmax,0,m))*rhs(xmax,0,m)/4.0;
		sum+=conj((*this)(xmax,ymax,m))*rhs(xmax,ymax,m)/4.0;
		sum+=conj((*this)(0,ymax,m))*rhs(0,ymax,m)/4.0;
	}


	sum*=(xmax*ymax*zmax*grid.getxIncrement()*grid.getyIncrement()*grid.getzIncrement());

	return sum;
}


template <typename T,typename gridtype>
void DiscreteFunction3D<T,gridtype>::smooth_mean(){
	long int xmax=grid.getSizeX()-1;
	long int ymax=grid.getSizeY()-1;
	long int zmax=grid.getSizeZ()-1;

	/*vertex sum*/
	(*this)(0,0,0)=((*this)(0,0,0)+(*this)(1,0,0)+(*this)(0,1,0)+(*this)(0,0,1))/4.0;
	(*this)(0,0,zmax)=((*this)(0,0,zmax)+(*this)(1,0,zmax)+(*this)(0,1,zmax)+(*this)(0,0,zmax-1))/4.0;
	(*this)(0,ymax,0)=((*this)(0,ymax,0)+(*this)(1,ymax,0)+(*this)(0,ymax-1,0)+(*this)(0,ymax,1))/4.0;
	(*this)(0,ymax,zmax)=((*this)(0,ymax,zmax)+(*this)(1,ymax,zmax)+(*this)(0,ymax-1,zmax)+(*this)(0,ymax,zmax-1))/4.0;
	(*this)(xmax,0,0)=((*this)(xmax,0,0)+(*this)(xmax-1,0,0)+(*this)(xmax,1,0)+(*this)(xmax,0,1))/4.0;
	(*this)(xmax,ymax,0)=((*this)(xmax,ymax,0)+(*this)(xmax-1,ymax,0)+(*this)(xmax,ymax-1,0)+(*this)(xmax,ymax,1))/4.0;
	(*this)(xmax,ymax,zmax)=((*this)(xmax,ymax,zmax)+(*this)(xmax-1,ymax,zmax)+(*this)(xmax,ymax-1,zmax)+(*this)(xmax,ymax,zmax-1))/4.0;
	(*this)(xmax,0,zmax)=((*this)(xmax,0,zmax)+(*this)(xmax-1,0,zmax)+(*this)(xmax,1,zmax)+(*this)(xmax,0,zmax-1))/4.0;




	#pragma omp for
	for(long int k=1; k<(xmax); k++){
		for(long int l=1; l<(ymax); l++){
			for(long int m=1; m<(zmax); m++){
				(*this)(k,l,m)=((*this)(k,l,m)+(*this)(k+1,l,m)+(*this)(k-1,l,m)+(*this)(k,l+1,m)+(*this)(k,l-1,m)+(*this)(k,l,m+1)+(*this)(k,l,m-1))/7.0;
			}
		}
	}


	/*faces sum*/
	for(long int l=1; l<(ymax); l++){
		for(long int m=1; m<(zmax); m++){
			(*this)(0,l,m)=((*this)(0,l,m)+(*this)(1,l,m)+(*this)(0,l+1,m)+(*this)(0,l-1,m)+(*this)(0,l,m+1)+(*this)(0,l,m-1))/6.0;
			(*this)(xmax,l,m)=((*this)(xmax,l,m)+(*this)(xmax-1,l,m)+(*this)(xmax,l+1,m)+(*this)(xmax,l-1,m)+(*this)(xmax,l,m+1)+(*this)(xmax,l,m-1))/6.0;
		}
	}

	for(long int k=1; k<(xmax); k++){
		for(long int m=1; m<(zmax); m++){
			(*this)(k,0,m)=((*this)(k,0,m)+(*this)(k,1,m)+(*this)(k+1,0,m)+(*this)(k-1,0,m)+(*this)(k,0,m+1)+(*this)(k,0,m-1))/6.0;
			(*this)(k,ymax,m)=((*this)(k,ymax,m)+(*this)(k+1,ymax,m)+(*this)(k,ymax-1,m)+(*this)(k-1,ymax,m)+(*this)(k,ymax,m+1)+(*this)(k,ymax,m-1))/6.0;

		}
	}

	for(long int k=1; k<(xmax); k++){
		for(long int l=1; l<(ymax); l++){
			(*this)(k,l,0)=((*this)(k,l,0)+(*this)(k,l,1)+(*this)(k-1,l,0)+(*this)(k+1,l,0)+(*this)(k,l-1,0)+(*this)(k,l+1,0))/6.0;
			(*this)(k,l,zmax)=((*this)(k,l,zmax)+(*this)(k,l,zmax-1)+(*this)(k-1,l,zmax)+(*this)(k+1,l,zmax)+(*this)(k,l-1,zmax)+(*this)(k,l+1,zmax))/6.0;
		}
	}


	/*edges sum*/
	for(long int l=1; l<(ymax); l++){
		(*this)(0,l,0)=((*this)(0,l,0)+(*this)(0,l+1,0)+(*this)(0,l-1,0)+(*this)(1,l,0)+(*this)(0,l,1))/5.0;
		(*this)(0,l,zmax)=((*this)(0,l,zmax)+(*this)(0,l+1,zmax)+(*this)(0,l-1,zmax)+(*this)(1,l,zmax)+(*this)(0,l,zmax-1))/5.0;
		(*this)(xmax,l,0)=((*this)(xmax,l,0)+(*this)(xmax,l+1,0)+(*this)(xmax,l-1,0)+(*this)(xmax-1,l,0)+(*this)(xmax,l,1))/5.0;
		(*this)(xmax,l,zmax)=((*this)(xmax,l,zmax)+(*this)(xmax,l+1,zmax)+(*this)(xmax,l-1,zmax)+(*this)(xmax-1,l,zmax)+(*this)(xmax,l,zmax-1))/5.0;
	}

	for(long int m=1; m<(zmax); m++){
		(*this)(0,0,m)=((*this)(0,0,m)+(*this)(0,0,m-1)+(*this)(0,0,m+1)+(*this)(1,0,m)+(*this)(0,1,m))/5.0;
		(*this)(xmax,0,m)=((*this)(xmax,0,m)+(*this)(xmax,0,m-1)+(*this)(xmax,0,m+1)+(*this)(xmax-1,0,m)+(*this)(xmax,1,m))/5.0;
		(*this)(xmax,ymax,m)=((*this)(xmax,ymax,m)+(*this)(xmax,ymax,m-1)+(*this)(xmax,ymax,m+1)+(*this)(xmax-1,ymax,m)+(*this)(xmax,ymax-1,m))/5.0;
		(*this)(0,ymax,m)=((*this)(0,ymax,m)+(*this)(0,ymax,m-1)+(*this)(0,ymax,m+1)+(*this)(1,ymax,m)+(*this)(0,ymax-1,m))/5.0;
	}

	for(long int k=1; k<(xmax); k++){
		(*this)(k,0,0)=((*this)(k,0,0)+(*this)(k-1,0,0)+(*this)(k+1,0,0)+(*this)(k,0,1)+(*this)(k,1,0))/5.0;
		(*this)(k,ymax,0)=((*this)(k,ymax,0)+(*this)(k-1,ymax,0)+(*this)(k+1,ymax,0)+(*this)(k,ymax,1)+(*this)(k,ymax-1,0))/5.0;
		(*this)(k,ymax,zmax)=((*this)(k,ymax,zmax)+(*this)(k-1,ymax,zmax)+(*this)(k+1,ymax,zmax)+(*this)(k,ymax,zmax-1)+(*this)(k,ymax-1,zmax))/5.0;
		(*this)(k,0,zmax)=((*this)(k,0,zmax)+(*this)(k-1,0,zmax)+(*this)(k+1,0,zmax)+(*this)(k,0,zmax-1)+(*this)(k,1,zmax))/5.0;
	}

}

template <typename T,typename gridtype>
T DiscreteFunction3D<T,gridtype>::Dot(DiscreteFunction3D<T,gridtype>&& rhs){
	long int xmax=grid.getSizeX()-1;
	long int ymax=grid.getSizeY()-1;
	long int zmax=grid.getSizeZ()-1;

	T sum=0;

	/*vertex sum*/
	sum+=conj((*this)(0,0,0))*rhs(0,0,0)/8.0;
	sum+=conj((*this)(0,ymax,0))*rhs(0,ymax,0)/8.0;
	sum+=conj((*this)(0,ymax,zmax))*rhs(0,ymax,zmax)/8.0;
	sum+=conj((*this)(0,0,zmax))*rhs(0,0,zmax)/8.0;
	sum+=conj((*this)(xmax,0,0))*rhs(xmax,0,0)/8.0;
	sum+=conj((*this)(xmax,ymax,0))*rhs(xmax,ymax,0)/8.0;
	sum+=conj((*this)(xmax,ymax,zmax))*rhs(xmax,ymax,zmax)/8.0;
	sum+=conj((*this)(xmax,0,zmax))*rhs(xmax,0,zmax)/8.0;


	/*inside + edges sum*/
	#pragma omp parallel shared(sum)
	{
	T partialsum = 0.0;

	#pragma omp for
	for(long int k=1; k<(xmax); k++){
		partialsum+=conj((*this)(k,0,0))*rhs(k,0,0)/4.0;
		partialsum+=conj((*this)(k,ymax,0))*rhs(k,ymax,0)/4.0;
		partialsum+=conj((*this)(k,ymax,zmax))*rhs(k,ymax,zmax)/4.0;
		partialsum+=conj((*this)(k,0,zmax))*rhs(k,0,zmax)/4.0;
		for(long int l=1; l<(ymax); l++){
			for(long int m=1; m<(zmax); m++){
				partialsum+=conj((*this)(k,l,m))*rhs(k,l,m);
			}
		}
	}

	#pragma omp critical
	{
		sum+=partialsum;
	}

	}

	/*faces sum*/
	for(long int l=1; l<(ymax); l++){
		for(long int m=1; m<(zmax); m++){
			sum+=conj((*this)(0,l,m))*rhs(0,l,m)/2.0;
			sum+=conj((*this)(xmax,l,m))*rhs(xmax,l,m)/2.0;
		}
	}

	for(long int k=1; k<(xmax); k++){
		for(long int m=1; m<(zmax); m++){
			sum+=conj((*this)(k,ymax,m))*rhs(k,ymax,m)/2.0;
			sum+=conj((*this)(k,0,m))*rhs(k,0,m)/2.0;
		}
	}

	for(long int k=1; k<(xmax); k++){
		for(long int l=1; l<(zmax); l++){
			sum+=conj((*this)(k,l,0))*rhs(k,l,0)/2.0;
			sum+=conj((*this)(k,l,zmax))*rhs(k,l,zmax)/2.0;
		}
	}


	/*edges sum*/
	for(long int l=1; l<(ymax); l++){
		sum+=conj((*this)(0,l,0))*rhs(0,l,0)/4.0;
		sum+=conj((*this)(0,l,zmax))*rhs(0,l,zmax)/4.0;
		sum+=conj((*this)(xmax,l,0))*rhs(xmax,l,0)/4.0;
		sum+=conj((*this)(xmax,l,zmax))*rhs(xmax,l,zmax)/4.0;
	}
	for(long int m=1; m<(zmax); m++){
		sum+=conj((*this)(0,0,m))*rhs(0,0,m)/4.0;
		sum+=conj((*this)(xmax,0,m))*rhs(xmax,0,m)/4.0;
		sum+=conj((*this)(xmax,ymax,m))*rhs(xmax,ymax,m)/4.0;
		sum+=conj((*this)(0,ymax,m))*rhs(0,ymax,m)/4.0;
	}


	sum*=(xmax*ymax*zmax*grid.getxIncrement()*grid.getyIncrement()*grid.getzIncrement());


	return sum;
}



template<typename ayType, typename axType>
const DiscreteFunction3D<ayType,axType> operator+(const DiscreteFunction3D<ayType,axType>& lhs,const DiscreteFunction3D<ayType,axType>& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize() || lhs.getGrid().getSizeX()!=rhs.getGrid().getSizeX() || lhs.getGrid().getSizeY()!=rhs.getGrid().getSizeY()){
		throw std::runtime_error("Incompatible grid at + operation in DiscreteFunction3D");
	}

	DiscreteFunction3D<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result[k]+=rhs(k);
	}

	return result;

}

template<typename ayType, typename axType>
const DiscreteFunction3D<ayType,axType> operator+(DiscreteFunction3D<ayType,axType>&& lhs,DiscreteFunction3D<ayType,axType>&& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize() || lhs.getGrid().getSizeX()!=rhs.getGrid().getSizeX() || lhs.getGrid().getSize()!=rhs.getGrid().getSize() || lhs.getGrid().getSizeX()!=rhs.getGrid().getSizeX() || lhs.getGrid().getSizeY()!=rhs.getGrid().getSizeY()){
		throw std::runtime_error("Incompatible grid at + operation in DiscreteFunction3D");
	}

	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs[k]+=rhs(k);
	}

	return lhs;
}

template<typename ayType, typename axType>
const DiscreteFunction3D<ayType,axType> operator-(const DiscreteFunction3D<ayType,axType>& lhs,const DiscreteFunction3D<ayType,axType>& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize() || lhs.getGrid().getSizeX()!=rhs.getGrid().getSizeX() || lhs.getGrid().getSizeY()!=rhs.getGrid().getSizeY()){
		throw std::runtime_error("Incompatible grid at - operation in DiscreteFunction3D");
	}

	DiscreteFunction3D<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result[k]-=rhs(k);
	}

	return result;

}
template<typename ayType, typename axType>
const DiscreteFunction3D<ayType,axType> operator-(DiscreteFunction3D<ayType,axType>&& lhs, DiscreteFunction3D<ayType,axType>&& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize() || lhs.getGrid().getSizeX()!=rhs.getGrid().getSizeX() || lhs.getGrid().getSizeY()!=rhs.getGrid().getSizeY()){
		throw std::runtime_error("Incompatible grid at - operation in DiscreteFunction3D");
	}

	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs[k]-=rhs(k);
	}

	return lhs;
}

template<typename ayType, typename axType>
const DiscreteFunction3D<ayType,axType> operator*(const DiscreteFunction3D<ayType,axType>& lhs,const DiscreteFunction3D<ayType,axType>& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize() || lhs.getGrid().getSizeX()!=rhs.getGrid().getSizeX() || lhs.getGrid().getSizeY()!=rhs.getGrid().getSizeY()){
		throw std::runtime_error("Incompatible grid at * operation in DiscreteFunction3D");
	}

	DiscreteFunction3D<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result[k]*=rhs(k);
	}

	return result;

}
template<typename ayType, typename axType>
const DiscreteFunction3D<ayType,axType> operator*(DiscreteFunction3D<ayType,axType>&& lhs,DiscreteFunction3D<ayType,axType>&& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize() || lhs.getGrid().getSizeX()!=rhs.getGrid().getSizeX() || lhs.getGrid().getSizeY()!=rhs.getGrid().getSizeY()){
		throw std::runtime_error("Incompatible grid at * operation in DiscreteFunction3D");
	}

	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs[k]*=rhs(k);
	}

	return lhs;
}

template<typename ayType, typename axType>
const DiscreteFunction3D<ayType,axType> operator/(const DiscreteFunction3D<ayType,axType>& lhs,const DiscreteFunction3D<ayType,axType>& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize() || lhs.getGrid().getSizeX()!=rhs.getGrid().getSizeX() || lhs.getGrid().getSizeY()!=rhs.getGrid().getSizeY()){
		throw std::runtime_error("Incompatible grid at / operation in DiscreteFunction3D");
	}

	DiscreteFunction3D<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result[k]/=rhs(k);
	}

	return result;

}
template<typename ayType, typename axType>
const DiscreteFunction3D<ayType,axType> operator/(DiscreteFunction3D<ayType,axType>&& lhs, DiscreteFunction3D<ayType,axType>&& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at / operation in DiscreteFunction3D");
	}

	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs[k]/=rhs(k);
	}

	return lhs;
}

template<typename T, typename gridtype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator+=(DiscreteFunction3D<T,gridtype>& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at += operation in DiscreteFunction3D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)+=func.data[k];
	}

	return *this;
}


template<typename T, typename gridtype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator+=( DiscreteFunction3D<T,gridtype>&& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at += operation in DiscreteFunction3D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)+=func.data[k];
	}

	return *this;
}

template<typename T, typename gridtype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator-=( DiscreteFunction3D<T,gridtype>& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at -= operation in DiscreteFunction3D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)-=func.data[k];
	}

	return *this;
}
template<typename T, typename gridtype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator-=(DiscreteFunction3D<T,gridtype>&& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at -= operation in DiscreteFunction3D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)-=func.data[k];
	}

	return *this;
}


template<typename T, typename gridtype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator-=(const DiscreteFunction3D<T,gridtype>& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at -= operation in DiscreteFunction3D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)-=func.data[k];
	}

	return *this;
}
template<typename T, typename gridtype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator-=(const DiscreteFunction3D<T,gridtype>&& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at -= operation in DiscreteFunction3D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)-=func.data[k];
	}

	return *this;
}

template<typename T, typename gridtype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator*=( DiscreteFunction3D<T,gridtype>& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at *= operation in DiscreteFunction3D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)*=func.data[k];
	}

	return *this;
}
template<typename T, typename gridtype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator*=( DiscreteFunction3D<T,gridtype>&& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at *= operation in DiscreteFunction3D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)*=func.data[k];
	}

	return *this;
}

template<typename T, typename gridtype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator/=( DiscreteFunction3D<T,gridtype>& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at /= operation in DiscreteFunction3D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)/=func.data[k];
	}

	return *this;
}
template<typename T, typename gridtype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator/=( DiscreteFunction3D<T,gridtype>&& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at /= operation in DiscreteFunction3D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)/=func.data[k];
	}

	return *this;
}

template<typename T, typename gridtype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator/=(const DiscreteFunction3D<T,gridtype>& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at /= operation in DiscreteFunction3D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)/=func.data[k];
	}

	return *this;
}

template<typename T, typename gridtype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator/=(const DiscreteFunction3D<T,gridtype>&& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at /= operation in DiscreteFunction3D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)/=func.data[k];
	}

	return *this;
}

template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction3D<ayType,axType> operator+(const DiscreteFunction3D<ayType,axType>& lhs,const numerictype& rhs){
	DiscreteFunction3D<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result.data[k]+=static_cast<ayType>(rhs);
	}
	return result;
}

template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction3D<ayType,axType> operator+(DiscreteFunction3D<ayType,axType>&& lhs,numerictype&& rhs){

	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs.data[k]+=static_cast<ayType>(rhs);
	}
	return lhs;
}

template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction3D<ayType,axType> operator-(const DiscreteFunction3D<ayType,axType>& lhs,const numerictype& rhs){
	DiscreteFunction3D<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result.data[k]-=static_cast<ayType>(rhs);
	}

	return result;
}
template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction3D<ayType,axType> operator-( DiscreteFunction3D<ayType,axType>&& lhs,numerictype&& rhs){
	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs.data[k]+=static_cast<ayType>(rhs);
	}
	return lhs;
}


template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction3D<ayType,axType> operator*(const DiscreteFunction3D<ayType,axType>& lhs,const numerictype& rhs){
	DiscreteFunction3D<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result.data[k]*=static_cast<ayType>(rhs);
	}

	return result;
}

template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction3D<ayType,axType> operator*(DiscreteFunction3D<ayType,axType>&& lhs,numerictype&& rhs){
	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs.data[k]+=static_cast<ayType>(rhs);
	}
	return lhs;
}


template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction3D<ayType,axType> operator/(const DiscreteFunction3D<ayType,axType>& lhs,const numerictype& rhs){
	DiscreteFunction3D<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result.data[k]/=static_cast<ayType>(rhs);
	}

	return result;
}
template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction3D<ayType,axType> operator/(DiscreteFunction3D<ayType,axType>&& lhs,numerictype&& rhs){
	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs.data[k]+=static_cast<ayType>(rhs);
	}
	return lhs;
}

template<typename T, typename gridtype>
template <typename numerictype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator+=(const numerictype& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)+=static_cast<T>(func);
	}

	return *this;
}
template<typename T, typename gridtype>
template <typename numerictype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator+=( numerictype&& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)+=static_cast<T>(func);
	}

	return *this;
}

template<typename T, typename gridtype>
template <typename numerictype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator-=(const numerictype& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)-=static_cast<T>(func);
	}

	return *this;
}


template<typename T, typename gridtype>
template <typename numerictype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator-=(numerictype&& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)-=static_cast<T>(func);
	}

	return *this;
}

template<typename T, typename gridtype>
template <typename numerictype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator*=(const numerictype& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)*=static_cast<T>(func);
	}

	return *this;
}
template<typename T, typename gridtype>
template <typename numerictype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator*=(numerictype&& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)*=static_cast<T>(func);
	}

	return *this;
}

template<typename T, typename gridtype>
template <typename numerictype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator/=(const numerictype& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)/=static_cast<T>(func);
	}

	return *this;
}

template<typename T, typename gridtype>
template <typename numerictype>
DiscreteFunction3D<T,gridtype>& DiscreteFunction3D<T,gridtype>::operator/=( numerictype&& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)/=static_cast<T>(func);
	}

	return *this;
}





} /* namespace epital */

#endif
