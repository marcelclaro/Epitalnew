/*
 * DFunction2D.cpp
 *
 *  Created on: Nov 7, 2014
 *      Author: marcel
 */

#include "DFunction2D.hpp"
#include <fftw3.h>

#ifndef DFUNCTION2D_CPP_
#define DFUNCTION2D_CPP_

namespace epital {


template <typename T,typename gridtype>
DiscreteFunction2D<T,gridtype>::DiscreteFunction2D(Grid2D<gridtype> grid): grid(grid)  {
	data = new T[grid.getSize()];
}


template <typename T,typename gridtype>
DiscreteFunction2D<T,gridtype>::~DiscreteFunction2D() {
	delete [] data;
}

template <typename T,typename gridtype>
DiscreteFunction2D<T,gridtype>::DiscreteFunction2D(const DiscreteFunction2D<T,gridtype>& tocopy): grid(tocopy.getGrid()){
	//Create function data
	data = new T[grid.getSize()];

	for(long int size=0; size<grid.getSize(); ++size){
		data[size]= tocopy.data[size];
	}
}

template <typename T,typename gridtype>
DiscreteFunction2D<T,gridtype>::DiscreteFunction2D(DiscreteFunction2D<T,gridtype>& tocopy): grid(tocopy.getGrid()){
	//Create function data
	data = new T[grid.getSize()];
	T* tocopydata = tocopy.getdata();

	for(long int size=0; size<grid.getSize(); ++size){
		data[size]= tocopydata[size];
	}
}

template <typename T,typename gridtype>
DiscreteFunction2D<T,gridtype>::DiscreteFunction2D(std::shared_ptr<DiscreteFunction2D<T,gridtype>> tocopy): grid(tocopy->getGrid()){
	//Create function data
	data = new T[grid->getSize()];
	T* tocopydata = tocopy->getdata();

	for(long int size=0; size<grid.getSize(); ++size){
		data[size]= tocopydata[size];
	}
}

template <typename T,typename gridtype>
DiscreteFunction2D<T,gridtype>::DiscreteFunction2D(DiscreteFunction2D<T,gridtype>&& tomove): grid(tomove.getGrid()){
	data=tomove.getdata();
	tomove.data=nullptr;
}

template <typename T,typename gridtype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator=(DiscreteFunction2D<T,gridtype>&& tomove){
	if(this == &tomove)
		return *this;
	if(grid.getSize()!=tomove.getGrid().getSize() || grid.getSizeX()!=tomove.getGrid().getSizeX() || grid.getSizeY()!=tomove.getGrid().getSizeY() ){
		throw std::runtime_error("Incompatible grid at = operation in DiscreteFunction2D");
	};

	delete [] data;
	data = nullptr;
	data = tomove.getdata();
	tomove.data=nullptr;
	return *this;
}

template <typename T,typename gridtype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator=(DiscreteFunction2D<T,gridtype>& tocopy){
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
inline T& DiscreteFunction2D<T,gridtype>::operator()(Grid::Position2D position) {
	return data[position[0]+position[1]*grid.getSizeX()];
}

template <typename T,typename gridtype>
inline T& DiscreteFunction2D<T,gridtype>::operator()(long int x, long int y) {
	return data[x+y*grid.getSizeX()];
}



template <typename T,typename gridtype>
inline T* DiscreteFunction2D<T,gridtype>::getdata(){
	return data;
}


template <typename T,typename gridtype>
const Grid2D<gridtype>& DiscreteFunction2D<T,gridtype>::getGrid() const{
    return grid;
}

template <typename T,typename gridtype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::FFT(){
	fftw_complex *in;
	fftw_plan p;
	fftw_init_threads();

	fftw_plan_with_nthreads(4);
	in = reinterpret_cast<fftw_complex*>(this->getdata());
    p = fftw_plan_dft_2D(grid.getSizeX(),grid.getSizeY(), in, in, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p); /* repeat as needed */

    fftw_destroy_plan(p);
    fftw_cleanup_threads();
    return *this;
}


template <typename T,typename gridtype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::iFFT(){
	fftw_complex *in;
	fftw_plan p;
	fftw_init_threads();

	fftw_plan_with_nthreads(4);
	in = reinterpret_cast<fftw_complex*>(this->getdata());
    p = fftw_plan_dft_2D(grid.getSizeX(),grid.getSizeY(), in, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(p); /* repeat as needed */

    fftw_destroy_plan(p);
    fftw_cleanup_threads();
    return *this;
}

template <typename T,typename gridtype>
T DiscreteFunction2D<T,gridtype>::Modulus(){
	long int xmax=grid.getSizeX()-1;
	long int ymax=grid.getSizeY()-1;

	T sum=0;

	/*vertex sum*/
	sum+=(*this)(0,0)*conj((*this)(0,0))/4.0;
	sum+=(*this)(0,ymax)*conj((*this)(0,ymax))/4.0;
	sum+=(*this)(xmax,0)*conj((*this)(xmax,0))/4.0;
	sum+=(*this)(xmax,ymax)*conj((*this)(xmax,ymax))/4.0;

	/*inside + edges sum*/
	#pragma omp parallel shared(sum)
	{
	T partialsum = 0.0;

	#pragma omp for
	for(long int k=1; k<(xmax); k++){
		partialsum+=(*this)(k,0)*conj((*this)(k,0))/2.0;
		partialsum+=(*this)(k,ymax)*conj((*this)(k,ymax))/2.0;
		for(long int l=1; l<(ymax); l++){
				partialsum+=(*this)(k,l)*conj((*this)(k,l));
		}
	}

	#pragma omp critical
	{
		sum+=partialsum;
	}

	}

	/*faces sum*/

	for(long int l=1; l<(ymax); l++){
		sum+=(*this)(0,l)*conj((*this)(0,l))/2.0;
		sum+=(*this)(xmax,l)*conj((*this)(xmax,l))/2.0;
	}



	sum*=(xmax*ymax*grid.getxIncrement()*grid.getyIncrement());
	sum=sqrt(sum);

	return sum;

}


template <typename T,typename gridtype>
void DiscreteFunction2D<T,gridtype>::saveHDF5(std::string filename){
	long int xmax=grid.getSizeX();
	long int ymax=grid.getSizeY();

	arma::Cube<T> datatosave(xmax,ymax);

	#pragma omp for
	for(long int k=0; k<xmax; k++){
		for(long int l=0; l<ymax; l++){
				datatosave(k,l)=(*this)(k,l);
		}
	}

	datatosave.save(filename,arma::hdf5_binary);
}

template <typename T,typename gridtype>
void DiscreteFunction2D<T,gridtype>::fillGaussian(T sigma, T norm){
	for(long int k=0; k<grid.getSizeX();k++)
		for(long int l=0; l<grid.getSizeY();l++){
				T arg(-((k-grid.getSizeX()/2)*(k-grid.getSizeX()/2.0)+(l-grid.getSizeY()/2)*(l-grid.getSizeY()/2.0))/(2.0*sigma));
				(*this)(k,l)=(exp(arg)*norm/2.0)+(exp(arg)*std::complex<double>(0.0,1.0)*norm/2.0);
			}
}

template <typename T,typename gridtype>
void DiscreteFunction2D<T,gridtype>::fillGaussian2(T sigma, T norm){
	for(long int k=0; k<grid.getSizeX();k++)
		for(long int l=0; l<grid.getSizeY();l++){
				T arg(-((k-grid.getSizeX()/2)*(k-grid.getSizeX()/2.0)+(l-grid.getSizeY()/2)*(l-grid.getSizeY()/2.0))/(2.0*sigma));
				(*this)(k,l)=(k-grid.getSizeX()/2.0)*(l-grid.getSizeY()/2.0)*(exp(arg)*norm/2.0)+(exp(arg)*std::complex<double>(0.0,1.0)*norm/2.0);
			}
}


template <typename T,typename gridtype>
void DiscreteFunction2D<T,gridtype>::normalize(T modulus){
	T actualmod = Modulus();
	(*this)/=(actualmod/modulus);
}



template <typename T,typename gridtype>
void DiscreteFunction2D<T,gridtype>::Gradient(DiscreteFunction2D<T,gridtype>& x, DiscreteFunction2D<T,gridtype>& y){
	long int xmax=grid.getSizeX();
	long int ymax=grid.getSizeY();


	for(long int l=0;l<ymax;++l){
			//edge: order h and h²
			x(0,l) = ((*this)(1,l)-(*this)(0,l))/grid.getxIncrement();
			x(1,l) = ((*this)(2,l)-(*this)(0,l))/(2*grid.getxIncrement());

			//central difference order h⁴
			#pragma omp parallel for
			for(long int k = 2; k < (xmax-2); k++){
				x(k,l)=((*this)(k+1,l)*8.0-(*this)(k+2,l)-(*this)(k-1,l)*8.0+(*this)(k-2,l))/(12.0*grid.getxIncrement());
			}

			//edge: order h and h²
			x(xmax-2,l)= ((*this)(xmax-1,l)-(*this)(xmax-3,l))/(2.0*grid.getxIncrement());
			x(xmax-1,l)= ((*this)(xmax-1,l)-(*this)(xmax-2,l))/grid.getxIncrement();
	}

	for(long int k=0;k<xmax;++k){
			//edge: order h and h²
			y(k,0) = ((*this)(k,1)-(*this)(k,0))/grid.getyIncrement();
			y(k,1) = ((*this)(k,2)-(*this)(k,0))/(2*grid.getyIncrement());

			//central difference order h⁴
			#pragma omp parallel for
			for(long int l = 2; l < (ymax-2); l++){
				y(k,l)=((*this)(k,l+1)*8.0-(*this)(k,l+2)-(*this)(k,l-1)*8.0+(*this)(k,l-2))/(12.0*grid.getyIncrement());
			}

			//edge: order h and h²
			y(k,ymax-2)= ((*this)(k,ymax-1)-(*this)(k,ymax-3))/(2.0*grid.getyIncrement());
			y(k,ymax-1)= ((*this)(k,ymax-1)-(*this)(k,ymax-2))/grid.getyIncrement();

	}


}



template <>
template <typename TC>
DiscreteFunction2D<double,double> DiscreteFunction2D<std::complex<double>,double>::getDensity(){
	DiscreteFunction2D<double,double> density(grid);
	for(long int k=0; k<grid.getSizeX();k++)
		for(long int l=0; l<grid.getSizeY();l++){
				density(k,l)=std::real((*this)(k,l)*conj((*this)(k,l)));
		}
	return density;
}

template <typename T,typename gridtype>
T DiscreteFunction2D<T,gridtype>::Dot(DiscreteFunction2D<T,gridtype>& rhs){
	long int xmax=grid.getSizeX()-1;
	long int ymax=grid.getSizeY()-1;

	T sum=0;

	/*vertex sum*/
	sum+=rhs(0,0)*conj((*this)(0,0))/4.0;
	sum+=rhs(0,ymax)*conj((*this)(0,ymax))/4.0;
	sum+=rhs(xmax,0)*conj((*this)(xmax,0))/4.0;
	sum+=rhs(xmax,ymax)*conj((*this)(xmax,ymax))/4.0;

	/*inside + edges sum*/
	#pragma omp parallel shared(sum)
	{
	T partialsum = 0.0;

	#pragma omp for
	for(long int k=1; k<(xmax); k++){
		partialsum+=rhs(k,0)*conj((*this)(k,0))/2.0;
		partialsum+=rhs(k,ymax)*conj((*this)(k,ymax))/2.0;
		for(long int l=1; l<(ymax); l++){
				partialsum+=rhs(k,l)*conj((*this)(k,l));
		}
	}

	#pragma omp critical
	{
		sum+=partialsum;
	}

	}

	/*faces sum*/

	for(long int l=1; l<(ymax); l++){
		sum+=rhs(0,l)*conj((*this)(0,l))/2.0;
		sum+=rhs(xmax,l)*conj((*this)(xmax,l))/2.0;
	}



	sum*=(xmax*ymax*grid.getxIncrement()*grid.getyIncrement());
	sum=sqrt(sum);

	return sum;

}


template <typename T,typename gridtype>
T DiscreteFunction2D<T,gridtype>::Dot(DiscreteFunction2D<T,gridtype>&& rhs){
	long int xmax=grid.getSizeX()-1;
	long int ymax=grid.getSizeY()-1;

	T sum=0;

	/*vertex sum*/
	sum+=rhs(0,0)*conj((*this)(0,0))/4.0;
	sum+=rhs(0,ymax)*conj((*this)(0,ymax))/4.0;
	sum+=rhs(xmax,0)*conj((*this)(xmax,0))/4.0;
	sum+=rhs(xmax,ymax)*conj((*this)(xmax,ymax))/4.0;

	/*inside + edges sum*/
	#pragma omp parallel shared(sum)
	{
	T partialsum = 0.0;

	#pragma omp for
	for(long int k=1; k<(xmax); k++){
		partialsum+=rhs(k,0)*conj((*this)(k,0))/2.0;
		partialsum+=rhs(k,ymax)*conj((*this)(k,ymax))/2.0;
		for(long int l=1; l<(ymax); l++){
				partialsum+=rhs(k,l)*conj((*this)(k,l));
		}
	}

	#pragma omp critical
	{
		sum+=partialsum;
	}

	}

	/*faces sum*/

	for(long int l=1; l<(ymax); l++){
		sum+=rhs(0,l)*conj((*this)(0,l))/2.0;
		sum+=rhs(xmax,l)*conj((*this)(xmax,l))/2.0;
	}



	sum*=(xmax*ymax*grid.getxIncrement()*grid.getyIncrement());
	sum=sqrt(sum);

	return sum;
}



template<typename ayType, typename axType>
const DiscreteFunction2D<ayType,axType> operator+(const DiscreteFunction2D<ayType,axType>& lhs,const DiscreteFunction2D<ayType,axType>& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize() || lhs.getGrid().getSizeX()!=rhs.getGrid().getSizeX() || lhs.getGrid().getSizeY()!=rhs.getGrid().getSizeY()){
		throw std::runtime_error("Incompatible grid at + operation in DiscreteFunction2D");
	}

	DiscreteFunction2D<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result[k]+=rhs(k);
	}

	return result;

}

template<typename ayType, typename axType>
const DiscreteFunction2D<ayType,axType> operator+(DiscreteFunction2D<ayType,axType>&& lhs,DiscreteFunction2D<ayType,axType>&& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize() || lhs.getGrid().getSizeX()!=rhs.getGrid().getSizeX() || lhs.getGrid().getSize()!=rhs.getGrid().getSize() || lhs.getGrid().getSizeX()!=rhs.getGrid().getSizeX() || lhs.getGrid().getSizeY()!=rhs.getGrid().getSizeY()){
		throw std::runtime_error("Incompatible grid at + operation in DiscreteFunction2D");
	}

	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs[k]+=rhs(k);
	}

	return lhs;
}

template<typename ayType, typename axType>
const DiscreteFunction2D<ayType,axType> operator-(const DiscreteFunction2D<ayType,axType>& lhs,const DiscreteFunction2D<ayType,axType>& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize() || lhs.getGrid().getSizeX()!=rhs.getGrid().getSizeX() || lhs.getGrid().getSizeY()!=rhs.getGrid().getSizeY()){
		throw std::runtime_error("Incompatible grid at - operation in DiscreteFunction2D");
	}

	DiscreteFunction2D<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result[k]-=rhs(k);
	}

	return result;

}
template<typename ayType, typename axType>
const DiscreteFunction2D<ayType,axType> operator-(DiscreteFunction2D<ayType,axType>&& lhs, DiscreteFunction2D<ayType,axType>&& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize() || lhs.getGrid().getSizeX()!=rhs.getGrid().getSizeX() || lhs.getGrid().getSizeY()!=rhs.getGrid().getSizeY()){
		throw std::runtime_error("Incompatible grid at - operation in DiscreteFunction2D");
	}

	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs[k]-=rhs(k);
	}

	return lhs;
}

template<typename ayType, typename axType>
const DiscreteFunction2D<ayType,axType> operator*(const DiscreteFunction2D<ayType,axType>& lhs,const DiscreteFunction2D<ayType,axType>& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize() || lhs.getGrid().getSizeX()!=rhs.getGrid().getSizeX() || lhs.getGrid().getSizeY()!=rhs.getGrid().getSizeY()){
		throw std::runtime_error("Incompatible grid at * operation in DiscreteFunction2D");
	}

	DiscreteFunction2D<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result[k]*=rhs(k);
	}

	return result;

}
template<typename ayType, typename axType>
const DiscreteFunction2D<ayType,axType> operator*(DiscreteFunction2D<ayType,axType>&& lhs,DiscreteFunction2D<ayType,axType>&& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize() || lhs.getGrid().getSizeX()!=rhs.getGrid().getSizeX() || lhs.getGrid().getSizeY()!=rhs.getGrid().getSizeY()){
		throw std::runtime_error("Incompatible grid at * operation in DiscreteFunction2D");
	}

	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs[k]*=rhs(k);
	}

	return lhs;
}

template<typename ayType, typename axType>
const DiscreteFunction2D<ayType,axType> operator/(const DiscreteFunction2D<ayType,axType>& lhs,const DiscreteFunction2D<ayType,axType>& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize() || lhs.getGrid().getSizeX()!=rhs.getGrid().getSizeX() || lhs.getGrid().getSizeY()!=rhs.getGrid().getSizeY()){
		throw std::runtime_error("Incompatible grid at / operation in DiscreteFunction2D");
	}

	DiscreteFunction2D<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result[k]/=rhs(k);
	}

	return result;

}
template<typename ayType, typename axType>
const DiscreteFunction2D<ayType,axType> operator/(DiscreteFunction2D<ayType,axType>&& lhs, DiscreteFunction2D<ayType,axType>&& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at / operation in DiscreteFunction2D");
	}

	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs[k]/=rhs(k);
	}

	return lhs;
}

template<typename T, typename gridtype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator+=(DiscreteFunction2D<T,gridtype>& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at += operation in DiscreteFunction2D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)+=func.data[k];
	}

	return *this;
}


template<typename T, typename gridtype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator+=( DiscreteFunction2D<T,gridtype>&& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at += operation in DiscreteFunction2D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)+=func.data[k];
	}

	return *this;
}

template<typename T, typename gridtype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator-=( DiscreteFunction2D<T,gridtype>& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at -= operation in DiscreteFunction2D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)-=func.data[k];
	}

	return *this;
}
template<typename T, typename gridtype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator-=(DiscreteFunction2D<T,gridtype>&& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at -= operation in DiscreteFunction2D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)-=func.data[k];
	}

	return *this;
}


template<typename T, typename gridtype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator-=(const DiscreteFunction2D<T,gridtype>& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at -= operation in DiscreteFunction2D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)-=func.data[k];
	}

	return *this;
}
template<typename T, typename gridtype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator-=(const DiscreteFunction2D<T,gridtype>&& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at -= operation in DiscreteFunction2D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)-=func.data[k];
	}

	return *this;
}

template<typename T, typename gridtype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator*=( DiscreteFunction2D<T,gridtype>& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at *= operation in DiscreteFunction2D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)*=func.data[k];
	}

	return *this;
}
template<typename T, typename gridtype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator*=( DiscreteFunction2D<T,gridtype>&& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at *= operation in DiscreteFunction2D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)*=func.data[k];
	}

	return *this;
}

template<typename T, typename gridtype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator/=( DiscreteFunction2D<T,gridtype>& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at /= operation in DiscreteFunction2D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)/=func.data[k];
	}

	return *this;
}
template<typename T, typename gridtype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator/=( DiscreteFunction2D<T,gridtype>&& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at /= operation in DiscreteFunction2D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)/=func.data[k];
	}

	return *this;
}

template<typename T, typename gridtype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator/=(const DiscreteFunction2D<T,gridtype>& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at /= operation in DiscreteFunction2D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)/=func.data[k];
	}

	return *this;
}

template<typename T, typename gridtype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator/=(const DiscreteFunction2D<T,gridtype>&& func){
	if(grid.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at /= operation in DiscreteFunction2D");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)/=func.data[k];
	}

	return *this;
}

template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction2D<ayType,axType> operator+(const DiscreteFunction2D<ayType,axType>& lhs,const numerictype& rhs){
	DiscreteFunction2D<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result.data[k]+=static_cast<ayType>(rhs);
	}
	return result;
}

template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction2D<ayType,axType> operator+(DiscreteFunction2D<ayType,axType>&& lhs,numerictype&& rhs){

	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs.data[k]+=static_cast<ayType>(rhs);
	}
	return lhs;
}

template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction2D<ayType,axType> operator-(const DiscreteFunction2D<ayType,axType>& lhs,const numerictype& rhs){
	DiscreteFunction2D<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result.data[k]-=static_cast<ayType>(rhs);
	}

	return result;
}
template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction2D<ayType,axType> operator-( DiscreteFunction2D<ayType,axType>&& lhs,numerictype&& rhs){
	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs.data[k]+=static_cast<ayType>(rhs);
	}
	return lhs;
}


template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction2D<ayType,axType> operator*(const DiscreteFunction2D<ayType,axType>& lhs,const numerictype& rhs){
	DiscreteFunction2D<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result.data[k]*=static_cast<ayType>(rhs);
	}

	return result;
}

template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction2D<ayType,axType> operator*(DiscreteFunction2D<ayType,axType>&& lhs,numerictype&& rhs){
	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs.data[k]+=static_cast<ayType>(rhs);
	}
	return lhs;
}


template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction2D<ayType,axType> operator/(const DiscreteFunction2D<ayType,axType>& lhs,const numerictype& rhs){
	DiscreteFunction2D<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result.data[k]/=static_cast<ayType>(rhs);
	}

	return result;
}
template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction2D<ayType,axType> operator/(DiscreteFunction2D<ayType,axType>&& lhs,numerictype&& rhs){
	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs.data[k]+=static_cast<ayType>(rhs);
	}
	return lhs;
}

template<typename T, typename gridtype>
template <typename numerictype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator+=(const numerictype& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)+=static_cast<T>(func);
	}

	return *this;
}
template<typename T, typename gridtype>
template <typename numerictype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator+=( numerictype&& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)+=static_cast<T>(func);
	}

	return *this;
}

template<typename T, typename gridtype>
template <typename numerictype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator-=(const numerictype& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)-=static_cast<T>(func);
	}

	return *this;
}


template<typename T, typename gridtype>
template <typename numerictype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator-=(numerictype&& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)-=static_cast<T>(func);
	}

	return *this;
}

template<typename T, typename gridtype>
template <typename numerictype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator*=(const numerictype& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)*=static_cast<T>(func);
	}

	return *this;
}
template<typename T, typename gridtype>
template <typename numerictype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator*=(numerictype&& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)*=static_cast<T>(func);
	}

	return *this;
}

template<typename T, typename gridtype>
template <typename numerictype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator/=(const numerictype& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)/=static_cast<T>(func);
	}

	return *this;
}

template<typename T, typename gridtype>
template <typename numerictype>
DiscreteFunction2D<T,gridtype>& DiscreteFunction2D<T,gridtype>::operator/=( numerictype&& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid.getSize(); k++){
		*(data+k)/=static_cast<T>(func);
	}

	return *this;
}





} /* namespace epital */

#endif
