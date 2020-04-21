/*
 * DFunction.cpp
 *
 *  Created on: Jul 4, 2012
 *      Author: marcel
 */

/**
 * @file DFunction.cpp
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

#include "Grid.hpp"
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <functional>
#include <complex>
#include <fftw3.h>
#include <fstream>
#include "DFunction.hpp"

#ifndef DFUNCTION_CPP_
#define DFUNCTION_CPP_

using namespace std;

namespace epital{

template<typename yType, typename xType>
DiscreteFunction<yType,xType>::DiscreteFunction(const Grid1D<xType> grid): grid_(grid){
	data = new yType[grid_.getSize()];
}

template<typename yType, typename xType>
DiscreteFunction<yType,xType>::DiscreteFunction(const Grid1D<xType> grid,std::function<yType(xType)> basefunction): grid_(grid){
	//Exception if is a invalid function pointer
	if(!basefunction)
		throw std::invalid_argument("Null Function Pointer to create DiscreteFunction");

	//Create function data
	data = new yType[grid_.getSize()];

	//Initialize data with basefunction values;
	double step = grid_.getIncrement();
	double begin = grid_.getInferiorLimit();
	//double position = grid_.getInferiorLimit();
	#pragma omp parallel for
	for(long int pointnumber = 0 ; pointnumber < grid_.getSize();pointnumber++){
		*(data+pointnumber)=basefunction(begin+pointnumber*step);
	}
}

/**
 * Copy constructor
 * @param tocopy
 */
template<typename yType, typename xType>
DiscreteFunction<yType,xType>::DiscreteFunction(const DiscreteFunction<yType,xType>& tocopy): grid_(tocopy.getGrid()){

	//Create function data
	data = new yType[tocopy.getGrid().getSize()];

	long int pos = 0;
	for(auto& i : tocopy){
		*(data+pos)= i;
		++pos;
	}
}

/**
 * Copy constructor
 * @param tocopy
 */
template<typename yType, typename xType>
DiscreteFunction<yType,xType>::DiscreteFunction(DiscreteFunction<yType,xType>& tocopy): grid_(tocopy.getGrid()){

		//Create function data
		data = new yType[tocopy.getGrid().getSize()];

		long int pos = 0;
		for(auto& i : tocopy){
			*(data+pos)= i;
			++pos;
		}
}

/**
 * Copy constructor
 * @param tocopy
 */
template<typename yType, typename xType>
DiscreteFunction<yType,xType>::DiscreteFunction(std::shared_ptr<DiscreteFunction<yType,xType>> tocopy): grid_(tocopy->getGrid()){

		//Create function data
		data = new yType[tocopy->getGrid().getSize()];

		long int pos = 0;
		for(auto& i : *tocopy){
			*(data+pos)= i;
			++pos;
		}
}



/**
 * Move constructor
 * @param tomove
 */
template<typename yType, typename xType>
DiscreteFunction<yType,xType>::DiscreteFunction(DiscreteFunction<yType,xType>&& tomove): grid_(tomove.getGrid()){
	data=tomove.getdata();
	tomove.data=nullptr;
}

template<typename yType, typename xType>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator=(DiscreteFunction<yType,xType>&& tomove){
	if(this == &tomove)
		return *this;
	if(grid_.getSize()!=tomove.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at = operation in DiscreteFunction");
	};

	delete [] data;
	data = nullptr;
	data = tomove.getdata();
	tomove.data=nullptr;
	return *this;
}

template<typename yType, typename xType>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator=(DiscreteFunction<yType,xType>& tocopy){
	if(this == &tocopy)
		return *this;

	if(grid_.getSize()!=tocopy.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at = operation in DiscreteFunction");
	};

	long int pos = 0;
	for(auto& i : tocopy){
		*(data+pos)= i;
		++pos;
	}

	grid_=tocopy.getGrid();
	return *this;
}

template<typename yType, typename xType>
DiscreteFunction<yType,xType>::~DiscreteFunction() {
	delete [] data;
}

template<typename yType, typename xType>
typename DiscreteFunction<yType,xType>::iterator DiscreteFunction<yType,xType>::begin() const{
	return iterator(data);
}

template<typename yType, typename xType>
typename DiscreteFunction<yType,xType>::iterator DiscreteFunction<yType,xType>::end() const{
	iterator tmp;
	tmp = iterator(data+grid_.getSize());
	return tmp;
}

template<typename yType, typename xType>
inline yType DiscreteFunction<yType,xType>::getValue(Grid::Position position) const{
	return *(data+position);
}

template<typename yType, typename xType>
inline void DiscreteFunction<yType,xType>::setValue(Grid::Position position,yType value){
	*(data+position)=value;
}

template<typename yType, typename xType>
inline yType& DiscreteFunction<yType,xType>::operator[](Grid::Position position){
	return *(data+position);
}

template<typename yType, typename xType>
inline const yType& DiscreteFunction<yType,xType>::operator()(Grid::Position position) const{
	return *(data+position);
}

template<typename yType, typename xType>
inline yType* DiscreteFunction<yType,xType>::getdata(){
    return data;
}

template<typename yType, typename xType>
const Grid1D<xType>& DiscreteFunction<yType,xType>::getGrid() const{
    return grid_;
}

template<typename yType, typename xType>
inline bool DiscreteFunction<yType,xType>::isComplex() const{
	return false;
}


template<>
inline bool DiscreteFunction<std::complex<double>, double>::isComplex() const{
	return true;
}

template<>
inline bool DiscreteFunction<std::complex<long double>, double>::isComplex() const{
	return true;
}

//


template<>
inline bool DiscreteFunction<std::complex<double>, float>::isComplex() const{
	return true;
}


template<>
inline bool DiscreteFunction<std::complex<long double>, float>::isComplex() const{
	return true;
}

//

template<>
inline bool DiscreteFunction<std::complex<double>,long double>::isComplex() const{
	return true;
}

template<>
inline bool DiscreteFunction<std::complex<long double>,long double>::isComplex() const{
	return true;
}


//General data(long double,doubles, int...)
template<typename yType, typename xType>
template <typename T>
DiscreteFunction<T,xType> DiscreteFunction<yType,xType>::Real(){
	DiscreteFunction<T,xType> newfunc = DiscreteFunction<T,xType>(grid_);

	#pragma omp parallel for shared(newfunc)
	for(long int pos = 0; pos < grid_.getSize(); pos++){
		newfunc[pos]=*(data+pos);
	}

	return newfunc;
}

//Complex functions:

template<>
template <typename T>
DiscreteFunction<T,std::complex<double>> DiscreteFunction<std::complex<double>,std::complex<double>>::Real(){
	DiscreteFunction<T,std::complex<double>> newfunc = DiscreteFunction<T,std::complex<double>>(grid_);

	#pragma omp parallel for shared(newfunc)
	for(long int pos = 0; pos < grid_.getSize(); pos++){
		newfunc[pos]=(data+pos)->real();
	}

	return newfunc;
}

template<>
template <typename T>
DiscreteFunction<T,double> DiscreteFunction<std::complex<double>,double>::Real(){
	DiscreteFunction<T,double> newfunc = DiscreteFunction<T,double>(grid_);

	#pragma omp parallel for shared(newfunc)
	for(long int pos = 0; pos < grid_.getSize(); pos++){
		newfunc[pos]=(data+pos)->real();
	}

	return newfunc;
}

template<>
template <typename T>
DiscreteFunction<T,float> DiscreteFunction<std::complex<double>,float>::Real(){
	DiscreteFunction<T,float> newfunc = DiscreteFunction<T,float>(grid_);

	#pragma omp parallel for shared(newfunc)
	for(long int pos = 0; pos < grid_.getSize(); pos++){
		newfunc[pos]=(data+pos)->real();
	}

	return newfunc;
}

template<>
template <typename T>
DiscreteFunction<T,long double> DiscreteFunction<std::complex<double>,long double>::Real(){
	DiscreteFunction<T,long double> newfunc = DiscreteFunction<T,long double>(grid_);

	#pragma omp parallel for shared(newfunc)
	for(long int pos = 0; pos < grid_.getSize(); pos++){
		newfunc[pos]=(data+pos)->real();
	}

	return newfunc;
}

template<>
template <typename T>
DiscreteFunction<long double,std::complex<double>> DiscreteFunction<std::complex<long double>,std::complex<double>>::Real(){
	DiscreteFunction<long double,std::complex<double>> newfunc = DiscreteFunction<long double,std::complex<double>>(grid_);

	#pragma omp parallel for shared(newfunc)
	for(long int pos = 0; pos < grid_.getSize(); pos++){
		newfunc[pos]=(data+pos)->real();
	}

	return newfunc;
}

template<>
template <typename T>
DiscreteFunction<long double,double> DiscreteFunction<std::complex<long double>,double>::Real(){
	DiscreteFunction<long double,double> newfunc = DiscreteFunction<long double,double>(grid_);

	#pragma omp parallel for shared(newfunc)
	for(long int pos = 0; pos < grid_.getSize(); pos++){
		newfunc[pos]=(data+pos)->real();
	}

	return newfunc;
}

template<>
template <typename T>
DiscreteFunction<long double,long double> DiscreteFunction<std::complex<long double>,long double>::Real(){
	DiscreteFunction<long double,long double> newfunc = DiscreteFunction<long double,long double>(grid_);

	#pragma omp parallel for shared(newfunc)
	for(long int pos = 0; pos < grid_.getSize(); pos++){
		newfunc[pos]=(data+pos)->real();
	}

	return newfunc;
}

//General data(long double,doubles, int...)
template<typename yType, typename xType>
template <typename T>
DiscreteFunction<T,xType> DiscreteFunction<yType,xType>::Imaginary(){
	DiscreteFunction<T,xType> newfunc =DiscreteFunction<T,xType>(grid_);

	#pragma omp parallel for shared(newfunc)
	for(long int pos = 0; pos < grid_.getSize(); pos++){
		newfunc[pos]=0.0;
	}

	return newfunc;
}

//Complex functions:

template<>
template <typename T>
DiscreteFunction<T,std::complex<double>> DiscreteFunction<std::complex<double>,std::complex<double>>::Imaginary(){
	DiscreteFunction<T,std::complex<double>> newfunc = DiscreteFunction<T,std::complex<double>>(grid_);

	#pragma omp parallel for shared(newfunc)
	for(long int pos = 0; pos < grid_.getSize(); pos++){
		newfunc[pos]=(data+pos)->imag();
	}

	return newfunc;
}

template<>
template <typename T>
DiscreteFunction<T,double> DiscreteFunction<std::complex<double>,double>::Imaginary(){
	DiscreteFunction<T,double> newfunc = DiscreteFunction<T,double>(grid_);

	#pragma omp parallel for shared(newfunc)
	for(long int pos = 0; pos < grid_.getSize(); pos++){
		newfunc[pos]=(data+pos)->imag();
	}

	return newfunc;
}

template<>
template <typename T>
DiscreteFunction<T,long double> DiscreteFunction<std::complex<double>,long double>::Imaginary(){
	DiscreteFunction<T,long double> newfunc = DiscreteFunction<T,long double>(grid_);

	#pragma omp parallel for shared(newfunc)
	for(long int pos = 0; pos < grid_.getSize(); pos++){
		newfunc[pos]=(data+pos)->imag();
	}

	return newfunc;
}

template<>
template <typename T>
DiscreteFunction<long double,std::complex<double>> DiscreteFunction<std::complex<long double>,std::complex<double>>::Imaginary(){
	DiscreteFunction<long double,std::complex<double>> newfunc = DiscreteFunction<long double,std::complex<double>>(grid_);


	#pragma omp parallel for shared(newfunc)
	for(long int pos = 0; pos < grid_.getSize(); pos++){
		newfunc[pos]=(data+pos)->imag();
	}


	return newfunc;
}

template<>
template <typename T>
DiscreteFunction<long double,double> DiscreteFunction<std::complex<long double>,double>::Imaginary(){
	DiscreteFunction<long double,double> newfunc = DiscreteFunction<long double,double>(grid_);


	#pragma omp parallel for shared(newfunc)
	for(long int pos = 0; pos < grid_.getSize(); pos++){
		newfunc[pos]=(data+pos)->imag();
	}


	return newfunc;
}

template<>
template <typename T>
DiscreteFunction<long double,float> DiscreteFunction<std::complex<long double>,float>::Imaginary(){
	DiscreteFunction<long double,float> newfunc = DiscreteFunction<long double,float>(grid_);


	#pragma omp parallel for shared(newfunc)
	for(long int pos = 0; pos < grid_.getSize(); pos++){
		newfunc[pos]=(data+pos)->imag();
	}


	return newfunc;
}

template<>
template <typename T>
DiscreteFunction<long double,long double> DiscreteFunction<std::complex<long double>,long double>::Imaginary(){
	DiscreteFunction<long double,long double> newfunc = DiscreteFunction<long double,long double>(grid_);


	#pragma omp parallel for shared(newfunc)
	for(long int pos = 0; pos < grid_.getSize(); pos++){
		newfunc[pos]=(data+pos)->imag();
	}


	return newfunc;
}



template<typename ayType, typename axType>
const DiscreteFunction<ayType,axType> operator+(const DiscreteFunction<ayType,axType>& lhs,const DiscreteFunction<ayType,axType>& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at + operation in DiscreteFunction");
	}

	DiscreteFunction<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result[k]+=rhs(k);
	}

	return result;

}

template<typename ayType, typename axType>
const DiscreteFunction<ayType,axType> operator+(DiscreteFunction<ayType,axType>&& lhs,DiscreteFunction<ayType,axType>&& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at + operation in DiscreteFunction");
	}

	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs[k]+=rhs(k);
	}

	return lhs;
}

template<typename ayType, typename axType>
const DiscreteFunction<ayType,axType> operator-(const DiscreteFunction<ayType,axType>& lhs,const DiscreteFunction<ayType,axType>& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at - operation in DiscreteFunction");
	}

	DiscreteFunction<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result[k]-=rhs(k);
	}

	return result;

}
template<typename ayType, typename axType>
const DiscreteFunction<ayType,axType> operator-(DiscreteFunction<ayType,axType>&& lhs, DiscreteFunction<ayType,axType>&& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at - operation in DiscreteFunction");
	}

	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs[k]-=rhs(k);
	}

	return lhs;
}

template<typename ayType, typename axType>
const DiscreteFunction<ayType,axType> operator*(const DiscreteFunction<ayType,axType>& lhs,const DiscreteFunction<ayType,axType>& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at * operation in DiscreteFunction");
	}

	DiscreteFunction<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result[k]*=rhs(k);
	}

	return result;

}
template<typename ayType, typename axType>
const DiscreteFunction<ayType,axType> operator*(DiscreteFunction<ayType,axType>&& lhs,DiscreteFunction<ayType,axType>&& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at * operation in DiscreteFunction");
	}

	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs[k]*=rhs(k);
	}

	return lhs;
}

template<typename ayType, typename axType>
const DiscreteFunction<ayType,axType> operator/(const DiscreteFunction<ayType,axType>& lhs,const DiscreteFunction<ayType,axType>& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at / operation in DiscreteFunction");
	}

	DiscreteFunction<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result[k]/=rhs(k);
	}

	return result;

}
template<typename ayType, typename axType>
const DiscreteFunction<ayType,axType> operator/(DiscreteFunction<ayType,axType>&& lhs, DiscreteFunction<ayType,axType>&& rhs){
	if(lhs.getGrid().getSize()!=rhs.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at / operation in DiscreteFunction");
	}

	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs[k]/=rhs(k);
	}

	return lhs;
}

template<typename yType, typename xType>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator+=(DiscreteFunction<yType,xType>& func){
	if(grid_.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at += operation in DiscreteFunction");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid_.getSize(); k++){
		*(data+k)+=func(k);
	}

	return *this;
}


template<typename yType, typename xType>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator+=( DiscreteFunction<yType,xType>&& func){
	if(grid_.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at += operation in DiscreteFunction");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid_.getSize(); k++){
		*(data+k)+=func(k);
	}

	return *this;
}

template<typename yType, typename xType>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator-=( DiscreteFunction<yType,xType>& func){
	if(grid_.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at -= operation in DiscreteFunction");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid_.getSize(); k++){
		*(data+k)-=func(k);
	}

	return *this;
}
template<typename yType, typename xType>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator-=(DiscreteFunction<yType,xType>&& func){
	if(grid_.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at -= operation in DiscreteFunction");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid_.getSize(); k++){
		*(data+k)-=func(k);
	}

	return *this;
}

template<typename yType, typename xType>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator*=( DiscreteFunction<yType,xType>& func){
	if(grid_.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at *= operation in DiscreteFunction");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid_.getSize(); k++){
		*(data+k)*=func(k);
	}

	return *this;
}
template<typename yType, typename xType>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator*=( DiscreteFunction<yType,xType>&& func){
	if(grid_.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at *= operation in DiscreteFunction");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid_.getSize(); k++){
		*(data+k)*=func(k);
	}

	return *this;
}

template<typename yType, typename xType>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator/=( DiscreteFunction<yType,xType>& func){
	if(grid_.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at /= operation in DiscreteFunction");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid_.getSize(); k++){
		*(data+k)/=func(k);
	}

	return *this;
}
template<typename yType, typename xType>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator/=( DiscreteFunction<yType,xType>&& func){
	if(grid_.getSize()!=func.getGrid().getSize()){
		throw std::runtime_error("Incompatible grid at /= operation in DiscreteFunction");
	}

	#pragma omp parallel for
	for(long int k = 0; k < grid_.getSize(); k++){
		*(data+k)/=func(k);
	}

	return *this;
}

template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction<ayType,axType> operator+(const DiscreteFunction<ayType,axType>& lhs,const numerictype& rhs){
	DiscreteFunction<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result[k]+=static_cast<ayType>(rhs);
	}
	return result;
}

template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction<ayType,axType> operator+(DiscreteFunction<ayType,axType>&& lhs,numerictype&& rhs){

	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs[k]+=static_cast<ayType>(rhs);
	}
	return lhs;
}

template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction<ayType,axType> operator-(const DiscreteFunction<ayType,axType>& lhs,const numerictype& rhs){
	DiscreteFunction<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result[k]-=static_cast<ayType>(rhs);
	}

	return result;
}
template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction<ayType,axType> operator-( DiscreteFunction<ayType,axType>&& lhs,numerictype&& rhs){
	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs[k]+=static_cast<ayType>(rhs);
	}
	return lhs;
}


template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction<ayType,axType> operator*(const DiscreteFunction<ayType,axType>& lhs,const numerictype& rhs){
	DiscreteFunction<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result[k]*=static_cast<ayType>(rhs);
	}

	return result;
}
template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction<ayType,axType> operator*(DiscreteFunction<ayType,axType>&& lhs,numerictype&& rhs){
	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs[k]+=static_cast<ayType>(rhs);
	}
	return lhs;
}


template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction<ayType,axType> operator/(const DiscreteFunction<ayType,axType>& lhs,const numerictype& rhs){
	DiscreteFunction<ayType,axType> result(lhs);

	#pragma omp parallel for shared(result)
	for(long int k = 0; k < result.getGrid().getSize(); k++){
		result[k]/=static_cast<ayType>(rhs);
	}

	return result;
}
template <typename ayType, typename axType,typename numerictype>
const DiscreteFunction<ayType,axType> operator/(DiscreteFunction<ayType,axType>&& lhs,numerictype&& rhs){
	#pragma omp parallel for shared(lhs)
	for(long int k = 0; k < lhs.getGrid().getSize(); k++){
		lhs[k]+=static_cast<ayType>(rhs);
	}
	return lhs;
}

template<typename yType, typename xType>
template <typename numerictype>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator+=(const numerictype& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid_.getSize(); k++){
		*(data+k)+=static_cast<yType>(func);
	}

	return *this;
}
template<typename yType, typename xType>
template <typename numerictype>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator+=( numerictype&& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid_.getSize(); k++){
		*(data+k)+=static_cast<yType>(func);
	}

	return *this;
}

template<typename yType, typename xType>
template <typename numerictype>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator-=(const numerictype& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid_.getSize(); k++){
		*(data+k)-=static_cast<yType>(func);
	}

	return *this;
}
template<typename yType, typename xType>
template <typename numerictype>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator-=(numerictype&& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid_.getSize(); k++){
		*(data+k)-=static_cast<yType>(func);
	}

	return *this;
}

template<typename yType, typename xType>
template <typename numerictype>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator*=(const numerictype& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid_.getSize(); k++){
		*(data+k)*=static_cast<yType>(func);
	}

	return *this;
}
template<typename yType, typename xType>
template <typename numerictype>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator*=(numerictype&& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid_.getSize(); k++){
		*(data+k)*=static_cast<yType>(func);
	}

	return *this;
}

template<typename yType, typename xType>
template <typename numerictype>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator/=(const numerictype& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid_.getSize(); k++){
		*(data+k)/=static_cast<yType>(func);
	}

	return *this;
}
template<typename yType, typename xType>
template <typename numerictype>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::operator/=( numerictype&& func){
	#pragma omp parallel for
	for(long int k = 0; k < grid_.getSize(); k++){
		*(data+k)/=static_cast<yType>(func);
	}

	return *this;
}

template<typename yType, typename xType>
yType DiscreteFunction<yType,xType>::Integral() const{
	yType value = 0.0;
	yType gridspacing = grid_.getIncrement();

	#pragma omp parallel shared(value)
	{

	yType partialsum = 0.0;

	//Trapezoidal approximation
	#pragma omp for
	for(long int pos = 0; pos < grid_.getSize()-1; pos++){
		partialsum+=gridspacing*(*(data+pos)+*(data+pos+1))/2.0;
	}

	#pragma omp critical
	{
	value+=partialsum;
	}

	}
	return value;
}

template<typename yType, typename xType>
yType DiscreteFunction<yType,xType>::Modulus() const{
	yType value = 0.0;
	yType gridspacing = grid_.getIncrement();

	#pragma omp parallel shared(value)
	{

	yType partialsum = 0.0;

	//Trapezoidal approximation
	#pragma omp for
	for(long int pos = 0; pos < grid_.getSize()-1; pos++){
		partialsum+= gridspacing*( (*(data+pos))*(*(data+pos))+ (*(data+pos+1))*(*(data+pos+1)) )/2.0;
	}


	#pragma omp critical
	{
	value+=partialsum;
	}

	}

	return sqrt(value);
}

//TODO stardardization with complex and other types;
template<typename yType, typename xType>
yType DiscreteFunction<yType,xType>::Dot(DiscreteFunction<yType,xType>& rhs) const{

	yType dotproduct=0.0;

	yType gridspacing = grid_.getIncrement();

	#pragma omp parallel shared(dotproduct)
	{

	yType threadsum = 0.0;

	//Trapezoidal approximation
	#pragma omp for
	for(long int k = 0; k < grid_.getSize()-1; k++){
		threadsum=threadsum + gridspacing*(std::conj(*(data+k)) * rhs(k) + std::conj(*(data+k+1)) * rhs(k+1) )/2.0;
	}

	#pragma omp critical
	{
	dotproduct+=threadsum;
	}


	}
	return dotproduct;

}

//TODO stardardization with complex and other types;
template<typename yType, typename xType>
yType DiscreteFunction<yType,xType>::Dot(DiscreteFunction<yType,xType>&& rhs) const{

	yType dotproduct=0.0;

	yType gridspacing = grid_.getIncrement();

	#pragma omp parallel shared(dotproduct)
	{

	yType threadsum = 0.0;

	//Trapezoidal approximation
	#pragma omp for
	for(long int k = 0; k < grid_.getSize()-1; k++){
		threadsum=threadsum + gridspacing*(std::conj(*(data+k)) * rhs(k) + std::conj(*(data+k+1)) * rhs(k+1) )/2.0;
	}

	#pragma omp critical
	{
	dotproduct+=threadsum;
	}


	}
	return dotproduct;

}

template<typename yType, typename xType>
DiscreteFunction<yType,xType> DiscreteFunction<yType,xType>::getDifferential(){

	DiscreteFunction<yType,xType> newfunc(grid_);

	yType* newdata = newfunc.getdata();

	//Central difference approximation

	//edge: order h and h²
	*newdata = (*(data+1)-*data)/grid_.getIncrement();
	*(newdata+1) = (*(data+2)-*(data))/(2*grid_.getIncrement());

	//central difference order h⁴
	#pragma omp parallel for
	for(long int k = 2; k < (grid_.getSize()-2); k++){
		*(newdata+k)=(*(data+k+1)*8.0-*(data+k+2)-*(data+k-1)*8.0+*(data+k-2))/(12.0*grid_.getIncrement());
	}

	//edge: order h and h²
	*(newdata+grid_.getSize()-2)= (*(data+grid_.getSize()-1)-*(data+grid_.getSize()-3))/(2.0*grid_.getIncrement());
	*(newdata+grid_.getSize()-1)= (*(data+grid_.getSize()-1)-*(data+grid_.getSize()-2))/grid_.getIncrement();

	return newfunc;
}

template<typename yType, typename xType>
DiscreteFunction<yType,xType> DiscreteFunction<yType,xType>::getFFT(){

	DiscreteFunction<yType,xType> result(grid_);

	fftw_complex *in,*out;
	fftw_plan p;

	in = reinterpret_cast<fftw_complex*>(data);
    out = reinterpret_cast<fftw_complex*>(result.getdata());
    p = fftw_plan_dft_1d(grid_.getSize(), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p); /* repeat as needed */

    fftw_destroy_plan(p);

    return result;
}

template<typename yType, typename xType>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::Differential(){

	DiscreteFunction<yType,xType> newfunc(grid_);

	yType* newdata = newfunc.getdata();

	//Central difference approximation

	//edge: order h and h²
	*newdata = (*(data+1)-*data)/grid_.getIncrement();
	*(newdata+1) = (*(data+2)-*(data))/(2*grid_.getIncrement());

	//central difference order h⁴
	#pragma omp parallel for
	for(long int k = 2; k < (grid_.getSize()-2); k++){
		*(newdata+k)=(*(data+k+1)*8.0-*(data+k+2)-*(data+k-1)*8.0+*(data+k-2))/(12.0*grid_.getIncrement());
	}

	//edge: order h and h²
	*(newdata+grid_.getSize()-2)= (*(data+grid_.getSize()-1)-*(data+grid_.getSize()-3))/(2.0*grid_.getIncrement());
	*(newdata+grid_.getSize()-1)= (*(data+grid_.getSize()-1)-*(data+grid_.getSize()-2))/grid_.getIncrement();

	*this=std::move(newfunc);

	return *this;
}

template<typename yType, typename xType>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::FFT(){


	fftw_complex *in;
	fftw_plan p;

	in = reinterpret_cast<fftw_complex*>(data);
    p = fftw_plan_dft_1d(grid_.getSize(), in, in, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p); /* repeat as needed */

    fftw_destroy_plan(p);

    return *this;
}

template<typename yType, typename xType>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::normalize(yType norm){

	yType modulii = this->Modulus();
	*this/=(modulii/norm);

	return *this;

}

template<typename yType, typename xType>
DiscreteFunction<yType,xType>& DiscreteFunction<yType,xType>::smooth(long int passes){

	for(long int k = 0; k < passes ; k++){
		for(long int pos = 1; pos < grid_.getSize()-2; pos++){
			data[pos]=(data[pos-1]+data[pos]+data[pos+1])/3.0;
		}
		for(long int pos = grid_.getSize()-2; pos > 0; pos--){
			data[pos]=(data[pos-1]+data[pos]+data[pos+1])/3.0;
		}
	}

	return *this;
}

template<typename yType, typename xType>
DiscreteFunction<yType,xType> DiscreteFunction<yType,xType>::probabilitydensity(){

	DiscreteFunction<yType,xType> result(grid_);

	#pragma omp parallel for
	for(long int k = 0; k < grid_.getSize(); k++){
		result[k]=std::conj(*(data+k))*(*(data+k));
	}

	return result;
}

template<typename yType, typename xType>
xType DiscreteFunction<yType,xType>::meanLocalization(){

	xType result=0;

	#pragma omp parallel shared(result)
	{

	xType partialsum = 0.0;

		//Trapezoidal approximation
	#pragma omp for
	for(long int k = 0; k < (grid_.getSize()-1); k++){
		partialsum+=0.5*grid_.getIncrement()*(std::conj(*(data+k))*(grid_.getIncrement()*k+grid_.getInferiorLimit())*(*(data+k))+std::conj(*(data+k+1))*(grid_.getIncrement()*(k+1)+grid_.getInferiorLimit())*(*(data+k+1))).real();
	}

	#pragma omp critical
	{
		result+=partialsum;
	}

	}
	return result;
}

template<typename yType, typename xType>
bool DiscreteFunction<yType,xType>::savetoFile(std::string filename){
	file=filename;
	std::ofstream myfile;
	myfile.open(file, ios::out | ios::trunc | ios::binary);
	size_t xtype = typeid(xType).hash_code();
	size_t ytype = typeid(yType).hash_code();
	if(myfile.is_open()){
		myfile.write( reinterpret_cast<const char*>(&xtype),sizeof(size_t));
		myfile.write( reinterpret_cast<const char*>(&ytype),sizeof(size_t));
		myfile.write(reinterpret_cast<const char*>(&grid_),sizeof(const Grid1D<xType>));
		myfile.write(reinterpret_cast<const char*>(data),grid_.getSize()*sizeof(yType));
		myfile.close();
		return true;
	}
	else{
		return false;
	}

}

template<typename yType, typename xType>
bool DiscreteFunction<yType,xType>::savetoFileASCII(std::string filename){
	file=filename;
	std::ofstream myfile;
	myfile.open(file, ios::out | ios::trunc);
	xType position = grid_.getInferiorLimit();
	xType increment = grid_.getIncrement();
	if(myfile.is_open()){
		for(int i=0; i < grid_.getSize();i++){
			myfile << std::setprecision(4) << std::setw(10) << position << "\t" << std::setw(10) << data[i] << "\n";
			position += increment;
		}
		myfile.close();
		return true;
	}
	else{
		return false;
	}

}

template<typename yType, typename xType>
DiscreteFunction<yType,xType>::DiscreteFunction(std::string filename): grid_(0,0,1){
	file=filename;
	std::ifstream myfile;
	myfile.open(file, ios::in | ios::binary);
	size_t xtype = typeid(xType).hash_code();
	size_t ytype = typeid(yType).hash_code();
	size_t filextype;
	size_t fileytype;
	if(myfile.is_open()){
		myfile.read(reinterpret_cast<char*>(&filextype),sizeof(size_t));
		myfile.read(reinterpret_cast<char*>(&fileytype),sizeof(size_t));
		if(xtype!=filextype||ytype!=fileytype)
			throw std::runtime_error("Incompatible xType and/or yType file");
		myfile.read(reinterpret_cast<char*>(const_cast<Grid1D<xType>*>(&grid_)),sizeof(Grid1D<xType>));
		data = new yType[grid_.getSize()];
		myfile.read(reinterpret_cast<char*>(data),grid_.getSize()*sizeof(yType));
		myfile.close();
	}
	else
		throw std::runtime_error("Function File error");


}

template<typename yType, typename xType>
DiscreteFunction<yType,xType>::DiscreteFunction(yType* data, Grid1D<xType> grid): data(data), grid_(grid){
}



}

#endif
