/*
 * PWFunction.cpp
 *
 *  Created on: Jun 20, 2013
 *      Author: marcel
 */


/**
 * @file PWFunction.cpp
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
#include "PWFunction.hpp"
#include <vector>
#include <utility>
#include <iostream>


#ifndef PWFUNCTION_CPP_
#define PWFUNCTION_CPP_

using namespace std;

namespace epital {

template<typename complextype, typename xType>
PlaneWavesFunction<complextype,xType>::PlaneWavesFunction(long int layers, vector<pair<xType,xType> > limitslst): layers(layers), periodic(false), limitslst(limitslst) {
	vector<pair<std::complex<complextype>,std::complex<complextype>>> zerodata;
	data.resize(layers,zerodata);
}

template<typename complextype, typename xType>
PlaneWavesFunction<complextype,xType>::PlaneWavesFunction(const PlaneWavesFunction<complextype,xType>& tocopy) {
	layers=tocopy.layers;
	limitslst=tocopy.limitslst;
	data=tocopy.data;
	periodic=tocopy.periodic;

}

template<typename complextype, typename xType>
PlaneWavesFunction<complextype,xType>::PlaneWavesFunction(PlaneWavesFunction<complextype,xType>& tocopy) {
	layers=tocopy.layers;
	limitslst=tocopy.limitslst;
	data=tocopy.data;
	periodic=tocopy.periodic;
}



template<typename complextype, typename xType>
PlaneWavesFunction<complextype,xType>::PlaneWavesFunction(PlaneWavesFunction<complextype,xType>&& tomove) {
	layers=tomove.layers;
	limitslst=std::move(tomove.limitslst);
	data=std::move(tomove.data);
	periodic=tomove.periodic;
}


template<typename complextype, typename xType>
PlaneWavesFunction<complextype,xType>::~PlaneWavesFunction() {

}

template<typename complextype, typename xType>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator=(PlaneWavesFunction<complextype,xType>&& tomove){
	layers=tomove.layers;
	limitslst=move(tomove.limitslst);
	data=move(tomove.data);
	return *this;
}

template<typename complextype, typename xType>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator=(PlaneWavesFunction<complextype,xType>& tocopy){
	layers=tocopy.layers;
	limitslst=tocopy.limitslst;
	data=tocopy.data;
	return *this;
}

template<typename acomplextype, typename axType>
const PlaneWavesFunction<acomplextype,axType> operator+(const PlaneWavesFunction<acomplextype,axType>& lhs,const PlaneWavesFunction<acomplextype,axType>& rhs){
	if(lhs.layers!=rhs.layers){
		throw std::runtime_error("Incompatible layer numbers at PlaneWavesFunction(PWFunction) ");
	}

	PlaneWavesFunction<acomplextype,axType> result(lhs.layers,lhs.limitslst);
	result.data = lhs.data;
	for(long int i = 0; i< result.layers; i++ ){
		for(pair<std::complex<acomplextype>,std::complex<acomplextype>>& j: rhs.data[i])
			result.data[1].push_back(j);
	}

	return result;

}

template<typename acomplextype, typename axType>
const PlaneWavesFunction<acomplextype,axType> operator+(PlaneWavesFunction<acomplextype,axType>&& lhs,PlaneWavesFunction<acomplextype,axType>&& rhs){
	if(lhs.layers!=rhs.layers){
		throw std::runtime_error("Incompatible layer numbers at PlaneWavesFunction(PWFunction) ");
	}

	for(long int i = 0; i< lhs.layers; i++ ){
	for(pair<std::complex<acomplextype>,std::complex<acomplextype>>& j: rhs.data[i])
		lhs.data[i].push_back(j);
	}

	return lhs;
}

template<typename acomplextype, typename axType>
const PlaneWavesFunction<acomplextype,axType> operator-(const PlaneWavesFunction<acomplextype,axType>& lhs,const PlaneWavesFunction<acomplextype,axType>& rhs){
	if(lhs.layers!=rhs.layers){
		throw std::runtime_error("Incompatible layer numbers at PlaneWavesFunction(PWFunction) ");
	}

	PlaneWavesFunction<acomplextype,axType> result(lhs.layers,lhs.limitslst);
	result.data = lhs.data;
	for(long int i = 0; i< result.layers; i++ ){
		for(pair<std::complex<acomplextype>,std::complex<acomplextype>>& j: rhs.data[i]){
			result.data[i].push_back(make_pair(-j.first,j.second));
		}
	}

	return result;

}
template<typename acomplextype, typename axType>
const PlaneWavesFunction<acomplextype,axType> operator-(PlaneWavesFunction<acomplextype,axType>&& lhs, PlaneWavesFunction<acomplextype,axType>&& rhs){
	if(lhs.layers!=rhs.layers){
		throw std::runtime_error("Incompatible layer numbers at PlaneWavesFunction(PWFunction) ");
	}

	for(long int i = 0; i< lhs.layers; i++ ){
		for(pair<std::complex<acomplextype>,std::complex<acomplextype>>& j: rhs.data[i]){
			lhs.data[i].push_back(make_pair(-j.first,j.second));
		}
	}

	return lhs;
}

template<typename acomplextype, typename axType>
const PlaneWavesFunction<acomplextype,axType> operator*(const PlaneWavesFunction<acomplextype,axType>& lhs,const PlaneWavesFunction<acomplextype,axType>& rhs){
	if(lhs.layers!=rhs.layers){
		throw std::runtime_error("Incompatible layer numbers at PlaneWavesFunction(PWFunction) ");
	}

	PlaneWavesFunction<acomplextype,axType> result(lhs.layers,lhs.limitslst);
	for(long int n = 0; n< lhs.layers; n++ ){
		for(pair<std::complex<acomplextype>,std::complex<acomplextype>>& i: lhs.data[n])
			for(pair<std::complex<acomplextype>,std::complex<acomplextype>>& j: rhs.data[n]){
				result.data[n].push_back(make_pair(i.first*j.first,i.second+j.second));
			}
	}

	return result;

}
template<typename acomplextype, typename axType>
const PlaneWavesFunction<acomplextype,axType> operator*(PlaneWavesFunction<acomplextype,axType>&& lhs,PlaneWavesFunction<acomplextype,axType>&& rhs){
	if(lhs.layers!=rhs.layers){
		throw std::runtime_error("Incompatible layer numbers at PlaneWavesFunction(PWFunction) ");
	}

	PlaneWavesFunction<acomplextype,axType> result(lhs.layers,lhs.limitslst);
	for(long int n = 0; n< lhs.layers; n++ ){
		for(pair<std::complex<acomplextype>,std::complex<acomplextype>>& i: lhs.data[n])
			for(pair<std::complex<acomplextype>,std::complex<acomplextype>>& j: rhs.data[n]){
				result.data[n].push_back(make_pair(i.first*j.first,i.second+j.second));
			}
	}

	return result;
}

template<typename acomplextype, typename axType>
const PlaneWavesFunction<acomplextype,axType> operator/(const PlaneWavesFunction<acomplextype,axType>& lhs,const PlaneWavesFunction<acomplextype,axType>& rhs){
	if(lhs.layers!=rhs.layers){
		throw std::runtime_error("Incompatible layer numbers at PlaneWavesFunction(PWFunction) ");
	}

	PlaneWavesFunction<acomplextype,axType> result(lhs.layers,lhs.limitslst);
	for(long int n = 0; n< lhs.layers; n++ ){
		for(pair<std::complex<acomplextype>,std::complex<acomplextype>>& i: lhs.data[n])
			for(pair<std::complex<acomplextype>,std::complex<acomplextype>>& j: rhs.data[n]){
				result.data.push_back(make_pair(i.first/j.first,i.second-j.second));
			}
	}

	return result;

}
template<typename acomplextype, typename axType>
const PlaneWavesFunction<acomplextype,axType> operator/(PlaneWavesFunction<acomplextype,axType>&& lhs, PlaneWavesFunction<acomplextype,axType>&& rhs){
	if(lhs.layers!=rhs.layers){
		throw std::runtime_error("Incompatible layer numbers at PlaneWavesFunction(PWFunction) ");
	}

	PlaneWavesFunction<acomplextype,axType> result(lhs.layers,lhs.limitslst);
	for(long int n = 0; n< lhs.layers; n++ ){
		for(pair<std::complex<acomplextype>,std::complex<acomplextype>>& i: lhs.data[n])
			for(pair<std::complex<acomplextype>,std::complex<acomplextype>>& j: rhs.data[n]){
				result.data[n].push_back(make_pair(i.first/j.first,i.second-j.second));
			}
	}

	return result;
}

template<typename complextype, typename xType>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator+=(PlaneWavesFunction<complextype,xType>& func){
	if(this->layers!=func.layers){
		throw std::runtime_error("Incompatible layer numbers at PlaneWavesFunction(PWFunction) ");
	}

	for(long int n = 0; n< layers; n++ ){
		for(pair<std::complex<complextype>,std::complex<complextype>>& j: func.data[n])
			data[n].push_back(j);
	}

	return *this;
}


template<typename complextype, typename xType>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator+=( PlaneWavesFunction<complextype,xType>&& func){
	if(this->layers!=func.layers){
		throw std::runtime_error("Incompatible layer numbers at PlaneWavesFunction(PWFunction) ");
	}

	for(long int n = 0; n< layers; n++ ){
		for(pair<std::complex<complextype>,std::complex<complextype>>& j: func.data[n])
			data[n].push_back(j);
	}


	return *this;
}

template<typename complextype, typename xType>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator-=( PlaneWavesFunction<complextype,xType>& func){
	if(this->layers!=func.layers){
		throw std::runtime_error("Incompatible layer numbers at PlaneWavesFunction(PWFunction) ");
	}

	for(long int n = 0; n< layers; n++ ){
		for(pair<std::complex<complextype>,std::complex<complextype>>& j: func.data[n]){
			data[n].push_back(make_pair(-j.first,j.second));
		}
	}

	return *this;
}
template<typename complextype, typename xType>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator-=(PlaneWavesFunction<complextype,xType>&& func){
	if(this->layers!=func.layers){
		throw std::runtime_error("Incompatible layer numbers at PlaneWavesFunction(PWFunction) ");
	}

	for(long int n = 0; n< layers; n++ ){
		for(pair<std::complex<complextype>,std::complex<complextype>>& j: func.data[n]){
			data[n].push_back(make_pair(-j.first,j.second));
		}
	}

	return *this;
}

template<typename complextype, typename xType>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator*=( PlaneWavesFunction<complextype,xType>& func){
	if(this->layers!=func.layers){
		throw std::runtime_error("Incompatible layer numbers at PlaneWavesFunction(PWFunction) ");
	}

	vector<vector<pair<std::complex<complextype>,std::complex<complextype>>>> copy=data;
	data.clear();
	data.resize(layers);
	for(long int n = 0; n< layers; n++ ){
	for(pair<std::complex<complextype>,std::complex<complextype>>& i: copy[n])
		for(pair<std::complex<complextype>,std::complex<complextype>>& j: func.data[n]){
			data[n].push_back(make_pair(i.first*j.first,i.second+j.second));
		}
	}

	return *this;
}
template<typename complextype, typename xType>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator*=( PlaneWavesFunction<complextype,xType>&& func){
	if(this->layers!=func.layers){
		throw std::runtime_error("Incompatible layer numbers at PlaneWavesFunction(PWFunction) ");
	}

	vector<vector<pair<std::complex<complextype>,std::complex<complextype>>>> copy=data;
	data.clear();
	data.resize(layers);
	for(long int n = 0; n< layers; n++ ){
	for(pair<std::complex<complextype>,std::complex<complextype>>& i: copy[n])
		for(pair<std::complex<complextype>,std::complex<complextype>>& j: func.data[n]){
			data[n].push_back(make_pair(i.first*j.first,i.second+j.second));
		}
	}

	return *this;
}

template<typename complextype, typename xType>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator/=( PlaneWavesFunction<complextype,xType>& func){
	if(this->layers!=func.layers){
		throw std::runtime_error("Incompatible layer numbers at PlaneWavesFunction(PWFunction) ");
	}

	vector<vector<pair<std::complex<complextype>,std::complex<complextype>>>> copy=data;
	data.clear();
	data.resize(layers);
	for(long int n = 0; n< layers; n++ ){
	for(pair<std::complex<complextype>,std::complex<complextype>> &i: copy[n])
		for(pair<std::complex<complextype>,std::complex<complextype>>& j: func[n]){
			data[n].push_back(make_pair(i.first/j.first,i.second-j.second));
		}
	}

	return *this;
}
template<typename complextype, typename xType>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator/=( PlaneWavesFunction<complextype,xType>&& func){
	if(this->layers!=func.layers){
		throw std::runtime_error("Incompatible layer numbers at PlaneWavesFunction(PWFunction) ");
	}

	vector<vector<pair<std::complex<complextype>,std::complex<complextype>>>> copy=data;
	data.clear();
	data.resize(layers);
	for(long int n = 0; n< layers; n++ ){
	for(pair<std::complex<complextype>,std::complex<complextype>>& i: copy[n])
		for(pair<std::complex<complextype>,std::complex<complextype>>& j: func[n]){
			data[n].push_back(make_pair(i.first/j.first,i.second-j.second));
		}
	}

	return *this;
}

template <typename acomplextype, typename axType,typename numerictype>
const PlaneWavesFunction<acomplextype,axType> operator+(const PlaneWavesFunction<acomplextype,axType>& lhs,const numerictype& rhs){
	PlaneWavesFunction<acomplextype,axType> result(lhs);

	for(long int n = 0; n< result.layers; n++ ){
		result.data[n].push_back(make_pair(std::complex<acomplextype>(rhs),std::complex<acomplextype>(0.0)));
	}
	return result;
}

template <typename acomplextype, typename axType,typename numerictype>
const PlaneWavesFunction<acomplextype,axType> operator+(PlaneWavesFunction<acomplextype,axType>&& lhs,numerictype&& rhs){
	for(long int n = 0; n< lhs.layers; n++ ){
		lhs.data[n].push_back(make_pair(std::complex<acomplextype>(rhs),std::complex<acomplextype>(0.0)));
	}

	return lhs;
}

template <typename acomplextype, typename axType,typename numerictype>
const PlaneWavesFunction<acomplextype,axType> operator-(const PlaneWavesFunction<acomplextype,axType>& lhs,const numerictype& rhs){
	PlaneWavesFunction<acomplextype,axType> result(lhs);

	for(long int n = 0; n< result.layers; n++ ){
		result.data[n].push_back(make_pair(std::complex<acomplextype>(-rhs),std::complex<acomplextype>(0.0)));
	}
	return result;
}

template <typename acomplextype, typename axType,typename numerictype>
const PlaneWavesFunction<acomplextype,axType> operator-(PlaneWavesFunction<acomplextype,axType>&& lhs,numerictype&& rhs){
	for(long int n = 0; n< lhs.layers; n++ ){
		lhs.data[n].push_back(make_pair(std::complex<acomplextype>(-rhs),std::complex<acomplextype>(0.0)));
	}

	return lhs;
}



template <typename acomplextype, typename axType,typename numerictype>
const PlaneWavesFunction<acomplextype,axType> operator*(const PlaneWavesFunction<acomplextype,axType>& lhs,const numerictype& rhs){
	PlaneWavesFunction<acomplextype,axType> result(lhs.layers,lhs.limitslst);

	result.data.resize(lhs.layers);
	for(long int n = 0; n< result.layers; n++ ){
		for(pair<std::complex<acomplextype>,std::complex<acomplextype>>& j: lhs.data[n]){
			result.data[n].push_back(make_pair(static_cast<acomplextype>(rhs)*j.first,j.second));
		}
	}
	return result;
}

template <typename acomplextype, typename axType,typename numerictype>
const PlaneWavesFunction<acomplextype,axType> operator*(PlaneWavesFunction<acomplextype,axType>&& lhs,numerictype&& rhs){
	for(long int n = 0; n< lhs.layers; n++ ){
		for(pair<std::complex<acomplextype>,std::complex<acomplextype>>& j: lhs.data[n]){
			j.first=static_cast<acomplextype>(rhs)*j.first;
		}
	}

	return lhs;
}


template <typename acomplextype, typename axType,typename numerictype>
const PlaneWavesFunction<acomplextype,axType> operator/(const PlaneWavesFunction<acomplextype,axType>& lhs,const numerictype& rhs){
	PlaneWavesFunction<acomplextype,axType> result(lhs.layers,lhs.limitslst);

	result.data.resize(lhs.layers);
	for(long int n = 0; n< result.layers; n++ ){
		for(pair<std::complex<acomplextype>,std::complex<acomplextype>>& j: lhs.data[n]){
			result.data[n].push_back(make_pair(j.first/static_cast<acomplextype>(rhs),j.second));
		}
	}
	return result;
}
template <typename acomplextype, typename axType,typename numerictype>
const PlaneWavesFunction<acomplextype,axType> operator/(PlaneWavesFunction<acomplextype,axType>&& lhs,numerictype&& rhs){
	for(long int n = 0; n< lhs.layers; n++ ){
		for(pair<std::complex<acomplextype>,std::complex<acomplextype>>& j: lhs.data[n]){
			j.first=j.first/static_cast<acomplextype>(rhs);
		}
	}

	return lhs;
}

template<typename complextype, typename xType>
template <typename numerictype>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator+=(const numerictype& func){
	for(long int n = 0; n< layers; n++ ){
		data[n].push_back(make_pair(std::complex<complextype>(func),std::complex<complextype>(0.0)));
	}

	return *this;
}
template<typename complextype, typename xType>
template <typename numerictype>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator+=(numerictype&& func){
	for(long int n = 0; n< layers; n++ ){
		data[n].push_back(make_pair(std::complex<complextype>(func),std::complex<complextype>(0.0)));
	}

	return *this;
}

template<typename complextype, typename xType>
template <typename numerictype>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator-=(const numerictype& func){
	for(long int n = 0; n< layers; n++ ){
		data[n].push_back(make_pair(std::complex<complextype>(-func),std::complex<complextype>(0.0)));
	}

	return *this;
}
template<typename complextype, typename xType>
template <typename numerictype>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator-=(numerictype&& func){
	for(long int n = 0; n< layers; n++ ){
		data[n].push_back(make_pair(std::complex<complextype>(-func),std::complex<complextype>(0.0)));
	}

	return *this;
}


template<typename complextype, typename xType>
template <typename numerictype>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator*=(const numerictype& func){
	for(long int n = 0; n< layers; n++ ){
		for(pair<std::complex<complextype>,std::complex<complextype>>& j: data[n]){
			j.first=std::complex<complextype>(func)*j.first;
		}
	}

	return *this;
}
template<typename complextype, typename xType>
template <typename numerictype>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator*=(numerictype&& func){
	for(long int n = 0; n< layers; n++ ){
		for(pair<std::complex<complextype>,std::complex<complextype>>& j: data[n]){
			j.first=std::complex<complextype>(func)*j.first;
		}
	}

	return *this;
}

template<typename complextype, typename xType>
template <typename numerictype>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator/=(const numerictype& func){
	for(long int n = 0; n< layers; n++ ){
		for(pair<std::complex<complextype>,std::complex<complextype>>& j: data[n]){
			j.first=(j.first/std::complex<complextype>(func));
		}
	}

	return *this;
}
template<typename complextype, typename xType>
template <typename numerictype>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator/=( numerictype&& func){
	for(long int n = 0; n< layers; n++ ){
		for(pair<std::complex<complextype>,std::complex<complextype>>& j: data[n]){
			j.first=(j.first/std::complex<complextype>(func));
		}
	}

	return *this;
}

template<typename complextype, typename xType>
inline std::complex<complextype> PlaneWavesFunction<complextype,xType>::operator()(xType position) const{
	std::complex<complextype> value(0.0);
	//cout << "position=" << position << endl;
	if(periodic){
		if(position<limitslst[0].first||position>limitslst[layers-1].second){
			xType period=limitslst[layers-1].second-limitslst[0].first;
			//cout << "period=" << period << endl;
			if(position<0.0)
				position=(period+fmod(position,period));
			else
				position = fmod(position,period);
			//cout << "new position=" << position << endl;
		}
		bool vfound=false;
		if(position>=limitslst[0].first && position<=limitslst[0].second){
			for(pair<std::complex<complextype>,std::complex<complextype>> j: data[0]){
				value+=j.first*std::exp(j.second*std::complex<complextype>(position));
			}
			vfound=true;
		}
		for(long int n = 1; n<layers&&(!vfound); n++ ){
			if(position>limitslst[n].first && position<=limitslst[n].second){
				for(pair<std::complex<complextype>,std::complex<complextype>> j: data[n]){
					value+=j.first*std::exp(j.second*std::complex<complextype>(position));
				}
				vfound=true;
			}
		}
		if(vfound==false)
			throw std::runtime_error("error in operator() of PlaneWaveFunction: position are not delimited in defined layers (periodic)");

	}
	else{

		if(position<limitslst[0].first||position>limitslst[layers-1].second)
			throw std::runtime_error("error in operator() of PlaneWaveFunction: out of defined range in PlaneWaveFunction");

		bool vfound=false;
		if(position>=limitslst[0].first && position<=limitslst[0].second){
			for(pair<std::complex<complextype>,std::complex<complextype>> j: data[0]){
				value+=j.first*std::exp(j.second*std::complex<complextype>(position));
			}
			vfound=true;
		}
		for(long int n = 1; n<layers&&(!vfound); n++ ){
			if(position>limitslst[n].first && position<=limitslst[n].second){
				for(pair<std::complex<complextype>,std::complex<complextype>> j: data[n]){
					value+=j.first*std::exp(j.second*std::complex<complextype>(position));
				}
				vfound=true;
			}
		}
		if(vfound==false)
			throw std::runtime_error("error in operator() of PlaneWaveFunction: position are not delimited in defined layers (non periodic)");
	}

	return std::complex<complextype>(static_cast<complextype>(value.real()),static_cast<complextype>(value.imag()));
}

template<typename complextype, typename xType>
inline std::complex<complextype> PlaneWavesFunction<complextype,xType>::operator[](xType position) const{
	std::complex<complextype> value(0.0);
	//cout << "position=" << position << endl;
	if(periodic){
		if(position<limitslst[0].first||position>limitslst[layers-1].second){
			xType period=limitslst[layers-1].second-limitslst[0].first;
			//cout << "period=" << period << endl;
			if(position<0.0)
				position=(period+fmod(position,period));
			else
				position = fmod(position,period);
			//cout << "new position=" << position << endl;
		}
		bool vfound=false;
		if(position>=limitslst[0].first && position<=limitslst[0].second){
			for(pair<std::complex<complextype>,std::complex<complextype>> j: data[0]){
				value+=j.first*std::exp(j.second*std::complex<complextype>(position));
			}
			vfound=true;
		}
		for(long int n = 1; n<layers&&(!vfound); n++ ){
			if(position>limitslst[n].first && position<=limitslst[n].second){
				for(pair<std::complex<complextype>,std::complex<complextype>> j: data[n]){
					value+=j.first*std::exp(j.second*std::complex<complextype>(position));
				}
				vfound=true;
			}
		}
		if(vfound==false)
			throw std::runtime_error("error in operator() of PlaneWaveFunction: position are not delimited in defined layers (periodic)");

	}
	else{

		if(position<limitslst[0].first||position>limitslst[layers-1].second)
			throw std::runtime_error("error in operator() of PlaneWaveFunction: out of defined range in PlaneWaveFunction");

		bool vfound=false;
		if(position>=limitslst[0].first && position<=limitslst[0].second){
			for(pair<std::complex<complextype>,std::complex<complextype>> j: data[0]){
				value+=j.first*std::exp(j.second*std::complex<complextype>(position));
			}
			vfound=true;
		}
		for(long int n = 1; n<layers&&(!vfound); n++ ){
			if(position>limitslst[n].first && position<=limitslst[n].second){
				for(pair<std::complex<complextype>,std::complex<complextype>> j: data[n]){
					value+=j.first*std::exp(j.second*std::complex<complextype>(position));
				}
				vfound=true;
			}
		}
		if(vfound==false)
			throw std::runtime_error("error in operator() of PlaneWaveFunction: position are not delimited in defined layers (non periodic)");
	}

	return value;
}


template<typename complextype, typename xType>
long int PlaneWavesFunction<complextype,xType>::getlayers() const{
	return layers;
}

template<typename complextype, typename xType>
pair<xType,xType> PlaneWavesFunction<complextype,xType>::getBounds() const{
	return make_pair(limitslst.front().first,limitslst.back().second);
}


/*
template<typename complextype, typename xType>
template <typename numerictype>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::operator/=( numerictype&& func){
	for(long int n = 0; n< layers; n++ ){
		for(pair<std::complex<complextype>,std::complex<complextype>> j: data[n]){
			j.first=j.first/static_cast<complextype>(func);
		}
	}

	return *this;
}
*/

/**
 * Return the integral of the function
 * @return integral value;
 */
template<typename complextype, typename xType>
std::complex<complextype> PlaneWavesFunction<complextype,xType>::Integral() const{
	std::complex<complextype> value(0.0);
	for(long int n = 0; n< layers; n++ ){
		for(pair<std::complex<complextype>,std::complex<complextype>> j: data[n]){
			if(j.second==std::complex<complextype>(0.0))
				value+=(j.first*(limitslst[n].second-limitslst[n].first));
			else
				value+=(j.first*(std::exp(j.second*limitslst[n].second)-std::exp(j.second*limitslst[n].first))/j.second);
		}
	}

	return value;
}

/**
 * Return the modulus of function (<f*|f>)^1/2
 * @return modulus value;
 */
template<typename complextype, typename xType>
complextype PlaneWavesFunction<complextype,xType>::Modulus() const{
	vector<vector<pair<std::complex<complextype>,std::complex<complextype>>>> temp;
	temp.resize(layers);
	for(long int n = 0; n< layers; n++ ){
	for(pair<std::complex<complextype>,std::complex<complextype>> i: data[n])
		for(pair<std::complex<complextype>,std::complex<complextype>> j: data[n]){
			temp[n].push_back(make_pair(std::conj(i.first)*j.first,std::conj(i.second)+j.second));
		}
	}
	std::complex<complextype> value(0.0);
	for(long int n = 0; n< layers; n++ ){
		for(pair<std::complex<complextype>,std::complex<complextype>> j: temp[n]){
			if(j.second==std::complex<complextype>(0.0))
				value+=(j.first*std::complex<complextype>(limitslst[n].second-limitslst[n].first));
			else
				value+=(j.first*(std::exp(j.second*std::complex<complextype>(limitslst[n].second))-std::exp(j.second*std::complex<complextype>(limitslst[n].first)))/j.second);
		}
	}

	return static_cast<complextype>((std::sqrt(value)).real());

}

template<typename complextype, typename xType>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::normalize(complextype norm){

	complextype modulii = this->Modulus();
	*this/=(modulii/norm);

	return *this;

}

template<typename complextype, typename xType>
void  PlaneWavesFunction<complextype,xType>::setPeriodic(){
	periodic=true;
}

template<typename complextype, typename xType>
xType PlaneWavesFunction<complextype,xType>::meanLocalization(){
	vector<vector<pair<std::complex<complextype>,std::complex<complextype>>>> temp;
	temp.resize(layers);
	for(long int n = 0; n< layers; n++ ){
	for(pair<std::complex<complextype>,std::complex<complextype>> i: data[n])
		for(pair<std::complex<complextype>,std::complex<complextype>> j: data[n]){
			temp[n].push_back(make_pair(std::conj(i.first)*j.first,std::conj(i.second)+j.second));
		}
	}
	std::complex<complextype> value(0.0);
	for(long int n = 0; n< layers; n++ ){
		for(pair<std::complex<complextype>,std::complex<complextype>> k: temp[n]){
			if(k.second==std::complex<complextype>(0.0))
				value+=(0.5*k.first*(limitslst[n].second*limitslst[n].second-limitslst[n].first*limitslst[n].first));
			else
				value+=(k.first/(k.second*k.second))*( (k.second*limitslst[n].second-1.0)*std::exp(k.second*limitslst[n].second) - (k.second*limitslst[n].first-1.0)*std::exp(k.second*limitslst[n].first)  );
		}
	}

	return static_cast<xType>(value.real());
}

template<typename complextype, typename xType>
PlaneWavesFunction<complextype,xType>& PlaneWavesFunction<complextype,xType>::addPhase(complextype phase){
	for(long int n = 0; n< layers; n++ ){
		for(pair<std::complex<complextype>,std::complex<complextype>>& j: data[n]){
			j.second=(j.second+std::complex<complextype>(0.0l,phase));
		}
	}

	return *this;
}

template<typename complextype, typename xType>
std::shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> PlaneWavesFunction<complextype,xType>::getDiscrete(Grid1D<xType> grid){

	std::shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> discrete = make_shared<DiscreteFunction<std::complex<complextype>,xType>>(grid);
	//Initialize data with basefunction values;
	double step = grid.getIncrement();
	double begin = grid.getInferiorLimit();
	//double position = grid_.getInferiorLimit();
	for(long int pointnumber = 0 ; pointnumber < grid.getSize();pointnumber++){
		(*discrete)[pointnumber]=operator()(begin+pointnumber*step);
	}

	return discrete;

}


} /* namespace epital */
#endif
