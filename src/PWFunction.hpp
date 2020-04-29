/*
 * PWFunction.hpp
 *
 *  Created on: Jun 20, 2013
 *      Author: marcel
 */


/**
 * @file PWFunction.hpp
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

#ifndef PWFUNCTION_HPP_
#define PWFUNCTION_HPP_

#include <utility>
#include <vector>
#include <stdexcept>
#include <functional>
#include <complex>
#include <memory>
#include <string>
#include "DFunction.hpp"
#include "Grid.hpp"

using namespace std;

namespace epital {

template<typename complextype, typename xType>
class SolverTM;

/*
template<typename complextype, typename xType>
class DiscreteFunction;


template<typename xType>
class Grid1D;
*/

template<typename complextype, typename xType = double>
class PlaneWavesFunction {
public:
	PlaneWavesFunction() = delete;
	PlaneWavesFunction(long int layers,vector<pair<xType,xType> > limitslst);
	PlaneWavesFunction(const PlaneWavesFunction<complextype,xType>& tocopy);
	PlaneWavesFunction(PlaneWavesFunction<complextype,xType>& tocopy);
	PlaneWavesFunction(PlaneWavesFunction<complextype,xType>&& tomove);

	/**
	 * Equal move operator
	 * @param tomove
	 * @return created discrete function
	 */
	PlaneWavesFunction<complextype,xType>& operator=(PlaneWavesFunction<complextype,xType>&& tomove);

	/**
	 * Equal copy operator
	 * @param tomove
	 * @return created discrete function
	 */
	PlaneWavesFunction<complextype,xType>& operator=(PlaneWavesFunction<complextype,xType>& tocopy);

	virtual ~PlaneWavesFunction();

	/**
	 * Operator for fast value read
	 * @param position Position read
	 * @return Reference to value
	 */
	inline std::complex<complextype> operator()(xType position) const;


	/**
	 * Operator for fast value read
	 * @param position Position read
	 * @return Reference to value
	 */
	inline std::complex<complextype> operator[](xType position) const;

	/**
	 * Get the number of layers
	 * @return number of layers
	 */
	long int getlayers() const;

	/**
	 * Get boundary of definition of function
	 * @return a pair with boundaries
	 */
	pair<xType,xType> getBounds() const;



	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template<typename acomplextype, typename axType> friend const PlaneWavesFunction<acomplextype,axType> operator+(const PlaneWavesFunction<acomplextype,axType>& lhs,const PlaneWavesFunction<acomplextype,axType>& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template<typename acomplextype, typename axType> friend const PlaneWavesFunction<acomplextype,axType> operator+(PlaneWavesFunction<acomplextype,axType>&& lhs,PlaneWavesFunction<acomplextype,axType>&& rhs);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template<typename acomplextype, typename axType> friend const PlaneWavesFunction<acomplextype,axType> operator-(const PlaneWavesFunction<acomplextype,axType>& lhs,const PlaneWavesFunction<acomplextype,axType>& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template<typename acomplextype, typename axType> friend const PlaneWavesFunction<acomplextype,axType> operator-(PlaneWavesFunction<acomplextype,axType>&& lhs, PlaneWavesFunction<acomplextype,axType>&& rhs);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template<typename acomplextype, typename axType> friend const PlaneWavesFunction<acomplextype,axType> operator*(const PlaneWavesFunction<acomplextype,axType>& lhs,const PlaneWavesFunction<acomplextype,axType>& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template<typename acomplextype, typename axType> friend const PlaneWavesFunction<acomplextype,axType> operator*(PlaneWavesFunction<acomplextype,axType>&& lhs,PlaneWavesFunction<acomplextype,axType>&& rhs);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template<typename acomplextype, typename axType> friend const PlaneWavesFunction<acomplextype,axType> operator/(const PlaneWavesFunction<acomplextype,axType>& lhs,const PlaneWavesFunction<acomplextype,axType>& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template<typename acomplextype, typename axType> friend const PlaneWavesFunction<acomplextype,axType> operator/(PlaneWavesFunction<acomplextype,axType>&& lhs, PlaneWavesFunction<acomplextype,axType>&& rhs);


	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	PlaneWavesFunction<complextype,xType>& operator+=(PlaneWavesFunction<complextype,xType>& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	PlaneWavesFunction<complextype,xType>& operator+=(PlaneWavesFunction<complextype,xType>&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	PlaneWavesFunction<complextype,xType>& operator-=(PlaneWavesFunction<complextype,xType>& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	PlaneWavesFunction<complextype,xType>& operator-=(PlaneWavesFunction<complextype,xType>&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	PlaneWavesFunction<complextype,xType>& operator*=(PlaneWavesFunction<complextype,xType>& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	PlaneWavesFunction<complextype,xType>& operator*=(PlaneWavesFunction<complextype,xType>&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	PlaneWavesFunction<complextype,xType>& operator/=(PlaneWavesFunction<complextype,xType>& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	PlaneWavesFunction<complextype,xType>& operator/=(PlaneWavesFunction<complextype,xType>&& func);


	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template <typename acomplextype, typename axType,typename numerictype> friend const PlaneWavesFunction<acomplextype,axType> operator+(const PlaneWavesFunction<acomplextype,axType>& lhs,const numerictype& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template <typename acomplextype, typename axType,typename numerictype> friend const PlaneWavesFunction<acomplextype,axType> operator+( PlaneWavesFunction<acomplextype,axType>&& lhs,numerictype&& rhs);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template <typename acomplextype, typename axType,typename numerictype> friend const PlaneWavesFunction<acomplextype,axType> operator-(const PlaneWavesFunction<acomplextype,axType>& lhs,const numerictype& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template <typename acomplextype, typename axType,typename numerictype> friend const PlaneWavesFunction<acomplextype,axType> operator-( PlaneWavesFunction<acomplextype,axType>&& lhs,numerictype&& rhs);


	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template <typename acomplextype, typename axType,typename numerictype> friend const PlaneWavesFunction<acomplextype,axType> operator*(const PlaneWavesFunction<acomplextype,axType>& lhs,const numerictype& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template <typename acomplextype, typename axType,typename numerictype> friend const PlaneWavesFunction<acomplextype,axType> operator*( PlaneWavesFunction<acomplextype,axType>&& lhs,numerictype&& rhs);


	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template <typename acomplextype, typename axType,typename numerictype> friend const PlaneWavesFunction<acomplextype,axType> operator/(const PlaneWavesFunction<acomplextype,axType>& lhs,const numerictype& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template <typename acomplextype, typename axType,typename numerictype> friend const PlaneWavesFunction<acomplextype,axType> operator/( PlaneWavesFunction<acomplextype,axType>&& lhs,numerictype&& rhs);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	PlaneWavesFunction<complextype,xType>& operator+=(const numerictype& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	PlaneWavesFunction<complextype,xType>& operator+=(numerictype&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	PlaneWavesFunction<complextype,xType>& operator-=(const numerictype& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	PlaneWavesFunction<complextype,xType>& operator-=(numerictype&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	PlaneWavesFunction<complextype,xType>& operator*=(const numerictype& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	PlaneWavesFunction<complextype,xType>& operator*=(numerictype&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	PlaneWavesFunction<complextype,xType>& operator/=(const numerictype& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	PlaneWavesFunction<complextype,xType>& operator/=( numerictype&& func);

	/**
	 * Normalize function
	 * @param norm norm of the function (default = 1.0)
	 * @return reference to initial/transformed function
	 */
	PlaneWavesFunction<complextype,xType>& normalize(complextype norm=1.0);


	/**
	 * Return the integral of the function
	 * @return integral value;
	 */
	std::complex<complextype> Integral() const;

	/**
	 * Return the modulus of function (<f*|f>)^1/2
	 * @return modulus value;
	 */
	complextype Modulus() const;


	void setPeriodic();

	PlaneWavesFunction<complextype,xType>& addPhase(complextype phase);


	/**
    * Give the mean position of the wavefunction f ( <f*|z|f> )
	* @return probability density function
	 */
	xType meanLocalization();

	/**
	 * Dot product with another function
	 * @param rhs another function
	 * @return dot product value
	 */
	std::complex<complextype> Dot(PlaneWavesFunction<complextype,xType>& rhs) const;

	/**
	 * Dot product with another function(move semantic version)
	 * @param rhs another function
	 * @return dot product value
	 */
	std::complex<complextype> Dot(PlaneWavesFunction<complextype,xType>&& rhs) const;

	std::shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> getDiscrete(Grid1D<xType> grid);

private:
	long int layers;
	bool periodic;
	vector<pair<xType,xType> > limitslst;
	vector<vector<pair<std::complex<complextype>,std::complex<complextype>>>> data;
	friend SolverTM<complextype,xType>;

};

} /* namespace epital */

#include "PWFunction.cpp"

#endif /* PWFUNCTION_HPP_ */
