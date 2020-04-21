/*
 * Function.hpp
 *
 *  Created on: Jul 4, 2012
 *      Author: marcel
 */

/**
 * @file DFunction.hpp
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
#include <functional>
#include <complex>
#include <memory>
#include <string>

#ifndef DFUNCTION_HPP_
#define DFUNCTION_HPP_

namespace epital{

/**
 * @class DiscreteFunction
 * @brief Define a function in discrete space based on a grid
 */
template<typename yType, typename xType = double>
class DiscreteFunction {
public:

	/**
	 * Create a discrete function above one grid
	 * @param grid Base Grid
	 */
	DiscreteFunction(const Grid1D<xType> grid);

	/**
    * Create a discrete function from a continuous function
	* @param grid Base Grid
	* @param function
	 */
	DiscreteFunction(const Grid1D<xType> grid,std::function<yType(xType)> basefunction);

	/**
	 * Copy constructor
	 * @param tocopy
	 */
	DiscreteFunction(const DiscreteFunction<yType,xType>& tocopy);

	/**
	 * Copy constructor
	 * @param tocopy
	 */
	DiscreteFunction(std::shared_ptr<DiscreteFunction<yType,xType>> tocopy);

	/**
	 * Copy constructor
	 * @param tocopy
	 */
	DiscreteFunction(DiscreteFunction<yType,xType>& tocopy);


	/**
	 * Move constructor
	 * @param tomove
	 */
	DiscreteFunction(DiscreteFunction<yType,xType>&& tomove);

	/**
	 * Function from file
	 * @param filename
	 */
	DiscreteFunction(std::string filename);


	DiscreteFunction(yType* data, Grid1D<xType> grid);


	/**
	 * Equal move operator
	 * @param tomove
	 * @return created discrete function
	 */
	DiscreteFunction<yType,xType>& operator=(DiscreteFunction<yType,xType>&& tomove);

	/**
	 * Equal copy operator
	 * @param tomove
	 * @return created discrete function
	 */
	DiscreteFunction<yType,xType>& operator=(DiscreteFunction<yType,xType>& tocopy);

	/**
	 * Destructor
	 */
	virtual ~DiscreteFunction();




	/**
	 * Return function value at a position
	 * @param position  Position read
	 * @return Function value
	 */
	inline yType getValue(Grid::Position position) const;


	/**
	 * Set value of function
	 * @param position Position to write
	 * @param value Value to write
	 */
	inline void setValue(Grid::Position position,yType value);

	/**
	 * Operator for fast value read
	 * @param position Position read
	 * @return Reference to value
	 */
	inline const yType& operator()(Grid::Position position) const;

	/**
	 * Operator for fast value read and write
	 * @param position Position read or write
	 * @return Reference to value
	 */
	inline yType& operator[](Grid::Position position);

	/**
	 * Access to raw data
	 * @return Pointer to data
	 */
	inline yType* getdata();


	/**
	 * Say if is a complex number function
	 * @return True if is complex
	 */
	inline bool isComplex() const;

	/**
	 * Return real value of a function
	 * Return function itself if is a real or other undefined type
	 * @return the complex value of a function
	 */
	template <typename T = double>
	DiscreteFunction<T,xType> Real();

	/**
	 * Return imaginary value of a function
	 * Return  0 function itself if is a real type or other undefined type
	 * @return the complex value of a function
	 */
	template <typename T = double>
	DiscreteFunction<T,xType> Imaginary();

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template<typename ayType, typename axType> friend const DiscreteFunction<ayType,axType> operator+(const DiscreteFunction<ayType,axType>& lhs,const DiscreteFunction<ayType,axType>& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template<typename ayType, typename axType> friend const DiscreteFunction<ayType,axType> operator+(DiscreteFunction<ayType,axType>&& lhs,DiscreteFunction<ayType,axType>&& rhs);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template<typename ayType, typename axType> friend const DiscreteFunction<ayType,axType> operator-(const DiscreteFunction<ayType,axType>& lhs,const DiscreteFunction<ayType,axType>& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template<typename ayType, typename axType> friend const DiscreteFunction<ayType,axType> operator-(DiscreteFunction<ayType,axType>&& lhs, DiscreteFunction<ayType,axType>&& rhs);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template<typename ayType, typename axType> friend const DiscreteFunction<ayType,axType> operator*(const DiscreteFunction<ayType,axType>& lhs,const DiscreteFunction<ayType,axType>& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template<typename ayType, typename axType> friend const DiscreteFunction<ayType,axType> operator*(DiscreteFunction<ayType,axType>&& lhs,DiscreteFunction<ayType,axType>&& rhs);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template<typename ayType, typename axType> friend const DiscreteFunction<ayType,axType> operator/(const DiscreteFunction<ayType,axType>& lhs,const DiscreteFunction<ayType,axType>& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template<typename ayType, typename axType> friend const DiscreteFunction<ayType,axType> operator/(DiscreteFunction<ayType,axType>&& lhs, DiscreteFunction<ayType,axType>&& rhs);


	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	DiscreteFunction<yType,xType>& operator+=(DiscreteFunction<yType,xType>& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	DiscreteFunction<yType,xType>& operator+=(DiscreteFunction<yType,xType>&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	DiscreteFunction<yType,xType>& operator-=(DiscreteFunction<yType,xType>& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	DiscreteFunction<yType,xType>& operator-=(DiscreteFunction<yType,xType>&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	DiscreteFunction<yType,xType>& operator*=(DiscreteFunction<yType,xType>& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	DiscreteFunction<yType,xType>& operator*=(DiscreteFunction<yType,xType>&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	DiscreteFunction<yType,xType>& operator/=(DiscreteFunction<yType,xType>& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	DiscreteFunction<yType,xType>& operator/=(DiscreteFunction<yType,xType>&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template <typename ayType, typename axType,typename numerictype> friend const DiscreteFunction<ayType,axType> operator+(const DiscreteFunction<ayType,axType>& lhs,const numerictype& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template <typename ayType, typename axType,typename numerictype> friend const DiscreteFunction<ayType,axType> operator+( DiscreteFunction<ayType,axType>&& lhs,numerictype&& rhs);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template <typename ayType, typename axType,typename numerictype> friend const DiscreteFunction<ayType,axType> operator-(const DiscreteFunction<ayType,axType>& lhs,const numerictype& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template <typename ayType, typename axType,typename numerictype> friend const DiscreteFunction<ayType,axType> operator-( DiscreteFunction<ayType,axType>&& lhs,numerictype&& rhs);


	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template <typename ayType, typename axType,typename numerictype> friend const DiscreteFunction<ayType,axType> operator*(const DiscreteFunction<ayType,axType>& lhs,const numerictype& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template <typename ayType, typename axType,typename numerictype> friend const DiscreteFunction<ayType,axType> operator*( DiscreteFunction<ayType,axType>&& lhs,numerictype&& rhs);


	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template <typename ayType, typename axType,typename numerictype> friend const DiscreteFunction<ayType,axType> operator/(const DiscreteFunction<ayType,axType>& lhs,const numerictype& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template <typename ayType, typename axType,typename numerictype> friend const DiscreteFunction<ayType,axType> operator/( DiscreteFunction<ayType,axType>&& lhs,numerictype&& rhs);


	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	DiscreteFunction<yType,xType>& operator+=(const numerictype& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	DiscreteFunction<yType,xType>& operator+=( numerictype&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	DiscreteFunction<yType,xType>& operator-=(const numerictype& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	DiscreteFunction<yType,xType>& operator-=(numerictype&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	DiscreteFunction<yType,xType>& operator*=(const numerictype& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	DiscreteFunction<yType,xType>& operator*=(numerictype&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	DiscreteFunction<yType,xType>& operator/=(const numerictype& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	DiscreteFunction<yType,xType>& operator/=( numerictype&& func);


	/**
	 * Return the integral of the function
	 * @return integral value;
	 */
	yType Integral() const;

	/**
	 * Return the modulus of function (<f*|f>)^1/2
	 * @return modulus valeu;
	 */
	yType Modulus() const;

	/**
	 * Dot product with another function
	 * @param rhs another function
	 * @return dot product value
	 */
	yType Dot(DiscreteFunction<yType,xType>& rhs) const;

	/**
	 * Dot product with another function(move semantic version)
	 * @param rhs another function
	 * @return dot product value
	 */
	yType Dot(DiscreteFunction<yType,xType>&& rhs) const;

	/**
	 * Get fast fourier transform (FFT) of the function
	 * TODO need review
	 * @return FFT function of the function
	 */
	DiscreteFunction<yType,xType> getFFT();

	/**
	 * Get differential of function
	 *
	 * @return the result function from differentiation
	 */
	DiscreteFunction<yType,xType> getDifferential();

	/**
	 * Differentiation of the function
	 * @return reference to differentiated function
	 */
	DiscreteFunction<yType,xType>& Differential();

	/**
	 * Fast Fourier Transform of the function
	 * TODO need review
	 * @return reference to initial/transformed function
	 */
	DiscreteFunction<yType,xType>& FFT();

	/**
	 * Normalize function
	 * @param norm norm of the function (default = 1.0)
	 * @return reference to initial/transformed function
	 */
	DiscreteFunction<yType,xType>& normalize(yType norm=1.0);

	/**
	 * Smooth function
	 * @param passes number of smooth passes
	 * @return reference to initial/transformed function
	 */
	DiscreteFunction<yType,xType>& smooth(long int passes);

	/**
	 * Give the probability density of the wavefunction f ( <f*|f> )
	 * @return probability density function
	 */
	DiscreteFunction<yType,xType> probabilitydensity();


	/**
    * Give the mean position of the wavefunction f ( <f*|z|f> )
	* @return probability density function
	 */
	xType meanLocalization();


	/**
	 * Save to a file
	 * @param filename
	 * @return true if succeed false otherwise
	 */
	bool savetoFile(std::string filename);


	/**
	 * Save to a file
	 * @param filename
	 * @return true if succeed false otherwise
	 */
	bool savetoFileASCII(std::string filename);


	/**
	 * Return base grid
	 * @return Base grid
	 */
	inline const Grid1D<xType>& getGrid() const;

	/**
	 * @class iterator iterator for DiscretFunction
	 * @brief DiscreteFunction iterator
	 * This is a random access iterator for function
	 */
	class iterator:
			public std::iterator<std::random_access_iterator_tag, yType> {
	   public:

	      iterator(yType* p) : pointer_(p) {}

	      iterator() : pointer_(nullptr) {}

	      ~iterator() {}

	      // The assignment and relational operators are straightforward
	      inline iterator& operator=(const iterator& other)
	      {
	    	 pointer_ = other.pointer_;
	         return *this;
	      }

	      inline bool operator==(const iterator& other)
	      {
	         return (pointer_ == other.pointer_);
	      }

	      inline bool operator!=(const iterator& other)
	      {
	         return (pointer_ != other.pointer_);
	      }

	      // Update my state such that I refer to the next element in the
	      // SQueue.
	      inline iterator& operator++(){
	    	  pointer_+=1;
	    	  return (*this);
	      }


	      inline iterator operator++(int){
	    	  iterator tmp(*this);
	    	  ++(*this);
	    	  return tmp;
	      }

	      inline iterator& operator--(){
	      	   pointer_-=1;
	      	   return (*this);
	      }


	      inline iterator operator--(int){
	      	   iterator tmp(*this);
	      	   --(*this);
	      	   return tmp;
	      }

	      inline iterator& operator+=(long int n){
	      	   pointer_+=n;
	      	   return *this;
	      }

	      inline iterator& operator-=(long int n){
	      	   pointer_-=n;
	      	   return *this;
	      }

	      inline iterator operator+(long int n){
	      	  iterator tmp(pointer_+n);
	      	  return tmp;
	      }

	      inline iterator operator-(long int n){
	    	  iterator tmp(pointer_-n);
	    	  return tmp;
	      }

	      inline long int operator-(const iterator& iter){
	    	  long int tmp;
	    	  tmp = iter.pointer_-pointer_;
	    	  return tmp;
	      }

	      inline bool operator<(const iterator& other)
	      {
	         return (pointer_ < other.pointer_);
	      }

	      inline bool operator>(const iterator& other)
	      {
	         return (pointer_ > other.pointer_);
	      }

	      inline bool operator<=(const iterator& other)
	      {
	         return (pointer_ <= other.pointer_);
	      }

	      inline bool operator>=(const iterator& other)
	      {
	         return (pointer_ >= other.pointer_);
	      }

	      inline yType& operator[](long int n){
	    	  return *(pointer_+n);
	      }

	      // Return a reference to the value in the node.  I do this instead
	      // of returning by value so a caller can update the value in the
	      // node directly.
	      inline yType& operator*()
	      {
	         return(*pointer_);
	      }

	      // Return the address of the value referred to.
	      inline yType* operator->()
	      {
	         return(&*(DiscreteFunction<yType,xType>::iterator)*this);
	      }

	      friend class DiscreteFunction<yType,xType>::const_iterator;

	   private:
	      yType* pointer_;
	   };


	/**
	 * iterator to begin of function
	 * @return iterator
	 */
	DiscreteFunction<yType,xType>::iterator begin() const;

    /**
     * iterator to past end of function
     * @return iterator
     */
	DiscreteFunction<yType,xType>::iterator end() const;

	/**
	 * @class const_iterator constant iterator for DiscretFunction
	 * @brief DiscreteFunction constant iterator
	 * This is a random access iterator for function
	 */
	class const_iterator:
			public std::iterator<std::random_access_iterator_tag, yType> {
	   public:

		  const_iterator(yType* p) : pointer_(p) {}

		  const_iterator(const iterator& other) : pointer_(other.pointer_) {}

		  const_iterator() : pointer_(nullptr) {}

	      ~const_iterator() {}

	      // The assignment and relational operators are straightforward
	      inline const_iterator& operator=(const iterator& other)
	      {
	    	 pointer_ = other.pointer_;
	         return *this;
	      }

	      // The assignment and relational operators are straightforward
	      inline const_iterator& operator=(const const_iterator& other)
	      {
	    	 pointer_ = other.pointer_;
	         return *this;
	      }


	      inline bool operator==(const iterator& other)
	      {
	         return (pointer_ == other.pointer_);
	      }

	      inline bool operator!=(const iterator& other)
	      {
	         return (pointer_ != other.pointer_);
	      }

	      // Update my state such that I refer to the next element in the
	      // SQueue.
	      inline const_iterator& operator++(){
	    	  pointer_+=1;
	    	  return (*this);
	      }


	      inline const_iterator operator++(int){
	    	  const_iterator tmp(*this);
	    	  ++(*this);
	    	  return tmp;
	      }

	      inline const_iterator& operator--(){
	      	   pointer_-=1;
	      	   return (*this);
	      }


	      inline const_iterator operator--(int){
	    	  const_iterator tmp(*this);
	      	   --(*this);
	      	   return tmp;
	      }

	      inline const_iterator& operator+=(long int n){
	      	   pointer_+=n;
	      	   return *this;
	      }

	      inline const_iterator& operator-=(long int n){
	      	   pointer_-=n;
	      	   return *this;
	      }

	      inline const_iterator operator+(long int n){
	    	  const_iterator tmp(pointer_+n);
	      	  return tmp;
	      }

	      inline const_iterator operator-(long int n){
	    	  const_iterator tmp(pointer_-n);
	    	  return tmp;
	      }

	      inline long int operator-(const const_iterator& iter){
	    	  long int tmp;
	    	  tmp = iter.pointer_-pointer_;
	    	  return tmp;
	      }

	      inline bool operator<(const const_iterator& other)
	      {
	         return (pointer_ < other.pointer_);
	      }

	      inline bool operator>(const const_iterator& other)
	      {
	         return (pointer_ > other.pointer_);
	      }

	      inline bool operator<=(const const_iterator& other)
	      {
	         return (pointer_ <= other.pointer_);
	      }

	      inline bool operator>=(const const_iterator& other)
	      {
	         return (pointer_ >= other.pointer_);
	      }

	      inline const yType& operator[](long int n){
	    	  return *(pointer_+n);
	      }

	      // Return a reference to the value in the node.  I do this instead
	      // of returning by value so a caller can update the value in the
	      // node directly.
	      inline const yType& operator*()
	      {
	         return(*pointer_);
	      }

	      // Return the address of the value referred to.
	      inline yType* operator->()
	      {
	         return(&*(DiscreteFunction<yType,xType>::iterator)*this);
	      }

	   private:
	      yType* pointer_;
	   };

private:
	yType* data; ///< The raw data
	const Grid1D<xType> grid_; ///< Base grid
	std::string file;
};


}

#include "DFunction.cpp"



#endif /* FUNCTION_HPP_ */
