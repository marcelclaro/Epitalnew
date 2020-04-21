/*
 * DFunction2D.hpp
 *
 *  Created on: Nov 7, 2014
 *      Author: marcel
 */


#include <armadillo>
#include <armadillo_bits/config.hpp>
#include "Grid.hpp"
#include <memory>

#ifndef DFUNCTION2D_HPP_
#define DFUNCTION2D_HPP_

namespace epital {

template <typename T,typename gridtype=double>
class DiscreteFunction2D {
public:

	DiscreteFunction2D(Grid2D<gridtype> grid);

	/**
	 * Copy constructor
	 * @param tocopy
	 */
	DiscreteFunction2D(const DiscreteFunction2D<T,gridtype>& tocopy);

	/**
	 * Copy constructor
	 * @param tocopy
	 */
	DiscreteFunction2D(std::shared_ptr<DiscreteFunction2D<T,gridtype>> tocopy);

	/**
	 * Copy constructor
	 * @param tocopy
	 */
	DiscreteFunction2D(DiscreteFunction2D<T,gridtype>& tocopy);


	/**
	 * Move constructor
	 * @param tomove
	 */
	DiscreteFunction2D(DiscreteFunction2D<T,gridtype>&& tomove);



	/**
	 * Equal move operator
	 * @param tomove
	 * @return created discrete function
	 */
	DiscreteFunction2D<T,gridtype>& operator=(DiscreteFunction2D<T,gridtype>&& tomove);

	/**
	 * Equal copy operator
	 * @param tomove
	 * @return created discrete function
	 */
	DiscreteFunction2D<T,gridtype>& operator=(DiscreteFunction2D<T,gridtype>& tocopy);

	/**
	 * Operator for fast value read and write
	 * @param position Position read and write
	 * @return Reference to value
	 */
	inline T& operator()(Grid::Position2D position);

	/**
	 * Operator for fast value read and write
	 * @param position Position read and write
	 * @return Reference to value
	 */
	inline T& operator()(long int x, long int y);


	/**
	 * Access to raw data
	 * @return Pointer to data
	 */
	inline T* getdata();



	virtual ~DiscreteFunction2D();

	/**
	 * Return base grid
	 * @return Base grid
	 */
	inline const Grid2D<gridtype>& getGrid() const;



	/**
	 * Fast Fourier Transform of the function
	 * @return reference to initial/transformed function
	 */
	DiscreteFunction2D<T,gridtype>& FFT();


	/**
	 * Fast Fourier Transform of the function
	 * @return reference to initial/transformed function
	 */
	DiscreteFunction2D<T,gridtype>& iFFT();

	/**
	 * Dot product with another function
	 * @param rhs another function
	 * @return dot product value
	 */
	T Dot(DiscreteFunction2D<T,gridtype>& rhs) ;

	/**
	 * Dot product with another function(move semantic version)
	 * @param rhs another function
	 * @return dot product value
	 */
	T Dot(DiscreteFunction2D<T,gridtype>&& rhs) ;


	/**
	 * Return the modulus of function (<f*|f>)^1/2
	 * @return modulus valeu;
	 */
	T Modulus();

	void saveHDF5(std::string filename);

	void normalize(T modulus);

	void Gradient(DiscreteFunction2D<T,gridtype>& x, DiscreteFunction2D<T,gridtype>& y);

	template <typename TC = double>
	DiscreteFunction2D<TC,gridtype> getDensity();


	void fillGaussian(T sigma, T norm);
	void fillGaussian2(T sigma, T norm);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template<typename ayType, typename axType> friend const DiscreteFunction2D<ayType,axType> operator+(const DiscreteFunction2D<ayType,axType>& lhs,const DiscreteFunction2D<ayType,axType>& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template<typename ayType, typename axType> friend const DiscreteFunction2D<ayType,axType> operator+(DiscreteFunction2D<ayType,axType>&& lhs,DiscreteFunction2D<ayType,axType>&& rhs);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template<typename ayType, typename axType> friend const DiscreteFunction2D<ayType,axType> operator-(const DiscreteFunction2D<ayType,axType>& lhs,const DiscreteFunction2D<ayType,axType>& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template<typename ayType, typename axType> friend const DiscreteFunction2D<ayType,axType> operator-(DiscreteFunction2D<ayType,axType>&& lhs, DiscreteFunction2D<ayType,axType>&& rhs);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template<typename ayType, typename axType> friend const DiscreteFunction2D<ayType,axType> operator*(const DiscreteFunction2D<ayType,axType>& lhs,const DiscreteFunction2D<ayType,axType>& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template<typename ayType, typename axType> friend const DiscreteFunction2D<ayType,axType> operator*(DiscreteFunction2D<ayType,axType>&& lhs,DiscreteFunction2D<ayType,axType>&& rhs);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template<typename ayType, typename axType> friend const DiscreteFunction2D<ayType,axType> operator/(const DiscreteFunction2D<ayType,axType>& lhs,const DiscreteFunction2D<ayType,axType>& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template<typename ayType, typename axType> friend const DiscreteFunction2D<ayType,axType> operator/(DiscreteFunction2D<ayType,axType>&& lhs, DiscreteFunction2D<ayType,axType>&& rhs);


	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	DiscreteFunction2D<T,gridtype>& operator+=(DiscreteFunction2D<T,gridtype>& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	DiscreteFunction2D<T,gridtype>& operator+=(DiscreteFunction2D<T,gridtype>&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	DiscreteFunction2D<T,gridtype>& operator-=(DiscreteFunction2D<T,gridtype>& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	DiscreteFunction2D<T,gridtype>& operator-=(DiscreteFunction2D<T,gridtype>&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	DiscreteFunction2D<T,gridtype>& operator*=(DiscreteFunction2D<T,gridtype>& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	DiscreteFunction2D<T,gridtype>& operator*=(DiscreteFunction2D<T,gridtype>&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	DiscreteFunction2D<T,gridtype>& operator/=(DiscreteFunction2D<T,gridtype>& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	DiscreteFunction2D<T,gridtype>& operator/=(DiscreteFunction2D<T,gridtype>&& func);


	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	DiscreteFunction2D<T,gridtype>& operator+=(const DiscreteFunction2D<T,gridtype>& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	DiscreteFunction2D<T,gridtype>& operator+=(const DiscreteFunction2D<T,gridtype>&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	DiscreteFunction2D<T,gridtype>& operator-=(const DiscreteFunction2D<T,gridtype>& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	DiscreteFunction2D<T,gridtype>& operator-=(const DiscreteFunction2D<T,gridtype>&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	DiscreteFunction2D<T,gridtype>& operator*=(const DiscreteFunction2D<T,gridtype>& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	DiscreteFunction2D<T,gridtype>& operator*=(const DiscreteFunction2D<T,gridtype>&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	DiscreteFunction2D<T,gridtype>& operator/=(const DiscreteFunction2D<T,gridtype>& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	DiscreteFunction2D<T,gridtype>& operator/=(const DiscreteFunction2D<T,gridtype>&& func);


	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template <typename ayType, typename axType,typename numerictype> friend const DiscreteFunction2D<ayType,axType> operator+(const DiscreteFunction2D<ayType,axType>& lhs,const numerictype& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template <typename ayType, typename axType,typename numerictype> friend const DiscreteFunction2D<ayType,axType> operator+( DiscreteFunction2D<ayType,axType>&& lhs,numerictype&& rhs);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template <typename ayType, typename axType,typename numerictype> friend const DiscreteFunction2D<ayType,axType> operator-(const DiscreteFunction2D<ayType,axType>& lhs,const numerictype& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template <typename ayType, typename axType,typename numerictype> friend const DiscreteFunction2D<ayType,axType> operator-( DiscreteFunction2D<ayType,axType>&& lhs,numerictype&& rhs);


	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template <typename ayType, typename axType,typename numerictype> friend const DiscreteFunction2D<ayType,axType> operator*(const DiscreteFunction2D<ayType,axType>& lhs,const numerictype& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template <typename ayType, typename axType,typename numerictype> friend const DiscreteFunction2D<ayType,axType> operator*( DiscreteFunction2D<ayType,axType>&& lhs,numerictype&& rhs);


	/**
	 * Operator for arithmetic operations with discrete function
	 * @return resultant function
	 */
	template <typename ayType, typename axType,typename numerictype> friend const DiscreteFunction2D<ayType,axType> operator/(const DiscreteFunction2D<ayType,axType>& lhs,const numerictype& rhs);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return resultant function
	 */
	template <typename ayType, typename axType,typename numerictype> friend const DiscreteFunction2D<ayType,axType> operator/( DiscreteFunction2D<ayType,axType>&& lhs,numerictype&& rhs);


	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	DiscreteFunction2D<T,gridtype>& operator+=(const numerictype& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	DiscreteFunction2D<T,gridtype>& operator+=( numerictype&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	DiscreteFunction2D<T,gridtype>& operator-=(const numerictype& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	DiscreteFunction2D<T,gridtype>& operator-=(numerictype&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	DiscreteFunction2D<T,gridtype>& operator*=(const numerictype& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	DiscreteFunction2D<T,gridtype>& operator*=(numerictype&& func);

	/**
	 * Operator for arithmetic operations with discrete function
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	DiscreteFunction2D<T,gridtype>& operator/=(const numerictype& func);
	/**
	 * Operator for arithmetic operations with discrete function(move semantic version)
	 * @return reference to the initial/new function
	 */
	template <typename numerictype>
	DiscreteFunction2D<T,gridtype>& operator/=( numerictype&& func);



private:
	Grid2D<gridtype> grid;
	std::string filename;
	T* data;
};

} /* namespace epital */

#include "DFunction2D.cpp"

#endif /* DFUNCTION2D_HPP_ */
