/*
 * Grid.h
 *
 *  Created on: Jul 5, 2012
 *      Author: marcel
 */

/**
 * @file Grid.hpp
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

#include <array>
#include <stdexcept>

#ifndef GRID_HPP_
#define GRID_HPP_

namespace epital{

/**
 * Define types for position in a multi(uni)-dimensional grid
 * @brief types for position in a multi(uni)-dimensional grid
 */
namespace Grid{
	typedef long int Position; ///< 1D position
	typedef std::array<long int,2> Position2D; ///< 2D position
	typedef std::array<long int,3> Position3D; ///< 3D position
}

/**
 * @brief This class define a 1D grid/mesh where a function or a space can be defined
 */
template <typename T = double>
class Grid1D {
public:
	/**
	 * Create a 1D regular grid.
	 * @param inferiorlimit value in inferior limit of the grid
	 * @param superiorlimit value in superior limit of the grid
	 * @param gridincrement distance between two point in grid
	 */
	Grid1D(T inferiorlimit, T superiorlimit, long int pointnumbers);
	/**
	 * standard destructor
	 */
	virtual ~Grid1D();

	/**
	 * Return number of points in this grid
	 * @return Number of points
	 */
	long int getSize() const;

	/**
	 * Value of superior limit of the grid
	 * @return superior limit
	 */
	T getSuperiorLimit() const;

	/**
	 * Value of inferior limit of the grid
	 * @return inferior limit
	 */
	T getInferiorLimit() const;

	/**
	 * Distance between two point in grid
	 * @return distance
	 */
	T getIncrement() const;

	/**
	 * Return a position on the grid from a true coordinate of space
	 * @param coordinate coordinate in real space
	 * @return Corresponding coordinate
	 */
	Grid::Position getPosition(T coordinate,bool exact=true) const  ;

	/**
	 * Return true if is the same grid
	 * @param tocompare
	 * @return true or false
	 */
	template <typename J>
	friend bool	operator==(const Grid1D<J>& lhs,const Grid1D<J>& rhs);

	/**
	* Return true if is not the same grid
	* @param tocompare
	* @return true or false
	*/
	template <typename J>
	friend bool	operator!=(const Grid1D<J>& lhs,const Grid1D<J>& rhs);

private:
	long int gridsize; ///< Total grid size
	T inflimit; ///< inferior limit of the grid
	T suplimit; ///< superior limit of the grid
	T increment; ///< Distance between two point in grid
};




/**
 * @brief This class define a 3D grid/mesh where a function or a space can be defined
 */
template <typename T = double>
class Grid3D {
public:
	/**
	 * Create a 3D regular grid
	 * @param xinferiorlimit value in inferior limit of the grid in the x coordinate
	 * @param xsuperiorlimit value in superior limit of the grid in the x coordinate
	 * @param xgridincrement distance between two point in grid in the x coordinate
	 * @param yinferiorlimit value in inferior limit of the grid in the y coordinate
	 * @param ysuperiorlimit value in superior limit of the grid in the y coordinate
	 * @param ygridincrement distance between two point in grid in the y coordinate
	 */
	Grid3D(T xinferiorlimit, T xsuperiorlimit, T xgridincrement,
			T yinferiorlimit, T ysuperiorlimit, T ygridincrement,
			T zinferiorlimit, T zsuperiorlimit, T zgridincrement);

	/**
	 * Standard destructor
	 */
	virtual ~Grid3D();

	/**
	 * Total number of points in the grid
	 * @return number of points x*y*z
	 */
	long int getSize() const;

	/**
	 * Number of points in the x coordinate
	 * @return Number of points in x
	 */
	long int getSizeX() const;

	/**
	 * Number of points in the y coordinate
	 * @return Number of points in y
	 */
	long int getSizeY() const;

	/**
	 * Number of points in the y coordinate
	 * @return Number of points in z
	 */
	long int getSizeZ() const;



	/**
	 * Value of superior limit of x coordinate in the grid
	 * @return superior limit in x
	 */
	T getxSuperiorLimit() const;

	/**
	* Value of inferior limit of x coordinate in the grid
	* @return inferior limit in x
	*/
	T getxInferiorLimit() const;

	/**
	 * Value of superior limit of y coordinate in the grid
	 * @return superior limit in y
	 */
	T getySuperiorLimit() const;

	/**
	* Value of inferior limit of y coordinate in the grid
	* @return inferior limit in y
	*/
	T getyInferiorLimit() const;

	/**
	 * Value of superior limit of y coordinate in the grid
	 * @return superior limit in y
	 */
	T getzSuperiorLimit() const;

	/**
	* Value of inferior limit of y coordinate in the grid
	* @return inferior limit in y
	*/
	T getzInferiorLimit() const;

	/**
	 * Distance between two point on x coordinate of the grid
	 * @return distance in x
	 */
	T getxIncrement() const;

	/**
	 * Distance between two point on y coordinate of the grid
	 * @return distance in y
	 */
	T getyIncrement() const;

	/**
	 * Distance between two point on y coordinate of the grid
	 * @return distance in y
	 */
	T getzIncrement() const;

	/**
	 * Return a position on the grid from a true coordinate of space
	 * @param coordinatex x coordinate in real space
	 * @param coordinatey y coordinate in real space
	 * @return position
	 */
	Grid::Position3D getPosition(T coordinatex, T coordinatey, T coordinatez, bool exact=true) const  ;
private:
	long int gridsize; ///<Total grid size
	long int xgridsize; ///< grid size in x coordinate
	long int ygridsize; ///< grid size in y coordinate
	long int zgridsize; ///< grid size in z coordinate
	double xinflimit; ///< inferior limit of the grid in x coordinate
	double xsuplimit; ///< superior limit of the grid in x coordinate
	double xincrement; ///< Distance between two point in x coordinate of the grid
	double yinflimit; ///< inferior limit of the grid in y coordinate
	double ysuplimit; ///< superior limit of the grid in y coordinate
	double yincrement; ///< Distance between two point in y coordinate of the grid
	double zinflimit; ///< inferior limit of the grid in z coordinate
	double zsuplimit; ///< superior limit of the grid in z coordinate
	double zincrement; ///< Distance between two point in z coordinate of the grid
};

/**
 * @brief This class define a 2D grid/mesh where a function or a space can be defined
 */
template <typename T = double>
class Grid2D {
public:
	/**
	 * Create a 2D regular grid
	 * @param xinferiorlimit value in inferior limit of the grid in the x coordinate
	 * @param xsuperiorlimit value in superior limit of the grid in the x coordinate
	 * @param xgridincrement distance between two point in grid in the x coordinate
	 * @param yinferiorlimit value in inferior limit of the grid in the y coordinate
	 * @param ysuperiorlimit value in superior limit of the grid in the y coordinate
	 * @param ygridincrement distance between two point in grid in the y coordinate
	 */
	Grid2D(T xinferiorlimit, T xsuperiorlimit, T xgridincrement,
			T yinferiorlimit, T ysuperiorlimit, T ygridincrement);

	/**
	 * Standard destructor
	 */
	virtual ~Grid2D();

	/**
	 * Total number of points in the grid
	 * @return number of points x*y*z
	 */
	long int getSize() const;

	/**
	 * Number of points in the x coordinate
	 * @return Number of points in x
	 */
	long int getSizeX() const;

	/**
	 * Number of points in the y coordinate
	 * @return Number of points in y
	 */
	long int getSizeY() const;



	/**
	 * Value of superior limit of x coordinate in the grid
	 * @return superior limit in x
	 */
	T getxSuperiorLimit() const;

	/**
	* Value of inferior limit of x coordinate in the grid
	* @return inferior limit in x
	*/
	T getxInferiorLimit() const;

	/**
	 * Value of superior limit of y coordinate in the grid
	 * @return superior limit in y
	 */
	T getySuperiorLimit() const;

	/**
	* Value of inferior limit of y coordinate in the grid
	* @return inferior limit in y
	*/
	T getyInferiorLimit() const;


	/**
	 * Distance between two point on x coordinate of the grid
	 * @return distance in x
	 */
	T getxIncrement() const;

	/**
	 * Distance between two point on y coordinate of the grid
	 * @return distance in y
	 */
	T getyIncrement() const;


	/**
	 * Return a position on the grid from a true coordinate of space
	 * @param coordinatex x coordinate in real space
	 * @param coordinatey y coordinate in real space
	 * @return position
	 */
	Grid::Position2D getPosition(T coordinatex, T coordinatey, bool exact=true) const  ;
private:
	long int gridsize; ///<Total grid size
	long int xgridsize; ///< grid size in x coordinate
	long int ygridsize; ///< grid size in y coordinate
	double xinflimit; ///< inferior limit of the grid in x coordinate
	double xsuplimit; ///< superior limit of the grid in x coordinate
	double xincrement; ///< Distance between two point in x coordinate of the grid
	double yinflimit; ///< inferior limit of the grid in y coordinate
	double ysuplimit; ///< superior limit of the grid in y coordinate
	double yincrement; ///< Distance between two point in y coordinate of the grid
};


/*Implementation: */

template <typename J>
bool	operator==(const Grid1D<J>& lhs,const Grid1D<J>& rhs){
	return (lhs.gridsize==rhs.gridsize&&lhs.inflimit==rhs.inflimit&&lhs.suplimit==rhs.suplimit&&lhs.increment==rhs.increment);
}

template <typename J>
bool	operator==(const Grid3D<J>& lhs,const Grid3D<J>& rhs){
	return (lhs.gridsize==rhs.gridsize&&lhs.xinflimit==rhs.xinflimit&&lhs.xsuplimit==rhs.xsuplimit&&lhs.xincrement==rhs.xincrement
			&&lhs.yinflimit==rhs.yinflimit&&lhs.ysuplimit==rhs.ysuplimit&&lhs.yincrement==rhs.yincrement
			&&lhs.zinflimit==rhs.zinflimit&&lhs.zsuplimit==rhs.zsuplimit&&lhs.zincrement==rhs.zincrement);
}

template <typename J>
bool	operator==(const Grid2D<J>& lhs,const Grid2D<J>& rhs){
	return (lhs.gridsize==rhs.gridsize&&lhs.xinflimit==rhs.xinflimit&&lhs.xsuplimit==rhs.xsuplimit&&lhs.xincrement==rhs.xincrement
			&&lhs.yinflimit==rhs.yinflimit&&lhs.ysuplimit==rhs.ysuplimit&&lhs.yincrement==rhs.yincrement);
}

template <typename J>
bool	operator!=(const Grid1D<J>& lhs,const Grid1D<J>& rhs){
	return (lhs.gridsize!=rhs.gridsize&&lhs.inflimit!=rhs.inflimit&&lhs.suplimit!=rhs.suplimit&&lhs.increment!=rhs.increment);
}

template <typename J>
bool	operator!=(const Grid3D<J>& lhs,const Grid3D<J>& rhs){
	return (lhs.gridsize!=rhs.gridsize
			&&lhs.xinflimit!=rhs.xinflimit&&lhs.xsuplimit!=rhs.xsuplimit&&lhs.xincrement!=rhs.xincrement
			&&lhs.yinflimit!=rhs.yinflimit&&lhs.ysuplimit!=rhs.ysuplimit&&lhs.yincrement!=rhs.yincrement
			&&lhs.zinflimit!=rhs.zinflimit&&lhs.zsuplimit!=rhs.zsuplimit&&lhs.zincrement!=rhs.zincrement);
}


template <typename J>
bool	operator!=(const Grid2D<J>& lhs,const Grid2D<J>& rhs){
	return (lhs.gridsize!=rhs.gridsize&&
			lhs.xinflimit!=rhs.xinflimit&&lhs.xsuplimit!=rhs.xsuplimit&&lhs.xincrement!=rhs.xincrement&&
			lhs.yinflimit!=rhs.yinflimit&&lhs.ysuplimit!=rhs.ysuplimit&&lhs.yincrement!=rhs.yincrement);
}

template <typename T>
Grid1D<T>::Grid1D(T inferiorlimit, T superiorlimit, long int pointnumbers): gridsize(pointnumbers),inflimit(inferiorlimit),suplimit(superiorlimit){
	increment=(suplimit - inflimit)/static_cast<double>(pointnumbers);
}

template <typename T>
Grid3D<T>::Grid3D(T xinferiorlimit, T xsuperiorlimit, T xgridincrement,
		T yinferiorlimit, T ysuperiorlimit, T ygridincrement,
		T zinferiorlimit, T zsuperiorlimit, T zgridincrement)
: xinflimit(xinferiorlimit),xsuplimit(xsuperiorlimit),xincrement(xgridincrement),
  yinflimit(yinferiorlimit),ysuplimit(ysuperiorlimit),yincrement(ygridincrement),
  zinflimit(zinferiorlimit),zsuplimit(zsuperiorlimit),zincrement(zgridincrement){
	xgridsize=((xsuplimit-xinflimit)/xincrement)+1;
	ygridsize=((ysuplimit-yinflimit)/yincrement)+1;
	zgridsize=((zsuplimit-zinflimit)/zincrement)+1;
	gridsize=xgridsize*ygridsize*zgridsize;
}

template <typename T>
Grid2D<T>::Grid2D(T xinferiorlimit, T xsuperiorlimit, T xgridincrement,
		T yinferiorlimit, T ysuperiorlimit, T ygridincrement)
: xinflimit(xinferiorlimit),xsuplimit(xsuperiorlimit),xincrement(xgridincrement),
  yinflimit(yinferiorlimit),ysuplimit(ysuperiorlimit),yincrement(ygridincrement){
	xgridsize=((xsuplimit-xinflimit)/xincrement)+1;
	ygridsize=((ysuplimit-yinflimit)/yincrement)+1;
	gridsize=xgridsize*ygridsize;
}

template <typename T>
Grid1D<T>::~Grid1D() {

}

template <typename T>
Grid3D<T>::~Grid3D() {

}

template <typename T>
Grid2D<T>::~Grid2D() {

}

template <typename T>
long int Grid1D<T>::getSize() const{
	return gridsize;
}

template <typename T>
long int Grid3D<T>::getSize() const{
	return gridsize;
}

template <typename T>
long int Grid2D<T>::getSize() const{
	return gridsize;
}

template <typename T>
long int Grid3D<T>::getSizeX() const{
	return xgridsize;
}

template <typename T>
long int Grid2D<T>::getSizeX() const{
	return xgridsize;
}


template <typename T>
long int Grid3D<T>::getSizeY() const{
	return ygridsize;
}

template <typename T>
long int Grid2D<T>::getSizeY() const{
	return ygridsize;
}



template <typename T>
long int Grid3D<T>::getSizeZ() const{
	return zgridsize;
}



template <typename T>
T Grid1D<T>::getSuperiorLimit() const{
	return suplimit;
}

template <typename T>
T Grid1D<T>::getInferiorLimit() const{
	return inflimit;
}

template <typename T>
T Grid1D<T>::getIncrement() const{
	return increment;
}

/*/////////////////////////////////////*/

template <typename T>
T Grid3D<T>::getxSuperiorLimit() const{
	return xsuplimit;
}

template <typename T>
T Grid3D<T>::getxInferiorLimit() const{
	return xinflimit;
}

template <typename T>
T Grid2D<T>::getxSuperiorLimit() const{
	return xsuplimit;
}

template <typename T>
T Grid2D<T>::getxInferiorLimit() const{
	return xinflimit;
}

template <typename T>
T Grid3D<T>::getxIncrement() const{
	return xincrement;
}

template <typename T>
T Grid2D<T>::getxIncrement() const{
	return xincrement;
}

template <typename T>
T Grid3D<T>::getySuperiorLimit() const{
	return ysuplimit;
}

template <typename T>
T Grid3D<T>::getyInferiorLimit() const{
	return yinflimit;
}

template <typename T>
T Grid3D<T>::getyIncrement() const{
	return yincrement;
}

template <typename T>
T Grid2D<T>::getySuperiorLimit() const{
	return ysuplimit;
}

template <typename T>
T Grid2D<T>::getyInferiorLimit() const{
	return yinflimit;
}

template <typename T>
T Grid2D<T>::getyIncrement() const{
	return yincrement;
}

template <typename T>
T Grid3D<T>::getzSuperiorLimit() const{
	return zsuplimit;
}

template <typename T>
T Grid3D<T>::getzInferiorLimit() const{
	return zinflimit;
}

template <typename T>
T Grid3D<T>::getzIncrement() const{
	return zincrement;
}


template <typename T>
Grid::Position Grid1D<T>::getPosition(T coordinate,bool exact) const  {
	Grid::Position position;
	T temp = (coordinate-inflimit);

	//Test if coordinate is valid (resides on the grid)
	if(exact){
		T test = increment-fmod(temp,increment);
		T tolerance=increment*0.01;
		if(test>tolerance&&test<(increment-tolerance))
			throw std::invalid_argument("Invalid coordinate: Must be on a grid point!");
	}

	//calculate position
	position = static_cast<long int>(floor(temp/increment+0.5));
	return position;
}

template <typename T>
Grid::Position2D Grid2D<T>::getPosition(T coordinatex,T coordinatey,bool exact) const  {
	Grid::Position2D position;
	T tempx = (coordinatex-xinflimit);
	T tempy = (coordinatey-yinflimit);


	//Test if coordinate is valid (resides on the grid)
	if(exact){
		T testx = xincrement-fmod(tempx,xincrement);
		T testy = yincrement-fmod(tempy,yincrement);
		T tolerancex=xincrement*0.01;
		T tolerancey=yincrement*0.01;

		if(testx>tolerancex&&testx<(xincrement-tolerancex))
			throw std::invalid_argument("Invalid coordinate X: Must be on a grid point!");
		if(testy>tolerancey&&testy<(yincrement-tolerancey))
					throw std::invalid_argument("Invalid coordinate Y: Must be on a grid point!");

	}

	//calculate position
	position[0] = static_cast<long int>(floor(tempx/xincrement+0.5));
	position[1] = static_cast<long int>(floor(tempy/yincrement+0.5));
	return position;
}

template <typename T>
Grid::Position3D Grid3D<T>::getPosition(T coordinatex,T coordinatey,T coordinatez,bool exact) const  {
	Grid::Position3D position;
	T tempx = (coordinatex-xinflimit);
	T tempy = (coordinatey-yinflimit);
	T tempz = (coordinatez-zinflimit);

	//Test if coordinate is valid (resides on the grid)
	if(exact){
		T testx = xincrement-fmod(tempx,xincrement);
		T testy = yincrement-fmod(tempy,yincrement);
		T testz = zincrement-fmod(tempz,zincrement);
		T tolerancex=xincrement*0.01;
		T tolerancey=yincrement*0.01;
		T tolerancez=zincrement*0.01;
		if(testx>tolerancex&&testx<(xincrement-tolerancex))
			throw std::invalid_argument("Invalid coordinate X: Must be on a grid point!");
		if(testy>tolerancey&&testy<(yincrement-tolerancey))
					throw std::invalid_argument("Invalid coordinate Y: Must be on a grid point!");
		if(testz>tolerancez&&testz<(zincrement-tolerancez))
					throw std::invalid_argument("Invalid coordinate Z: Must be on a grid point!");
	}

	//calculate position
	position[0] = static_cast<long int>(floor(tempx/xincrement+0.5));
	position[1] = static_cast<long int>(floor(tempy/yincrement+0.5));
	position[2] = static_cast<long int>(floor(tempz/zincrement+0.5));
	return position;
}

}

#endif /* GRID_HPP_ */
