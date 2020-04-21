/*
 * Graphics.hpp
 *
 *  Created on: Jul 27, 2012
 *      Author: marcel
 */

/**
 * @file Graphics.hpp
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

#include <string>
#include <stdexcept>
#include <functional>
#include <unistd.h>
#include <string>
#include <cstring>
#include <cstdio>
#include <utility>
#include <initializer_list>

#include "DFunction.hpp"
#include "PWFunction.hpp"


#ifndef GRAPHICS_HPP_
#define GRAPHICS_HPP_



using namespace std;


namespace epital{

/**
 * @class Graphics
 * @brief Object that plot graphics of functions and data.
 */
class Graphics {
public:

	/**
	 * Cronstructor: Open a process of a Gnuplot a create a Unix pipe to communicate with them.
	 */
	Graphics();

	/**
	 * Create a graphic from others graphics:
	 *
	 * Example: Graphics multiplot({&Graph1, &Graph2, ...});
	 * @param graphlist A list of Graphics addresses.
	 */
	Graphics(initializer_list<const Graphics*> graphlist);

	enum Complexplot{
		REAL,
		IMAGINARY,
		SQUARED
	};

	/**
	 * Create a Graphics from a function string (of gnuplot).
	 *
	 * Example: Graphics plot("cos(x)*sin(x)");
	 * @param stringfunction string of function.
	 */
	Graphics(string stringfunction);

	/**
	 * Create a Graphics from a standard  c++ function/method.
	 * @param thefunction Function to be ploted.
	 * @param rangelow Lower "x" point of function.
	 * @param rangehigh Higher "x" point of function.
	 * @param numberofpoints Number of points ploted.
	 */
	template<typename T>
	Graphics(std::function<T(double)> thefunction, double rangelow, double rangehigh, long int numberofpoints=200);


	/**
	 * Create a Graphics from a discrete function.
	 * @param thefunction Discrete function to be ploted
	 */
	template<typename yType, typename xType>
	Graphics(const DiscreteFunction<yType,xType> &thefunction);

	/**
	 * Create a Graphics from a plane wavefunction.
	 * @param thefunction planewavefunction to be ploted
	 */
	template<typename yType, typename xType>
	Graphics(const PlaneWavesFunction<yType,xType> &thefunction, long int numberofpoints, Graphics::Complexplot plottype, yType yoffset=0);


	/**
	 * Standard destructor
	 */
	virtual ~Graphics();

	/**
	 * Enumeration for line plot style.
	 */
	enum class Plotstyle {Lines, Points, Linepoints, Dots};

	/**
	 * Set data title.
	 * @param str string with data title.
	 */
	void setTitle(string str);

	/**
	 * Set plot X range
	 * @param lower lower limit
	 * @param higher higher limit
	 */
	void setXrange(double lower, double higher);

	/**
	 * Set plot Y range
	 * @param lower lower limit
	 * @param higher higher limit
	 */
	void setYrange(double lower, double higher);

	/**
	 * Set line style.
	 * @param styletype The style
	 */
	void setStyle(Plotstyle styletype);

	/**
	 * Send a command to gnuplot.
	 * @param str string of command.
	 */
	void sendcommands(string str);

	/**
	 * Set logaritmic scale in x axis
	 * @param tof true of false
	 */
	void setXlogscale(bool tof=true);

	/**
	* Set logaritmic scale in y axis
	* @param tof true of false
	*/
	void setYlogscale(bool tof=true);

	/**
	 * Plot Graphics.
	 */
	void Plot();

	/**
	 * Convert data (except Graphics from string) of graphic using a conversion function
	 * @param convertionfunction the conversion function
	 */
	void DataConvertionY(function<double(double)> convertionfunction);

	/**
	 * Convert data (except Graphics from string) of graphic using a conversion function
	 * @param convertionfunction the conversion function
	 */
	void DataConvertionX(function<double(double)> convertionfunction);

	/**
	 * @class GraphData
	 * @brief store x-y data to plot
	 */
	class GraphData{
	public:
		GraphData(): points_(nullptr), numpoints(0){};
		GraphData(pair<double,double> *points, long int pointsnumber): points_(points), numpoints(pointsnumber){};
		GraphData(long int pointsnumber): numpoints(pointsnumber){ points_ = new pair<double,double>[numpoints]; };
		~GraphData(){delete [] points_;};
	protected:
		pair<double,double> *points_; ///<Raw Plot data
		long int numpoints;  ///< Number of points stored

		friend class Graphics;
	};

protected:
	FILE* gnuplot;  ///< Unix-Pipe to gnuplot.
	string ini_command; ///< string to initial command of plot ( "plot" );
	string range; ///< string with range command.
	string xrange; ///< string with x range to compose range string.
	string yrange; ///< string with y range to compose range string.
	string title; ///< string with data title.
	string style; ///< string with style command.
	string dataformat;  ///< string with data format or function string.
	string allcommands; ///< string with all plot commands ( dataformat + style + title ).
	bool logscalex;  ///<Indicate logaritmic scale in x axis
	bool logscaley; ///<Indicate logaritmic scale in y axis
	GraphData* data; ///< Data to be ploted

	string generatedCommands();

};


template<typename T>
Graphics::Graphics(std::function<T(double)> thefunction, double rangelow, double rangehigh, long int numberofpoints): Graphics(){

	//Allocate memory to data
	pair<double, double>* tmpdata = new pair<double, double>[numberofpoints];

	//Extract data from function
	double point = rangelow;
	double step = (rangehigh-rangelow)/numberofpoints;
	for(long int i = 0;i<numberofpoints;i++){
		tmpdata[i]=make_pair(point,static_cast<double>(thefunction(point)));
		point+=step;
	}
	//Create data object
	data = new GraphData(tmpdata,numberofpoints);

	//Configure associated commands
	dataformat = string("'-' using 1:2 binary format='%double%double' ") + " record=" + to_string(data->numpoints) + "  ";
	allcommands = dataformat + style + title;
}

template<typename yType, typename xType>
Graphics::Graphics(const DiscreteFunction<yType,xType> &thefunction): Graphics(){

	//Allocate memory to data
	long int numberofpoints = thefunction.getGrid().getSize();
	pair<double, double>* tmpdata = new pair<double, double>[numberofpoints];

	//Extract data from function
	double point = thefunction.getGrid().getInferiorLimit();
	double step = thefunction.getGrid().getIncrement();
	for(long int i = 0;i<numberofpoints;i++){
		tmpdata[i]=make_pair(point,static_cast<double>(const_cast<DiscreteFunction<yType,xType> &>(thefunction)(i)));
		point+=step;
	}

	//Create data object
	data = new GraphData(tmpdata,numberofpoints);

	//Configure associated commands
	dataformat = string("'-' using 1:2 binary format='%double%double' ") + " record=" + to_string(data->numpoints) + "  ";
	allcommands = dataformat + style + title;

}

template<typename complextype, typename xType>
Graphics::Graphics(const PlaneWavesFunction<complextype,xType> &thefunction, long int numberofpoints, Graphics::Complexplot plottype, complextype yoffset): Graphics(){

	//Allocate memory to data
	pair<double, double>* tmpdata = new pair<double, double>[numberofpoints];

	//Extract data from function
	double point = thefunction.getBounds().first;
	double step = (thefunction.getBounds().second-thefunction.getBounds().first)/(numberofpoints+1);
	for(long int i = 0;i<numberofpoints;i++){
		if(plottype==Complexplot::REAL)
			tmpdata[i]=make_pair(point,static_cast<double>(std::real(thefunction(point)))+yoffset);
		if(plottype==Complexplot::IMAGINARY)
			tmpdata[i]=make_pair(point,static_cast<double>(std::imag(thefunction(point)))+yoffset);
		if(plottype==Complexplot::SQUARED)
			tmpdata[i]=make_pair(point,static_cast<double>(std::real((std::conj(thefunction(point))*thefunction(point))))+yoffset);
		point+=step;
	}

	//Create data object
	data = new GraphData(tmpdata,numberofpoints);

	//Configure associated commands
	dataformat = string("'-' using 1:2 binary format='%double%double' ") + " record=" + to_string(data->numpoints) + "  ";
	allcommands = dataformat + style + title;

}


}
#endif /* GRAPHICS_HPP_ */
