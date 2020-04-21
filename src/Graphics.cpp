/*
 * Graphics.cpp
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

#include "Graphics.hpp"
#include <string>
#include <cstring>
#include <cstdio>
#include <unistd.h>
#include <initializer_list>

namespace epital{

Graphics::Graphics(): logscalex(false), logscaley(false), data(nullptr) {
	//Create pipe to gnuplot
	if(!(gnuplot = popen("gnuplot -persist", "w"))){
		throw std::runtime_error("Gnuplot initialization error");
	}

	//Initialize with default commands
	ini_command = "plot ";
	style = "with lines ";
	title = "notitle ";
	xrange = "[]";
	yrange= "[] ";
}

Graphics::~Graphics() {
	//Close pipe and gnuplot
	if(pclose(gnuplot)){
		printf("Gnuplot finalization error");
	}
}

void Graphics::Plot(){

	if(logscalex){
		if(logscaley)
			fputs("set logscale xy\n",gnuplot);
		else{
			fputs("set nologscale y\n",gnuplot);
			fputs("set logscale x\n",gnuplot);
		}

	}
	else{
		fputs("set nologscale x\n",gnuplot);
		if(logscaley)
			fputs("set logscale y\n",gnuplot);
	}

	string sendstring;  //String with all commands together

	//Send plot commands and data
	if(data==nullptr){		//If is a string function only
		sendstring = ini_command + range + allcommands + "\n";
		fputs(sendstring.c_str(),gnuplot);
	}
	else{ //If is a plot from data
     	sendstring = ini_command + range + allcommands + "\n";
		fputs(sendstring.c_str(),gnuplot);
		fwrite(data->points_,sizeof(pair<double,double>),data->numpoints,gnuplot); //Send points as binary
	}

}

Graphics::Graphics(string stringfunction): Graphics(){
	dataformat = stringfunction;
	allcommands = dataformat + style + title;
}

Graphics::Graphics(initializer_list<const Graphics*> graphlist): logscalex(false), logscaley(false), data(nullptr){
	//Create pipe to gnuplot
	if(!(gnuplot = popen("gnuplot -persist", "w"))){
		throw std::runtime_error("Gnuplot initialization error");
	}

	//Initialize with default commands
	ini_command = "plot ";
	xrange = "[]";
	yrange= "[] ";

	//Create commands for multiple plot
	for(auto& lst : graphlist){
		allcommands += (lst->allcommands + ", ");
	}
	allcommands.erase(allcommands.end()-2,allcommands.end()); //Erase ", " from last appended string

	//Count number of points from all Graphics
	long int numberofpoints = 0;
	for(auto lst : graphlist){
		if(lst->data != nullptr)
			numberofpoints += lst->data->numpoints;
	}

	//Allocate memory to data
	pair<double, double>* tmpdata = new pair<double, double>[numberofpoints];

	//Extract data from all plots
	long int currpos=0;
	for(auto lst : graphlist){
		if(lst->data != nullptr){
			for(long int localpos = 0; localpos < lst->data->numpoints; localpos++){
				tmpdata[currpos]=lst->data->points_[localpos];
				++currpos;
			}
		}
	}

	//Create data object
	data = new GraphData(tmpdata,numberofpoints);

}

void Graphics::setXrange(double lower, double higher){
	xrange = "[" + to_string(lower) + ":" + to_string(higher) + "]";
	range = xrange + yrange;
}
void Graphics::setYrange(double lower, double higher){
	yrange = "[" + to_string(lower) + ":" + to_string(higher) + "] ";
	range = xrange + yrange;
}

void Graphics::setTitle(string str){
	title = "title '" + str + "' ";
	allcommands = dataformat + style + title;
}
void Graphics::setStyle(Plotstyle styletype){
	switch(styletype){
	case Plotstyle::Points:
		style = "with points ";
		break;
	case Plotstyle::Linepoints:
		style = "with linespoints ";
		break;
	case Plotstyle::Dots:
		style = "with dots ";
		break;
	case Plotstyle::Lines:
		style = "with lines ";
		break;
	default:
		style = "with lines ";
		break;
	}
	allcommands = dataformat + style + title;
}

void Graphics::sendcommands(string str){
	str += "\n";
	fputs(str.c_str(),gnuplot);
}

void Graphics::setXlogscale(bool tof){
	logscalex = tof;
}

void Graphics::setYlogscale(bool tof){
	logscaley = tof;
}

void Graphics::DataConvertionY(function<double(double)> convertionfunction){
	if(data!=nullptr){
		for(long int i = 0;i<data->numpoints;i++){
			data->points_[i]=make_pair(data->points_[i].first,convertionfunction(data->points_[i].second));
		}
	}
	else
		throw runtime_error("Attempt conversion in null data set");
}

void Graphics::DataConvertionX(function<double(double)> convertionfunction){
	if(data!=nullptr){
		for(long int i = 0;i<data->numpoints;i++){
			data->points_[i]=make_pair(convertionfunction(data->points_[i].first),data->points_[i].second);
		}
	}
	else
		throw runtime_error("Attempt conversion in null data set");
}

}
