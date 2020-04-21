/*
 * WannierFunction.cpp
 *
 *  Created on: Aug 2, 2013
 *      Author: marcel
 */

#include "WannierFunction.hpp"

#ifndef WANNIERFUNCTION_CPP_
#define WANNIERFUNCTION_CPP_

namespace epital {

template<typename complextype, typename xType>
WannierFunction<complextype,xType>::WannierFunction(xType periodsize, xType vectorstep): levels(0), periodsize(periodsize), vectorstep(vectorstep){
}

template<typename complextype, typename xType>
WannierFunction<complextype,xType>::~WannierFunction(){

}


template<typename complextype, typename xType>
void WannierFunction<complextype,xType>::addFunction(shared_ptr<PlaneWavesFunction<complextype,xType>> wavefunction,std::complex<complextype> energy){


	waves.push_back(wavefunction);
	energies.push_back(energy);

	levels++;

}

template<typename complextype, typename xType>
shared_ptr<PlaneWavesFunction<complextype,xType>> WannierFunction<complextype,xType>::getWavefunction(long int number){
	if(number >= levels)
		throw std::invalid_argument("Wrong number of solutions in WannierFunction<complextype,xType>::getWavefunction");

	return waves[number];
}

template<typename complextype, typename xType>
std::complex<complextype> WannierFunction<complextype,xType>::getEnergy(long int number){
	if(number > levels)
		throw std::invalid_argument("Wrong number of solutions in WannierFunction<complextype,xType>::getEnergy");

	return energies[number];
}

template<typename complextype, typename xType>
long int WannierFunction<complextype,xType>::getLevels(){
	return levels;
}

template<typename complextype, typename xType>
xType WannierFunction<complextype,xType>::getPeriod(){
	return periodsize;
}

template<typename complextype, typename xType>
xType WannierFunction<complextype,xType>::getBlochVectorStep(){
	return vectorstep;
}

template<typename complextype, typename xType>
std::complex<complextype> WannierFunction<complextype,xType>::getValue(long int period, xType position){
	position=position-period*periodsize;
	std::complex<complextype> value(0.0,0.0);
	for(long int i = 0; i < levels-1; i++){
		value+=vectorstep*((*waves[i])(position)*std::exp(std::complex<complextype>(0.0,position*(i*vectorstep-Constant::pi/periodsize)))+(*waves[i+1])(position)*std::exp(std::complex<complextype>(0.0,position*((i+1)*vectorstep-Constant::pi/periodsize))))/2.0;
	}

	return sqrt(periodsize/(Constant::pi*2))*value;
}
/*
template<typename complextype, typename xType>
std::complex<complextype> WannierFunction<complextype,xType>::getValue_efloat(long int period, xType position){
	position=position-period*periodsize;
	std::complex<complextype> value(0.0);
	for(long int i = 0; i < levels-1; i++){
		value+=vectorstep*((*waves[i])[position]*std::exp(std::complex<complextype>(0.0,position*(i*vectorstep-Constant::pi/periodsize)))+(*waves[i+1])[position]*std::exp(std::complex<complextype>(0.0,position*((i+1)*vectorstep-Constant::pi/periodsize))))/2.0;
	}

		return sqrt((periodsize/(Constant::pi*2.0)))*value;
}
*/

template<typename complextype, typename xType>
DiscreteFunction<std::complex<complextype>,xType> WannierFunction<complextype,xType>::getDiscrete(long int period, long int periodsback, long int periodsfar, long int points, bool normalize){

	std::function<std::complex<complextype>(xType)> funcpointer = std::bind(std::mem_fn(&WannierFunction<complextype,xType>::getValue),this, period, std::placeholders::_1);

	Grid1D<xType> grid(period*periodsize-periodsback*periodsize,period*periodsize+(periodsfar+1)*periodsize,points);

	if(normalize)
		return DiscreteFunction<std::complex<complextype>,xType>(grid,funcpointer).normalize(1.0);
	else
		return DiscreteFunction<std::complex<complextype>,xType>(grid,funcpointer);
}

/*
template<typename complextype, typename xType>
bool Solution<complextype,xType>::saveToFile(std::string filename){
	file=filename;
	std::ofstream myfile;
	myfile.open(file, ios::out | ios::trunc | ios::binary);
	size_t xtype = typeid(xType).hash_code();
	size_t ytype = typeid(complextype).hash_code();
	if(myfile.is_open()){
		myfile.write(reinterpret_cast<const char*>(&levels),sizeof(long int));
		myfile.write( reinterpret_cast<const char*>(&xtype),sizeof(size_t));
		myfile.write( reinterpret_cast<const char*>(&ytype),sizeof(size_t));
		Grid1D<xType> grid(0,0,1);
		std::complex<complextype> energy;
		std::complex<complextype>* data;
		auto enitr = energies.begin();
		auto wavitr = waves.begin();
		for(long int i=0 ; i < levels ; i++){
			grid=(*wavitr)->getGrid();
			energy=*enitr;
			data=(*wavitr)->getdata();
			myfile.write(reinterpret_cast<const char*>(&energy),sizeof(std::complex<complextype>));
			myfile.write(reinterpret_cast<const char*>(&grid),sizeof(const Grid1D<xType>));
			myfile.write(reinterpret_cast<const char*>(data),grid.getSize()*sizeof(std::complex<complextype>));
			enitr++;
			wavitr++;
		}
		myfile.close();
		return true;
	}
	else{
		return false;
	}

}

template<typename complextype, typename xType>
Solution<complextype,xType>::Solution(std::string filename){
	file=filename;
	levels=0;
	std::ifstream myfile;
	myfile.open(file, ios::in | ios::binary);
	size_t xtype = typeid(xType).hash_code();
	size_t ytype = typeid(complextype).hash_code();
	size_t filextype;
	size_t fileytype;
	Grid1D<xType> grid(0,0,1);
	std::complex<complextype> energy;
	std::complex<complextype>* data;
	long int levelnumber = 0;
	if(myfile.is_open()){
		myfile.read(reinterpret_cast<char*>(&levelnumber),sizeof(long int));
		myfile.read(reinterpret_cast<char*>(&filextype),sizeof(size_t));
		myfile.read(reinterpret_cast<char*>(&fileytype),sizeof(size_t));
		if(xtype!=filextype||ytype!=fileytype)
			throw std::runtime_error("Incompatible xType and/or yType file");

		for(long int i=0 ; i < levelnumber ; i++){
			myfile.read(reinterpret_cast<char*>(&energy),sizeof(std::complex<complextype>));
			myfile.read(reinterpret_cast<char*>(const_cast<Grid1D<xType>*>(&grid)),sizeof(Grid1D<xType>));
			data = new std::complex<complextype>[grid.getSize()];
			myfile.read(reinterpret_cast<char*>(data),grid.getSize()*sizeof(std::complex<complextype>));
			shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunction =
					std::make_shared<DiscreteFunction<std::complex<complextype>,xType>>(data,grid);
			addLevel(wavefunction,energy);
		}


		myfile.close();
	}
	else
		throw std::runtime_error("Function File error");


}
*/

} /* namespace epital */

#endif
