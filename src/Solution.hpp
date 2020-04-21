/*
 * Solution.hpp
 *
 *  Created on: Nov 5, 2012
 *      Author: marcel
 */

#include <thread>
#include <iostream>
#include <limits>
#include <string>
#include <cstring>
#include <cstdio>
#include <list>
#include <complex>
#include <chrono>
#include <ratio>
#include <fstream>

#include "CustomTypes.hpp"
#include "CMaterial.hpp"
#include "Materials.hpp"
#include "Units.hpp"
#include "Alloy.hpp"
#include "Grid.hpp"
#include "DFunction.hpp"
#include "Graphics.hpp"
#include "Heterostructure.hpp"
#include "WaveTimeEvolution.hpp"

#ifndef SOLUTION_HPP_
#define SOLUTION_HPP_

namespace epital{
/**
 * \brief This class are used as result of a solve method, it has the wavefunctions and energies.
 */
template<typename complextype, typename xType>
class Solution{
public:
	Solution();
	virtual ~Solution();

	Solution(string filename);

	/**
	 * Add a level to solution, append in order of energy
	 * @param wavefunction wavefunction to add
	 * @param energy corresponding energy
	 */
	void addLevel(shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunction,std::complex<complextype> energy);

	/**
	 * Get wavefunction
	 * @param number number of level 0=Ground State
	 * @return Function pointer
	 */
	shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> getWavefunction(long int number);

	/**
	 * Get energy of level
	 * @param number number of the level 0=Ground state
	 * @return
	 */
	std::complex<complextype> getEnergy(long int number);


	/**
	 * Get number of levels.
	 * @return return number of levels.
	 */
	long int getLevels();

	bool saveToFile(string filename);

protected:
	long int levels; ///< Number of levels.
	string file;

	list<shared_ptr<DiscreteFunction<std::complex<complextype>,xType>>> waves; ///<Pointers to wavefunctions
	list<std::complex<complextype>> energies; ///< Energies of founded wavefunctions

};

template<typename complextype, typename xType>
Solution<complextype,xType>::Solution(): levels(0){
}

template<typename complextype, typename xType>
Solution<complextype,xType>::~Solution(){

}


template<typename complextype, typename xType>
void Solution<complextype,xType>::addLevel(shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> wavefunction,std::complex<complextype> energy){
	if( (energy*conj(energy)).real() > (energies.back()*conj(energies.back())).real()  ){
		energies.push_back(energy);
		waves.push_back(wavefunction);
	}
	else{
		auto i = energies.begin();
		auto j = waves.begin();
		while( (energy*conj(energy)).real() > ((*i)*conj(*i)).real() ){
			i++;
			j++;
		}
		energies.insert(i,energy);
		waves.insert(j,wavefunction);

	}

	levels++;

}

template<typename complextype, typename xType>
shared_ptr<DiscreteFunction<std::complex<complextype>,xType>> Solution<complextype,xType>::getWavefunction(long int number){
	if(number >= levels)
		throw std::invalid_argument("Wrong number of solutions");

	auto j = waves.begin();
	for(long int n = 1; n<=number; n++){
		j++;
	}
	return *j;
}

template<typename complextype, typename xType>
std::complex<complextype> Solution<complextype,xType>::getEnergy(long int number){
	if(number > levels)
		throw std::invalid_argument("Wrong number of solutions");

	auto i = energies.begin();
		for(long int n = 1; n<=number; n++){
			i++;
		}
		return *i;
}

template<typename complextype, typename xType>
long int Solution<complextype,xType>::getLevels(){
	return levels;
}

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


}

#endif /* SOLUTION_HPP_ */
