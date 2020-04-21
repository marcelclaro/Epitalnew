/*
 * lab3D.cpp
 *
 *  Created on: Nov 6, 2014
 *      Author: marcel
 */

//============================================================================
// Name        : Epital.cpp
// Author      : Marcel Claro
// Version     :
// Copyright   : Your copyright notice
// Description :
//============================================================================



#include <thread>
#include <iostream>
#include <limits>
#include <string>
#include <cstring>
#include <cstdio>
#include <vector>
#include <complex>
#include <chrono>
#include <ratio>
#include <memory>

#include "../CustomTypes.hpp"
#include "../CMaterial.hpp"
#include "../Materials.hpp"
#include "../Units.hpp"
#include "../Alloy.hpp"
#include "../Grid.hpp"
#include "../DFunction.hpp"
#include "../Graphics.hpp"
#include "../Heterostructure.hpp"
#include "../Heterostructure3D.hpp"
#include "../Solution.hpp"
#include "../WaveTimeEvolution.hpp"
#include "../SolverSO.hpp"
#include "../SolverTM.hpp"
//#include "PhononScattering.hpp"
//#include "CarrierStatistics.hpp"
//#include "QCD.hpp"
#include "../PerturbationMethod.hpp"
#include "../DFunction3D.hpp"
#include "../WaveTimeEvolution3D.hpp"
//#include "strainCalculator.hpp"

//#include <o2scl/funct.h>
//#include <o2scl/root_brent_gsl.h>
//#include <o2scl/root_cern.h>

//#include <Eigen/Dense>

using namespace std;
using namespace epital;


int main(){

	cout.precision(15);


	std::cout << "Loading geometry..." << std::endl;
	shared_ptr<GaAs> subs = make_shared<GaAs>(77);
	//Heterostructure3D<double>* sample("./config.hetero3d",subs,2);
	Heterostructure3D<double>* sample = new Heterostructure3D<double>("./config.hetero3d",subs,77);
	//sample.printPotential("potential.vti");

	std::cout << "Getting Potential..." << std::endl;
	DiscreteFunction3D<double,double> potential = sample->ConductionBand_gamma();
	DiscreteFunction3D<double,double> effective_mass = sample->effectiveMass_gamma();



/*
	DiscreteFunction3D<double,double> potential = sample->ValenceBand_HH();
	DiscreteFunction3D<double,double> effective_masst = sample->effectiveMass_HHt();
	DiscreteFunction3D<double,double> effective_massl = sample->effectiveMass_HHl();
*/

	potential.saveHDF5("potential.hdf");
	//effective_masst.saveHDF5("mass.hdf");

	std::cout << "Filling initial function..." << std::endl;
	std::shared_ptr<DiscreteFunction3D<std::complex<double>,double>> initial = make_shared<DiscreteFunction3D<std::complex<double>,double>>(sample->getGrid());
	//std::shared_ptr<DiscreteFunction3D<std::complex<double>,double>> initial2 = make_shared<DiscreteFunction3D<std::complex<double>,double>>(sample.getGrid());
	std::shared_ptr<DiscreteFunction3D<std::complex<double>,double>> initial3 = make_shared<DiscreteFunction3D<std::complex<double>,double>>(sample->getGrid());
	//std::shared_ptr<DiscreteFunction3D<std::complex<double>,double>> initial4 = make_shared<DiscreteFunction3D<std::complex<double>,double>>(sample.getGrid());
	initial->fillGaussian(1000,1);
	//initial2->fillGaussian2(1000,1);
	initial3->fillGaussian2(1000,1);
	//initial4->fillGaussian2(1000,1);

	delete sample;

	std::vector<std::shared_ptr<DiscreteFunction3D<std::complex<double>,double>>> waves;

	WaveTimeEvolution3D<std::complex<double>,double> evoluter(WaveTimeEvolution3D<std::complex<double>,double>::Carrier::ELECTRON);
	//WaveTimeEvolution3D<std::complex<double>,double> evoluter(WaveTimeEvolution3D<std::complex<double>,double>::Carrier::HOLE);

/*
	DiscreteFunction3D<double,double> ux(sample.getGrid());
	DiscreteFunction3D<double,double> uy(sample.getGrid());
	DiscreteFunction3D<double,double> uz(sample.getGrid());

	strainCalculator<double>::solveDisplace(sample,ux,uy,uz);

	DiscreteFunction3D<double,double> exx(sample.getGrid());
	DiscreteFunction3D<double,double> exy(sample.getGrid());
	DiscreteFunction3D<double,double> exz(sample.getGrid());
	DiscreteFunction3D<double,double> ezx(sample.getGrid());
	DiscreteFunction3D<double,double> ezy(sample.getGrid());
	DiscreteFunction3D<double,double> ezz(sample.getGrid());


	ux.Divergent(exx,exy,exz);
	uz.Divergent(ezx,ezy,ezz);

	exx.saveHDF5("exx.hdf");
	exy.saveHDF5("exy.hdf");
	exy.saveHDF5("ezz.hdf");
*/

	std::cout << "calculating..." << std::endl;
/*
	for(long int steps = 0; steps< 81; ++steps){
	evoluter.groundlevel_finite_valence(std::complex<double>(0,-1.0_fs), 12, *initial, potential, effective_masst,effective_massl);
	cout << "Nivel 1 - finite: "<< energy_to_eV(real(evoluter.Energy(*initial,potential,effective_masst,effective_massl))) << std::endl;
	}
*/

	for(long int steps = 0; steps< 31; ++steps){
	evoluter.groundlevel_finite(std::complex<double>(0,-1.0_fs), 5, *initial, potential, effective_mass);
	cout << "Nivel 1 - finite: "<< energy_to_eV(real(evoluter.Energy(*initial,potential,effective_mass))) << std::endl;
	}


/*
	for(long int steps = 0; steps< 15; ++steps){
	evoluter.groundlevel_periodicFFT(std::complex<double>(0,-1.0_fs), 12, *initial, potential, InAs(77).effectiveMass_gamma());
	cout << "Nivel 1 - finite: "<< energy_to_eV(real(evoluter.Energy(*initial,potential,effective_mass))) << std::endl;
	}
*/
/*
	waves.push_back(initial);
	for(long int steps = 0; steps< 81; ++steps){
	evoluter.newlevel_finite(std::complex<double>(0,-1.0_fs), 5, *initial3, potential,waves, effective_mass);
	cout << "Nivel 2 - finite: "<< energy_to_eV(real(evoluter.Energy(*initial3,potential,effective_mass))) << std::endl;
	}
*/
	/*
	waves.push_back(initial);
	evoluter.newlevel_periodicFFT(std::complex<double>(0,-25.0_fs), 15*25, *initial2, potential,waves, InAs(2).effectiveMass_gamma());
	cout << "Nivel 2: "<< energy_to_eV(real(evoluter.Energy(*initial2,potential,effective_mass))) << std::endl;
*/

	/*waves.push_back(initial2);
	evoluter.newlevel_periodicFFT(std::complex<double>(0,-200.0_fs), 250, *initial3, potential,waves, GaAs(77).effectiveMass_gamma());
	cout << "Nivel 3: "<< energy_to_eV(real(evoluter.Energy(*initial3,potential,effective_mass))) << std::endl;

	/*waves.push_back(initial3);
	evoluter.newlevel_periodicFFT(std::complex<double>(0,-200.0_fs), 250, *initial4, potential,waves, InAs(77).effectiveMass_gamma());
	 */


	DiscreteFunction3D<double,double> density = initial->getDensity();
	//DiscreteFunction3D<double,double> density2 = initial2->getDensity();
	DiscreteFunction3D<double,double> density3 = initial3->getDensity();
	/*DiscreteFunction3D<double,double> density4 = initial4->getDensity();*/

	std::cout << "Saving wavefunction..." << std::endl;
	density.saveHDF5("wavefunctiondens.hdf");
	//density2.saveHDF5("wavefunctiondens2.hdf");
	//density3.saveHDF5("wavefunctiondens3.hdf");
	/*density4.saveHDF5("wavefunctiondens4.hdf");*/
/*
	ux.saveHDF5("ux.hdf");
	uy.saveHDF5("uy.hdf");
	uz.saveHDF5("uz.hdf");
*/
	return 0;
}
