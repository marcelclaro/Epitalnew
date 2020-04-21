/*
 * PhononTest.cpp
 *
 *  Created on: Dec 6, 2012
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
#include "Solution.hpp"
#include "SolverSO.hpp"
#include "PhononScattering.hpp"
#include "CarrierStatistics.hpp"


using namespace std;
using namespace epital;

#include <unistd.h>
inline void mysleep(unsigned millis) {
	::usleep(millis * 1000);
}

int main(){


	shared_ptr<GaAs> subs = make_shared<GaAs>(77);

	vector<Heterostructure<double>::Epilayer> layers;
	vector<Heterostructure<double>::Epilayer> contact;
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),11.3_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),68.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),56.5_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),20.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),39.55_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),23.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),31.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),28.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),31.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),34.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),31.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),38.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),31.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),48.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),11.3_Angs));

	Heterostructure<double> sample(layers,contact,make_shared<GaAs>(2));
	Grid1D<double> basegrid(sample.Begin(),sample.End(),2000);  //Create a Grid

	//Function o potential
	DiscreteFunction<double,double> potential = sample.ConductionBand_gamma(basegrid);

	Solution<double,double> solution = SolverSO<double,double>::SolvePeriodic(&sample,Band::Conduction,9,100.0_fs);

	//Create graphics
	Graphics pt(solution.getWavefunction(0)->probabilitydensity().Real()+solution.getEnergy(0).real());
	Graphics pt2(solution.getWavefunction(1)->probabilitydensity().Real()+solution.getEnergy(1).real());
	Graphics pt3(solution.getWavefunction(2)->probabilitydensity().Real()+solution.getEnergy(2).real());
	Graphics pt4(solution.getWavefunction(3)->probabilitydensity().Real()+solution.getEnergy(3).real());
	Graphics pt5(solution.getWavefunction(4)->probabilitydensity().Real()+solution.getEnergy(4).real());
	Graphics pt6(solution.getWavefunction(5)->probabilitydensity().Real()+solution.getEnergy(5).real());
	Graphics pt7(solution.getWavefunction(6)->probabilitydensity().Real()+solution.getEnergy(6).real());
	Graphics pt8(solution.getWavefunction(7)->probabilitydensity().Real()+solution.getEnergy(7).real());
	Graphics pt9(solution.getWavefunction(8)->probabilitydensity().Real()+solution.getEnergy(8).real());

	Graphics gama(potential);
	Graphics all({&gama,&pt,&pt2,&pt3,&pt4,&pt5,&pt6,&pt7,&pt8,&pt9});

	for(int i=0; i< 9;i++)
		cout << "Energy-Level " << i << " :" << energy_to_eV(solution.getEnergy(i).real()) << endl;

	all.setXrange(sample.Begin(),sample.End());
	all.DataConvertion(energy_to_eV);
	all.Plot();


	double emmsionrate = PhononScattering<double>::AbsRateLO(80,solution.getWavefunction(1),solution.getWavefunction(0),solution.getEnergy(1).real(),solution.getEnergy(0).real(),std::sqrt(2*0.067*0.0_meV),subs,Band::Conduction);

	cout << "Phonon scattering rate = "<<  frequency_to_SI(emmsionrate) << endl;

	double fermilevel;

	fermilevel = CarrierStatistics<double>::FermiEnergyWells(solution,1.0e+18_under_cubic_cm*58.0_Angs,make_shared<GaAs>(77),77,Band::Conduction);

	cout << "Fermilevel = " << energy_to_eV(fermilevel) << endl;

	double meanemmsionrate = PhononScattering<double>::MeanAbsRateLO(80,solution.getWavefunction(1),solution.getWavefunction(0),solution.getEnergy(1).real(),solution.getEnergy(0).real(),fermilevel,fermilevel,subs,Band::Conduction);

	cout << "Scattering 1-2 (m-2s-1) ="<<  frequency_to_SI(CarrierStatistics<double>::CarriersinSubband(80,solution.getEnergy(0).real(),fermilevel,subs,Band::Conduction)*meanemmsionrate*3.57116e+20) << endl;

	meanemmsionrate = PhononScattering<double>::MeanAbsRateLO(80,solution.getWavefunction(2),solution.getWavefunction(0),solution.getEnergy(2).real(),solution.getEnergy(0).real(),fermilevel,fermilevel,subs,Band::Conduction);

	cout << "Scattering 1-3 (m-2s-1) ="<<  frequency_to_SI(CarrierStatistics<double>::CarriersinSubband(80,solution.getEnergy(0).real(),fermilevel,subs,Band::Conduction)*meanemmsionrate*3.57116e+20) << endl;

	meanemmsionrate = PhononScattering<double>::MeanAbsRateLO(80,solution.getWavefunction(3),solution.getWavefunction(0),solution.getEnergy(3).real(),solution.getEnergy(0).real(),fermilevel,fermilevel,subs,Band::Conduction);

	cout << "Scattering 1-4 (m-2s-1) ="<<  frequency_to_SI(CarrierStatistics<double>::CarriersinSubband(80,solution.getEnergy(0).real(),fermilevel,subs,Band::Conduction)*meanemmsionrate*3.57116e+20) << endl;

	meanemmsionrate = PhononScattering<double>::MeanAbsRateLO(80,solution.getWavefunction(2),solution.getWavefunction(1),solution.getEnergy(2).real(),solution.getEnergy(1).real(),fermilevel,fermilevel,subs,Band::Conduction);

	cout << "Scattering 2-3 (m-2s-1) ="<<  frequency_to_SI(CarrierStatistics<double>::CarriersinSubband(80,solution.getEnergy(1).real(),fermilevel,subs,Band::Conduction)*meanemmsionrate*3.57116e+20) << endl;

	meanemmsionrate = PhononScattering<double>::MeanAbsRateLO(80,solution.getWavefunction(3),solution.getWavefunction(1),solution.getEnergy(3).real(),solution.getEnergy(1).real(),fermilevel,fermilevel,subs,Band::Conduction);

	cout << "Scattering 2-4 (m-2s-1) ="<<  frequency_to_SI(CarrierStatistics<double>::CarriersinSubband(80,solution.getEnergy(1).real(),fermilevel,subs,Band::Conduction)*meanemmsionrate*3.57116e+20) << endl;

	meanemmsionrate = PhononScattering<double>::MeanAbsRateLO(80,solution.getWavefunction(4),solution.getWavefunction(1),solution.getEnergy(4).real(),solution.getEnergy(1).real(),fermilevel,fermilevel,subs,Band::Conduction);

	cout << "Scattering 2-5 (m-2s-1) ="<<  frequency_to_SI(CarrierStatistics<double>::CarriersinSubband(80,solution.getEnergy(1).real(),fermilevel,subs,Band::Conduction)*meanemmsionrate*3.57116e+20) << endl;

	meanemmsionrate = PhononScattering<double>::MeanAbsRateLO(80,solution.getWavefunction(3),solution.getWavefunction(2),solution.getEnergy(3).real(),solution.getEnergy(2).real(),fermilevel,fermilevel,subs,Band::Conduction);

	cout << "Scattering 3-4 (m-2s-1) ="<<  frequency_to_SI(CarrierStatistics<double>::CarriersinSubband(80,solution.getEnergy(2).real(),fermilevel,subs,Band::Conduction)*meanemmsionrate*3.57116e+20) << endl;

	meanemmsionrate = PhononScattering<double>::MeanAbsRateLO(80,solution.getWavefunction(4),solution.getWavefunction(2),solution.getEnergy(4).real(),solution.getEnergy(2).real(),fermilevel,fermilevel,subs,Band::Conduction);

	cout << "Scattering 3-5 (m-2s-1) ="<<  frequency_to_SI(CarrierStatistics<double>::CarriersinSubband(80,solution.getEnergy(2).real(),fermilevel,subs,Band::Conduction)*meanemmsionrate*3.57116e+20) << endl;

	meanemmsionrate = PhononScattering<double>::MeanAbsRateLO(80,solution.getWavefunction(4),solution.getWavefunction(3),solution.getEnergy(4).real(),solution.getEnergy(3).real(),fermilevel,fermilevel,subs,Band::Conduction);

	cout << "Scattering 4-5 (m-2s-1) ="<<  frequency_to_SI(CarrierStatistics<double>::CarriersinSubband(80,solution.getEnergy(3).real(),fermilevel,subs,Band::Conduction)*meanemmsionrate*3.57116e+20) << endl;

	meanemmsionrate = PhononScattering<double>::MeanAbsRateLO(80,solution.getWavefunction(5),solution.getWavefunction(4),solution.getEnergy(5).real(),solution.getEnergy(4).real(),fermilevel,fermilevel,subs,Band::Conduction);

	cout << "Scattering 5-6 (m-2s-1) ="<<  frequency_to_SI(CarrierStatistics<double>::CarriersinSubband(80,solution.getEnergy(4).real(),fermilevel,subs,Band::Conduction)*meanemmsionrate*3.57116e+20) << endl;


/*
Energy-Level 0 :1.57473
Energy-Level 1 :1.60191
Energy-Level 2 :1.62478
Energy-Level 3 :1.64255
Energy-Level 4 :1.66283
Energy-Level 5 :1.68991
Energy-Level 6 :1.71708
Energy-Level 7 :1.72509
Energy-Level 8 :1.8309
Phonon scattering rate = 2.13243e+08
Fermilevel = 1.59338
Scattering 1-2 (m-2s-1) =8.77606e+23
Scattering 1-3 (m-2s-1) =1.84156e+24
Scattering 1-4 (m-2s-1) =4.92232e+22
Scattering 2-3 (m-2s-1) =3.38206e+24
Scattering 2-4 (m-2s-1) =6.67699e+23
Scattering 2-5 (m-2s-1) =1.5558e+21
Scattering 3-4 (m-2s-1) =2.81195e+23
Scattering 3-5 (m-2s-1) =4.09824e+22
Scattering 4-5 (m-2s-1) =9.06194e+21
Scattering 5-6 (m-2s-1) =9.74734e+20
 */


	return 0;
}
