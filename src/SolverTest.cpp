/*
 * SolverTest.cpp
 *
 *  Created on: Dec 3, 2012
 *      Author: marcel
 */


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

#include <CustomTypes.hpp>
#include <CMaterial.hpp>
#include <Materials.hpp>
#include <Units.hpp>
#include <Alloy.hpp>
#include <Grid.hpp>
#include <DFunction.hpp>
#include <Graphics.hpp>
#include <Heterostructure.hpp>
#include <WaveTimeEvolution.hpp>
#include <Solution.hpp>
#include <SolverSO.hpp>
//#include <Eigen/Dense>

using namespace std;

#include <unistd.h>
inline void mysleep(unsigned millis) {
	::usleep(millis * 1000);
}

int main(){

	using namespace epital;

	GaAs subs(77);

	vector<Heterostructure<double>::Epilayer> layers;
	vector<Heterostructure<double>::Epilayer> contact;
	/*layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.44),11.3_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),100.0_Angs));
	/*layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),21.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.14),6.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),48.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.14),6.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),21.0_Angs));
	/*98A GaAs*//*
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.44),56.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.14),21.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.44),39.55_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.14),30.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.44),31.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.14),43.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.44),31.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.14),62.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.44),25.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),31.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.44),11.3_Angs));*/
	//contact.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.99),300.0_Angs));
	//contact.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.99),300.0_Angs));


	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),11.3_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),68.0_Angs,1.0e+18_under_cubic_cm));
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

	Solution<double,double> solution = SolverSO<double,double>::SolvePeriodic(&sample,Band::Conduction,9,0.0,200.0_fs,0.01_fs);

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

	for(int i=0; i< 8;i++)
		cout << "Energy-Level " << i << " :" << energy_to_eV(solution.getEnergy(i).real()) << endl;



	all.setXrange(sample.Begin(),sample.End());
	all.DataConvertion(energy_to_eV);
	all.Plot();


	return 0;
}
