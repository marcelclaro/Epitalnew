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
#include "SolverTM.hpp"
#include "PhononScattering.hpp"
#include "CarrierStatistics.hpp"
#include "QCD.hpp"
#include "PerturbationMethod.hpp"

#include <o2scl/funct.h>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/root_cern.h>

//#include <Eigen/Dense>

using namespace std;
using namespace epital;

#include <unistd.h>
inline void mysleep(unsigned millis) {
	::usleep(millis * 1000);
}

int main(){

	cout.precision(15);

	shared_ptr<GaAs> subs = make_shared<GaAs>(77);

	vector<Heterostructure<double>::Epilayer> layers;
	vector<Heterostructure<double>::Epilayer> contact;

	/*layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),11.3_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),70.0_Angs,1.0e+18_under_cubic_cm));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),56.5_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.10),37.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),39.55_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.10),47.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),31.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.10),70.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),31.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),34.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),31.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),38.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),31.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),48.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),11.3_Angs));*/

	/*layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.47),9.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaAs>(2),30.0_Angs,3e+18_under_cubic_cm));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.47),33.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.26),20.5_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.47),27.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.26),25.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.47),24.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.26),41.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.47),20.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaAs>(2),14.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.47),28.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaAs>(2),17.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.47),39.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaAs>(2),24.0_Angs,3e+18_under_cubic_cm));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.47),9.0_Angs));*/
/*
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.44),11.3_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),5.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),88.0_Angs,1.0e+18_under_cubic_cm));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),5.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.44),56.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.14),21.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.44),39.6_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.14),30.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.44),31.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.14),43.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.44),31.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.14),62.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.44),25.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),31.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.44),11.3_Angs));
*/
/*
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),12.5_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),65.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.33),12.5_Angs));
*/

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
	/*layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.20),50.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.0),50.0_Angs))*/;
	//layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(2,0.20),25.0_Angs));

	Heterostructure<double> sample(layers,contact,subs);


	//Solution<double,double> solution = SolverSO<double,double>::SolvePeriodic(&sample,Band::Conduction,9,100.0_fs);

	QCD<double,double> qcd(layers,contact,subs);

	Grid1D<double> basegrid(sample.Begin(),sample.End(),2000);  //Create a Grid
	DiscreteFunction<double,double> potential = sample.ConductionBand_gamma(basegrid);
	DiscreteFunction<double,double> mass = sample.effectiveMass_gamma(basegrid);

	/*potential.savetoFile("Solution/potential.func");

	DiscreteFunction<double,double> potential2("Solution/potential.func");

	Graphics test(potential2);
	test.Plot();*/

	Solution<double,double> solutionfst = SolverSO<double,double>::SolvePeriodic(sample,Band::Conduction,9,200.0_fs,0.05_fs,600);
	/*solutionfst.saveToFile("solution3620.sol");*/

	/*Solution<double,double> solution("solution3620.sol");*/

	//Create graphics
	Graphics pt(solutionfst.getWavefunction(0)->probabilitydensity().Real()+solutionfst.getEnergy(0).real());
	Graphics pt2(solutionfst.getWavefunction(1)->probabilitydensity().Real()+solutionfst.getEnergy(1).real());
	Graphics pt3(solutionfst.getWavefunction(2)->probabilitydensity().Real()+solutionfst.getEnergy(2).real());
	Graphics pt4(solutionfst.getWavefunction(3)->probabilitydensity().Real()+solutionfst.getEnergy(3).real());
	Graphics pt5(solutionfst.getWavefunction(4)->probabilitydensity().Real()+solutionfst.getEnergy(4).real());
	Graphics pt6(solutionfst.getWavefunction(5)->probabilitydensity().Real()+solutionfst.getEnergy(5).real());
	Graphics pt7(solutionfst.getWavefunction(6)->probabilitydensity().Real()+solutionfst.getEnergy(6).real());
	Graphics pt8(solutionfst.getWavefunction(7)->probabilitydensity().Real()+solutionfst.getEnergy(7).real());
	Graphics pt9(solutionfst.getWavefunction(8)->probabilitydensity().Real()+solutionfst.getEnergy(8).real());
	//Graphics pt10(solution.getWavefunction(9)->probabilitydensity().Real()+solution.getEnergy(9).real());
	//Graphics pt11(solution.getWavefunction(10)->probabilitydensity().Real()+solution.getEnergy(10).real());
	/*Graphics pt12(solution.getWavefunction(11)->probabilitydensity().Real()+solution.getEnergy(11).real());
	Graphics pt13(solution.getWavefunction(12)->probabilitydensity().Real()+solution.getEnergy(12).real());
	Graphics pt14(solution.getWavefunction(13)->probabilitydensity().Real()+solution.getEnergy(13).real());
	Graphics pt15(solution.getWavefunction(14)->probabilitydensity().Real()+solution.getEnergy(14).real());
	Graphics pt16(solution.getWavefunction(15)->probabilitydensity().Real()+solution.getEnergy(15).real());
	Graphics pt17(solution.getWavefunction(16)->probabilitydensity().Real()+solution.getEnergy(16).real());
	Graphics pt18(solution.getWavefunction(17)->probabilitydensity().Real()+solution.getEnergy(17).real());
	Graphics pt19(solution.getWavefunction(18)->probabilitydensity().Real()+solution.getEnergy(18).real());
	Graphics pt20(solution.getWavefunction(19)->probabilitydensity().Real()+solution.getEnergy(19).real());*/

	Graphics gama(potential);
	Graphics all({&gama,&pt,&pt2,&pt3,&pt4,&pt5,&pt6,&pt7,&pt8,&pt9/*,&pt10/*,&pt11,&pt12,&pt13,&pt14,&pt15,&pt16,&pt17,&pt18,&pt19,&pt20*/});

	for(int i=0; i< 9;i++)
		cout << "Energy-Level " << i << " :" << energy_to_eV(solutionfst.getEnergy(i).real()) << endl;

	all.setXrange(sample.Begin(),sample.End());
	all.DataConvertionY(energy_to_eV);
	all.Plot();


	//cout << "100/0.01	200/0.01	300/0.01	400/0.01	200/0.005	200/0.001	200/0.0005" << endl;

	/*Solution<double,double> solution = SolverSO<double,double>::SolvePeriodic(&sample,Band::Conduction,8,0.0,200.0_fs,0.01_fs,1690);/
	double fermilevel = CarrierStatistics<double>::FermiEnergyWells(solution,sample.getTotalDoping(),sample.getSubstrate(),77,Band::Conduction);
	cout << "Fermi:" <<energy_to_eV(fermilevel) << endl;


	/*for(double j = 0.0; j < 2; j+=0.1){
		Solution<double,double> solutione = SolverSO<double,double>::SolvePeriodic(&sample,Band::Conduction,1,j*Constant::pi/(sample.getTotalWidth()),400.0_fs,0.005_fs,100);
		cout << energy_to_eV(solutione.getEnergy(0).real()) << endl;

	}
	*/

	//PerturbationMethod<std::complex<double>,double> method(PerturbationMethod<std::complex<double>,double>::Carrier::ELECTRON);
	//cout << energy_to_eV(std::real(method.EnergyBrilluoin(0,solution,potential,mass,0.2*(Constant::pi/501.65_Angs)))) << endl;
	//cout << energy_to_eV(std::real(method.EnergyBrilluoin(0,solution,potential,mass,0.0*(Constant::pi/501.65_Angs)))) << endl;


	//Solution<double,double> solutione = SolverSO<double,double>::SolvePeriodic(&sample,Band::Conduction,20,200.0_fs,0.01_fs,1000);
	//solutione.saveToFile("solution.sol");

	/*Solution<double,double> solutionee("solution.sol");

	cout << energy_to_eV(solutionee.getEnergy(8).real()) << endl;

	for(double j = 0.0; j < 1; j+=0.05){
		cout << energy_to_eV(std::real(method.EnergyBrilluoin(6,solutionee,mass,j*(Constant::pi/sample.getTotalWidth())))) << endl;

	}


	cout << endl << "2º" << endl;
	cout << energy_to_eV(solutionee.getEnergy(9).real()) << endl;

	for(double j = 0.0; j < 1; j+=0.05){
		cout << energy_to_eV(std::real(method.EnergyBrilluoin(7,solutionee,mass,j*(Constant::pi/sample.getTotalWidth())))) << endl;

	}
	 */

	//qcd.setSolveParams(Band::Conduction,10);

	//for(Temperature temp = 80; temp<160; temp+=20){
		//double resistance = qcd.ResistanceArea(temp);
		//cout << "Resistance at " << temp << "oC = "<< resistance*7.19015e-11 << endl;
	//}

	/*SolverTM<double,double> solver(sample,Band::Conduction);

	auto func = std::mem_fn(&SolverTM<double,double>::TraceTM);
	std::function<double(double)> funcpointer = std::bind(func, solver, std::placeholders::_1);

	/*Graphics plot(funcpointer,energy_from_eV(1.512),energy_from_eV(1.85),10000000);

	plot.setYrange(-1.0,1.0);
	plot.DataConvertionX(energy_to_eV);
	plot.Plot();*/

	/*double root=energy_from_eV(1.798);
	double rootmax=energy_from_eV(1.802);
	o2scl::gsl_root_brent<> solverr;
	o2scl::funct_mfptr<SolverTM<double,double>> tosolve(&solver,&SolverTM<double,double>::TraceTM);
	solverr.solve_bkt(root,rootmax,tosolve);

	cout << energy_to_eV(root) << " (root0)= " << funcpointer(root)<< endl;

	double root2=energy_from_eV(1.798);
	o2scl:: gsl_root_brent<> solverb;
	o2scl::funct_mfptr<SolverTM<double,double>> tosolve2(&solver,&SolverTM<double,double>::TraceTMm);
	solverb.solve_bkt(root2,rootmax,tosolve2);

	cout << energy_to_eV(root2) << " (rootb)= " << funcpointer(root2)<< endl;

	double root3=energy_from_eV(1.798);
	o2scl::gsl_root_brent<> solvert;
	o2scl::funct_mfptr<SolverTM<double,double>> tosolve3(&solver,&SolverTM<double,double>::TraceTMp);
	solvert.solve_bkt(root3,rootmax,tosolve3);

	cout << energy_to_eV(root3) << " (roott)= " << funcpointer(root3)<< endl;

	cout << energy_to_eV(root2)-energy_to_eV(root3) << " (miniband)= " << funcpointer(root3)<< endl;*/



	//cout << funcpointer(GaAs(2).gap_gamma()+energy_from_eV(0.2));

	return 0;
}
