/*
 * lab.cpp
 *
 *  Created on: Jul 24, 2013
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
#include "SolutionWannier.hpp"
#include "SolverSO.hpp"
#include "SolverTM.hpp"
//#include "PhononScattering.hpp"
//#include "PhotonScattering.hpp"
#include "CarrierStatistics.hpp"
//#include "QCD.hpp"
#include "PerturbationMethod.hpp"
#include "BrentsRootFinder.hpp"
/*
#include <o2scl/funct.h>
#include <o2scl/gsl_root_brent.h>
#include <o2scl/cern_mroot_root.h>
*/
//#include <Eigen/Dense>

using namespace std;
using namespace epital;

#include <unistd.h>
inline void mysleep(unsigned millis) {
	::usleep(millis * 1000);
}

int main(){

	cout.precision(15);

	shared_ptr<InGaAs> subs = make_shared<InGaAs>(77,0.53);

	vector<Heterostructure<double>::Epilayer> layers;
	vector<Heterostructure<double>::Epilayer> contact;

	//layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(300,0.53),100.0_Angs));
	/*layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgSe>(300),60.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnSe>(300),5.9_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<CdSe>(300),11.7_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgSe>(300),5.9_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<CdSe>(300),11.7_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnSe>(300),5.9_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgSe>(300),60.0_Angs));*/
	//layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(300,0.53),100.0_Angs));

/*
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),100.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.51),21.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),100.0_Angs));

*//*

	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(300,0.44,0.15),100.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(300,0.51),41.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(300,0.44,0.15),100.0_Angs));
*/

	/*

	layers.push_back(Heterostructure<double>::Epilayer(make_shared<>(300,0.51),41.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<>(300,0.44,0.15),100.0_Angs));


	*/
/*
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<InAlGaAs>(77,0.53,0.12),150.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<InGaAs>(77,0.60),95.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<InAlGaAs>(77,0.53,0.12),150.0_Angs));
*/


	double thicknesserror = 1.00;

	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.307,0.36),thicknesserror*100.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*23.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.307,0.36),thicknesserror*100.0_Angs));

/*


	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*100.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*14.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*38.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*28.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*11.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*28.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*11.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*28.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*13.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*30.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*16.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*26.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*19.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*28.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*23.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*24.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*29.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*14.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*100.0_Angs));
	*/


	Heterostructure<double> sample(layers,contact,subs,0);


	//Solution<double,double> solution = SolverSO<double,double>::SolvePeriodic(&sample,Band::Conduction,9,100.0_fs);

	//QCD<double,double> qcd(layers,contact,subs);

	Grid1D<double> basegrid(sample.Begin(),sample.End(),2000);  //Create a Grid
	DiscreteFunction<double,double> potential = sample.ConductionBand_gamma(basegrid);
	DiscreteFunction<double,double> heavyhole = sample.ValenceBand_HH(basegrid);/*
	DiscreteFunction<double,double> lighthole = sample.ValenceBand_lH(basegrid);
	DiscreteFunction<double,double> mass = sample.effectiveMass_gamma(basegrid);*/

	/*potential.savetoFile("Solution/potential.func");

	DiscreteFunction<double,double> potential2("Solution/potential.func");

	Graphics test(potential2);
	test.Plot();*/

	//Solution<double,double> solutionfst = SolverSO<double,double>::SolvePeriodicFFT(sample,Band::Conduction,3,100.0_fs,0.01_fs,2000);
	/*solutionfst.saveToFile("solution3620.sol");*/
	//Graphics pt(solutionfst.getWavefunction(0)->probabilitydensity().Real()+solutionfst.getEnergy(0).real());
	//Graphics pt2(solutionfst.getWavefunction(1)->probabilitydensity().Real()+solutionfst.getEnergy(1).real());
	//Graphics pt3(solutionfst.getWavefunction(2)->probabilitydensity().Real()+solutionfst.getEnergy(2).real());
	/*Graphics pt4(solutionfst.getWavefunction(3)->probabilitydensity().Real()+solutionfst.getEnergy(3).real());
	Graphics pt5(solutionfst.getWavefunction(4)->probabilitydensity().Real()+solutionfst.getEnergy(4).real());
	Graphics pt6(solutionfst.getWavefunction(5)->probabilitydensity().Real()+solutionfst.getEnergy(5).real());
	Graphics pt7(solutionfst.getWavefunction(6)->probabilitydensity().Real()+solutionfst.getEnergy(6).real());
	Graphics pt8(solutionfst.getWavefunction(7)->probabilitydensity().Real()+solutionfst.getEnergy(7).real());
	Graphics pt9(solutionfst.getWavefunction(8)->probabilitydensity().Real()+solutionfst.getEnergy(8).real());*/
	//Graphics gama(potential);
	//Graphics all({&gama,&pt,&pt2,&pt3/*,&pt4,&pt5,&pt6,&pt7,&pt8,&pt9/*,&pt10/*,&pt11,&pt12,&pt13,&pt14,&pt15,&pt16,&pt17,&pt18,&pt19,&pt20*/});
/*
	for(int i=0; i< 1;i++)
		cout << "Energy-Level " << i << " :" << energy_to_eV(solutionfst.getEnergy(i).real()) << endl;

	all.setXrange(sample.Begin(),sample.End());
	all.DataConvertionY(energy_to_eV);
	all.Plot();
*/
	//Solution<double,double> solution("solution3620.sol");


    SolutionPW<double,double> solution = SolverTM<double,double>::SolvePeriodic(sample,Band::Conduction, 3, 0.001_eV,0.96_eV,1.84_eV);
    //SolutionPW<double,double> solutionholes = SolverTM<double,double>::SolvePeriodic(sample,Band::HeavyHole, 9, 0.001_eV,-0.97_eV,-1.14_eV);
    cout << "Niveis=" << solution.getLevels() << " \n";
	Graphics ptb(*solution.getWavefunction(0),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(0).real());
	Graphics pt2b(*solution.getWavefunction(1),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(1).real());
	Graphics pt3b(*solution.getWavefunction(2),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(2).real());
	/*Graphics pt4b(*solution.getWavefunction(3),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(3).real());
	Graphics pt5b(*solution.getWavefunction(4),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(4).real());
	Graphics pt6b(*solution.getWavefunction(5),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(5).real());
	Graphics pt7b(*solution.getWavefunction(6),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(6).real());
	Graphics pt8b(*solution.getWavefunction(7),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(7).real());
	Graphics pt9b(*solution.getWavefunction(8),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(8).real());
/*	Graphics ptbh(*solutionholes.getWavefunction(0),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(0).real());
	Graphics pt2bh(*solutionholes.getWavefunction(1),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(1).real());
	Graphics pt3bh(*solutionholes.getWavefunction(2),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(2).real());
	Graphics pt4bh(*solutionholes.getWavefunction(3),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(3).real());
	Graphics pt5bh(*solutionholes.getWavefunction(4),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(4).real());
	Graphics pt6bh(*solutionholes.getWavefunction(5),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(5).real());
	Graphics pt7bh(*solutionholes.getWavefunction(6),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(6).real());
	Graphics pt8bh(*solutionholes.getWavefunction(7),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(7).real());
	Graphics pt9bh(*solutionholes.getWavefunction(8),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(8).real());*/

	/*Graphics pt10b(*solution.getWavefunction(9),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(9).real());
	Graphics pt11b(*solution.getWavefunction(10),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(10).real());*/
	Graphics gamab(potential);
	Graphics hh(heavyhole);
	//potential.savetoFileASCII("hibridelectron.txt");
	//heavyhole.savetoFileASCII("hybridholes.txt");
	//Graphics lh(lighthole);*/
	Graphics allb({&gamab,&hh/*,&lh*/,&ptb,&pt2b,&pt3b/*,&pt4b,&pt5b,&pt6b,&pt7b,&pt8b,&pt9b/*,&ptbh,&pt2bh,&pt3bh,&pt4bh,&pt5bh,&pt6bh,&pt7bh,&pt8bh,&pt9bh/*,&pt10b,&pt11b*/});
	allb.setXrange(sample.Begin(),sample.End());
	allb.DataConvertionY(energy_to_eV);
	allb.Plot();

	cout << energy_to_eV(solution.getEnergy(0).real()) << "\n" ;
	cout << energy_to_eV(solution.getEnergy(1).real()) << "\n" ;
	cout << energy_to_eV(solution.getEnergy(2).real()) << "\n" ;

	/*
	cout << "PL: \n";
	cout << energy_to_eV(solution.getEnergy(0).real()-solutionholes.getEnergy(0).real()) << "\n" ;
	cout << energy_to_eV(solution.getEnergy(1).real()-solutionholes.getEnergy(1).real()) << "\n" ;
	cout << energy_to_eV(solution.getEnergy(2).real()-solutionholes.getEnergy(2).real()) << "\n" ;
	cout << energy_to_eV(solution.getEnergy(3).real()-solutionholes.getEnergy(3).real()) << "\n" ;
	cout << energy_to_eV(solution.getEnergy(4).real()-solutionholes.getEnergy(4).real()) << "\n" ;
	cout << energy_to_eV(solution.getEnergy(5).real()-solutionholes.getEnergy(5).real()) << "\n" ;
	cout << energy_to_eV(solution.getEnergy(6).real()-solutionholes.getEnergy(6).real()) << "\n" ;
	cout << energy_to_eV(solution.getEnergy(7).real()-solutionholes.getEnergy(7).real()) << "\n" ;*/



	/*cout << "value=" << (*solution.getWavefunction(0))(0.0) << "/" << (*solution.getWavefunction(0))(sample.getTotalWidth()) << "-%="  << 100.0*std::real((*solution.getWavefunction(0))(0.0)-(*solution.getWavefunction(0))(sample.getTotalWidth()))/(*solution.getWavefunction(0))(0.0) << endl;
	cout << "value=" <<(*solution.getWavefunction(1))(0.0) << "/" << (*solution.getWavefunction(1))(sample.getTotalWidth()) << "-%="  << 100.0*std::real((*solution.getWavefunction(1))(0.0)-(*solution.getWavefunction(1))(sample.getTotalWidth()))/(*solution.getWavefunction(1))(0.0) << endl;
	cout << "value=" <<(*solution.getWavefunction(2))(0.0) << "/" << (*solution.getWavefunction(2))(sample.getTotalWidth()) << "-%="  << 100.0*std::real((*solution.getWavefunction(2))(0.0)-(*solution.getWavefunction(2))(sample.getTotalWidth()))/(*solution.getWavefunction(2))(0.0) << endl;
	cout << "value=" <<(*solution.getWavefunction(3))(0.0) << "/" << (*solution.getWavefunction(3))(sample.getTotalWidth()) << "-%="  << 100.0*std::real((*solution.getWavefunction(3))(0.0)-(*solution.getWavefunction(3))(sample.getTotalWidth()))/(*solution.getWavefunction(3))(0.0) << endl;
	cout << "value=" <<(*solution.getWavefunction(4))(0.0) << "/" << (*solution.getWavefunction(4))(sample.getTotalWidth()) << "-%="  << 100.0*std::real((*solution.getWavefunction(4))(0.0)-(*solution.getWavefunction(4))(sample.getTotalWidth()))/(*solution.getWavefunction(4))(0.0) << endl;
	cout << "value=" <<(*solution.getWavefunction(5))(0.0) << "/" << (*solution.getWavefunction(5))(sample.getTotalWidth()) << "-%="  << 100.0*std::real((*solution.getWavefunction(5))(0.0)-(*solution.getWavefunction(5))(sample.getTotalWidth()))/(*solution.getWavefunction(5))(0.0) << endl;
	cout << "value=" <<(*solution.getWavefunction(6))(0.0) << "/" << (*solution.getWavefunction(6))(sample.getTotalWidth()) <<  "-%=" << 100.0*std::real((*solution.getWavefunction(6))(0.0)-(*solution.getWavefunction(6))(sample.getTotalWidth()))/(*solution.getWavefunction(6))(0.0) << endl;
	cout << "value=" <<(*solution.getWavefunction(7))(0.0) << "/" << (*solution.getWavefunction(7))(sample.getTotalWidth()) <<  "-%=" << 100.0*std::real((*solution.getWavefunction(7))(0.0)-(*solution.getWavefunction(7))(sample.getTotalWidth()))/(*solution.getWavefunction(7))(0.0) << endl;
*/
	//double fermilevel = CarrierStatistics<double>::FermiEnergyWells(solution,sample.getTotalDoping(),make_shared<InGaAs>(2,0.14),77,Band::Conduction);
	//cout << "Fermi:" <<energy_to_eV(fermilevel) << endl;

	//vector<double> widths = SolverTM<double,double>::getminibandwidth(sample,Band::Conduction,solution);

/*
	for(long int p=1; p < 4480 ; p+=10)
		cout << real(conj(solutionb.getWavefunction(7)->getValue(0,p))*solutionb.getWavefunction(7)->getValue(0,p)) << endl;
*/
/*
	SolutionWannier<double,double> solutionb = SolverTM<double,double>::getWannierFunctions(sample,Band::Conduction,solution,100);

	DiscreteFunction<std::complex<double>,double> wannier = solutionb.getWavefunction(0)->getDiscrete(0,3,3,3000);
	Graphics pt(wannier.probabilitydensity().Real()+solution.getEnergy(0).real());
	DiscreteFunction<std::complex<double>,double> wannier2 = solutionb.getWavefunction(1)->getDiscrete(0,3,3,3000);
	Graphics pt2(wannier2.probabilitydensity().Real()+solution.getEnergy(1).real());
	DiscreteFunction<std::complex<double>,double> wannier3 = solutionb.getWavefunction(2)->getDiscrete(0,3,3,3000);
	Graphics pt3(wannier3.probabilitydensity().Real()+solution.getEnergy(2).real());
	DiscreteFunction<std::complex<double>,double> wannier4 = solutionb.getWavefunction(3)->getDiscrete(0,3,3,3000);
	Graphics pt4(wannier4.probabilitydensity().Real()+solution.getEnergy(3).real());
	DiscreteFunction<std::complex<double>,double> wannier5 = solutionb.getWavefunction(4)->getDiscrete(0,3,3,3000);
	Graphics pt5(wannier5.probabilitydensity().Real()+solution.getEnergy(4).real());
	DiscreteFunction<std::complex<double>,double> wannier6 = solutionb.getWavefunction(5)->getDiscrete(0,3,3,3000);
	Graphics pt6(wannier6.probabilitydensity().Real()+solution.getEnergy(5).real());
	DiscreteFunction<std::complex<double>,double> wannier7 = solutionb.getWavefunction(6)->getDiscrete(0,3,3,3000);
	Graphics pt7(wannier7.probabilitydensity().Real()+solution.getEnergy(6).real());
	DiscreteFunction<std::complex<double>,double> wannier8 = solutionb.getWavefunction(7)->getDiscrete(0,3,3,3000);
	Graphics pt8(wannier8.probabilitydensity().Real()+solution.getEnergy(7).real());
	Graphics gama(potential);

	Graphics all({&gama,&pt,&pt2,&pt3,&pt4,&pt5,&pt6,&pt7,&pt8/*,&pt9*//*});

	//all.setXrange(sample.Begin(),sample.End());
	all.DataConvertionY(energy_to_eV);
	all.Plot();
/*
	cout << "2d density=" << uni_density_to_SI(sample.getTotalDoping())*5.2917721092e-11*10e-4 << endl;
	cout << "Transition rate=" << frequency_to_SI(sample.getTotalDoping()*PhotonScattering<double>::MeanTransRate(77,shared_ptr<DiscreteFunction<std::complex<double>,double>>(new DiscreteFunction<std::complex<double>,double>(wannier)),shared_ptr<DiscreteFunction<std::complex<double>,double>>(new DiscreteFunction<std::complex<double>,double>(wannier8)),solution.getEnergy(7).real(),solution.getEnergy(0).real(),0.0,sample.getTotalDoping(),powerperarea_from_SI(cos(Constant::pi/4)*3.35e-8/(400e-6*400e-6)),sample.getSubstrate(),Band::Conduction)*3.57106e+20) << endl;
	cout << "Signal(g->1)=" << 40*10e4*1.609e-19*400e-6*400e-6*frequency_to_SI(sample.getTotalDoping()*PhotonScattering<double>::MeanTransRate(77,shared_ptr<DiscreteFunction<std::complex<double>,double>>(new DiscreteFunction<std::complex<double>,double>(wannier)),shared_ptr<DiscreteFunction<std::complex<double>,double>>(new DiscreteFunction<std::complex<double>,double>(wannier8)),solution.getEnergy(7).real(),solution.getEnergy(0).real(),0.0,sample.getTotalDoping(),powerperarea_from_SI(cos(Constant::pi/4)*3.35e-8/(400e-6*400e-6)),sample.getSubstrate(),Band::Conduction)*3.57106e+20) << endl;

	//cout << funcpointer(GaAs(2).gap_gamma()+energy_from_eV(0.2));*/

	return 0;
}



