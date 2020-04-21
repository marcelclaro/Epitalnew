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
#include "PhononScattering.hpp"
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

	shared_ptr<Material> subs = make_shared<InGaAs>(77,0.53);

	vector<Heterostructure<double>::Epilayer> layers;
	vector<Heterostructure<double>::Epilayer> contact;
	double thicknesserror = 1.00;
/*
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*14.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*38.0_Angs, 3e+18_under_cubic_cm));
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
*/


	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*10.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*42.0_Angs,3.0e+18_under_cubic_cm));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*33.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*13.5_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*22.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*14.4_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*22.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*16.3_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*22.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*18.4_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*22.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*21.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*22.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*24.2_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*20.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*28.1_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*20.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<ZnCdSe>(77,0.48),thicknesserror*33.5_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<MgZnCdSe>(77,0.29,0.45),thicknesserror*10.0_Angs));




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


    SolutionPW<double,double> solution = SolverTM<double,double>::SolvePeriodic(sample,Band::Conduction, 9, 0.001_eV,1.21_eV,1.84_eV);
    SolutionPW<double,double> solutionholes = SolverTM<double,double>::SolvePeriodic(sample,Band::HeavyHole, 9, 0.001_eV,-0.97_eV,-1.14_eV);
    cout << "Niveis=" << solution.getLevels() << " \n";
	Graphics ptb(*solution.getWavefunction(0),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(0).real());
	Graphics pt2b(*solution.getWavefunction(1),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(1).real());
	Graphics pt3b(*solution.getWavefunction(2),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(2).real());
	Graphics pt4b(*solution.getWavefunction(3),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(3).real());
	Graphics pt5b(*solution.getWavefunction(4),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(4).real());
	Graphics pt6b(*solution.getWavefunction(5),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(5).real());
	Graphics pt7b(*solution.getWavefunction(6),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(6).real());
	Graphics pt8b(*solution.getWavefunction(7),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(7).real());
	Graphics pt9b(*solution.getWavefunction(8),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(8).real());
	Graphics ptbh(*solutionholes.getWavefunction(0),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(0).real());
	Graphics pt2bh(*solutionholes.getWavefunction(1),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(1).real());
	Graphics pt3bh(*solutionholes.getWavefunction(2),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(2).real());
	Graphics pt4bh(*solutionholes.getWavefunction(3),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(3).real());
	Graphics pt5bh(*solutionholes.getWavefunction(4),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(4).real());
	Graphics pt6bh(*solutionholes.getWavefunction(5),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(5).real());
	Graphics pt7bh(*solutionholes.getWavefunction(6),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(6).real());
	Graphics pt8bh(*solutionholes.getWavefunction(7),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(7).real());
	Graphics pt9bh(*solutionholes.getWavefunction(8),1000,Graphics::Complexplot::SQUARED,solutionholes.getEnergy(8).real());

	/*Graphics pt10b(*solution.getWavefunction(9),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(9).real());
	Graphics pt11b(*solution.getWavefunction(10),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(10).real());*/
	Graphics gamab(potential);
	Graphics hh(heavyhole);
	//Graphics lh(lighthole);*/
	Graphics allb({&gamab/*,&hh/*,&lh*/,&ptb,&pt2b,&pt3b,&pt4b,&pt5b,&pt6b,&pt7b,&pt8b,&pt9b/*,&ptbh,&pt2bh,&pt3bh,&pt4bh,&pt5bh,&pt6bh,&pt7bh,&pt8bh,&pt9bh/*,&pt10b,&pt11b*/});
	allb.setXrange(sample.Begin(),sample.End());
	allb.DataConvertionY(energy_to_eV);
	allb.Plot();


	cout << "Levels e: \n";
	cout << energy_to_eV(solution.getEnergy(0).real()) << "\n" ;
	cout << energy_to_eV(solution.getEnergy(1).real()) << "\n" ;
	cout << energy_to_eV(solution.getEnergy(2).real()) << "\n" ;
	cout << energy_to_eV(solution.getEnergy(3).real()) << "\n" ;
	cout << energy_to_eV(solution.getEnergy(4).real()) << "\n" ;
	cout << energy_to_eV(solution.getEnergy(5).real()) << "\n" ;
	cout << energy_to_eV(solution.getEnergy(6).real()) << "\n" ;
	cout << energy_to_eV(solution.getEnergy(7).real()) << "\n" ;


	cout << "levels HH: \n";
	cout << energy_to_eV(solutionholes.getEnergy(0).real()) << "\n" ;
	cout << energy_to_eV(solutionholes.getEnergy(1).real()) << "\n" ;
	cout << energy_to_eV(solutionholes.getEnergy(2).real()) << "\n" ;
	cout << energy_to_eV(solutionholes.getEnergy(3).real()) << "\n" ;
	cout << energy_to_eV(solutionholes.getEnergy(4).real()) << "\n" ;
	cout << energy_to_eV(solutionholes.getEnergy(5).real()) << "\n" ;
	cout << energy_to_eV(solutionholes.getEnergy(6).real()) << "\n" ;
	cout << energy_to_eV(solutionholes.getEnergy(7).real()) << "\n" ;

	cout << "PL: \n";
	cout << energy_to_eV(solution.getEnergy(0).real()-solutionholes.getEnergy(0).real())-0.0208 << "\n" ;
	cout << energy_to_eV(solution.getEnergy(1).real()-solutionholes.getEnergy(1).real())-0.0208 << "\n" ;
	cout << energy_to_eV(solution.getEnergy(2).real()-solutionholes.getEnergy(2).real())-0.0208 << "\n" ;
	cout << energy_to_eV(solution.getEnergy(3).real()-solutionholes.getEnergy(3).real())-0.0208 << "\n" ;
	cout << energy_to_eV(solution.getEnergy(4).real()-solutionholes.getEnergy(4).real())-0.0208 << "\n" ;
	cout << energy_to_eV(solution.getEnergy(5).real()-solutionholes.getEnergy(5).real())-0.0208 << "\n" ;
	cout << energy_to_eV(solution.getEnergy(6).real()-solutionholes.getEnergy(6).real())-0.0208 << "\n" ;
	cout << energy_to_eV(solution.getEnergy(7).real()-solutionholes.getEnergy(7).real())-0.0208 << "\n" ;

	double fermilevel = CarrierStatistics<double>::FermiEnergyWells(solution,sample.getTotalDoping(),make_shared<ZnCdSe>(77,0.48),77,Band::Conduction);
	cout << "Fermi:" <<energy_to_eV(fermilevel) << endl;

	cout << "Transition rate: \n" << std::scientific;
	shared_ptr<DiscreteFunction<std::complex<double>,double>> level0 = solution.getWavefunction(0)->getDiscrete(basegrid);
	shared_ptr<DiscreteFunction<std::complex<double>,double>> level1 = solution.getWavefunction(1)->getDiscrete(basegrid);
	shared_ptr<DiscreteFunction<std::complex<double>,double>> level2 = solution.getWavefunction(0)->getDiscrete(basegrid);
	shared_ptr<DiscreteFunction<std::complex<double>,double>> level3 = solution.getWavefunction(1)->getDiscrete(basegrid);
	shared_ptr<DiscreteFunction<std::complex<double>,double>> level4 = solution.getWavefunction(0)->getDiscrete(basegrid);
	shared_ptr<DiscreteFunction<std::complex<double>,double>> level5 = solution.getWavefunction(1)->getDiscrete(basegrid);
	shared_ptr<DiscreteFunction<std::complex<double>,double>> level6 = solution.getWavefunction(0)->getDiscrete(basegrid);
	cout << "(1->0) "<< frequency_to_SI(PhononScattering<double>::MeanEmissionRateLO(77,level0,level1,solution.getEnergy(0).real(),solution.getEnergy(1).real(),fermilevel,fermilevel,subs,Band::Conduction)) << endl;
	cout << "(2->1) "<< frequency_to_SI(PhononScattering<double>::MeanEmissionRateLO(77,level1,level2,solution.getEnergy(1).real(),solution.getEnergy(2).real(),fermilevel,fermilevel,subs,Band::Conduction)) << endl;
	cout << "(3->2) "<< frequency_to_SI(PhononScattering<double>::MeanEmissionRateLO(77,level2,level3,solution.getEnergy(2).real(),solution.getEnergy(3).real(),fermilevel,fermilevel,subs,Band::Conduction)) << endl;
	cout << "(4->3) "<< frequency_to_SI(PhononScattering<double>::MeanEmissionRateLO(77,level3,level4,solution.getEnergy(3).real(),solution.getEnergy(4).real(),fermilevel,fermilevel,subs,Band::Conduction)) << endl;
	cout << "(5->4) "<< frequency_to_SI(PhononScattering<double>::MeanEmissionRateLO(77,level4,level5,solution.getEnergy(4).real(),solution.getEnergy(5).real(),fermilevel,fermilevel,subs,Band::Conduction)) << endl;
	cout << "(6->5) "<< frequency_to_SI(PhononScattering<double>::MeanEmissionRateLO(77,level5,level6,solution.getEnergy(5).real(),solution.getEnergy(6).real(),fermilevel,fermilevel,subs,Band::Conduction)) << endl;


/*
	cout << "2d density=" << uni_density_to_SI(sample.getTotalDoping())*5.2917721092e-11*10e-4 << endl;
	cout << "Transition rate=" << frequency_to_SI(sample.getTotalDoping()*PhotonScattering<double>::MeanTransRate(77,shared_ptr<DiscreteFunction<std::complex<double>,double>>(new DiscreteFunction<std::complex<double>,double>(wannier)),shared_ptr<DiscreteFunction<std::complex<double>,double>>(new DiscreteFunction<std::complex<double>,double>(wannier8)),solution.getEnergy(7).real(),solution.getEnergy(0).real(),0.0,sample.getTotalDoping(),powerperarea_from_SI(cos(Constant::pi/4)*3.35e-8/(400e-6*400e-6)),sample.getSubstrate(),Band::Conduction)*3.57106e+20) << endl;
	cout << "Signal(g->1)=" << 40*10e4*1.609e-19*400e-6*400e-6*frequency_to_SI(sample.getTotalDoping()*PhotonScattering<double>::MeanTransRate(77,shared_ptr<DiscreteFunction<std::complex<double>,double>>(new DiscreteFunction<std::complex<double>,double>(wannier)),shared_ptr<DiscreteFunction<std::complex<double>,double>>(new DiscreteFunction<std::complex<double>,double>(wannier8)),solution.getEnergy(7).real(),solution.getEnergy(0).real(),0.0,sample.getTotalDoping(),powerperarea_from_SI(cos(Constant::pi/4)*3.35e-8/(400e-6*400e-6)),sample.getSubstrate(),Band::Conduction)*3.57106e+20) << endl;

	//cout << funcpointer(GaAs(2).gap_gamma()+energy_from_eV(0.2));*/

	return 0;
}



