/*
 * Lab26.cpp
 *
 *  Created on: Nov 20, 2017
 *      Author: marcel
 */

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
//#include "CarrierStatistics.hpp"
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

	shared_ptr<InGaAs> subs = make_shared<InGaAs>(77,0.52);

	vector<Heterostructure<double>::Epilayer> layers;
	vector<Heterostructure<double>::Epilayer> contact;

	vector<double> widths;
	vector<double> widthslh;
	vector<double> widthshh;
/*
	layers.clear();
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<InAlAs>(77,0.48),150.0_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<InGaAs>(77,0.0.47),15.7_Angs));
	layers.push_back(Heterostructure<double>::Epilayer(make_shared<InAlAs>(77,0.48),150.0_Angs));

	Heterostructure<double> *sample = new Heterostructure<double>(layers,contact,subs,0);

	Grid1D<double> basegrid(sample->Begin(),sample->End(),2000);  //Create a Grid
	DiscreteFunction<double,double> potential = sample->ConductionBand_gamma(basegrid);
	DiscreteFunction<double,double> potentiallight = sample->ValenceBand_lH(basegrid);
	DiscreteFunction<double,double> potentialheavy = sample->ValenceBand_HH(basegrid);


    SolutionPW<double,double> solution = SolverTM<double,double>::SolvePeriodic(*sample,Band::Conduction, 2, 0.001_eV,0.9_eV,3.0_eV);
    //SolutionPW<double,double> solutionlight = SolverTM<double,double>::SolvePeriodic(*sample,Band::LightHole, 1, 0.001_eV,0.0_eV,-2.0_eV);
    //SolutionPW<double,double> solutionheavy = SolverTM<double,double>::SolvePeriodic(*sample,Band::HeavyHole, 1, 0.001_eV,0.0_eV,-2.0_eV);


	widths.clear();
	//widths = SolverTM<double,double>::getminibandwidth(*sample,Band::Conduction,solution,1);

	Graphics ptb(*solution.getWavefunction(0),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(0).real());
	//Graphics pt2b(*solutionheavy.getWavefunction(0),1000,Graphics::Complexplot::SQUARED,solutionheavy.getEnergy(0).real());
	//Graphics pt3b(*solutionlight.getWavefunction(0),1000,Graphics::Complexplot::SQUARED,solutionlight.getEnergy(0).real());
	Graphics pt4b(*solution.getWavefunction(1),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(1).real());
	Graphics gamab(potential);
	Graphics lh(potentiallight);
	Graphics hh(potentialheavy);
	//Graphics allb({&gamab,/*&hh,&lh,&ptb/*,&pt2b,&pt3b,&pt4b/*,&pt5b,&pt6b,&pt7b,&pt8b,&pt9b,&pt10b,&pt11b*//*});
	allb.setXrange(sample->Begin(),sample->End());
	allb.DataConvertionY(energy_to_eV);
	allb.Plot();

*/
	for(double wellwidth = 35.0_Angs; wellwidth <= 100.0_Angs ; wellwidth += 5.0_Angs){
		layers.clear();
		layers.push_back(Heterostructure<double>::Epilayer(make_shared<InAlAs>(77,0.52),100.0_Angs));
		//layers.push_back(Heterostructure<double>::Epilayer(make_shared<InAlAs>(77,0.47),45.0_Angs));
		layers.push_back(Heterostructure<double>::Epilayer(make_shared<InGaAs>(77,0.53),wellwidth));
		//layers.push_back(Heterostructure<double>::Epilayer(make_shared<InAlAs>(77,0.47),45.0_Angs));
		layers.push_back(Heterostructure<double>::Epilayer(make_shared<InAlAs>(77,0.52),100.0_Angs));

		Heterostructure<double> *sample = new Heterostructure<double>(layers,contact,subs,0);

		Grid1D<double> basegrid(sample->Begin(),sample->End(),2000);  //Create a Grid
		DiscreteFunction<double,double> potential = sample->ConductionBand_gamma(basegrid);
		//DiscreteFunction<double,double> potentiallight = sample->ValenceBand_lH(basegrid);
		//DiscreteFunction<double,double> potentialheavy = sample->ValenceBand_HH(basegrid);


    	SolutionPW<double,double> solution = SolverTM<double,double>::SolvePeriodic(*sample,Band::Conduction, 2, 0.001_eV,0.8_eV,3.0_eV);
    	//SolutionPW<double,double> solutionlight = SolverTM<double,double>::SolvePeriodic(*sample,Band::LightHole, 1, 0.001_eV,0.0_eV,-2.0_eV);
    	//SolutionPW<double,double> solutionheavy = SolverTM<double,double>::SolvePeriodic(*sample,Band::HeavyHole, 1, 0.001_eV,0.0_eV,-2.0_eV);
    	//Solution<double,double> solutionfst = SolverSO<double,double>::SolvePeriodicFFT(*sample,Band::Conduction,2,200.0_fs,0.005_fs,2000);

    	Graphics pt(*solution.getWavefunction(0),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(0).real());
    	Graphics pt2(*solution.getWavefunction(1),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(1).real());
    	//Graphics pt3(solutionfst.getWavefunction(2)->probabilitydensity().Real()+solutionfst.getEnergy(2).real());
    	/*Graphics pt4(solutionfst.getWavefunction(3)->probabilitydensity().Real()+solutionfst.getEnergy(3).real());
    	Graphics pt5(solutionfst.getWavefunction(4)->probabilitydensity().Real()+solutionfst.getEnergy(4).real());
    	Graphics pt6(solutionfst.getWavefunction(5)->probabilitydensity().Real()+solutionfst.getEnergy(5).real());
    	Graphics pt7(solutionfst.getWavefunction(6)->probabilitydensity().Real()+solutionfst.getEnergy(6).real());
    	Graphics pt8(solutionfst.getWavefunction(7)->probabilitydensity().Real()+solutionfst.getEnergy(7).real());
    	Graphics pt9(solutionfst.getWavefunction(8)->probabilitydensity().Real()+solutionfst.getEnergy(8).real());*/
    	Graphics gama(potential);
    	Graphics all({&gama,&pt,&pt2/*,&pt3/*,&pt4,&pt5,&pt6,&pt7,&pt8,&pt9/*,&pt10/*,&pt11,&pt12,&pt13,&pt14,&pt15,&pt16,&pt17,&pt18,&pt19,&pt20*/});
    	all.setXrange(sample->Begin(),sample->End());
    	all.DataConvertionY(energy_to_eV);
    	all.Plot();



		widths.clear();
		/*widthslh.clear();
		widthshh.clear();
		widths = SolverTM<double,double>::getminibandwidth(*sample,Band::Conduction,solution,1);
		widthslh = SolverTM<double,double>::getminibandwidth(*sample,Band::LightHole,solutionlight,1);
		widthshh = SolverTM<double,double>::getminibandwidth(*sample,Band::HeavyHole,solutionheavy,1);
*/
		cout << length_to_SI(wellwidth) << '\t' << energy_to_eV(solution.getEnergy(0).real()) << '\t' << energy_to_eV(solution.getEnergy(1).real()) /*<< '\t' << energy_to_eV(solutionfst.getEnergy(0).real()) << '\t' << energy_to_eV(solutionfst.getEnergy(1).real()) /*<< '\t' << energy_to_eV(widths[0]) << '\t' << energy_to_eV(widthslh[0]) << '\t' << energy_to_eV(widthshh[0])*/ << '\n';

	}


	return 0;
}






