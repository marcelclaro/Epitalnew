
/*Basic libraries*/
#include <thread>
#include <iostream>
#include <limits>
#include <string>
#include <cstring>
#include <cstdio>
#include <vector>
#include <complex>
#include <ratio>

/*Epital libraries*/
#include "../CustomTypes.hpp"
#include "../CMaterial.hpp"
#include "../Materials.hpp"
#include "../Units.hpp"
#include "../Alloy.hpp"
#include "../Graphics.hpp"
#include "../Heterostructure.hpp"
#include "../Solution.hpp"
#include "../SolverSO.hpp"


using namespace std;
using namespace epital; //All the functions are in the namespace "epital"

int main(){

    /*Create a shared pointer of a Material to be used as a substrate*/
    /*This is the standard lattice constant of the whole structure, different lattice constants will be strained*/
    /* In this case it is GaAs at 77 K*/
    shared_ptr<Material> subs = make_shared<GaAs>(77);

    /*Creat a list of layers in the active region*/
    vector<Heterostructure<double>::Epilayer> layers;
    /*Creat a list of layers also in the contacts*/
    /*It is not used in most of the functions*/
    vector<Heterostructure<double>::Epilayer> contact;

    /*Add layers to the heterostructure*/
    /*In this case it will be a quantum well*/
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(77,0.4),50.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaAs>(77),20.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(77,0.4),50.0_Angs));

    /*Creates the heterostructure using active region, contacts (empty in this case) and substrate at some temperature (77K)*/
    Heterostructure<double> sample(layers,contact,subs,0);

    /*It finds the solution using plane waves base (transfer matrix method) in the conduction band*/
    Solution<double,double> solution = SolverSO<double,double>::SolvePeriodic(sample,Band::Conduction,1,100.0_fs);


    /*Print the number of levels found*/
    cout << "Levels=" << solution.getLevels() << " \n";

    /*Print the energy*/
    for(int i=0;i<solution.getLevels();i++)
      cout << "Level " << i << " - Energy:" << energy_to_eV(solution.getEnergy(i).real()) << "\n";

    /*Plot functions using gnuplot:*/

    /*Create the potential function (discrete with 1000 points)*/
    Grid1D<double> basegrid(sample.Begin(),sample.End(),1000);  //Create a Grid
    DiscreteFunction<double,double> potential = sample.ConductionBand_gamma(basegrid);
    /*simple discrete function plot creation*/
    Graphics electronpot(potential);

    /*create plot from plane waves*/
    Graphics firstlevel(solution.getWavefunction(0)->probabilitydensity().Real()+solution.getEnergy(0).real());


    /*Join plots*/
    Graphics allgraphs({&electronpot,&firstlevel});

    /*Configure axis and plot*/
    allgraphs.setXrange(sample.Begin(),sample.End());
    allgraphs.DataConvertionY(energy_to_eV); /*convert from default atomic units to eV*/
    allgraphs.Plot();


}
