
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
#include "../SolverTM.hpp"


using namespace std;
using namespace epital; //All the functions are in the namespace "epital"

int main(){

    /*Create a shared pointer of a Material to be used as a substrate*/
    /*This is the standard lattice constant of the whole structure, different lattice constants will be strained*/
    /* In this case it is GaAs at 77 K*/
    shared_ptr<Material> subs = make_shared<InGaSe>(15,0.5);

    /*Creat a list of layers in the active region*/
    vector<Heterostructure<double>::Epilayer> layers;
    /*Creat a list of layers also in the contacts (could be empty)*/
    /*It is not used in most of the functions*/
    vector<Heterostructure<double>::Epilayer> contact;

    /*Add layers to the heterostructure*/
    /*In this case it will be a quantum well*/
    /*layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),4*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),4*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),4*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),4*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),4*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),4*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),4*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),4*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),4*8.0_Angs));


    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),5*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),5*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),5*8.0_Angs));

    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),4*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),4*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),4*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),4*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),4*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),4*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),4*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),4*8.0_Angs));*/



    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),2.5*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),2.5*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),2.5*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),2.5*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),2.5*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),2.5*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),2.5*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),2.5*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),2.5*8.0_Angs));


    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),4*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),4*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),4*8.0_Angs));

    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),2.5*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),2.5*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),2.5*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),2.5*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),2.5*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),2.5*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaSe>(15),2.5*8.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<InSe>(15),2.5*8.0_Angs));




    /*Creates the heterostructure using active region, contacts (empty in this case) and substrate at some temperature (77K)*/
    Heterostructure<double> sample(layers,contact,subs,0);

    /*It finds the solution using plane waves base (transfer matrix method) in the conduction band*/
    /*It looks for 4 levels from the band edge (1.42eV) to the top of the barrier, setting minimum distance of 0.001eV between the levels (degenerency prevention).*/
    SolutionPW<double,double> solution = SolverTM<double,double>::SolvePeriodic(sample,Band::Conduction, 8, 0.00001_eV,0.001_eV,1.6_eV);
    SolutionPW<double,double> solutionhole = SolverTM<double,double>::SolvePeriodic(sample,Band::HeavyHole, 8, 0.00001_eV,-0.7_eV,-1.9_eV);

    /*Print the number of levels found*/
    cout << "Levels=" << solution.getLevels()+ solutionhole.getLevels() << " \n";

    /*Print the energy*/
    for(int i=0;i<solution.getLevels();i++)
      cout << "Level (electron)" << i << " - Energy:" << energy_to_eV(solution.getEnergy(i).real()) << "\n";

    for(int i=0;i<solutionhole.getLevels();i++)
      cout << "Level (hole)" << i << " - Energy:" << energy_to_eV(solutionhole.getEnergy(i).real()) << "\n";

    /*Plot functions using gnuplot:*/

    /*Create the potential function (discrete with 1000 points)*/
    Grid1D<double> basegrid(sample.Begin(),sample.End(),1000);  //Create a Grid
    DiscreteFunction<double,double> potential = sample.ConductionBand_gamma(basegrid);
    DiscreteFunction<double,double> potentialhole = sample.ValenceBand_HH(basegrid);
    /*simple discrete function plot creation*/
    Graphics electronpot(potential);
    Graphics holepot(potentialhole);
    //save as file:
    //potential.savetoFileASCII("ConductionBand_gamma.txt");
    //potentialhole.savetoFileASCII("ConductionBand_gamma.txt");

    /*create plot from plane waves*/
    Graphics e1(*solution.getWavefunction(0),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(0).real());
    Graphics e2(*solution.getWavefunction(1),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(1).real());
    Graphics e3(*solution.getWavefunction(2),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(2).real());
    Graphics e4(*solution.getWavefunction(3),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(3).real());
    Graphics e5(*solution.getWavefunction(4),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(4).real());
    Graphics e6(*solution.getWavefunction(5),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(5).real());
    Graphics e7(*solution.getWavefunction(6),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(6).real());
    Graphics e8(*solution.getWavefunction(7),1000,Graphics::Complexplot::SQUARED,solution.getEnergy(7).real());
    Graphics h1(*solutionhole.getWavefunction(0),1000,Graphics::Complexplot::SQUARED,solutionhole.getEnergy(0).real());
    Graphics h2(*solutionhole.getWavefunction(1),1000,Graphics::Complexplot::SQUARED,solutionhole.getEnergy(1).real());
    Graphics h3(*solutionhole.getWavefunction(2),1000,Graphics::Complexplot::SQUARED,solutionhole.getEnergy(2).real());
    Graphics h4(*solutionhole.getWavefunction(3),1000,Graphics::Complexplot::SQUARED,solutionhole.getEnergy(3).real());
    Graphics h5(*solutionhole.getWavefunction(4),1000,Graphics::Complexplot::SQUARED,solutionhole.getEnergy(4).real());
    Graphics h6(*solutionhole.getWavefunction(5),1000,Graphics::Complexplot::SQUARED,solutionhole.getEnergy(5).real());
    Graphics h7(*solutionhole.getWavefunction(6),1000,Graphics::Complexplot::SQUARED,solutionhole.getEnergy(6).real());
    Graphics h8(*solutionhole.getWavefunction(7),1000,Graphics::Complexplot::SQUARED,solutionhole.getEnergy(7).real());

    cout << "PL(eV/nm) = " << (energy_to_eV(solution.getEnergy(0).real())-energy_to_eV(solutionhole.getEnergy(0).real())) << " / " << 1239.8/(energy_to_eV(solution.getEnergy(0).real())-energy_to_eV(solutionhole.getEnergy(0).real())) << "\n";

    /*Join plots*/
    Graphics allgraphs({&electronpot,&holepot,&e1,&e2,&e3,&e4,&e5,&e6,&e7,&e8,&h1,&h2,&h3,&h4,&h5,&h6,&h7,&h8});

    /*Configure axis and plot*/
    allgraphs.setXrange(sample.Begin(),sample.End());
    allgraphs.DataConvertionY(energy_to_eV); /*convert from default atomic units to eV*/
    //allgraphs.DataConvertionX(length_to_SI); /*convert from default atomic units to eV*/
    allgraphs.Plot();


}
