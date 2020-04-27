#pragma once

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
#include "../src/CustomTypes.hpp"
#include "../src/CMaterial.hpp"
#include "../src/Materials.hpp"
#include "../src/Units.hpp"
#include "../src/Alloy.hpp"
#include "../src/Graphics.hpp"
#include "../src/Heterostructure.hpp"
#include "../src/Solution.hpp"
#include "../src/SolverSO.hpp"




using namespace std;
using namespace epital; //All the functions are in the namespace "epital"

double test_splitoperator(double wellwidth_angs){

    /*Create a shared pointer of a Material to be used as a substrate*/
    /*This is the standard lattice constant of the whole structure, different lattice constants will be strained*/
    /* In this case it is GaAs at 77 K*/
    shared_ptr<Material> subs = make_shared<GaAs>(77);

    /*Creat a list of layers in the active region*/
    vector<Heterostructure<double>::Epilayer> layers;
    /*Creat a list of layers also in the contacts*/
    vector<Heterostructure<double>::Epilayer> contact;

    /*Add layers to the heterostructure*/
    /*In this case it will be a quantum well*/
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(77,0.4),50.0_Angs));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<GaAs>(77),length_from_SI(wellwidth_angs*1E-10)));
    layers.push_back(Heterostructure<double>::Epilayer(make_shared<AlGaAs>(77,0.4),50.0_Angs));

    /*Creates the heterostructure using active region, contacts (empty in this case) and substrate at some temperature (77K)*/
    Heterostructure<double> sample(layers,contact,subs,0);

    /*It finds the solution using plane waves base (transfer matrix method) in the conduction band*/
    /*It looks for 4 levels from the band edge (1.42eV) to the top of the barrier, setting minimum distance of 0.001eV between the levels (degenerency prevention).*/
    Solution<double,double> solution = SolverSO<double,double>::SolvePeriodic(sample,Band::Conduction,1,100.0_fs);

    /*Print the energy*/
    return energy_to_eV(solution.getEnergy(0).real());

}
