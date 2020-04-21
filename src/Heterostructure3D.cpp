/*
 * Heterostructure3D.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: marcel
 */

#include "Heterostructure3D.hpp"
#include <iostream>

#ifndef Heterostructure3D_CPP_
#define Heterostructure3D_CPP_

namespace epital {

template <typename TYPE>
Heterostructure3D<TYPE>::Heterostructure3D(std::string geometryfile, std::shared_ptr<Material> background, Temperature temp): geometryfile(geometryfile), background(background), temp(temp){

	double dimensionx=0, dimensiony=0,dimensionz=0;
	TYPE spacingx=0,spacingy=0,spacingz=0;

	std::ifstream in;
	in.open(geometryfile);
	if ( ! in ) {
	      cout << "Error: Can't open the Geometry file .\n";
	      exit(1);
	}
	std::string readlines;
	std::getline(in,readlines);  // Get the frist line from the file, if any.
	if ( in ) {  // Continue if the line was sucessfully read.
		if(readlines=="#dimensions"){
			std::getline(in,readlines);
	        dimensionx=atof(readlines.c_str());
	        std::getline(in,readlines);
	        dimensiony=atof(readlines.c_str());
	        std::getline(in,readlines);
	        dimensionz=atof(readlines.c_str());
	      }
	      else
	    	  std::cout << "#dimension in Geometry file format error stupid!" << std::endl;

	}

	std::getline(in,readlines);  // Get the frist line from the file, if any.
	if ( in ) {  // Continue if the line was sucessfully read.
	    if(readlines=="#spacing"){
	        std::getline(in,readlines);
	        spacingx=atof(readlines.c_str());
	        std::getline(in,readlines);
	        spacingy=atof(readlines.c_str());
	        std::getline(in,readlines);
	        spacingz=atof(readlines.c_str());
	    }
	    else
	        std::cout << "#spacing in Geometry file format error stupid!" << std::endl;

	}

	basegrid= new Grid3D<TYPE>(0,length_from_SI(dimensionx*1e-9),length_from_SI(spacingx*1e-9),0,length_from_SI(dimensiony*1e-9),length_from_SI(spacingy*1e-9),0,length_from_SI(dimensionz*1e-9),length_from_SI(spacingz*1e-9));

	double vtkspacing[3];
	vtkspacing[0]=spacingx;
	vtkspacing[1]=spacingy;
	vtkspacing[2]=spacingz;

	vtk_conductiongama = vtkImageData::New();
	vtk_conductiongama->SetSpacing(vtkspacing);
	vtk_effectivemass_gamma = vtkImageData::New();
	vtk_effectivemass_gamma->SetSpacing(vtkspacing);
	vtk_lattice = vtkImageData::New();
	vtk_lattice->SetSpacing(vtkspacing);
	vtk_c11 = vtkImageData::New();
	vtk_c11->SetSpacing(vtkspacing);
	vtk_c12 = vtkImageData::New();
	vtk_c12->SetSpacing(vtkspacing);
	vtk_c44 = vtkImageData::New();
	vtk_c44->SetSpacing(vtkspacing);

	vtk_valencelh = vtkImageData::New();
	vtk_valencelh->SetSpacing(vtkspacing);
	vtk_valencehh = vtkImageData::New();
	vtk_valencehh->SetSpacing(vtkspacing);
	vtk_effectivemass_lht = vtkImageData::New();
	vtk_effectivemass_lht->SetSpacing(vtkspacing);
	vtk_effectivemass_lhl = vtkImageData::New();
	vtk_effectivemass_lhl->SetSpacing(vtkspacing);
	vtk_effectivemass_hht = vtkImageData::New();
	vtk_effectivemass_hht->SetSpacing(vtkspacing);
	vtk_effectivemass_hhl = vtkImageData::New();
	vtk_effectivemass_hhl->SetSpacing(vtkspacing);




	vtk_conductiongama->SetExtent(0, dimensionx/spacingx, 0, dimensiony/spacingy, 0, dimensionz/spacingz);
	vtk_effectivemass_gamma->SetExtent(0, dimensionx/spacingx, 0, dimensiony/spacingy, 0, dimensionz/spacingz);
	vtk_lattice->SetExtent(0, dimensionx/spacingx, 0, dimensiony/spacingy, 0, dimensionz/spacingz);
	vtk_c11->SetExtent(0, dimensionx/spacingx, 0, dimensiony/spacingy, 0, dimensionz/spacingz);
	vtk_c12->SetExtent(0, dimensionx/spacingx, 0, dimensiony/spacingy, 0, dimensionz/spacingz);
	vtk_c44->SetExtent(0, dimensionx/spacingx, 0, dimensiony/spacingy, 0, dimensionz/spacingz);

	vtk_valencelh->SetExtent(0, dimensionx/spacingx, 0, dimensiony/spacingy, 0, dimensionz/spacingz);
	vtk_valencehh->SetExtent(0, dimensionx/spacingx, 0, dimensiony/spacingy, 0, dimensionz/spacingz);
	vtk_effectivemass_lht->SetExtent(0, dimensionx/spacingx, 0, dimensiony/spacingy, 0, dimensionz/spacingz);
	vtk_effectivemass_lhl->SetExtent(0, dimensionx/spacingx, 0, dimensiony/spacingy, 0, dimensionz/spacingz);
	vtk_effectivemass_hht->SetExtent(0, dimensionx/spacingx, 0, dimensiony/spacingy, 0, dimensionz/spacingz);
	vtk_effectivemass_hhl->SetExtent(0, dimensionx/spacingx, 0, dimensiony/spacingy, 0, dimensionz/spacingz);

	vtkIdType count = vtk_conductiongama->GetNumberOfPoints();


	vtk_conductiongama->AllocateScalars(VTK_DOUBLE,1);
  // fill the image with foreground voxels:
	double backpotential = background->gap_gamma()+background->ValenceBandEnergy();
	for (vtkIdType i = 0; i < count; ++i){
		vtk_conductiongama->GetPointData()->GetScalars()->SetTuple1(i, backpotential);
	}

	vtk_effectivemass_gamma->AllocateScalars(VTK_DOUBLE,1);
  // fill the image with foreground voxels:
	double backmass_gamma = background->effectiveMass_gamma();
	for (vtkIdType i = 0; i < count; ++i){
		vtk_effectivemass_gamma->GetPointData()->GetScalars()->SetTuple1(i, backmass_gamma);
	}

	vtk_lattice->AllocateScalars(VTK_DOUBLE,1);
  // fill the image with foreground voxels:
	double backlattice = background->latticeparam();
	for (vtkIdType i = 0; i < count; ++i){
		vtk_lattice->GetPointData()->GetScalars()->SetTuple1(i, backlattice);
	}

	vtk_c11->AllocateScalars(VTK_DOUBLE,1);
  // fill the image with foreground voxels:
	double back_c11 = background->getC11();
	for (vtkIdType i = 0; i < count; ++i){
		vtk_c11->GetPointData()->GetScalars()->SetTuple1(i, back_c11);
	}

	vtk_c12->AllocateScalars(VTK_DOUBLE,1);
  // fill the image with foreground voxels:
	double back_c12 = background->getC12();
	for (vtkIdType i = 0; i < count; ++i){
		vtk_c12->GetPointData()->GetScalars()->SetTuple1(i, back_c12);
	}

	vtk_c44->AllocateScalars(VTK_DOUBLE,1);
  // fill the image with foreground voxels:
	double back_c44 = background->getC44();
	for (vtkIdType i = 0; i < count; ++i){
		vtk_c44->GetPointData()->GetScalars()->SetTuple1(i, back_c44);
	}

	vtk_valencelh->AllocateScalars(VTK_DOUBLE,1);
  // fill the image with foreground voxels:
	double back_valencelh = background->ValenceBandEnergy();
	for (vtkIdType i = 0; i < count; ++i){
		vtk_valencelh->GetPointData()->GetScalars()->SetTuple1(i, back_valencelh);
	}

	vtk_valencehh->AllocateScalars(VTK_DOUBLE,1);
  // fill the image with foreground voxels:
	double back_valencehh = background->ValenceBandEnergy();
	for (vtkIdType i = 0; i < count; ++i){
		vtk_valencehh->GetPointData()->GetScalars()->SetTuple1(i, back_valencehh);
	}

	vtk_effectivemass_lht->AllocateScalars(VTK_DOUBLE,1);
  // fill the image with foreground voxels:
	double back_effectivemass_lht = background->effectiveMass_lHt();
	for (vtkIdType i = 0; i < count; ++i){
		vtk_effectivemass_lht->GetPointData()->GetScalars()->SetTuple1(i, back_effectivemass_lht);
	}

	vtk_effectivemass_lhl->AllocateScalars(VTK_DOUBLE,1);
  // fill the image with foreground voxels:
	double back_effectivemass_lhl= background->effectiveMass_lHl();
	for (vtkIdType i = 0; i < count; ++i){
		vtk_effectivemass_lhl->GetPointData()->GetScalars()->SetTuple1(i, back_effectivemass_lhl);
	}

	vtk_effectivemass_hht->AllocateScalars(VTK_DOUBLE,1);
  // fill the image with foreground voxels:
	double back_effectivemass_hht = background->effectiveMass_HHt();
	for (vtkIdType i = 0; i < count; ++i){
		vtk_effectivemass_hht->GetPointData()->GetScalars()->SetTuple1(i, back_effectivemass_hht);
	}

	vtk_effectivemass_hhl->AllocateScalars(VTK_DOUBLE,1);
  // fill the image with foreground voxels:
	double back_effectivemass_hhl = background->effectiveMass_HHl();
	for (vtkIdType i = 0; i < count; ++i){
		vtk_effectivemass_hhl->GetPointData()->GetScalars()->SetTuple1(i, back_effectivemass_hhl);
	}

	std::getline(in,readlines);  // Get the frist line from the file, if any.
	double material_potential=backpotential;
	double material_valencepot_hh=back_valencehh;
	double material_valencepot_lh=back_valencelh;
	double material_massgamma=backmass_gamma;
	double material_massval_lht=back_effectivemass_lht;
	double material_massval_lhl=back_effectivemass_lhl;
	double material_massval_hht=back_effectivemass_hht;
	double material_massval_hhl=back_effectivemass_hhl;
	double material_lattice=backlattice;
	double material_c11=back_c11;
	double material_c12=back_c12;
	double material_c44=back_c44;
	if ( in ) {  // Continue if the line was sucessfully read.
		if(readlines=="#geometrie"){
	        bool haveothers=true;
	        std::string filename;
	        std::string materialname;
	        std::string material_addoon;
  	        std::getline(in,filename);
  	        std::cout << "reading... " << filename << std:: endl;
  	        std::getline(in,materialname);
  	        std::cout << "Material: " << materialname << std:: endl;
  	        std::getline(in,material_addoon);
  	        std::cout << "Add-on: " << material_addoon << std:: endl;
	        while(haveothers){
	        	std::string inputFilename = filename;

	            if(materialname=="GaAs"){
	                  Material* mat = new GaAs(temp);
	                  material_massgamma=mat->effectiveMass_gamma();
	                  material_massval_lht=mat->effectiveMass_lHt();
	                  material_massval_lhl=mat->effectiveMass_lHl();
	                  material_massval_hht=mat->effectiveMass_HHt();
	                  material_massval_hhl=mat->effectiveMass_HHl();
	                  material_lattice=mat->latticeparam();
	                  material_c11=mat->getC11();
	                  material_c12=mat->getC12();
	                  material_c44=mat->getC44();

	                  TYPE exx_yy;
	                  TYPE ezz;
	                  TYPE normalstrain;
	                  TYPE qfactorstrain;
	                  TYPE p;
	                  TYPE q;
	                  exx_yy=(background->latticeparam()-mat->latticeparam())/mat->latticeparam();
	                  ezz=-2*(mat->getC12()/mat->getC11())*exx_yy;
	                  normalstrain=(2.0*exx_yy+ezz);
	                  qfactorstrain=(2.0*exx_yy-2.0*ezz);
	                  p=(mat->getStrain_av()*normalstrain);
	                  q=-(mat->getStrain_b()*0.5*qfactorstrain);
	                  material_potential=mat->gap_gamma()+mat->ValenceBandEnergy()+mat->getStrain_ac()*normalstrain;
	                  material_valencepot_hh = mat->ValenceBandEnergy()-p-q;
	                  material_valencepot_lh = mat->ValenceBandEnergy()-p+q;
	                  delete mat;

	              }
	              else if(materialname=="AlAs"){
	                  Material* mat = new AlAs(temp);
	                  material_massgamma=mat->effectiveMass_gamma();
	                  material_lattice=mat->latticeparam();
	                  material_massval_lht=mat->effectiveMass_lHt();
	                  material_massval_lhl=mat->effectiveMass_lHl();
	                  material_massval_hht=mat->effectiveMass_HHt();
	                  material_massval_hhl=mat->effectiveMass_HHl();
	                  material_c11=mat->getC11();
	                  material_c12=mat->getC12();
	                  material_c44=mat->getC44();

	                  TYPE exx_yy;
	                  TYPE ezz;
	                  TYPE normalstrain;
	                  TYPE qfactorstrain;
	                  TYPE p;
	                  TYPE q;
	                  exx_yy=(background->latticeparam()-mat->latticeparam())/mat->latticeparam();
	                  ezz=-2*(mat->getC12()/mat->getC11())*exx_yy;
	                  normalstrain=(2.0*exx_yy+ezz);
	                  qfactorstrain=(2.0*exx_yy-2.0*ezz);
	                  p=(mat->getStrain_av()*normalstrain);
	                  q=-(mat->getStrain_b()*0.5*qfactorstrain);
	                  material_potential=mat->gap_gamma()+mat->ValenceBandEnergy()+mat->getStrain_ac()*normalstrain;
	                  material_valencepot_hh = mat->ValenceBandEnergy()-p-q;
	                  material_valencepot_lh = mat->ValenceBandEnergy()-p+q;
	                  delete mat;
	              }
	              else if(materialname=="InAs"){
	                  Material* mat = new InAs(temp);
	                  material_massgamma=mat->effectiveMass_gamma();
	                  material_massval_lht=mat->effectiveMass_lHt();
	                  material_massval_lhl=mat->effectiveMass_lHl();
	                  material_massval_hht=mat->effectiveMass_HHt();
	                  material_massval_hhl=mat->effectiveMass_HHl();
	                  material_lattice=mat->latticeparam();
	                  material_c11=mat->getC11();
	                  material_c12=mat->getC12();
	                  material_c44=mat->getC44();

	                  TYPE exx_yy;
	                  TYPE ezz;
	                  TYPE normalstrain;
	                  TYPE qfactorstrain;
	                  TYPE p;
	                  TYPE q;
	                  exx_yy=(background->latticeparam()-mat->latticeparam())/mat->latticeparam();
	                  ezz=-2*(mat->getC12()/mat->getC11())*exx_yy;
	                  normalstrain=(2.0*exx_yy+ezz);
	                  qfactorstrain=(2.0*exx_yy-2.0*ezz);
	                  p=(mat->getStrain_av()*normalstrain);
	                  q=-(mat->getStrain_b()*0.5*qfactorstrain);
	                  material_potential=mat->gap_gamma()+mat->ValenceBandEnergy()+mat->getStrain_ac()*normalstrain;
	                  material_valencepot_hh = mat->ValenceBandEnergy()-p-q;
	                  material_valencepot_lh = mat->ValenceBandEnergy()-p+q;
	                  delete mat;
	              }
	              else if(materialname=="InAlAs"){
	            	  float composition = atof(material_addoon.c_str());
	                  Material* mat = new InAlAs(temp,composition);
	                  material_massgamma=mat->effectiveMass_gamma();
	                  material_massval_lht=mat->effectiveMass_lHt();
	                  material_massval_lhl=mat->effectiveMass_lHl();
	                  material_massval_hht=mat->effectiveMass_HHt();
	                  material_massval_hhl=mat->effectiveMass_HHl();
	                  material_lattice=mat->latticeparam();
	                  material_c11=mat->getC11();
	                  material_c12=mat->getC12();
	                  material_c44=mat->getC44();

	                  TYPE exx_yy;
	                  TYPE ezz;
	                  TYPE normalstrain;
	                  TYPE qfactorstrain;
	                  TYPE p;
	                  TYPE q;
	                  exx_yy=(background->latticeparam()-mat->latticeparam())/mat->latticeparam();
	                  ezz=-2*(mat->getC12()/mat->getC11())*exx_yy;
	                  normalstrain=(2.0*exx_yy+ezz);
	                  qfactorstrain=(2.0*exx_yy-2.0*ezz);
	                  p=(mat->getStrain_av()*normalstrain);
	                  q=-(mat->getStrain_b()*0.5*qfactorstrain);
	                  material_potential=mat->gap_gamma()+mat->ValenceBandEnergy()+mat->getStrain_ac()*normalstrain;
	                  material_valencepot_hh = mat->ValenceBandEnergy()-p-q;
	                  material_valencepot_lh = mat->ValenceBandEnergy()-p+q;
	                  delete mat;
	              }
	              else if(materialname=="InGaAs"){
	            	  float composition = atof(material_addoon.c_str());
	                  Material* mat = new InGaAs(temp,composition);
	                  material_massgamma=mat->effectiveMass_gamma();
	                  material_massval_lht=mat->effectiveMass_lHt();
	                  material_massval_lhl=mat->effectiveMass_lHl();
	                  material_massval_hht=mat->effectiveMass_HHt();
	                  material_massval_hhl=mat->effectiveMass_HHl();
	                  material_lattice=mat->latticeparam();
	                  material_c11=mat->getC11();
	                  material_c12=mat->getC12();
	                  material_c44=mat->getC44();

	                  TYPE exx_yy;
	                  TYPE ezz;
	                  TYPE normalstrain;
	                  TYPE qfactorstrain;
	                  TYPE p;
	                  TYPE q;
	                  exx_yy=(background->latticeparam()-mat->latticeparam())/mat->latticeparam();
	                  ezz=-2*(mat->getC12()/mat->getC11())*exx_yy;
	                  normalstrain=(2.0*exx_yy+ezz);
	                  qfactorstrain=(2.0*exx_yy-2.0*ezz);
	                  p=(mat->getStrain_av()*normalstrain);
	                  q=-(mat->getStrain_b()*0.5*qfactorstrain);
	                  material_potential=mat->gap_gamma()+mat->ValenceBandEnergy()+mat->getStrain_ac()*normalstrain;
	                  material_valencepot_hh = mat->ValenceBandEnergy()-p-q;
	                  material_valencepot_lh = mat->ValenceBandEnergy()-p+q;
	                  delete mat;
	              }
	              else if(materialname=="AlGaAs"){
	            	  float composition = atof(material_addoon.c_str());
	                  Material* mat = new AlGaAs(temp,composition);
	                  material_massgamma=mat->effectiveMass_gamma();
	                  material_massval_lht=mat->effectiveMass_lHt();
	                  material_massval_lhl=mat->effectiveMass_lHl();
	                  material_massval_hht=mat->effectiveMass_HHt();
	                  material_massval_hhl=mat->effectiveMass_HHl();
	                  material_lattice=mat->latticeparam();
	                  material_c11=mat->getC11();
	                  material_c12=mat->getC12();
	                  material_c44=mat->getC44();

	                  TYPE exx_yy;
	                  TYPE ezz;
	                  TYPE normalstrain;
	                  TYPE qfactorstrain;
	                  TYPE p;
	                  TYPE q;
	                  exx_yy=(background->latticeparam()-mat->latticeparam())/mat->latticeparam();
	                  ezz=-2*(mat->getC12()/mat->getC11())*exx_yy;
	                  normalstrain=(2.0*exx_yy+ezz);
	                  qfactorstrain=(2.0*exx_yy-2.0*ezz);
	                  p=(mat->getStrain_av()*normalstrain);
	                  q=-(mat->getStrain_b()*0.5*qfactorstrain);
	                  material_potential=mat->gap_gamma()+mat->ValenceBandEnergy()+mat->getStrain_ac()*normalstrain;
	                  material_valencepot_hh = mat->ValenceBandEnergy()-p-q;
	                  material_valencepot_lh = mat->ValenceBandEnergy()-p+q;
	                  delete mat;
	              }
	              else if(materialname=="InAlGaAs"){
	            	  float compositionin,compositional;
	            	  std::sscanf(material_addoon.c_str(),"%f %f",&compositionin,&compositional);
	            	  std::cout << "In%: " << compositionin << " Al%: " << compositional << std::endl;
	                  Material* mat = new InAlGaAs(temp,compositionin,compositional);
	                  material_massgamma=mat->effectiveMass_gamma();
	                  material_massval_lht=mat->effectiveMass_lHt();
	                  material_massval_lhl=mat->effectiveMass_lHl();
	                  material_massval_hht=mat->effectiveMass_HHt();
	                  material_massval_hhl=mat->effectiveMass_HHl();
	                  material_lattice=mat->latticeparam();
	                  material_c11=mat->getC11();
	                  material_c12=mat->getC12();
	                  material_c44=mat->getC44();

	                  TYPE exx_yy;
	                  TYPE ezz;
	                  TYPE normalstrain;
	                  TYPE qfactorstrain;
	                  TYPE p;
	                  TYPE q;
	                  exx_yy=(background->latticeparam()-mat->latticeparam())/mat->latticeparam();
	                  ezz=-2*(mat->getC12()/mat->getC11())*exx_yy;
	                  normalstrain=(2.0*exx_yy+ezz);
	                  qfactorstrain=(2.0*exx_yy-2.0*ezz);
	                  p=(mat->getStrain_av()*normalstrain);
	                  q=-(mat->getStrain_b()*0.5*qfactorstrain);
	                  material_potential=mat->gap_gamma()+mat->ValenceBandEnergy()+mat->getStrain_ac()*normalstrain;
	                  material_valencepot_hh = mat->ValenceBandEnergy()-p-q;
	                  material_valencepot_lh = mat->ValenceBandEnergy()-p+q;
	                  delete mat;
	              }
	              else
	              {
	            	  std::cout << "Invalid material in geometry file" << std::endl;
	              }
	              vtkSmartPointer<vtkImageData> imagedata2 = vtkImageData::New();
	              imagedata2->DeepCopy(vtk_conductiongama);
	              vtkSmartPointer<vtkImageData> imagedata3 = vtkImageData::New();
	              imagedata3->DeepCopy(vtk_effectivemass_gamma);
	              vtkSmartPointer<vtkImageData> imagedata4 = vtkImageData::New();
	              imagedata4->DeepCopy(vtk_lattice);
	              vtkSmartPointer<vtkImageData> imagedata5 = vtkImageData::New();
	              imagedata5->DeepCopy(vtk_c11);
	              vtkSmartPointer<vtkImageData> imagedata6 = vtkImageData::New();
	              imagedata6->DeepCopy(vtk_c12);
	              vtkSmartPointer<vtkImageData> imagedata7 = vtkImageData::New();
	              imagedata7->DeepCopy(vtk_c44);
	              vtkSmartPointer<vtkImageData> imagedata8 = vtkImageData::New();
	              imagedata8->DeepCopy(vtk_valencehh);
	              vtkSmartPointer<vtkImageData> imagedata9 = vtkImageData::New();
	              imagedata9->DeepCopy(vtk_valencelh);
	              vtkSmartPointer<vtkImageData> imagedata10 = vtkImageData::New();
	              imagedata10->DeepCopy(vtk_effectivemass_hht);
	              vtkSmartPointer<vtkImageData> imagedata11 = vtkImageData::New();
	              imagedata11->DeepCopy(vtk_effectivemass_hhl);
	              vtkSmartPointer<vtkImageData> imagedata12 = vtkImageData::New();
	              imagedata12->DeepCopy(vtk_effectivemass_lht);
	              vtkSmartPointer<vtkImageData> imagedata13 = vtkImageData::New();
	              imagedata13->DeepCopy(vtk_effectivemass_lhl);

	              vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
	              reader->SetFileName(inputFilename.c_str());
	              reader->Update();
	              // Visualize
	              vtkSmartPointer<vtkPolyData> pd = reader->GetOutput();
	              reader->Update();
	              // polygonal data --> image stencil:
	              vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc =
	              vtkSmartPointer<vtkPolyDataToImageStencil>::New();
	              vtkSmartPointer<vtkPolyDataToImageStencil> pol3stenc =
	              vtkSmartPointer<vtkPolyDataToImageStencil>::New();
	              vtkSmartPointer<vtkPolyDataToImageStencil> pol4stenc =
	              vtkSmartPointer<vtkPolyDataToImageStencil>::New();
	              vtkSmartPointer<vtkPolyDataToImageStencil> pol5stenc =
	              vtkSmartPointer<vtkPolyDataToImageStencil>::New();
	              vtkSmartPointer<vtkPolyDataToImageStencil> pol6stenc =
	              vtkSmartPointer<vtkPolyDataToImageStencil>::New();
	              vtkSmartPointer<vtkPolyDataToImageStencil> pol7stenc =
	              vtkSmartPointer<vtkPolyDataToImageStencil>::New();
	              vtkSmartPointer<vtkPolyDataToImageStencil> pol8stenc =
	              vtkSmartPointer<vtkPolyDataToImageStencil>::New();
	              vtkSmartPointer<vtkPolyDataToImageStencil> pol9stenc =
	              vtkSmartPointer<vtkPolyDataToImageStencil>::New();
	              vtkSmartPointer<vtkPolyDataToImageStencil> pol10stenc =
	              vtkSmartPointer<vtkPolyDataToImageStencil>::New();
	              vtkSmartPointer<vtkPolyDataToImageStencil> pol11stenc =
	              vtkSmartPointer<vtkPolyDataToImageStencil>::New();
	              vtkSmartPointer<vtkPolyDataToImageStencil> pol12stenc =
	              vtkSmartPointer<vtkPolyDataToImageStencil>::New();
	              vtkSmartPointer<vtkPolyDataToImageStencil> pol13stenc =
	              vtkSmartPointer<vtkPolyDataToImageStencil>::New();

	              pol2stenc->SetInputData(pd);
	              pol3stenc->SetInputData(pd);
	              pol4stenc->SetInputData(pd);
	              pol5stenc->SetInputData(pd);
	              pol6stenc->SetInputData(pd);
	              pol7stenc->SetInputData(pd);
	              pol8stenc->SetInputData(pd);
	              pol9stenc->SetInputData(pd);
	              pol10stenc->SetInputData(pd);
	              pol11stenc->SetInputData(pd);
	              pol12stenc->SetInputData(pd);
	              pol13stenc->SetInputData(pd);


	              pol2stenc->SetOutputSpacing(vtkspacing);
	              pol2stenc->SetOutputWholeExtent(imagedata2->GetExtent());
	              pol2stenc->Update();
	              pol3stenc->SetOutputSpacing(vtkspacing);
	              pol3stenc->SetOutputWholeExtent(imagedata3->GetExtent());
	              pol3stenc->Update();
	              pol4stenc->SetOutputSpacing(vtkspacing);
	              pol4stenc->SetOutputWholeExtent(imagedata4->GetExtent());
	              pol4stenc->Update();
	              pol5stenc->SetOutputSpacing(vtkspacing);
	              pol5stenc->SetOutputWholeExtent(imagedata5->GetExtent());
	              pol5stenc->Update();
	              pol6stenc->SetOutputSpacing(vtkspacing);
	              pol6stenc->SetOutputWholeExtent(imagedata6->GetExtent());
	              pol6stenc->Update();
	              pol7stenc->SetOutputSpacing(vtkspacing);
	              pol7stenc->SetOutputWholeExtent(imagedata7->GetExtent());
	              pol7stenc->Update();
	              pol7stenc->SetOutputSpacing(vtkspacing);
	              pol7stenc->SetOutputWholeExtent(imagedata7->GetExtent());
	              pol7stenc->Update();
	              pol8stenc->SetOutputSpacing(vtkspacing);
	              pol8stenc->SetOutputWholeExtent(imagedata7->GetExtent());
	              pol8stenc->Update();
	              pol9stenc->SetOutputSpacing(vtkspacing);
	              pol9stenc->SetOutputWholeExtent(imagedata7->GetExtent());
	              pol9stenc->Update();
	              pol10stenc->SetOutputSpacing(vtkspacing);
	              pol10stenc->SetOutputWholeExtent(imagedata7->GetExtent());
	              pol10stenc->Update();
	              pol11stenc->SetOutputSpacing(vtkspacing);
	              pol11stenc->SetOutputWholeExtent(imagedata7->GetExtent());
	              pol11stenc->Update();
	              pol12stenc->SetOutputSpacing(vtkspacing);
	              pol12stenc->SetOutputWholeExtent(imagedata7->GetExtent());
	              pol12stenc->Update();
	              pol13stenc->SetOutputSpacing(vtkspacing);
	              pol13stenc->SetOutputWholeExtent(imagedata7->GetExtent());
	              pol13stenc->Update();

	              // cut the corresponding white image and set the background:
	              vtkSmartPointer<vtkImageStencil> imgstenc =
	              vtkSmartPointer<vtkImageStencil>::New();
	              vtkSmartPointer<vtkImageStencil> imgstenc2 =
	              vtkSmartPointer<vtkImageStencil>::New();
	              vtkSmartPointer<vtkImageStencil> imgstenc3 =
	              vtkSmartPointer<vtkImageStencil>::New();
	              vtkSmartPointer<vtkImageStencil> imgstenc4 =
	              vtkSmartPointer<vtkImageStencil>::New();
	              vtkSmartPointer<vtkImageStencil> imgstenc5 =
	              vtkSmartPointer<vtkImageStencil>::New();
	              vtkSmartPointer<vtkImageStencil> imgstenc6 =
	              vtkSmartPointer<vtkImageStencil>::New();
	              vtkSmartPointer<vtkImageStencil> imgstenc7 =
	              vtkSmartPointer<vtkImageStencil>::New();
	              vtkSmartPointer<vtkImageStencil> imgstenc8 =
	              vtkSmartPointer<vtkImageStencil>::New();
	              vtkSmartPointer<vtkImageStencil> imgstenc9 =
	              vtkSmartPointer<vtkImageStencil>::New();
	              vtkSmartPointer<vtkImageStencil> imgstenc10 =
	              vtkSmartPointer<vtkImageStencil>::New();
	              vtkSmartPointer<vtkImageStencil> imgstenc11 =
	              vtkSmartPointer<vtkImageStencil>::New();
	              vtkSmartPointer<vtkImageStencil> imgstenc12 =
	              vtkSmartPointer<vtkImageStencil>::New();

	              imgstenc->SetInputData(imagedata2);
	              imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
	              imgstenc2->SetInputData(imagedata3);
	              imgstenc2->SetStencilConnection(pol3stenc->GetOutputPort());
	              imgstenc3->SetInputData(imagedata4);
	              imgstenc3->SetStencilConnection(pol4stenc->GetOutputPort());
	              imgstenc4->SetInputData(imagedata5);
	              imgstenc4->SetStencilConnection(pol5stenc->GetOutputPort());
	              imgstenc5->SetInputData(imagedata6);
	              imgstenc5->SetStencilConnection(pol6stenc->GetOutputPort());
	              imgstenc6->SetInputData(imagedata7);
	              imgstenc6->SetStencilConnection(pol7stenc->GetOutputPort());
	              imgstenc7->SetInputData(imagedata8);
	              imgstenc7->SetStencilConnection(pol8stenc->GetOutputPort());
	              imgstenc8->SetInputData(imagedata9);
	              imgstenc8->SetStencilConnection(pol9stenc->GetOutputPort());
	              imgstenc9->SetInputData(imagedata10);
	              imgstenc9->SetStencilConnection(pol10stenc->GetOutputPort());
	              imgstenc10->SetInputData(imagedata11);
	              imgstenc10->SetStencilConnection(pol11stenc->GetOutputPort());
	              imgstenc11->SetInputData(imagedata12);
	              imgstenc11->SetStencilConnection(pol12stenc->GetOutputPort());
	              imgstenc12->SetInputData(imagedata13);
	              imgstenc12->SetStencilConnection(pol13stenc->GetOutputPort());

	              imgstenc->ReverseStencilOn();
	              imgstenc->SetBackgroundValue(material_potential);
	              imgstenc->Update();
	              imgstenc2->ReverseStencilOn();
	              imgstenc2->SetBackgroundValue(material_massgamma);
	              imgstenc2->Update();
	              imgstenc3->ReverseStencilOn();
	              imgstenc3->SetBackgroundValue(material_lattice);
	              imgstenc3->Update();
	              imgstenc4->ReverseStencilOn();
	              imgstenc4->SetBackgroundValue(material_c11);
	              imgstenc4->Update();
	              imgstenc5->ReverseStencilOn();
	              imgstenc5->SetBackgroundValue(material_c12);
	              imgstenc5->Update();
	              imgstenc6->ReverseStencilOn();
	              imgstenc6->SetBackgroundValue(material_c44);
	              imgstenc6->Update();

	              imgstenc7->ReverseStencilOn();
	              imgstenc7->SetBackgroundValue(material_valencepot_hh);
	              imgstenc7->Update();
	              imgstenc8->ReverseStencilOn();
	              imgstenc8->SetBackgroundValue(material_valencepot_lh);
	              imgstenc8->Update();
	              imgstenc9->ReverseStencilOn();
	              imgstenc9->SetBackgroundValue(material_massval_hht);
	              imgstenc9->Update();
	              imgstenc10->ReverseStencilOn();
	              imgstenc10->SetBackgroundValue(material_massval_hhl);
	              imgstenc10->Update();
	              imgstenc11->ReverseStencilOn();
	              imgstenc11->SetBackgroundValue(material_massval_lht);
	              imgstenc11->Update();
	              imgstenc12->ReverseStencilOn();
	              imgstenc12->SetBackgroundValue(material_massval_lhl);
	              imgstenc12->Update();


	              vtk_conductiongama->DeepCopy(imgstenc->GetOutput());
	              vtk_effectivemass_gamma->DeepCopy(imgstenc2->GetOutput());
	              vtk_lattice->DeepCopy(imgstenc3->GetOutput());
	              vtk_c11->DeepCopy(imgstenc4->GetOutput());
	              vtk_c12->DeepCopy(imgstenc5->GetOutput());
	              vtk_c44->DeepCopy(imgstenc6->GetOutput());
	              vtk_valencehh->DeepCopy(imgstenc7->GetOutput());
	              vtk_valencelh->DeepCopy(imgstenc8->GetOutput());
	              vtk_effectivemass_hht->DeepCopy(imgstenc9->GetOutput());
	              vtk_effectivemass_hhl->DeepCopy(imgstenc10->GetOutput());
	              vtk_effectivemass_lht->DeepCopy(imgstenc11->GetOutput());
	              vtk_effectivemass_lhl->DeepCopy(imgstenc12->GetOutput());

	              std::getline(in,readlines);
	              if(readlines[0]=='#')
	            	  haveothers=false;
	              else{
	            	  filename=readlines;
	        	      std::cout << "reading... " << filename << std:: endl;
	            	  std::getline(in,materialname);
	            	  std::cout << "Material: " << materialname << std:: endl;
	            	  std::getline(in,material_addoon);
	            	  std::cout << "Add-on: " << material_addoon << std:: endl;
	              }
	          }
	      	}
	    	else
	    		std::cout << "#geometrie section in geometry file format error (stupid!)" << std::endl;

	  	}



}

template <typename TYPE>
Heterostructure3D<TYPE>::~Heterostructure3D() {
	delete basegrid;
}

template <typename TYPE>
void Heterostructure3D<TYPE>::printPotential(std::string filename){
	vtkSmartPointer<vtkXMLImageDataWriter> writer =
	vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(vtk_conductiongama);
	writer->Write();
}

template <typename TYPE>
std::shared_ptr<Material> Heterostructure3D<TYPE>::getBackground(){
	return background;
}

template <typename TYPE>
DiscreteFunction3D<TYPE,TYPE> Heterostructure3D<TYPE>::ConductionBand_gamma() const{
	DiscreteFunction3D<TYPE,TYPE> func(*basegrid);
	for(long int k=0;k<basegrid->getSizeX();k++){
		for(long int l=0;l<basegrid->getSizeY();l++){
			for(long int m=0;m<basegrid->getSizeZ();m++){
				func(k,l,m)=*(static_cast<TYPE*>(vtk_conductiongama->GetScalarPointer(k,l,m)));
			}
		}
	}
	return func;
}

template <typename TYPE>
DiscreteFunction3D<TYPE,TYPE> Heterostructure3D<TYPE>::effectiveMass_gamma() const{
	DiscreteFunction3D<TYPE,TYPE> func(*basegrid);
	for(long int k=0;k<basegrid->getSizeX();k++){
		for(long int l=0;l<basegrid->getSizeY();l++){
			for(long int m=0;m<basegrid->getSizeZ();m++){
				func(k,l,m)=*(static_cast<TYPE*>(vtk_effectivemass_gamma->GetScalarPointer(k,l,m)));
			}
		}
	}
	return func;
}

template <typename TYPE>
DiscreteFunction3D<TYPE,TYPE> Heterostructure3D<TYPE>::latticeparam() const{
	DiscreteFunction3D<TYPE,TYPE> func(*basegrid);
	for(long int k=0;k<basegrid->getSizeX();k++){
		for(long int l=0;l<basegrid->getSizeY();l++){
			for(long int m=0;m<basegrid->getSizeZ();m++){
				func(k,l,m)=*(static_cast<TYPE*>(vtk_lattice->GetScalarPointer(k,l,m)));
			}
		}
	}
	return func;
}

template <typename TYPE>
DiscreteFunction3D<TYPE,TYPE> Heterostructure3D<TYPE>::C11() const{
	DiscreteFunction3D<TYPE,TYPE> func(*basegrid);
	for(long int k=0;k<basegrid->getSizeX();k++){
		for(long int l=0;l<basegrid->getSizeY();l++){
			for(long int m=0;m<basegrid->getSizeZ();m++){
				func(k,l,m)=*(static_cast<TYPE*>(vtk_c11->GetScalarPointer(k,l,m)));
			}
		}
	}
	return func;
}

template <typename TYPE>
DiscreteFunction3D<TYPE,TYPE> Heterostructure3D<TYPE>::C12() const{
	DiscreteFunction3D<TYPE,TYPE> func(*basegrid);
	for(long int k=0;k<basegrid->getSizeX();k++){
		for(long int l=0;l<basegrid->getSizeY();l++){
			for(long int m=0;m<basegrid->getSizeZ();m++){
				func(k,l,m)=*(static_cast<TYPE*>(vtk_c12->GetScalarPointer(k,l,m)));
			}
		}
	}
	return func;
}

template <typename TYPE>
DiscreteFunction3D<TYPE,TYPE> Heterostructure3D<TYPE>::C44() const{
	DiscreteFunction3D<TYPE,TYPE> func(*basegrid);
	for(long int k=0;k<basegrid->getSizeX();k++){
		for(long int l=0;l<basegrid->getSizeY();l++){
			for(long int m=0;m<basegrid->getSizeZ();m++){
				func(k,l,m)=*(static_cast<TYPE*>(vtk_c11->GetScalarPointer(k,l,m)));
			}
		}
	}
	return func;
}

////

template <typename TYPE>
DiscreteFunction3D<TYPE,TYPE> Heterostructure3D<TYPE>::ValenceBand_HH() const{
	DiscreteFunction3D<TYPE,TYPE> func(*basegrid);
	for(long int k=0;k<basegrid->getSizeX();k++){
		for(long int l=0;l<basegrid->getSizeY();l++){
			for(long int m=0;m<basegrid->getSizeZ();m++){
				func(k,l,m)=*(static_cast<TYPE*>(vtk_valencehh->GetScalarPointer(k,l,m)));
			}
		}
	}
	return func;
}
template <typename TYPE>
DiscreteFunction3D<TYPE,TYPE> Heterostructure3D<TYPE>::ValenceBand_lH() const{
	DiscreteFunction3D<TYPE,TYPE> func(*basegrid);
	for(long int k=0;k<basegrid->getSizeX();k++){
		for(long int l=0;l<basegrid->getSizeY();l++){
			for(long int m=0;m<basegrid->getSizeZ();m++){
				func(k,l,m)=*(static_cast<TYPE*>(vtk_valencelh->GetScalarPointer(k,l,m)));
			}
		}
	}
	return func;
}
template <typename TYPE>
DiscreteFunction3D<TYPE,TYPE> Heterostructure3D<TYPE>::effectiveMass_HHt() const{
	DiscreteFunction3D<TYPE,TYPE> func(*basegrid);
	for(long int k=0;k<basegrid->getSizeX();k++){
		for(long int l=0;l<basegrid->getSizeY();l++){
			for(long int m=0;m<basegrid->getSizeZ();m++){
				func(k,l,m)=*(static_cast<TYPE*>(vtk_effectivemass_hht->GetScalarPointer(k,l,m)));
			}
		}
	}
	return func;
}
template <typename TYPE>
DiscreteFunction3D<TYPE,TYPE> Heterostructure3D<TYPE>::effectiveMass_HHl() const{
	DiscreteFunction3D<TYPE,TYPE> func(*basegrid);
	for(long int k=0;k<basegrid->getSizeX();k++){
		for(long int l=0;l<basegrid->getSizeY();l++){
			for(long int m=0;m<basegrid->getSizeZ();m++){
				func(k,l,m)=*(static_cast<TYPE*>(vtk_effectivemass_hhl->GetScalarPointer(k,l,m)));
			}
		}
	}
	return func;
}
template <typename TYPE>
DiscreteFunction3D<TYPE,TYPE> Heterostructure3D<TYPE>::effectiveMass_lHt() const{
	DiscreteFunction3D<TYPE,TYPE> func(*basegrid);
	for(long int k=0;k<basegrid->getSizeX();k++){
		for(long int l=0;l<basegrid->getSizeY();l++){
			for(long int m=0;m<basegrid->getSizeZ();m++){
				func(k,l,m)=*(static_cast<TYPE*>(vtk_effectivemass_lht->GetScalarPointer(k,l,m)));
			}
		}
	}
	return func;
}
template <typename TYPE>
DiscreteFunction3D<TYPE,TYPE> Heterostructure3D<TYPE>::effectiveMass_lHl() const{
	DiscreteFunction3D<TYPE,TYPE> func(*basegrid);
	for(long int k=0;k<basegrid->getSizeX();k++){
		for(long int l=0;l<basegrid->getSizeY();l++){
			for(long int m=0;m<basegrid->getSizeZ();m++){
				func(k,l,m)=*(static_cast<TYPE*>(vtk_effectivemass_lhl->GetScalarPointer(k,l,m)));
			}
		}
	}
	return func;
}



template <typename TYPE>
Grid3D<TYPE> Heterostructure3D<TYPE>::getGrid(){
	return *basegrid;
}



} /* namespace epital */

#endif
