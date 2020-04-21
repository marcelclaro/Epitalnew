/*
 * Heterostructure3D.hpp
 *
 *  Created on: Nov 5, 2014
 *      Author: marcel
 */



#ifndef HETEROSTRUCTURE3D_HPP_
#define HETEROSTRUCTURE3D_HPP_

#include <armadillo>
#include <string.h>
#include <memory>
#include "Grid.hpp"
#include "CustomTypes.hpp"
#include "CMaterial.hpp"
#include "Materials.hpp"
#include "Units.hpp"
#include "DFunction3D.hpp"
#include <vtkPolyData.h>
#include <vtkSTLReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVersion.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>


namespace epital {

template <typename TYPE>
class Heterostructure3D {
public:
	Grid3D<TYPE>* basegrid;
	std::string geometryfile;
	std::shared_ptr<Material> background;
	Temperature temp;

	vtkImageData* vtk_conductiongama;
	vtkImageData* vtk_effectivemass_gamma;
	vtkImageData* vtk_lattice;
	vtkImageData* vtk_valencelh;
	vtkImageData* vtk_valencehh;
	vtkImageData* vtk_effectivemass_lht;
	vtkImageData* vtk_effectivemass_lhl;
	vtkImageData* vtk_effectivemass_hht;
	vtkImageData* vtk_effectivemass_hhl;
	vtkImageData* vtk_c11;
	vtkImageData* vtk_c12;
	vtkImageData* vtk_c44;

	Heterostructure3D(std::string geometryfile,std::shared_ptr<Material> background, Temperature temp = 77.0);
	virtual ~Heterostructure3D();

	void printPotential(std::string filename);

	DiscreteFunction3D<TYPE,TYPE> ConductionBand_gamma() const;
	DiscreteFunction3D<TYPE,TYPE> effectiveMass_gamma() const;
	DiscreteFunction3D<TYPE,TYPE> effectiveMass_HHt() const;
	DiscreteFunction3D<TYPE,TYPE> effectiveMass_HHl() const;
	DiscreteFunction3D<TYPE,TYPE> effectiveMass_lHt() const;
	DiscreteFunction3D<TYPE,TYPE> effectiveMass_lHl() const;
	DiscreteFunction3D<TYPE,TYPE> ValenceBand_HH() const;
	DiscreteFunction3D<TYPE,TYPE> ValenceBand_lH() const;
	DiscreteFunction3D<TYPE,TYPE> latticeparam() const;
	DiscreteFunction3D<TYPE,TYPE> C11() const;
	DiscreteFunction3D<TYPE,TYPE> C12() const;
	DiscreteFunction3D<TYPE,TYPE> C44() const;
	std::shared_ptr<Material> getBackground();

	Grid3D<TYPE> getGrid();


};

} /* namespace epital */

#include "Heterostructure3D.cpp"

#endif /* HETEROSTRUCTURE3D_HPP_ */
