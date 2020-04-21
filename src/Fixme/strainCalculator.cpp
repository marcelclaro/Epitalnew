/*
 * strainCalculator.cpp
 *
 *  Created on: Dec 18, 2014
 *      Author: marcelclaro
 */

#include <strainCalculator.hpp>

#ifndef STRAINCALCULATOR_CPP_
#define STRAINCALCULATOR_CPP_

namespace epital {

template <typename TYPE>
strainCalculator<TYPE>::strainCalculator() {
	// TODO Auto-generated constructor stub

}

template <typename TYPE>
strainCalculator<TYPE>::~strainCalculator() {
	// TODO Auto-generated destructor stub
}

template <typename TYPE>
void strainCalculator<TYPE>::solveDisplace(Heterostructure3D<TYPE>& structure, DiscreteFunction3D<TYPE,TYPE>& ux,DiscreteFunction3D<TYPE,TYPE>& uy,DiscreteFunction3D<TYPE,TYPE>& uz){
	DiscreteFunction3D<TYPE,TYPE> lattice(structure.latticeparam());
	DiscreteFunction3D<TYPE,TYPE> C11(structure.C11());
	DiscreteFunction3D<TYPE,TYPE> C12(structure.C12());
	DiscreteFunction3D<TYPE,TYPE> C44(structure.C44());

	Grid3D<TYPE> basegrid = structure.getGrid();
	long int xsize = basegrid.getSizeX();
	long int ysize = basegrid.getSizeY();
	long int zsize = basegrid.getSizeZ();
	long int zpitch = xsize*ysize;

	long int fieldsize=xsize*ysize*zsize;

	long int xmax = xsize-1;
	long int ymax = ysize-1;
	long int zmax = zsize-1;

	TYPE dx = basegrid.getxIncrement();
	TYPE dy = basegrid.getyIncrement();
	TYPE dz = basegrid.getzIncrement();

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;

	tripletList.reserve(3*fieldsize*21);
	lattice.saveHDF5("lattice.hdf");

	std::shared_ptr<Material> substrate = structure.getBackground();

	DiscreteFunction3D<TYPE,TYPE> eps(structure.latticeparam());
	TYPE* data = eps.getdata();
	for(long int i=0; i< eps.getGrid().getSize();++i){
		*data=(*data-substrate->latticeparam())/substrate->latticeparam();
	}

	DiscreteFunction3D<TYPE,TYPE> epsdx(structure.getGrid());
	DiscreteFunction3D<TYPE,TYPE> epsdy(structure.getGrid());
	DiscreteFunction3D<TYPE,TYPE> epsdz(structure.getGrid());

	eps.Divergent(epsdx,epsdy,epsdz);

	C11/=4.0;
	C12/=4.0;
	C44/=4.0;


	C11.smooth_mean();
	C12.smooth_mean();
	C44.smooth_mean();
/*
	for(long int l=0; l< ysize;++l){
		for(long int k=0; k< xsize;++k){
			C11(k,l,0)=0;
			C12(k,l,0)=0;
			C44(k,l,0)=0;
			C11(k,l,zmax)=0;
			C12(k,l,zmax)=0;
			C44(k,l,zmax)=0;
		}
	}
	for(long int m=1; m< (zsize-1);++m){
		for(long int k=0; k< xsize;++k){
						C11(k,0,m)=0;
						C12(k,0,m)=0;
						C44(k,0,m)=0;
						C11(k,ymax,m)=0;
						C12(k,ymax,m)=0;
						C44(k,ymax,m)=0;
		}
		for(long int l=1; l< ymax;++l){
			for(long int k=0; k< xsize;k+=xmax){
				C11(0,l,m)=0;
				C12(0,l,m)=0;
				C44(0,l,m)=0;
				C11(xmax,l,m)=0;
				C12(xmax,l,m)=0;
				C44(xmax,l,m)=0;
			}
		}
	}

*/
	for(long int m=0; m< zsize;++m){
		for(long int l=0; l< ysize;++l){
			for(long int k=0; k< xsize;++k){
				long int pos=k+l*xsize+m*xsize*ysize;
				TYPE alpha=C11(k,l,m);
				TYPE beta=C44(k,l,m);
				TYPE gama=C12(k,l,m);

				/*ux*/
				tripletList.push_back(T(pos,pos,4.0*alpha/(dx*dx)+2.0*beta/(dy*dy)+2.0*beta/(dz*dz)));
				if(k>1)
					tripletList.push_back(T(pos,pos-1,(-4.0*alpha)/(dx*dx)));
				if(k>2)
					tripletList.push_back(T(pos,pos-2,(2.0*alpha)/(dx*dx)));
				if(k<xmax)
					tripletList.push_back(T(pos,pos+1,(-4.0*alpha)/(dx*dx)));
				if(k<(xmax-1))
					tripletList.push_back(T(pos,pos+2,(2.0*alpha)/(dx*dx)));

				if(l>1)
					tripletList.push_back(T(pos,pos-xsize,(-2.0*beta)/(dy*dy)));
				if(l>2)
					tripletList.push_back(T(pos,pos-2*xsize,(1.0*beta)/(dy*dy)));
				if(l<ymax)
					tripletList.push_back(T(pos,pos+xsize,(-2.0*beta)/(dy*dy)));
				if(l<(ymax-1))
					tripletList.push_back(T(pos,pos+2*xsize,(1.0*beta)/(dy*dy)));

				if(m>1)
					tripletList.push_back(T(pos,pos-zpitch,(-2.0*beta)/(dz*dz)));
				if(m>2)
					tripletList.push_back(T(pos,pos-2*zpitch,(1.0*beta)/(dz*dz)));
				if(m<zmax)
					tripletList.push_back(T(pos,pos+zpitch,(-2.0*beta)/(dz*dz)));
				if(m<(zmax-1))
					tripletList.push_back(T(pos,pos+2*zpitch,(1.0*beta)/(dz*dz)));

				if(k>1){
					if(l>1)
						tripletList.push_back(T(pos,pos+fieldsize-1-xsize,(+4.0*(beta+gama))/(4.0*dy*dx)));
					if(l<ymax)
						tripletList.push_back(T(pos,pos+fieldsize-1+xsize,(-4.0*(beta+gama))/(4.0*dy*dx)));
				}
				if(k<xmax){
					if(l>1)
						tripletList.push_back(T(pos,pos+fieldsize+1-xsize,(-4.0*(beta+gama))/(4.0*dy*dx)));
					if(l<ymax)
						tripletList.push_back(T(pos,pos+fieldsize+1+xsize,(+4.0*(beta+gama))/(4.0*dy*dx)));
				}


				if(k>1){
					if(m>1)
						tripletList.push_back(T(pos,pos+2*fieldsize-1-zpitch,(+4.0*(beta+gama))/(4.0*dz*dx)));
					if(m<zmax)
						tripletList.push_back(T(pos,pos+2*fieldsize-1+zpitch,(-4.0*(beta+gama))/(4.0*dz*dx)));
				}
				if(k<xmax){
					if(m>1)
						tripletList.push_back(T(pos,pos+2*fieldsize+1-zpitch,(-4.0*(beta+gama))/(4.0*dz*dx)));
					if(m<zmax)
						tripletList.push_back(T(pos,pos+2*fieldsize+1+zpitch,(+4.0*(beta+gama))/(4.0*dz*dx)));
				}

				long int posuy = pos+fieldsize;
				/*uy*/
				tripletList.push_back(T(posuy,posuy,2*beta/(dx*dx)+4*alpha/(dy*dy)+2*beta/(dz*dz)));
				if(k>1)
					tripletList.push_back(T(posuy,posuy-1,(-2.0*beta)/(dx*dx)));
				if(k>2)
					tripletList.push_back(T(posuy,posuy-2,(1.0*beta)/(dx*dx)));
				if(k<xmax)
					tripletList.push_back(T(posuy,posuy+1,(-2.0*beta)/(dx*dx)));
				if(k<(xmax-1))
					tripletList.push_back(T(posuy,posuy+2,(1.0*beta)/(dx*dx)));

				if(l>1)
					tripletList.push_back(T(posuy,posuy-xsize,(-2.0*alpha)/(dy*dy)));
				if(l>2)
					tripletList.push_back(T(posuy,posuy-2*xsize,(1.0*alpha)/(dy*dy)));
				if(l<ymax)
					tripletList.push_back(T(posuy,posuy+xsize,(-2.0*alpha)/(dy*dy)));
				if(l<(ymax-1))
					tripletList.push_back(T(posuy,posuy+2*xsize,(1.0*alpha)/(dy*dy)));

				if(m>1)
					tripletList.push_back(T(posuy,posuy-zpitch,(-2.0*beta)/(dz*dz)));
				if(m>2)
					tripletList.push_back(T(posuy,posuy-2*zpitch,(1.0*beta)/(dz*dz)));
				if(m<zmax)
					tripletList.push_back(T(posuy,posuy+zpitch,(-2.0*beta)/(dz*dz)));
				if(m<(zmax-1))
					tripletList.push_back(T(posuy,posuy+2*zpitch,(1.0*beta)/(dz*dz)));

				if(k>1){
					if(l>1)
						tripletList.push_back(T(posuy,posuy-fieldsize-1-xsize,(+4.0*(beta+gama))/(4.0*dy*dx)));
					if(l<ymax)
						tripletList.push_back(T(posuy,posuy-fieldsize-1+xsize,(-4.0*(beta+gama))/(4.0*dy*dx)));
				}
				if(k<xmax){
					if(l>1)
						tripletList.push_back(T(posuy,posuy-fieldsize+1-xsize,(-4.0*(beta+gama))/(4.0*dy*dx)));
					if(l<ymax)
						tripletList.push_back(T(posuy,posuy-fieldsize+1+xsize,(+4.0*(beta+gama))/(4.0*dy*dx)));
				}


				if(l>1){
					if(m>1)
						tripletList.push_back(T(posuy,posuy+fieldsize-xsize-zpitch,(+4.0*(beta+gama))/(4.0*dz*dy)));
					if(m<zmax)
						tripletList.push_back(T(posuy,posuy+fieldsize-xsize+zpitch,(-4.0*(beta+gama))/(4.0*dz*dy)));
				}
				if(l<ymax){
					if(m>1)
						tripletList.push_back(T(posuy,posuy+fieldsize+xsize-zpitch,(-4.0*(beta+gama))/(4.0*dz*dy)));
					if(m<zmax)
						tripletList.push_back(T(posuy,posuy+fieldsize+xsize+zpitch,(+4.0*(beta+gama))/(4.0*dz*dy)));
				}

				/*uz*/
				long int posuz = pos+2*fieldsize;
				tripletList.push_back(T(posuz,posuz,4*alpha/(dz*dz)+2*beta/(dy*dy)+2*beta/(dx*dx)));
				if(k>1)
					tripletList.push_back(T(posuz,posuz-1,(-2.0*beta)/(dx*dx)));
				if(k>2)
					tripletList.push_back(T(posuz,posuz-2,(1.0*beta)/(dx*dx)));
				if(k<xmax)
					tripletList.push_back(T(posuz,posuz+1,(-2.0*beta)/(dx*dx)));
				if(k<(xmax-1))
					tripletList.push_back(T(posuz,posuz+2,(1.0*beta)/(dx*dx)));

				if(l>1)
					tripletList.push_back(T(posuz,posuz-xsize,(-2.0*beta)/(dy*dy)));
				if(l>2)
					tripletList.push_back(T(posuz,posuz-2*xsize,(1.0*beta)/(dy*dy)));
				if(l<ymax)
					tripletList.push_back(T(posuz,posuz+xsize,(-2.0*beta)/(dy*dy)));
				if(l<(ymax-1))
					tripletList.push_back(T(posuz,posuz+2*xsize,(1.0*beta)/(dy*dy)));

				if(m>1)
					tripletList.push_back(T(posuz,posuz-zpitch,(-4.0*alpha)/(dz*dz)));
				if(m>2)
					tripletList.push_back(T(posuz,posuz-2*zpitch,(2.0*alpha)/(dz*dz)));
				if(m<zmax)
					tripletList.push_back(T(posuz,posuz+zpitch,(-4.0*alpha)/(dz*dz)));
				if(m<(zmax-1))
					tripletList.push_back(T(posuz,posuz+2*zpitch,(2.0*alpha)/(dz*dz)));

				if(m>1){
					if(l>1)
						tripletList.push_back(T(posuz,posuz-2*fieldsize-zpitch-xsize,(+4.0*(beta+gama))/(4.0*dy*dz)));
					if(l<ymax)
						tripletList.push_back(T(posuz,posuz-2*fieldsize-zpitch+xsize,(-4.0*(beta+gama))/(4.0*dy*dz)));
				}
				if(m<zmax){
					if(l>1)
						tripletList.push_back(T(posuz,posuz-2*fieldsize+zpitch-xsize,(-4.0*(beta+gama))/(4.0*dy*dz)));
					if(l<ymax)
						tripletList.push_back(T(posuz,posuz-2*fieldsize+zpitch+xsize,(+4.0*(beta+gama))/(4.0*dy*dz)));
				}


				if(k>1){
					if(m>1)
						tripletList.push_back(T(posuz,posuz-fieldsize-1-zpitch,(+4.0*(beta+gama))/(4.0*dz*dx)));
					if(m<zmax)
						tripletList.push_back(T(posuz,posuz-fieldsize-1+zpitch,(-4.0*(beta+gama))/(4.0*dz*dx)));
				}
				if(k<xmax){
					if(m>1)
						tripletList.push_back(T(posuz,posuz-fieldsize+1-zpitch,(-4.0*(beta+gama))/(4.0*dz*dx)));
					if(m<zmax)
						tripletList.push_back(T(posuz,posuz-fieldsize+1+zpitch,(+4.0*(beta+gama))/(4.0*dz*dx)));
				}


			}
		}
	}

/*
	Eigen::VectorXd b(3*fieldsize);



	for(long int m=0; m< zsize;++m){
		for(long int l=0; l< ysize;++l){
			for(long int k=0; k< xsize;++k){
				TYPE alpha=C11(k,l,m);
				TYPE gama=C12(k,l,m);
				b(k+l*xsize+m*zpitch)=4.0*(alpha+gama)*epsdx(k,l,m);
				b(k+l*xsize+m*zpitch+fieldsize)=4.0*(alpha+gama)*epsdy(k,l,m);
				b(k+l*xsize+m*zpitch+2*fieldsize)=4.0*(alpha+gama)*epsdz(k,l,m);
			}
		}
	}
*/
	Eigen::SparseMatrix<TYPE,Eigen::RowMajor> Matrix(3*fieldsize,3*fieldsize);
	Matrix.setFromTriplets(tripletList.begin(), tripletList.end());

	/*
	for (int k=0; k<Matrix.outerSize(); ++k)
	for (Eigen::SparseMatrix<double>::InnerIterator it(Matrix,k); it; ++it)
	{
	std::cout << "value:" << it.value() << std::endl;
	std::cout << "row:" << it.row() << std::endl; // row index
	std::cout << "col:" << it.col() << std::endl; // col index (here it is equal to k)
	std::cout << "index:" << it.index() << std::endl; // inner index, here it is equal to it.row()
	if(it.value()>0)
		std::getchar();
	}
*/
	//Eigen::SparseLU<Eigen::SparseMatrix<TYPE, Eigen::ColMajor>, Eigen::COLAMDOrdering<int>> solver;
	//Eigen::SuperLU<Eigen::SparseMatrix<TYPE, Eigen::ColMajor>> solver;
	//Eigen::UmfPackLU<Eigen::SparseMatrix<TYPE, Eigen::ColMajor>> solver;
	//Eigen::PastixLU<Eigen::SparseMatrix<TYPE, Eigen::ColMajor>> solver;

	// Compute the ordering permutation vector from the structural pattern of A
	//solver.compute(Matrix);
	//solver.analyzePattern(Matrix);
	// Compute the numerical factorization
	//solver.factorize(Matrix);

	//Eigen::SparseLU<Eigen::SparseMatrix<TYPE, Eigen::ColMajor>> solver;
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<TYPE, Eigen::ColMajor>,Eigen::Upper,Eigen::DiagonalPreconditioner<TYPE>> solver;
	//solver.compute(Matrix);

	//Eigen::BiCGSTAB<Eigen::SparseMatrix<TYPE, Eigen::ColMajor>> solver;
	//solver.compute(Matrix);
	//Use the factors to solve the linear system

	//SparseQR<Eigen::SparseMatrix<TYPE, Eigen::ColMajor>, COLAMDOrdering<int>> qrsolver;
	//qrsolver.compute(Matrix);
/*

	for(long int m=0; m< zsize;++m){
		for(long int l=0; l< ysize;++l){
			for(long int k=0; k< xsize;++k){
				x(k+l*xsize+m*zpitch)=(lattice(k,l,m)-substrate->latticeparam())/substrate->latticeparam();
				x(k+l*xsize+m*zpitch+fieldsize)=(lattice(k,l,m)-substrate->latticeparam())/substrate->latticeparam();
				x(k+l*xsize+m*zpitch+2*fieldsize)=(lattice(k,l,m)-substrate->latticeparam())/substrate->latticeparam();
			}
		}
	}


	x = solver.solveWithGuess(b,x);*/
	//Eigen::VectorXd x(3*fieldsize);
	//x = solver.solve(b);
	/*std::cout << "#iterations: " << solver.iterations() << std::endl;
	std::cout << "estimated error: " << solver.error() << std::endl;*/

/*
	viennacl::compressed_matrix<TYPE>  vcl_sparse_matrix(3*fieldsize,3*fieldsize);
	viennacl::vector<TYPE> vcl_b(3*fieldsize);
	viennacl::vector<TYPE> vcl_rec(3*fieldsize);
	viennacl::copy(Matrix, vcl_sparse_matrix);
	viennacl::copy(b, vcl_b);
	Eigen::VectorXd x(3*fieldsize);
	//viennacl::linalg::ilu0_tag ilu0_config;
	//viennacl::linalg::ilu0_precond< viennacl::compressed_matrix<TYPE> > vcl_ilut(vcl_sparse_matrix, ilu0_config);
	vcl_rec = viennacl::linalg::solve(vcl_sparse_matrix, vcl_b, viennacl::linalg::cg_tag());
	viennacl::copy(vcl_rec,x);
	*/
	/*Eigen::VectorXd x(3*fieldsize);
	viennacl::linalg::gmres_tag solvetag(1e-5,5000,200);
	x = viennacl::linalg::solve(Matrix,b, solvetag);
	cout << "iteractions: " << solvetag.iters() << std::endl;
	 */

	Eigen::saveMarket(Matrix,"matrix.temp.mtx");

	paralution::init_paralution();

	paralution::LocalMatrix<TYPE> mat;
	mat.AllocateCSR("csr_matrix",21*fieldsize,3*fieldsize,3*fieldsize);
	paralution::LocalVector<TYPE> vecb;
	vecb.Allocate("rhs_vector",3*fieldsize);
	paralution::LocalVector<TYPE> vecx;
	vecx.Allocate("x_vector",3*fieldsize);

	for(long int m=0; m< zsize;++m){
		for(long int l=0; l< ysize;++l){
			for(long int k=0; k< xsize;++k){
				vecx[k+l*xsize+m*zpitch]=eps(k,l,m);
				vecx[k+l*xsize+m*zpitch+fieldsize]=eps(k,l,m);
				vecx[k+l*xsize+m*zpitch+2*fieldsize]=eps(k,l,m);
			}
		}
	}

	for(long int m=0; m< zsize;++m){
		for(long int l=0; l< ysize;++l){
			for(long int k=0; k< xsize;++k){
				TYPE alpha=C11(k,l,m);
				TYPE gama=C12(k,l,m);
				vecb[k+l*xsize+m*zpitch]=4.0*(alpha+gama)*epsdx(k,l,m);
				vecb[k+l*xsize+m*zpitch+fieldsize]=4.0*(alpha+gama)*epsdy(k,l,m);
				vecb[k+l*xsize+m*zpitch+2*fieldsize]=4.0*(alpha+gama)*epsdz(k,l,m);
			}
		}
	}



	mat.ReadFileMTX("matrix.temp.mtx");

	paralution::GMRES<paralution::LocalMatrix<TYPE>, paralution::LocalVector<TYPE>, TYPE > ls ;
	//paralution::ILU<paralution::LocalMatrix<TYPE>, paralution::LocalVector<TYPE>, TYPE > p;
	ls.Init(1e-10,1e-8,1e8,5000);
	ls.SetBasisSize(30);
	ls.SetOperator(mat);
	//p.SetRelaxation(1.3);
	//ls.SetPreconditioner(p);


	mat.MoveToAccelerator();
	vecb.MoveToAccelerator();
	vecx.MoveToAccelerator();


	ls.Build();
	ls.Solve(vecb,&vecx);

	vecx.MoveToHost();

	Eigen::VectorXd x(3*fieldsize);
	for(long int i =0; i<3*fieldsize;++i)
		x[i]=vecx[i];

	std::cout << "Solver status=" << ls.GetSolverStatus() << std::endl;

	mat.Clear();
	vecb.Clear();
	vecx.Clear();

	paralution::stop_paralution();




	for(long int m=0; m< zsize;++m){
		for(long int l=0; l< ysize;++l){
			for(long int k=0; k< xsize;++k){
				ux(k,l,m)=x(k+l*xsize+m*zpitch);
				uy(k,l,m)=x(k+l*xsize+m*zpitch+fieldsize);
				uz(k,l,m)=x(k+l*xsize+m*zpitch+2*fieldsize);
			}
		}
	}



};


} /* namespace epital */

#endif /* STRAINCALCULATOR_CPP_ */
