/*
 * CUDAfunctions.hpp
 *
 *  Created on: Nov 27, 2014
 *      Author: marcel
 */

#include <complex>

#ifndef CUDAFUNCTIONS_HPP_
#define CUDAFUNCTIONS_HPP_

extern void groundlevel_periodicFFT_CUDA(long int timesteps, std::complex<double>* initial, std::complex<double>* potential,std::complex<double>* fftmulti,
		double effectivemass,long int dimx,long int dimy,long int dimz, double xstep, double ystep, double zstep);


#endif /* CUDAFUNCTIONS_HPP_ */
