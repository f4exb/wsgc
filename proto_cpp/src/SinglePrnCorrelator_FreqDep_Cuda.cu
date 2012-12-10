/*
     Copyright 2012 Edouard Griffiths <f4exb at free dot fr>
 
     This file is part of WSGC. A Weak Signal transmission mode using Gold Codes.
 
     This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
 
     You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 51 Franklin Street, Boston, MA  02110-1301  USA

     Static not real time prototype in C++
      
     SinglePrnCorrelator_FreqDep
       
     Does a IFFT batch correlation and averaging given one single PRN to look for over a range of frequencies.
     
*/
#include "SinglePrnCorrelator_FreqDep_Cuda.h"
#include "GoldCodeGenerator.h"
#include "WsgcTypes.h"
#include "WsgcUtils.h"
#include "WsgcException.h"
#include "LocalCodesFFT_Cuda.h"
#include "PilotCorrelationAnalyzer.h"

#include "Cuda_Averagers.h"
#include "CudaManager.h"

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <cutil_inline.h>    // includes cuda.h and cuda_runtime_api.h

#include <iostream>


SinglePrnCorrelator_FreqDep_Cuda::SinglePrnCorrelator_FreqDep_Cuda(
		GoldCodeGenerator& gc_generator,
		CodeModulator& code_modulator,
		wsgc_float f_sampling,
		wsgc_float f_chip,
		std::vector<unsigned int>& pilot_symbols,
		unsigned int nb_f_bins,
		unsigned int prn_per_block,
		unsigned int nb_batch_prns,
		unsigned int freq_step_division,
		unsigned int cuda_device) :
    SinglePrnCorrelator_FreqDep::SinglePrnCorrelator_FreqDep(gc_generator.get_nb_code_samples(f_sampling,f_chip), nb_f_bins, prn_per_block, nb_batch_prns, freq_step_division),
    _local_codes(code_modulator, gc_generator, f_sampling, f_chip, pilot_symbols),
    _d_ifft_in(2*nb_batch_prns*_fft_N*freq_step_division*nb_f_bins),
    _d_ifft_out(2*nb_batch_prns*_fft_N*freq_step_division*nb_f_bins),
    _nb_pilot_prns(_local_codes.get_nb_codes()),
    _cuda_device(cuda_device),
    _pilot_correlation_analyzer(0)
{
	// allocate CuBLAS handle
	cublasStatus_t stat = cublasCreate(&_cublas_handle);
	if (stat != CUBLAS_STATUS_SUCCESS)
	{
		throw WsgcException("CUBLAS Error: Failed to initialize library");
	}

    // copy local code of PRN to device (_d_local_code)
	//const wsgc_complex *local_code = _local_codes.get_local_code(_prn);
	//_d_local_code.assign(reinterpret_cast<const cuComplex *>(local_code), reinterpret_cast<const cuComplex *>(local_code+_fft_N));

	// Allocate FFT plan

    _n[0] = _fft_N;
    _inembed[0] = _fft_N;
    _onembed[0] = _fft_N;

    cufftResult_t fft_stat = cufftPlanMany(&_ifft_plan, 1, _n,
		_inembed, 2*_nb_batch_prns, 2*_nb_batch_prns*_fft_N,
		_onembed, 2*_nb_batch_prns, 2*_nb_batch_prns*_fft_N,
		CUFFT_C2C, _nb_f_bins*_freq_step_division);

    if (fft_stat != CUFFT_SUCCESS)
    {
    	std::ostringstream err_os;
    	err_os << "CUFFT Error: Unable to create plan for pilot IFFT RC=" << fft_stat;
    	throw WsgcException(err_os.str());
    }

    // Select CUDA device
    cutilSafeCall(cudaSetDevice(_cuda_device));
}


SinglePrnCorrelator_FreqDep_Cuda::~SinglePrnCorrelator_FreqDep_Cuda()
{
    // IFFT plan
    cufftDestroy(_ifft_plan);

    //cuBLAS library handle
    cublasDestroy(_cublas_handle);
}
        

void SinglePrnCorrelator_FreqDep_Cuda::multiply_and_ifft(const thrust::device_vector<cuComplex>& d_source_block, unsigned int pilot_prn_index, unsigned int prn_position)
{
	const thrust::device_vector<cuComplex>& d_local_codes = _local_codes.get_local_codes();

	// Multiplication. Linear index ranging [0..T*Isf*Ihf[ making all FFT samples for all harmonic frequencies and sub-frequencies

    thrust::for_each(
        thrust::make_zip_iterator(
            thrust::make_tuple(
                thrust::make_permutation_iterator(d_source_block.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_A(_fft_N, _freq_step_division))),
                thrust::make_permutation_iterator(d_local_codes.begin() + pilot_prn_index*_fft_N, thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_B(_fft_N, _freq_step_division, _nb_f_bins))),
                thrust::make_permutation_iterator(_d_ifft_in.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_C(_nb_batch_prns, prn_position)))
            )
        ),
        thrust::make_zip_iterator(
            thrust::make_tuple(
				thrust::make_permutation_iterator(d_source_block.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(_fft_N*_freq_step_division*_nb_f_bins), transpose_index_A(_fft_N, _freq_step_division))),
				thrust::make_permutation_iterator(d_local_codes.begin() + (pilot_prn_index+1)*_fft_N, thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(_fft_N*_freq_step_division*_nb_f_bins), transpose_index_B(_fft_N, _freq_step_division, _nb_f_bins))),
				thrust::make_permutation_iterator(_d_ifft_in.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(_fft_N*_freq_step_division*_nb_f_bins), transpose_index_C(_nb_batch_prns, prn_position)))
            )
        ),
        cmulc_functor()
    );

    /*
    if (cudaThreadSynchronize() != cudaSuccess)
    {
        throw WsgcException("CUDA Multiplication Error: unable to synchronize threads");
    }
    */

    // IFFT

    if (cufftExecC2C(_ifft_plan,
    		         thrust::raw_pointer_cast(&_d_ifft_in[prn_position]),
    		         thrust::raw_pointer_cast(&_d_ifft_out[prn_position]),
    		         CUFFT_INVERSE) != CUFFT_SUCCESS)
    {
    	throw WsgcException("CUFFT Error: Failed to execute IFFT");
    }

    // Debug
    /*
    int cublas_max_index;

	cublasStatus_t stat = cublasIcamax(_cublas_handle, (_fft_N*_freq_step_division*_nb_f_bins),
			thrust::raw_pointer_cast(&_d_ifft_out[prn_position]),
			2*_nb_batch_prns, &cublas_max_index);

	if (stat != CUBLAS_STATUS_SUCCESS)
	{
		std::ostringstream err_os;
		err_os << "CUBLAS Error: cublasIcamax failed with RC=" << stat;
		std::cout << err_os.str() << std::endl;
		throw WsgcException(err_os.str());
	}

    cuComplex z = _d_ifft_out[prn_position + cublas_max_index*2*_nb_batch_prns];
    std::cout << prn_position << ": IFFT max: " << cublas_max_index << " : " << mag_algebraic_functor()(z) << std::endl;
    */
}


void SinglePrnCorrelator_FreqDep_Cuda::execute_averaging(bool first_half)
{
	assert(_pilot_correlation_analyzer != 0);

	AveragingDimensions_t ad;
	ad._B = _nb_batch_prns;
	ad._T = _fft_N;
	ad._Ifs = _freq_step_division;
	ad._Ifh = _nb_f_bins;
	unsigned int b_shift = (first_half ? 0 : _nb_batch_prns);

    if (cudaThreadSynchronize() != cudaSuccess)
    {
        throw WsgcException("CUDA Averaging Error: unable to synchronize threads before averaging");
    }

    // Sum average
    
    _pilot_correlation_analyzer->_pilot_cuda_avg_times.push_back(PilotCorrelationAnalyzer::tmp_time);
    clock_gettime(PilotCorrelationAnalyzer::_time_option, &_pilot_correlation_analyzer->_pilot_cuda_avg_times.back()._beg);

	if (_prn_per_block == 4) // supported size
	{
		Averager<4> averager;
		averager.run(_d_ifft_out, ad, b_shift);
	}
	else if (_prn_per_block == 6) // supported size
	{
		Averager<6> averager;
		averager.run(_d_ifft_out, ad, b_shift);
	}
	else if (_prn_per_block == 8) // supported size
	{
		Averager<8> averager;
		averager.run(_d_ifft_out, ad, b_shift);
	}
	else
	{
		throw WsgcException("CUDA Averaging Error: averaging not supported for this average size");
	}

    if (cudaThreadSynchronize() != cudaSuccess)
    {
        throw WsgcException("CUDA Averaging Error: unable to synchronize threads after averaging");
    }

    clock_gettime(PilotCorrelationAnalyzer::_time_option, &_pilot_correlation_analyzer->_pilot_cuda_avg_times.back()._end);

	// Reductions to find best match (max magnitude peak)
    
    _pilot_correlation_analyzer->_pilot_cuda_reduc_times.push_back(PilotCorrelationAnalyzer::tmp_time);
    clock_gettime(PilotCorrelationAnalyzer::_time_option, &_pilot_correlation_analyzer->_pilot_cuda_reduc_times.back()._beg);

    _batch_max_magnitudes.assign(_nb_batch_prns, 0.0);
    _batch_max_composite_indexes.assign(_nb_batch_prns, 0);
    
    thrust::tuple<cuComplex, int> init(_d_ifft_out[0], 0);
    thrust::tuple<cuComplex, int> max_tuple;

    for (unsigned int pi=0; pi < _nb_batch_prns; pi++)
    {
    	int cublas_max_index;
    	unsigned int shift = pi+(first_half ? 0 : ad._B);

    	cublasStatus_t stat = cublasIcamax(_cublas_handle, (ad._T*ad._Ifs*ad._Ifh),
    			thrust::raw_pointer_cast(&_d_ifft_out[shift]),
    			2*ad._B, &cublas_max_index);

    	if (stat != CUBLAS_STATUS_SUCCESS)
    	{
    		std::ostringstream err_os;
    		err_os << "CUBLAS Error: cublasIcamax failed with RC=" << stat;
    		std::cout << err_os.str() << std::endl;
    		throw WsgcException(err_os.str());
    	}

    	_batch_max_composite_indexes[pi] = cublas_max_index-1; // CUBLAS indexes are 1-based
    	/*
    	_batch_max_magnitudes[pi] = mag_squared_functor()(_d_ifft_out[]);
    	cuComplex z = _d_ifft_out[absolute_max_index+shift];
    	_batch_complex_values_max[pi].real() = z.x;
    	_batch_complex_values_max[pi].real() = z.y;
    	*/

    	/*
        max_tuple = thrust::reduce(    
            thrust::make_zip_iterator(
                thrust::make_tuple(
                    thrust::make_permutation_iterator(_d_ifft_out.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_avgC(ad._B, ad._T, ad._Ifs, ad._Ifh, 0, first_half))),
                    thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_avgC(ad._B, ad._T, ad._Ifs, ad._Ifh, 0, first_half))
                )
            ),
            thrust::make_zip_iterator(
                thrust::make_tuple(
                    thrust::make_permutation_iterator(_d_ifft_out.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(ad._T*ad._Ifs*ad._Ifh*ad._B), transpose_index_avgC(ad._B, ad._T, ad._Ifs, ad._Ifh, 0, first_half))),
                    thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(ad._T*ad._Ifs*ad._Ifh*ad._B), transpose_index_avgC(ad._B, ad._T, ad._Ifs, ad._Ifh, 0, first_half))
                )
            ),
            init,
            bigger_mag_tuple()
        );
        unsigned int absolute_max_index = thrust::get<1>(max_tuple);
        
        _batch_max_composite_indexes[pi] = absolute_max_index / ad._B;
        _batch_max_magnitudes[pi] = mag_squared_functor()(thrust::get<0>(max_tuple));
        _batch_complex_values_max[pi].real() = thrust::get<0>(max_tuple).x;
        _batch_complex_values_max[pi].imag() = thrust::get<0>(max_tuple).y;
        */
    }

    if (cudaThreadSynchronize() != cudaSuccess)
    {
        throw WsgcException("CUDA Reduction Error: unable to synchronize threads");
    }

    for (unsigned int pi=0; pi < _nb_batch_prns; pi++)
    {
    	/*
    	unsigned int t_max_index = _batch_max_composite_indexes[pi] % ad._T;
    	unsigned int f_max_index = _batch_max_composite_indexes[pi] / ad._T;
    	unsigned int fs_max_index = f_max_index % ad._Ifh;
    	unsigned int fh_max_index = f_max_index / ad._Ifh;
    	*/
    	unsigned int shift = pi+(first_half ? 0 : ad._B);

    	//cuComplex z = _d_ifft_out[t_max_index + fs_max_index*ad._T + fh_max_index *ad._T *ad._Ifs + shift];
    	cuComplex z = _d_ifft_out[shift + _batch_max_composite_indexes[pi]*2*ad._B];
    	_batch_max_magnitudes[pi] = mag_algebraic_functor()(z);
    	//std::cout << "_batch_max_magnitude=" << _batch_max_magnitudes[pi] << std::endl;
    	_batch_complex_values_max[pi].real() = z.x;
    	_batch_complex_values_max[pi].imag() = z.y;
    }

    clock_gettime(PilotCorrelationAnalyzer::_time_option, &_pilot_correlation_analyzer->_pilot_cuda_reduc_times.back()._end);
}
