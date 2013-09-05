/*
     Copyright 2012 Edouard Griffiths <f4exb at free dot fr>
 
     This file is part of WSGC. A Weak Signal transmission mode using Gold Codes
 
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
      
     LocalCodes - CUDA implementation
      
     Creates a local copy of the complex conjugate of the FFT transform of all symbols gold codes
     This pre-calculates the code to be inserted in the final IFFT.
     Used with BPSK complex signals.
     
*/
#include "LocalCodesFFT_Cuda.h"
#include "GoldCodeGenerator.h"
#include "CodeModulator.h"
#include "WsgcException.h"
#include "Cuda_Operators.h"
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <cuComplex.h>
#include <cutil_inline.h>    // includes cuda.h and cuda_runtime_api.h
#include <cufft.h>
#include <string.h>
#include <assert.h>

LocalCodesFFT_Cuda::LocalCodesFFT_Cuda(CodeModulator& code_modulator, 
            const GoldCodeGenerator& gc_generator,
            wsgc_float f_sampling, 
            wsgc_float f_chip,
            std::vector<unsigned int>& symbols,
            unsigned int cuda_device) :
    CudaDeviceManager::CudaDeviceManager(cuda_device),
    LocalCodesFFT(code_modulator, gc_generator, f_sampling, f_chip, symbols),
    _nb_codes(symbols.size()),
    _d_code(_nb_code_samples*symbols.size())
{
    _h_fft_code_in = new wsgc_complex[_nb_code_samples];
    cufftResult_t fft_stat = cufftPlan1d(&_fft_plan, _nb_code_samples, CUFFT_C2C, _nb_codes);

    if (fft_stat != CUFFT_SUCCESS)
    {
    	std::ostringstream os;
    	os << "CUFFT error: Failed to allocate FFT plan rc=" << fft_stat;
    	throw WsgcException(os.str());
    }

    fill_codes_matrix();
}


LocalCodesFFT_Cuda::~LocalCodesFFT_Cuda()
{
    cufftDestroy(_fft_plan);
    delete[] _h_fft_code_in;
}


void LocalCodesFFT_Cuda::fill_codes_matrix()
{
    std::vector<char> code;
    thrust::device_vector<cuComplex> _d_fft_code_in(_nb_codes*_nb_code_samples);
    
    // Fill input codes samples
    
    std::vector<unsigned int>::iterator prni_it = _symbols.begin();
    const std::vector<unsigned int>::iterator prni_end = _symbols.end();
    unsigned int i=0;
    
    for (; prni_it != prni_end; ++prni_it, i++)
    {
        assert(*prni_it < _gc_generator.get_nb_codes());

        _gc_generator.make_code(code, *prni_it); // 0/1 bits
        _code_modulator.fill_code_samples(reinterpret_cast<wsgc_fftw_complex *>(_h_fft_code_in), code, _f_sampling, _f_chip);   // This is the modulation specific part
        
        // copy to device vector of codes
        thrust::copy(
            reinterpret_cast<const cuComplex *>(_h_fft_code_in), 
            reinterpret_cast<const cuComplex *>(_h_fft_code_in+_nb_code_samples), 
            _d_fft_code_in.begin() + (i*_nb_codes)
        );

        index_symbol(i, *prni_it);
    }
    
    // Do the FFT
    
    if (cufftExecC2C(_fft_plan, thrust::raw_pointer_cast(&_d_fft_code_in[0]), thrust::raw_pointer_cast(&_d_fft_code_in[0]), CUFFT_FORWARD) != CUFFT_SUCCESS)
    {
        throw WsgcException("CUFFT Error: Failed to do FFT of local code");
    }
        
    // Do the conjugate
    
    thrust::transform(_d_fft_code_in.begin(), _d_fft_code_in.end(), _d_code.begin(), conj_functor());
}
        
