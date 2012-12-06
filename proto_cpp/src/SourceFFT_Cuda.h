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
      
     Source FFT - CUDA version
      
     Multiplies input samples by each of the sub-frequency LOs
     Do the FFT of the result
     
*/
#ifndef __SOURCE_FFT_CUDA__
#define __SOURCE_FFT_CUDA__

#include "WsgcTypes.h"
#include "SourceFFT.h"
#include <vector>
#include <thrust/device_vector.h>
#include <cuComplex.h>
#include <cutil_inline.h>    // includes cuda.h and cuda_runtime_api.h
#include <cufft.h>

class CodeModulator;
class GoldCodeGenerator;
class ContinuousPhaseCarrier;

/**
 * \brief Local copy of codes to be used for frequency domain correlation
 *
   Creates a local copy of the complex conjugate of the FFT transform of all symbols gold codes
   This pre-calculates the code to be inserted in the final IFFT.
 *
 */
class SourceFFT_Cuda : public SourceFFT
{
public:
    /**
    * \param f_sampling Sampling frequency
    * \param f_chip Chip rate
    * \param fft_N FFT size
    * \param freq_step_division Number of frequency sub steps
    */
    SourceFFT_Cuda(wsgc_float f_sampling, 
                   wsgc_float f_chip,
                   unsigned int fft_N,
                   unsigned int freq_step_division);

    virtual ~SourceFFT_Cuda();

    /**
     * Do the FFT
     */
    const thrust::device_vector<cuComplex>& do_fft(wsgc_complex *source_block);

    /**
     * Get the vector of source sub-frequency multiplied FFT samples
     * \return Reference to the vector of source sub-frequency multiplied FFT samples
     */
    const thrust::device_vector<cuComplex>& get_last_fft_samples() const
    {
        return _d_fft_out;
    }
    

protected:
    std::vector<ContinuousPhaseCarrier*> _local_oscillators; //!< LOs for sub-frequency steps
    cufftHandle _fft_plan;                       //!< FFTW plan for signal FFT
    thrust::device_vector<cuComplex> _d_fft_in;  //!< Signal input samples for FFT
    thrust::device_vector<cuComplex> _d_fft_out; //!< Signal FFT output samples
    wsgc_complex *_h_fft_sample_in;              //!< Host input samples

};

#endif // __SOURCE_FFT_CUDA__
