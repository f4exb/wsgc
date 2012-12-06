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
     
*/
#ifndef __LOCAL_CODES_FFT_CUDA__
#define __LOCAL_CODES_FFT_CUDA__

#include "WsgcTypes.h"
#include "LocalCodesFFT.h"
#include <vector>
#include <cufft.h>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>

class CodeModulator;
class GoldCodeGenerator;

/**
 * \brief Local copy of codes to be used for frequency domain correlation
 *
   Creates a local copy of the complex conjugate of the FFT transform of all symbols gold codes
   This pre-calculates the code to be inserted in the final IFFT.
 *
 */
class LocalCodesFFT_Cuda : public LocalCodesFFT
{
    public:
        /**
        * \param code_modulator Modulator used to build the codes
        * \param gc_generator Gold Code generator used to build the codes
        * \param f_sampling Sampling frequency
        * \param f_chip Chip rate
        * \param symbols List of symbols to be processed        
        */
        LocalCodesFFT_Cuda(CodeModulator& code_modulator, 
            GoldCodeGenerator& gc_generator, 
            wsgc_float f_sampling, 
            wsgc_float f_chip,
            std::vector<unsigned int>& symbols);
            
        ~LocalCodesFFT_Cuda();
        
        /**
         * Get the vector of the local copy of the codes
         * \return Reference to the vector of the local copy of the codes
         */
        const thrust::device_vector<cuComplex>& get_local_codes() const
        {
            return _d_code;
        }

    protected:
        unsigned int _nb_codes;
        std::vector<thrust::device_vector<cuComplex> > _codes_matrix; //!< Matrix holding the local copy of the codes
        cufftHandle _fft_plan; //!< CUFFT transform plan for local code FFT computation
        wsgc_complex *_h_fft_code_in; //!< CUFFT input data. This is the modulated code
        thrust::device_vector<cuComplex> _d_code; //!< CUFFT input data. This is the modulated codes
        
        /**
         * Internal method to fill the matrix holding the local copy of the codes at construction time
         */
        void fill_codes_matrix();
};

#endif // __LOCAL_CODES_FFT_CUDA__
