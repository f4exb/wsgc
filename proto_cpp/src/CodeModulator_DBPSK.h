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
      
     CodeModulator for DBPSK signals
      
*/
#ifndef __CODE_MODULATOR_DBPSK_H__
#define __CODE_MODULATOR_DBPSK_H__

#include "CodeModulator.h"

/**
 * \brief Code modulator for Binary Phase Shift Keying (BPSK)
 */
class CodeModulator_DBPSK : public CodeModulator
{
    public:
        CodeModulator_DBPSK() : _seed_bit(1) {}
        virtual ~CodeModulator_DBPSK() {}

        /**
         * Fills samples in the given array according to given {0,1} code bit samples
         * \param fftw_code_in pointer to first element in the array. Array must be allocated with enough memory to fit all samples
         * \param code_samples samples of code bits taking binary {0,1} values
         */
        virtual void fill_code_samples(wsgc_fftw_complex *fftw_code_in, std::vector<char>& code_samples);

        /**
         * Modulates a waveform according to given {0,1} code bit samples
         * \param in pointer to first element in the input waveform array. It must be at least of the size of the code_bits vector.
         * \param out pointer to the first element in the output waveform array. It must be allocated with at least of the size of the code_bits vector.
         * \param code_samples samples of code bits taking binary {0,1} values
         */
        virtual void modulate(const wsgc_fftw_complex *in, wsgc_fftw_complex *out, std::vector<char>& code_samples);

        /**
         * Fills samples in the given array according to given {0,1} code bits
         * \param fftw_code_in pointer to first element in the array. Array must be allocated with enough memory to fit all samples
         * \param code_bits code bits taking binary {0,1} values
         * \param f_sampling Sampling frequency
         * \param f_chip Chip frequency
         */
        virtual void fill_code_samples(wsgc_fftw_complex *fftw_code_in, std::vector<char>& code_bits, wsgc_float f_sampling, wsgc_float f_chip);

        /**
         * Modulates a waveform according to given {0,1} code bits
         * \param in pointer to first element in the input waveform array. It must be at least of the size of the code_bits vector.
         * \param out pointer to the first element in the output waveform array. It must be allocated with at least of the size of the code_bits vector.
         * \param code_bits code bits taking binary {0,1} values
         * \param f_sampling Sampling frequency
         * \param f_chip Chip frequency
         */
        virtual void modulate(const wsgc_fftw_complex *in, wsgc_fftw_complex *out, std::vector<char>& code_bits, wsgc_float f_sampling, wsgc_float f_chip);

        /**
         * Set the seed bit
         * \param seed_bit Reference code bit value when starting
         */
        void set_seed_bit(char seed_bit)
        {
        	_seed_bit = seed_bit;
        }


    protected:
        char _seed_bit; //!< Reference code bit value when starting
};

#endif // __CODE_MODULATOR_BPSK_H__
