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
#include "SinglePrnCorrelator_FreqDep_Host.h"
#include "GoldCodeGenerator.h"
#include "WsgcUtils.h"
#include "LocalCodesFFT.h"
#include <iostream>

SinglePrnCorrelator_FreqDep_Host::SinglePrnCorrelator_FreqDep_Host(
			GoldCodeGenerator& gc_generator,
			CodeModulator& code_modulator,
			wsgc_float f_sampling,
			wsgc_float f_chip,
			std::vector<unsigned int>& pilot_symbols,
			unsigned int nb_f_bins,
			unsigned int prn_per_block,
			unsigned int nb_batch_prns,
			unsigned int freq_step_division) :
    SinglePrnCorrelator_FreqDep::SinglePrnCorrelator_FreqDep(gc_generator.get_nb_code_samples(f_sampling,f_chip), nb_f_bins, prn_per_block, nb_batch_prns, freq_step_division),
	_local_codes(code_modulator, gc_generator, f_sampling, f_chip, pilot_symbols)
{
    // IFFT items
    _ifft_code_out = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*_nb_f_bins*_freq_step_division*_storage_depth*sizeof(wsgc_fftw_complex));

    _ifft_code_in_tmp = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_fftw_complex));
    _ifft_code_out_tmp = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_fftw_complex));

    // IFFT organization
    // For IFFT elements place elements of successive batches next to each other so as to optimize average (sum) calculation later => as much as possible in the cache
    // <-- fpi ---------------><- fpi+1 --- ...   Elements for fpi then elements for fpi+1 etc... (nb_pilot_f_bins frequencies)
    // [......|......| ...    ][......|...
    //    ||      /\            |\        IFFT element j for batch #0 then IFFT element j for batch #1 etc... (_nb_batch_prns=3 => _storage_depth = 6)
    //    ||     /  \          |ti+1,0,0|...
    //    /\    /    \                     
    //   /  \   |ti,1,0|ti,1,1|...
    //  /    \
    // [ti,0,0|ti,0,1|ti,0,2||ti,0,3|ti,0,4|ti,0,5|
    // <- even -------------><- odd --------------><- even ...
    // 
    // ti,j,k with i = FFT sample index; j = batch index; k = frequency index
    // ti,j,k with i = frequency index; j = FFT sample index; k = batch index
    //

    // Do the IFFT in temporary buffer
    _ifft_plan = WSGC_FFTW_PLAN(_fft_N,
                                reinterpret_cast<wsgc_fftw_complex *>(_ifft_code_in_tmp),
                                reinterpret_cast<wsgc_fftw_complex *>(_ifft_code_out_tmp),
                                FFTW_BACKWARD, FFTW_ESTIMATE);
}


SinglePrnCorrelator_FreqDep_Host::~SinglePrnCorrelator_FreqDep_Host()
{
    // IFFT plan
    WSGC_FFTW_DESTROY_PLAN(_ifft_plan);
    
    // IFFT items
    //WSGC_FFTW_FREE(_ifft_code_in);
    WSGC_FFTW_FREE(_ifft_code_out);
}
        

void SinglePrnCorrelator_FreqDep_Host::multiply_and_ifft(const wsgc_complex *source_block, unsigned int pilot_prn_index, unsigned int prn_position)
{
    unsigned int fi;
    int fpi;

    for (fi = 0, fpi = -(_nb_f_bins/2); fi< _nb_f_bins; fi++, fpi++) // for each frequency step
    {
    	for (unsigned int fsi=0; fsi < _freq_step_division; fsi++)
    	{
			// multiply source block by local code conjugate FFT
			for (unsigned int ffti = 0; ffti < _fft_N; ffti++) // multiply with local code
			{
				int lci = (ffti + _fft_N + fpi) % _fft_N;
				_ifft_code_in_tmp[ffti] = source_block[fsi*_fft_N + ffti] * _local_codes.get_local_code(pilot_prn_index)[lci];
				//_ifft_code_in_tmp[ffti] /= _fft_N; // pre-scaling
			}

			// do one IFFT
			WSGC_FFTW_EXECUTE(_ifft_plan);

			// push back the result in the global IFFT output array
			for (unsigned int ffti = 0; ffti < _fft_N; ffti++)
			{
				_ifft_code_out[_fft_N*_storage_depth*_freq_step_division*fi + _fft_N*_storage_depth*fsi + ffti*_storage_depth + prn_position] = _ifft_code_out_tmp[ffti];
			}
    	}
    }
}


void SinglePrnCorrelator_FreqDep_Host::execute_averaging(bool first_half)
{
    unsigned int start_position = (first_half ? 0 : _nb_batch_prns);
    _batch_max_magnitudes.assign(_nb_batch_prns, 0.0);
    _batch_max_composite_indexes.assign(_nb_batch_prns, 0);

	for (unsigned int iffti = 0; iffti < _fft_N * _nb_f_bins * _freq_step_division; iffti++) // for each frequency step and each time bin
	{
		unsigned int pi, prn_position;

		for (pi = 0, prn_position = start_position; pi < _nb_batch_prns; pi++, prn_position++) // for all values in the batch
		{
			wsgc_float mag, sum_mag;

            /*
			// calculate average sum of modules
			WsgcUtils::magnitude_estimation(&_ifft_code_out[iffti*_storage_depth + prn_position], &sum_mag); // init with value at relative position 0 in the batch
			for (unsigned int ai = 1; ai < _prn_per_block; ai++) // sum with all next values in the batch within the averaging length (_prn_per_block)
			{
				// sum with next values in the batch
				WsgcUtils::magnitude_estimation(&_ifft_code_out[iffti*_storage_depth + ((prn_position+ai) % _storage_depth)], &mag);
				sum_mag += mag;
			}
            */
            
            // calculate average complex sum
			for (unsigned int ai = 1; ai < _prn_per_block; ai++) // sum with all next values in the batch within the averaging length (_prn_per_block)
			{
				// sum with next values in the batch
                _ifft_code_out[iffti*_storage_depth + prn_position] += _ifft_code_out[iffti*_storage_depth + ((prn_position+ai) % _storage_depth)];
            }

            // take magnitude of the complex sum
     		WsgcUtils::magnitude_estimation(&_ifft_code_out[iffti*_storage_depth + prn_position], &sum_mag);
     		//WsgcUtils::magnitude_algebraic(&_ifft_code_out[iffti*_storage_depth + prn_position], &sum_mag);

			if (sum_mag > _batch_max_magnitudes[pi])
			{
				_batch_max_magnitudes[pi] = sum_mag;
				_batch_max_composite_indexes[pi] = iffti;
				_batch_complex_values_max[pi] = _ifft_code_out[iffti*_storage_depth + prn_position];
			}
		}
    }
}


const wsgc_complex *SinglePrnCorrelator_FreqDep_Host::get_averaged_ifft_frames() const
{
    return _ifft_code_out;
}
