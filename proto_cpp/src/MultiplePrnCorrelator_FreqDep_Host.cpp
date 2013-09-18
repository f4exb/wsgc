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

     MultiplePrnCorrelator_FreqDep

     Does a IFFT batch correlation and averaging given a list of PRNs to look for over a range of frequencies.

*/
#include "MultiplePrnCorrelator_FreqDep_Host.h"
#include "GoldCodeGenerator.h"
#include "WsgcUtils.h"
#include "LocalCodesFFT.h"
#include <iostream>

MultiplePrnCorrelator_FreqDep_Host::MultiplePrnCorrelator_FreqDep_Host(
			const GoldCodeGenerator& gc_generator,
			CodeModulator& code_modulator,
			wsgc_float f_sampling,
			wsgc_float f_chip,
            unsigned int nb_message_symbols,
			std::vector<unsigned int>& message_symbols,
			unsigned int nb_f_bins,
			unsigned int prn_per_block,
			unsigned int nb_batch_prns,
			unsigned int freq_step_division) :
    MultiplePrnCorrelator_FreqDep::MultiplePrnCorrelator_FreqDep(gc_generator.get_nb_code_samples(f_sampling,f_chip), f_sampling, nb_message_symbols, nb_f_bins, prn_per_block, nb_batch_prns, freq_step_division),
    _source_fft(f_sampling, f_chip, _fft_N, freq_step_division),
	_local_codes(code_modulator, gc_generator, f_sampling, f_chip, message_symbols)
{
    _nb_codes = _local_codes.get_nb_codes();
    _batch_max_magnitudes.resize(_nb_codes);
    _batch_max_composite_indexes.resize(_nb_codes);
    _batch_complex_values_max.resize(_nb_codes);
    _ifft_codes_out.resize(_nb_codes);

    for (unsigned int prni=0; prni<_nb_codes; prni++)
    {
        // IFFT items
        _ifft_codes_out[prni] = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*_nb_f_bins*_freq_step_division*_storage_depth*sizeof(wsgc_fftw_complex));
        // Result storage
        _batch_max_magnitudes[prni].assign(_nb_batch_prns, 0.0);
        _batch_max_composite_indexes[prni].assign(_nb_batch_prns, 0);
        _batch_complex_values_max[prni].assign(_nb_batch_prns, cone);
    }

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


MultiplePrnCorrelator_FreqDep_Host::~MultiplePrnCorrelator_FreqDep_Host()
{
    // IFFT plan
    WSGC_FFTW_DESTROY_PLAN(_ifft_plan);

    // IFFT items
    for (unsigned int prni=0; prni<_nb_codes; prni++)
    {
        WSGC_FFTW_FREE(_ifft_codes_out[prni]);
    }
}


void MultiplePrnCorrelator_FreqDep_Host::multiply_and_ifft(const wsgc_complex *source_block, unsigned int prn_position)
{
    unsigned int fi;
    int fpi;

    const wsgc_complex *fft_ptr = _source_fft.get_fft_samples(source_block);

    for (unsigned int prn_i=0; prn_i<_nb_codes; prn_i++)
    {
        for (fi = 0, fpi = -(_nb_f_bins/2); fi< _nb_f_bins; fi++, fpi++) // for each frequency step
        {
            for (unsigned int fsi=0; fsi < _freq_step_division; fsi++)
            {
                // multiply source block by local code conjugate FFT
                for (unsigned int ffti = 0; ffti < _fft_N; ffti++) // multiply with local code
                {
                    int lci = (ffti + _fft_N + fpi) % _fft_N;
                    _ifft_code_in_tmp[ffti] = fft_ptr[fsi*_fft_N + ffti] * _local_codes.get_local_code_by_index(prn_i)[lci];
                    //_ifft_code_in_tmp[ffti] /= _fft_N; // pre-scaling
                }

                // do one IFFT
                WSGC_FFTW_EXECUTE(_ifft_plan);

                // push back the result in the global IFFT output array
                for (unsigned int ffti = 0; ffti < _fft_N; ffti++)
                {
                    _ifft_codes_out[prn_i][_fft_N*_storage_depth*_freq_step_division*fi + _fft_N*_storage_depth*fsi + ffti*_storage_depth + prn_position] = _ifft_code_out_tmp[ffti];
                }
            }
        }
    }
}


void MultiplePrnCorrelator_FreqDep_Host::execute_averaging(bool first_half)
{
    unsigned int start_position = (first_half ? 0 : _nb_batch_prns);

    for (unsigned int prn_i=0; prn_i<_nb_codes; prn_i++)
    {
        _batch_max_magnitudes[prn_i].assign(_nb_batch_prns, 0.0);
        _batch_max_composite_indexes[prn_i].assign(_nb_batch_prns, 0);

        for (unsigned int iffti = 0; iffti < _fft_N * _nb_f_bins * _freq_step_division; iffti++) // for each frequency step and each time bin
        {
            unsigned int pi, prn_position;

            for (pi = 0, prn_position = start_position; pi < _nb_batch_prns; pi++, prn_position++) // for all values in the batch
            {
                wsgc_float mag, sum_mag;

                // calculate average complex sum
                for (unsigned int ai = 1; ai < _prn_per_symbol; ai++) // sum with all next values in the batch within the averaging length (_prn_per_block)
                {
                    // sum with next values in the batch
                    _ifft_codes_out[prn_i][iffti*_storage_depth + prn_position] += _ifft_codes_out[prn_i][iffti*_storage_depth + ((prn_position+ai) % _storage_depth)];
                }

                // take magnitude of the complex sum
                WsgcUtils::magnitude_estimation(&_ifft_codes_out[prn_i][iffti*_storage_depth + prn_position], &sum_mag);

                if (sum_mag > _batch_max_magnitudes[prn_i][pi])
                {
                    _batch_max_magnitudes[prn_i][pi] = sum_mag;
                    _batch_max_composite_indexes[prn_i][pi] = iffti;
                    _batch_complex_values_max[prn_i][pi] = _ifft_codes_out[prn_i][iffti*_storage_depth + prn_position];
                }
            }
        }
    }

    _batch_index++;
}


const wsgc_complex *MultiplePrnCorrelator_FreqDep_Host::get_averaged_ifft_frames(unsigned int prn_i) const
{
    return _ifft_codes_out[prn_i];
}
