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
      
     PilotCorrelator - CUDA version
      
     Takes samples for the length of a PRN and processes it in order to:
     - process correlation over all frequency bins for one or two pilot code sequences
     - do peak estimation therefore delay and frequency tracking
     
*/
#include "PilotCorrelator_Cuda.h"
#include "PilotCorrelationRecord.h"
#include "PilotCorrelationAnalyzer.h"
#include "GoldCodeGenerator.h"
#include "LocalCodesFFT_Cuda.h"
#include "SinglePrnCorrelator_FreqDep_Cuda.h"
#include "WsgcUtils.h"
#include <string.h>
#include <assert.h>
#include <iostream>


PilotCorrelator_Cuda::PilotCorrelator_Cuda(
			const GoldCodeGenerator& gc_generator,
			CodeModulator& code_modulator,
			wsgc_float f_sampling,
			wsgc_float f_chip,
			std::vector<unsigned int>& pilot_symbols,
			unsigned int cuda_device,
			unsigned int prn_per_symbol,
			unsigned int nb_pilot_f_bins,
			unsigned int nb_batch_prns,
			unsigned int freq_step_division) :
	PilotCorrelator(f_sampling, f_chip, gc_generator.get_nb_code_samples(f_sampling,f_chip), pilot_symbols, prn_per_symbol, nb_pilot_f_bins, nb_batch_prns, freq_step_division),
	_source_fft(f_sampling, f_chip, _fft_N, freq_step_division, cuda_device),
	_ifft_correlator_pilot(gc_generator, code_modulator, f_sampling, f_chip, pilot_symbols, nb_pilot_f_bins, cuda_device, prn_per_symbol, nb_batch_prns, freq_step_division)
{
}


PilotCorrelator_Cuda::~PilotCorrelator_Cuda()
{
}


/*
void PilotCorrelator_Host::execute(PilotCorrelationAnalyzer& pilot_correlation_analyzer)
{
    unsigned int pilot_prn_index = pilot_correlation_analyzer.get_prn_index();
    unsigned int batch_number = pilot_prn_index / _nb_batch_prns;

    execute(pilot_correlation_analyzer, _ifft_correlator_pilot_1);

    _new_batch_processed = false; // a new cycle begins

    if (_ifft_correlator_pilot_2 != 0)
    {
        execute(pilot_correlation_analyzer, _ifft_correlator_pilot_2);
    }
    
    if (((pilot_prn_index % _nb_batch_prns) == (_nb_batch_prns - 1)) && (batch_number > 0))
    {
        _result_available = true; // results are valid from now on
        _new_batch_processed = true; // a new batch has been processed
    }
}
*/

void PilotCorrelator_Cuda::execute(PilotCorrelationAnalyzer& pilot_correlation_analyzer, unsigned int pilot_prn_code_index)
{
	unsigned int prn = _pilot_symbols[pilot_prn_code_index];
    unsigned int pilot_prn_index = pilot_correlation_analyzer.get_prn_index();
    unsigned int batch_number = pilot_prn_index / _nb_batch_prns;
    wsgc_complex pilot1_max, pilot2_max;
    unsigned int prn_position = pilot_prn_index % _storage_depth;

    _ifft_correlator_pilot.set_pilot_correlation_analyzer(&pilot_correlation_analyzer);

	// Do FFT of input samples
	const thrust::device_vector<cuComplex>& fft_in = _source_fft.do_fft(pilot_correlation_analyzer.get_last_samples());

    pilot_correlation_analyzer._pilot_mul_ifft_times.push_back(PilotCorrelationAnalyzer::tmp_time);
    clock_gettime(PilotCorrelationAnalyzer::_time_option, &pilot_correlation_analyzer._pilot_mul_ifft_times.back()._beg);

    _ifft_correlator_pilot.multiply_and_ifft(fft_in, pilot_prn_code_index, prn_position);

    clock_gettime(PilotCorrelationAnalyzer::_time_option, &pilot_correlation_analyzer._pilot_mul_ifft_times.back()._end);

    // if a batch is ready calculate average. There is one batch delay before average can be calculated so results of batch #0 should be skipped
    if (((pilot_prn_index % _nb_batch_prns) == (_nb_batch_prns - 1)) && (batch_number > 0))
    {
        bool even_batch = (batch_number % 2 == 0);

        // Process pilot 1
        pilot_correlation_analyzer._pilot_avg_times.push_back(PilotCorrelationAnalyzer::tmp_time);
        clock_gettime(PilotCorrelationAnalyzer::_time_option, &pilot_correlation_analyzer._pilot_avg_times.back()._beg);

        _ifft_correlator_pilot.execute_averaging(!even_batch); // because of the one batch delay, the first half is processed on odd batches

        clock_gettime(PilotCorrelationAnalyzer::_time_option, &pilot_correlation_analyzer._pilot_avg_times.back()._end);

        update_pilot_correlation_records(pilot_correlation_analyzer, prn, even_batch, pilot_prn_index, &_ifft_correlator_pilot);
    }

    // update process indicators

    _new_batch_processed = false; // a new cycle begins

    if (((pilot_prn_index % _nb_batch_prns) == (_nb_batch_prns - 1)) && (batch_number > 0))
    {
        _result_available = true; // results are valid from now on
        _new_batch_processed = true; // a new batch has been processed
    }
}
