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
      
     PilotCorrelator
      
     Takes samples for the length of a PRN and processes it in order to:
     - process correlation over all frequency bins for one or two pilot code sequences
     - do peak estimation therefore delay and frequency tracking
     
*/
#include "PilotCorrelator.h"
#include "PilotCorrelationRecord.h"
#include "PilotCorrelationAnalyzer.h"
#include "GoldCodeGenerator.h"
#include "LocalCodesFFT.h"
#include "SinglePrnCorrelator_FreqDep.h"
#include "WsgcUtils.h"
#include <string.h>
#include <assert.h>
#include <iostream>


PilotCorrelator::PilotCorrelator(
			GoldCodeGenerator& gc_generator,
			CodeModulator& code_modulator,
			wsgc_float f_sampling,
			wsgc_float f_chip,
			std::vector<unsigned int>& pilot_symbols,
			unsigned int prn_per_symbol,
			unsigned int nb_pilot_f_bins,
			unsigned int nb_batch_prns,
			unsigned int freq_step_division) :
    _gc_generator(gc_generator),
    _code_modulator(code_modulator),
    _f_sampling(f_sampling),
    _f_chip(f_chip),
    _pilot_symbols(pilot_symbols),
    _prn_per_symbol(prn_per_symbol),
    _freq_step_division(freq_step_division),
    _fft_N(gc_generator.get_nb_code_samples(f_sampling,f_chip)),
    _nb_batch_prns(nb_batch_prns),
    _result_available(false),
    _new_batch_processed(false),
    _storage_depth(2*nb_batch_prns),
    _batch_index(0)
{
    // ensure PRN batch factor is at least the number of PRNs in an average minus one
    assert(nb_batch_prns >= _prn_per_symbol - 1);

    // ensure an odd number of frequency bins. The center frequency being at the center of the central interval
    if (nb_pilot_f_bins < 2)
    {
        _nb_pilot_f_bins = 3; // at least 3 intervals
    }
    else if (nb_pilot_f_bins < _fft_N)
    {
        _nb_pilot_f_bins = 2*(nb_pilot_f_bins/2)+1; // number if odd else next odd number
    }
    else
    {
    	_nb_pilot_f_bins = _fft_N - 1;
    }
}


PilotCorrelator::~PilotCorrelator()
{
}


void PilotCorrelator::update_pilot_correlation_records(
		PilotCorrelationAnalyzer& pilot_correlation_analyzer,
        unsigned int pilot_prn,
		bool even_batch,
		unsigned int pilot_prn_index,
		SinglePrnCorrelator_FreqDep *ifft_correlator_pilot)
{
	unsigned int correlation_index = 0;
    const std::vector<wsgc_float>& batch_magnitudes_max = ifft_correlator_pilot->get_batch_magnitudes_max();
    const std::vector<unsigned int>& batch_composite_indexes_max = ifft_correlator_pilot->get_batch_composite_indexes_max();
    const std::vector<wsgc_complex> batch_complex_values = ifft_correlator_pilot->get_batch_complex_values_max();
    unsigned int batch_shift = (even_batch ? _nb_batch_prns : 0); // even batch number targets second half

    for (unsigned int pi=0; pi < _nb_batch_prns; pi++, correlation_index++)
    {
    	PilotCorrelationRecord& pilot_correlation_record = pilot_correlation_analyzer.new_pilot_correlation_record();
        unsigned int ifft_peak_index = batch_composite_indexes_max[pi];
        unsigned int global_prn_index = pilot_prn_index - get_pipeline_length() + pi;
        unsigned int global_averaging_block_index = global_prn_index / _prn_per_symbol;
        unsigned int f_index, fs_index, t_index;
        ifft_correlator_pilot->calculate_indexes(ifft_peak_index, f_index, fs_index, t_index);
        //std::cout << ifft_peak_index << ":" << f_index << ":" << fs_index << ":" << t_index << std::endl;
        
        pilot_correlation_record.block_count = global_averaging_block_index;
        pilot_correlation_record.prn_index = global_prn_index;
        pilot_correlation_record.pilot_index = pilot_prn;
        pilot_correlation_record.magnitude_max = batch_magnitudes_max[pi] / _fft_N; // divide by FFT scaling factor
        pilot_correlation_record.phase_at_max = atan2(batch_complex_values[pi].imag(), batch_complex_values[pi].real());
        pilot_correlation_record.f_index_max = f_index * _freq_step_division + fs_index;
        pilot_correlation_record.t_index_max = t_index;
        // frequency shift to get back to zero IF = relative step displacement / step size. Step size is the size of a FFT bin.
        wsgc_float f_bin_size = _f_sampling / _fft_N;
        pilot_correlation_record.delta_f = ((((wsgc_float) pilot_correlation_record.f_index_max) / _freq_step_division) - (_nb_pilot_f_bins/2)) * f_bin_size;
    }
}


unsigned int PilotCorrelator::get_pipeline_length() const
{
    return _storage_depth - 1; // as soon as the last PRN in a two batch cycle (for _storage_depth length) is processed the first result is available hence the "-1".
}
