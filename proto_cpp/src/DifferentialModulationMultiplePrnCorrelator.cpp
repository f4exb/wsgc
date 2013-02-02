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

     DifferentialModulationMultiplePrnCorrelator

     This flavour of correlator deals with PRNs encoded with differential modulation

*/

#include "DifferentialModulationMultiplePrnCorrelator.h"
#include "CorrelationRecord.h"
#include <iostream>
#include <assert.h>

class CorrelationRecord;
class TrainingCorrelationRecord;

DifferentialModulationMultiplePrnCorrelator::DifferentialModulationMultiplePrnCorrelator(
		wsgc_float f_sampling,
		wsgc_float f_chip,
		unsigned int prn_length,
        unsigned int prn_per_symbol,
		const std::vector<unsigned int>& prn_list,
		unsigned int symbol_window_size,
		std::vector<CorrelationRecord>& correlation_records,
		std::vector<TrainingCorrelationRecord>& training_correlation_records) :
		_f_sampling(f_sampling),
		_f_chip(f_chip),
		_prn_length(prn_length),
		_fft_N(int((f_sampling*prn_length)/f_chip)),
        _prn_per_symbol(prn_per_symbol),
        _global_prn_index(0),
		_prn_list(prn_list),
		_samples_length(0),
		_prns_length(0),
		_symbol_window_size(symbol_window_size),
		_correlation_records(correlation_records),
		_training_correlation_records(training_correlation_records)
{
    init_results();
}


DifferentialModulationMultiplePrnCorrelator::~DifferentialModulationMultiplePrnCorrelator()
{}

void DifferentialModulationMultiplePrnCorrelator::init_results()
{
    _max_sy_iffti.assign(_symbol_window_size, 0);
    _max_sy_prni.assign(_symbol_window_size, 0);
    _max_sy_mags.assign(_symbol_window_size, 0.0);
    _max_sy_mags_prni.assign(_symbol_window_size*_prn_list.size(), 0.0);
}


void DifferentialModulationMultiplePrnCorrelator::dump_correlation_records(std::ostringstream& os, wsgc_float mag_factor)
{
	CorrelationRecord::dump_banner(os);
	std::vector<CorrelationRecord>::const_iterator it = _correlation_records.begin();
	const std::vector<CorrelationRecord>::const_iterator it_end = _correlation_records.end();

	for (; it != it_end; ++it)
	{
		it->dump_line(mag_factor, os);
	}
}

