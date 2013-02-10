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

     This flavour of correlator deals with PRNs without use of a pilot sequence

*/

#include "UnpilotedTrainingMessageCorrelator.h"
#include "TrainingCorrelationRecord.h"
#include <iostream>
#include <assert.h>


//=================================================================================================
UnpilotedTrainingMessageCorrelator::UnpilotedTrainingMessageCorrelator(
		wsgc_float f_sampling,
		wsgc_float f_chip,
		unsigned int prn_length,
        unsigned int sequence_length,
        unsigned int averaging_length,
		const std::vector<unsigned int>& prn_list,
		std::vector<TrainingCorrelationRecord>& training_correlation_records) :
		_f_sampling(f_sampling),
		_f_chip(f_chip),
		_prn_length(prn_length),
		_fft_N(int((f_sampling*prn_length)/f_chip)),
		_training_time_analyzer(_fft_N),
        _sequence_length(sequence_length),
        _averaging_length(averaging_length),
        _global_prn_index(0),
		_prn_list(prn_list),
		_samples_length(0),
		_prn_in_avg_count(0),
		_prn_in_seq_count(0),
		_training_correlation_records(training_correlation_records)
{
}


//=================================================================================================
UnpilotedTrainingMessageCorrelator::~UnpilotedTrainingMessageCorrelator()
{}


//=================================================================================================
void UnpilotedTrainingMessageCorrelator::dump_correlation_records(std::ostringstream& os, wsgc_float mag_factor)
{
	TrainingCorrelationRecord::dump_banner(os);
	std::vector<TrainingCorrelationRecord>::const_iterator it = _training_correlation_records.begin();
	const std::vector<TrainingCorrelationRecord>::const_iterator it_end = _training_correlation_records.end();

	for (; it != it_end; ++it)
	{
		it->dump_line(os, mag_factor);
	}
}


//=================================================================================================
void UnpilotedTrainingMessageCorrelator::dump_time_analyzer_results(std::ostringstream& os)
{
	_training_time_analyzer.dump_histo_time_shift_occurences(os);
}
