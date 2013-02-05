/*
 * PilotedMultiplePrnCorrelator.cpp
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

     PilotedMultiplePrnCorrelator

     Processes the correlation of multiple possible PRNs using one (or two) pilot PRN(s)
     It works with a pilot correlator and a message correlator referenced from outside.
     - The pilot correlator looks for a correlation on pilot PRN(s) only and produces the information
       on the best possible frequency shift from zero IF and time delta of the PRN origin in the receiver time reference
     - The message correlator uses this information to look for the one among possible PRNs (the PRN "alphabet") that has a distinctive
       peak at the frequency shift and time delta.
 */

#include "PilotedMultiplePrnCorrelator.h"
#include "PilotCorrelationAnalyzer.h"
#include "PilotCorrelator.h"
#include "PilotedMessageCorrelator.h"
#include "PrnAutocorrelator.h"
#include "WsgcUtils.h"
#include <assert.h>
#include <iostream>


PilotedMultiplePrnCorrelator::PilotedMultiplePrnCorrelator(
        PilotCorrelationAnalyzer& pilot_correlation_analyzer,
		std::vector<CorrelationRecord>& correlation_records,
		PilotCorrelator& pilot_correlator, 
        PilotedMessageCorrelator& message_correlator,
        PrnAutocorrelator& prn_autocorrelator) :
    _pilot_correlation_analyzer(pilot_correlation_analyzer),
    _correlation_records(correlation_records),
    _pilot_correlator(pilot_correlator),
    _message_correlator(message_correlator),
    _prn_autocorrelator(prn_autocorrelator)
{}


PilotedMultiplePrnCorrelator::~PilotedMultiplePrnCorrelator()
{}


void PilotedMultiplePrnCorrelator::set_source_block(wsgc_fftw_complex *fftw_source_block)
{
	wsgc_complex *source_block = reinterpret_cast<wsgc_complex *>(fftw_source_block);
    _pilot_correlation_analyzer.store_prn_samples(source_block);
}


void PilotedMultiplePrnCorrelator::make_correlation(unsigned int pilot_prn_code_index)
{
    unsigned int prn_index = _pilot_correlation_analyzer.get_prn_index();

    // autocorrelation - removed not useful for piloted operation
    /*
    _prn_autocorrelator.set_source_block(_pilot_correlation_analyzer.get_last_samples(), prn_index);
    _prn_autocorrelator.make_correlation();
    */

    // pilot correlation
    _pilot_correlator.execute(_pilot_correlation_analyzer, pilot_prn_code_index);

    // message correlation
    if (_pilot_correlator.new_batch_processed())
    {
        for (int pi=_pilot_correlator.get_batch_size()-1; pi >= 0; pi--)
        {
            // analyze pilot correlation records with the analyzer and process a symbols batch when ready
            if (_pilot_correlation_analyzer.validate_pilot_correlation_records_back(pi))
            {
                unsigned int preferred_time_shift_start;
                unsigned int preferred_time_shift_length;
                
                if (!_pilot_correlation_analyzer.analyze(preferred_time_shift_start, preferred_time_shift_length))
                {
                    std::cout << "Preferred PRN time shift could not be identified - Analysis window is invalidated" << std::endl;
                }
                
                _pilot_correlation_analyzer._message_times.push_back(PilotCorrelationAnalyzer::tmp_time);
                clock_gettime(PilotCorrelationAnalyzer::_time_option, &_pilot_correlation_analyzer._message_times.back()._beg);

                _message_correlator.execute(_pilot_correlation_analyzer);

                clock_gettime(PilotCorrelationAnalyzer::_time_option, &_pilot_correlation_analyzer._message_times.back()._end);

                _pilot_correlation_analyzer.post_process_noise(); // average all valid noise samples in the analysis window
                _pilot_correlation_analyzer.reset_analysis();
            }
        }
    }
}


