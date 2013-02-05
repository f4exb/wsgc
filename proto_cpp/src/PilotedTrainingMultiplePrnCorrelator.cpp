/*
 * PilotedTrainingMultiplePrnCorrelator.cpp
     Copyright 2012-2013 Edouard Griffiths <f4exb at free dot fr>

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

     PilotedTrainingMultiplePrnCorrelator

     Processes the correlation of training PRNs using the training identifying pilot PRN
     It works with a pilot correlator and a training message correlator referenced from outside.
     - The pilot correlator looks for a correlation on pilot PRN(s) only and produces the information
       on the best possible frequency shift from zero IF and time delta of the PRN origin in the receiver time reference
     - The training message correlator uses this information to make a correlation at the time delay and center frequency
       given by the pilot for all possible PRNs. It accumulates the correlation peaks of all PRNs with appropriate 
       shift so that the sequence of training PRNs would make the peak accumulation on the same PRN bin. Normally one 
       would use an incremental sequence of PRNs (by index order) and shift results one place to the right each time.
 */

#include "PilotedTrainingMultiplePrnCorrelator.h"
#include "PilotCorrelationAnalyzer.h"
#include "PilotCorrelator.h"
#include "PilotedTrainingMessageCorrelator.h"
#include "WsgcUtils.h"
#include <assert.h>
#include <iostream>


PilotedTrainingMultiplePrnCorrelator::PilotedTrainingMultiplePrnCorrelator(
        PilotCorrelationAnalyzer& pilot_correlation_analyzer,
		PilotCorrelator& pilot_correlator, 
        PilotedTrainingMessageCorrelator& training_message_correlator) :
    _pilot_correlation_analyzer(pilot_correlation_analyzer),
    _pilot_correlator(pilot_correlator),
    _training_message_correlator(training_message_correlator)
{}


PilotedTrainingMultiplePrnCorrelator::~PilotedTrainingMultiplePrnCorrelator()
{}


void PilotedTrainingMultiplePrnCorrelator::set_source_block(wsgc_fftw_complex *fftw_source_block)
{
	wsgc_complex *source_block = reinterpret_cast<wsgc_complex *>(fftw_source_block);
    _pilot_correlation_analyzer.store_prn_samples(source_block);
}


void PilotedTrainingMultiplePrnCorrelator::make_correlation(unsigned int pilot_prn_code_index)
{
    unsigned int prn_index = _pilot_correlation_analyzer.get_prn_index();

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
                else
                {
					_pilot_correlation_analyzer._message_times.push_back(PilotCorrelationAnalyzer::tmp_time);
					clock_gettime(PilotCorrelationAnalyzer::_time_option, &_pilot_correlation_analyzer._message_times.back()._beg);

					_training_message_correlator.execute(_pilot_correlation_analyzer);

					clock_gettime(PilotCorrelationAnalyzer::_time_option, &_pilot_correlation_analyzer._message_times.back()._end);

					std::cout << "Training sequence correlation results:" << std::endl;
					std::cout << " Ai max = " << _training_message_correlator._prnai_max
							<< ", Max = " <<  _training_message_correlator._mag_max
							<< ", Max to average max = " << _training_message_correlator._maxtoavg_max
							<< " PRNi max = " << _training_message_correlator._prni_max
							<< std::endl;
                }

                _pilot_correlation_analyzer.reset_analysis();
            }
        }
    }
}


