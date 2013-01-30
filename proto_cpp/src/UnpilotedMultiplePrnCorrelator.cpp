/*
 * UnpilotedMultiplePrnCorrelator.h
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

     UnpilotedMultiplePrnCorrelator

     Processes the correlation of multiple possible PRNs using the frequency domain method
     - The message correlator attempts to correlate with all possible PRNs (the PRN "alphabet") to identify the one that has a distinctive
       peak. The index of the peak in the IFFT result identified the time delay of the start of the sent PRN sequence.

 */
#include "UnpilotedMultiplePrnCorrelator.h"
#include "CorrelationRecord.h"
#include "AutocorrelationRecord.h"

UnpilotedMultiplePrnCorrelator::UnpilotedMultiplePrnCorrelator(
    		std::vector<CorrelationRecord>& correlation_records,
            UnpilotedMessageCorrelator& message_correlator) :
    _correlation_records(correlation_records),
    _message_correlator(message_correlator),
    _global_prn_index(0)
{}


UnpilotedMultiplePrnCorrelator::~UnpilotedMultiplePrnCorrelator()
{}


void UnpilotedMultiplePrnCorrelator::set_source_block(wsgc_complex *source_block)
{
    _message_correlator.set_source_block(source_block);
    //_prn_autocorrelator.set_source_block(source_block, _global_prn_index); // This has proven useless
    _global_prn_index++;
}
    
    
void UnpilotedMultiplePrnCorrelator::make_correlation()
{
    if (_message_correlator.execute(_correlation_records)) 
    { // new results are available and put in the last elements of the correlation records vector
        // store shift maximums for process batch in shifts dictionnary
        for (unsigned int bi=_message_correlator.get_nb_batch_prns(); bi>0; bi--)
        {
            unsigned int shift_index_max = _correlation_records[_correlation_records.size() - bi].shift_index_max;
            std::map<unsigned int, unsigned int>::iterator shift_occurences_it;
            shift_occurences_it = _shift_occurences_dict.find(shift_index_max);
            
            if (shift_occurences_it == _shift_occurences_dict.end())
            {
                _shift_occurences_dict[shift_index_max] = 1;
            }
            else
            {
                _shift_occurences_dict[shift_index_max] += 1;
            }
        }
    }

    //_prn_autocorrelator.make_correlation(_autocorrelation_records); // This has proven useless
}
    

void UnpilotedMultiplePrnCorrelator::dump_message_correlation_records(wsgc_float mag_factor, std::ostringstream& os) const
{
    CorrelationRecord::dump_banner(os);
    
    std::vector<CorrelationRecord>::const_iterator corr_it = _correlation_records.begin();
    const std::vector<CorrelationRecord>::const_iterator corr_end = _correlation_records.end();
    
    for (;corr_it != corr_end; ++corr_it)
    {
        corr_it->dump_line(mag_factor, os);
    }
}
