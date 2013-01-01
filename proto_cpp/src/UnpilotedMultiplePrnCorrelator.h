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

#ifndef __UNPILOTED_MULTIPLE_PRN_CORRELATOR_H__
#define __UNPILOTED_MULTIPLE_PRN_CORRELATOR_H__

#include "WsgcTypes.h"
#include "UnpilotedMessageCorrelator.h"
#include "PrnAutocorrelator.h"
#include <vector>
#include <map>
#include <sstream>

class CorrelationRecord;
class AutocorrelationRecord;

/**
 * \brief Correlator engine to get the correlation estimation of PRNs sent in the message
 *
 *  It uses an "unpiloted" message correlator.  The message correlator attempts to correlate with all possible PRNs (the PRN "alphabet") 
 *  to identify the one that has a distinctive peak. The index of the peak in the IFFT result identified the time delay of the start of
 *  the sent PRN sequence.
 *
 *  Supports pipelined processing in order to better accomodate the CUDA implementation.
 *  - Pipeline length is set by the Pilot Correlator "PRN batch factor" which can be no less than the averaging length (in PRN numbers) minus one. The
 *    pipeline length is calculated by the Pilot Correlator and is available through a getter method
 *  - is_correlation_record_available method tells if a correlation record is available after execution
 *  - is_correlation_record_available_next method tells if a correlation record will be available next. i.e. is false if it is the last correlation record.
 *
 */
class UnpilotedMultiplePrnCorrelator
{
public:
    
    /**
    * Correlator engine constructor
    * \param correlation_records Reference to the (message) correlation records
    * \param autocorrelation_records Reference to the (PRN) autocorrelation records
    * \param message_correlator Reference to the message correlator
    * \param prn_autocorrelator Reference to the PRN autocorrelator
    */
    UnpilotedMultiplePrnCorrelator(
    		std::vector<CorrelationRecord>& correlation_records,
    		std::vector<AutocorrelationRecord>& autocorrelation_records,
            UnpilotedMessageCorrelator& message_correlator,
            PrnAutocorrelator& prn_autocorrelator);

    virtual ~UnpilotedMultiplePrnCorrelator();
    /**
     * Set the PRN source block pointer. It is assumed the samples are one PRN length.
     * \param fftw_source_block source samples
     */
    void set_source_block(wsgc_complex *source_block);
    
    /**
     * Execute one correlation operation with the current PRN samples. The Global PRN index is incremented after
     * each call to this method.
     */
    void make_correlation();
    
    /**
     * Get the dictionnary of PRN shift occurences
     * \return Reference  to the dictionnary of PRN shift occurences
     */
    const std::map<unsigned int, unsigned int>& get_shift_occurences() const
    {
        return _shift_occurences_dict;
    }
    
    /**
     * Print out the message correlation records
     * \param mag_factor Magnitude displays are divided by this constant
     * \param os The output string stream
     */
    void dump_message_correlation_records(wsgc_float mag_factor, std::ostringstream& os) const;
    
    /**
     * Print out the PRN autocorrelation records
     * \param mag_factor Magnitude displays are divided by this constant
     * \param os The output string stream
     */
    void dump_prn_autocorrelation_records(wsgc_float mag_factor, std::ostringstream& os) const;


protected:
    std::vector<CorrelationRecord>& _correlation_records; //!< Reference to the message correlation records
    std::vector<AutocorrelationRecord>& _autocorrelation_records; //!< Reference to the PRN autocorrelation records
    UnpilotedMessageCorrelator& _message_correlator; //!< Reference to the message correlator
    PrnAutocorrelator& _prn_autocorrelator; //!< Reference to the PRN autocorrelator
    unsigned int _global_prn_index; //!< Index of PRN in the whole transmission (PRN counter)
    std::map<unsigned int, unsigned int> _shift_occurences_dict; //!< Dictionnary of IFFT indexes with maximum magnitude
};

#endif // __PILOTED_MULTIPLE_PRN_CORRELATOR_H__
