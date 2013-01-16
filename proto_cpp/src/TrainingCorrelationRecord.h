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
      
     Correlation record for the training sequence

     Stores one correlation round data 
*/
#ifndef __TRAINING_CORRELATION_RECORD_H__
#define __TRAINING_CORRELATION_RECORD_H__

#include "WsgcTypes.h"
#include <sstream>

/**
 * \brief Correlation record for the training sequence.  Stores one correlation round data.
 */
class TrainingCorrelationRecord
{
public:
	/**
	 * Create a new Training Correlation Record
	 * \param sequence_length Length of the training sequence
	 * \param analysis_window_prns Length of analysis window in number of PRNs
	 */
	TrainingCorrelationRecord(unsigned int sequence_length, unsigned int analysis_window_prns);
	~TrainingCorrelationRecord();

	void reset();
	void dump(std::ostringstream& os, unsigned int magnitude_factor=1.0) const;
	void dump_line(std::ostringstream& os, wsgc_float magnitude_factor=1.0) const;
	static void dump_banner(std::ostringstream& os);

	unsigned int _sequence_length;      //!< Length of the training sequence
	unsigned int _analysis_window_prns; //!< Length of analysis window in number of PRNs
	unsigned int _global_prn_index;     //!< PRN index in the global sequence of received PRNs
	unsigned int _prn_index_max;        //!< Index of the PRN in GC sequences having maximum cumulative shifting correlation
	wsgc_float   _magnitude_max;        //!< Magnitude of maximum cumulative shifting correlation
	wsgc_float   _magnitude_avgsum;     //!< Sum of cumulative shifting correlation magnitudes
	bool         _max_selected;         //!< Maximum has been selected
	unsigned int _pilot_shift;          //!< Pilot sequence shift corresponding to this PRN (piloted operation only)
	bool         _selected;             //!< PRN sample has been selected as valid for message correlation by the pilot correlation analyser
};

#endif // __TRAINING_CORRELATION_RECORD_H__
