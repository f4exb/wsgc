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
      
     Correlation record

     Stores one correlation round data 
*/
#include "TrainingCorrelationRecord.h"
#include <iostream>
#include <iomanip>

TrainingCorrelationRecord::TrainingCorrelationRecord(unsigned int sequence_length, unsigned int analysis_window_prns) :
	_sequence_length(sequence_length),
	_analysis_window_prns(analysis_window_prns)
{
	reset();
}


TrainingCorrelationRecord::~TrainingCorrelationRecord()
{}


void TrainingCorrelationRecord::reset()
{
    _global_prn_index = 0;     // PRN index in the global sequence of received PRNs
    _prn_index_max = 0;        // Index of the PRN in GC sequences having maximum correlation
    _magnitude_max = 0.0;      // Magnitude of maximum correlation
    _magnitude_avgsum = 0.0;   // Shifted averaging sum of correlation magnitudes
    _max_selected = false;     // Maximum has been selected
    _pilot_shift = 0;          // Pilot sequence shift corresponding to this PRN (piloted operation only)
    _selected = false;         // PRN sample has been selected as valid for message correlation by the pilot correlation analyser
}
        
void TrainingCorrelationRecord::dump(std::ostringstream& os, unsigned int magnitude_factor) const
{
    os << std::setiosflags(std::ios_base::fixed);
    //"-- %2d:%2d>" % ((bn-1)/self.N,Ni), "Fdi:", fd_i, "Ni:", Ni, "PRN:", "%03d"% prn_index_max, "Max:", "%07.1f"% value_max, "Shift:", "%04d" % bin_index_max, "Phase @ max: %5.2f" % phase_max, "Fd: %.2f" % fd, "Locked:", frequency_locked
    os << "PNi: " << std::setw(3) << _global_prn_index
       << " Ai: " << std::setw(2) << _global_prn_index % _analysis_window_prns
       << " PRN: " << std::setw(3) << std::setfill('0') << _prn_index_max
       << " Max: " << std::setw(7) << std::setprecision(1) << _magnitude_max / magnitude_factor;

	if (_magnitude_avgsum != 0.0)
	{
	   os << " Max/Avg: " << std::setw(7) << std::setprecision(3) << (_magnitude_max / _magnitude_avgsum)*_sequence_length;
	}

	os << " MaxSelected: " << (_max_selected ? "yes" : "no");

    os << " Pilot shift: " << std::setw(4) << _pilot_shift
       << " Selected: " << (_selected ? "yes" : " no");
}


void TrainingCorrelationRecord::dump_line(std::ostringstream& os, wsgc_float magnitude_factor) const
{
    os << std::setiosflags(std::ios_base::fixed);
    os << std::setw(3) << _global_prn_index
       << " " << std::setw(2) << _global_prn_index % _analysis_window_prns
       << " " << std::setw(3) << std::setfill('0') << _prn_index_max
       << " " << std::setw(7) << std::setprecision(1) << _magnitude_max / magnitude_factor;

    if (_magnitude_avgsum == 0.0)
    {
    	os << " - N/A -";
    }
    else
    {
       os << " " << std::setw(7) << std::setprecision(3) << (_magnitude_max / _magnitude_avgsum)*_sequence_length;
    }

    os << " " << (_max_selected ? " Y" : " N");

    os << " " << std::setw(4) << _pilot_shift
       << " " << (_selected ? " Y" : " N")
       << std::endl;
}


void TrainingCorrelationRecord::dump_banner(std::ostringstream& os)
{
    os << "PNi Ai PN# Mag.Max Max/Avg MS Ti..  S" << std::endl;
}
