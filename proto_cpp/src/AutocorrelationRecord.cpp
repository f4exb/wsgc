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
#include "AutocorrelationRecord.h"
#include <iostream>
#include <iomanip>

AutocorrelationRecord::AutocorrelationRecord()
{
	reset();
}


AutocorrelationRecord::~AutocorrelationRecord()
{}


void AutocorrelationRecord::reset()
{
    global_prn_index = 0;     // PRN index in the global sequence of received PRNs
    prn_per_symbol_index = 0; // PRN index in the arbitrary PRN per symbol cycle
    magnitude = 0.0;         // Magnitude of maximum correlation
    magnitude_avg = 0.0;      // Average of correlation magnitudes
}
        
void AutocorrelationRecord::dump(unsigned int magnitude_factor, std::ostringstream& os) const
{
    os << std::setiosflags(std::ios_base::fixed);
    //"-- %2d:%2d>" % ((bn-1)/self.N,Ni), "Fdi:", fd_i, "Ni:", Ni, "PRN:", "%03d"% prn_index_max, "Max:", "%07.1f"% value_max, "Shift:", "%04d" % bin_index_max, "Phase @ max: %5.2f" % phase_max, "Fd: %.2f" % fd, "Locked:", frequency_locked
    os << "PNi: " << std::setw(3) << global_prn_index
       << " Ni: " << std::setw(2) << prn_per_symbol_index
       << " Mag: " << std::setw(7) << std::setprecision(1) << magnitude / magnitude_factor
       << " Avg: " << std::setw(7) << std::setprecision(1) << magnitude_avg / magnitude_factor;

	if (magnitude_avg != 0.0)
	{
	   os << " Mag/Avg: " << std::setw(7) << std::setprecision(3) << magnitude / magnitude_avg;
	}
}


void AutocorrelationRecord::dump_line(wsgc_float magnitude_factor, std::ostringstream& os) const
{
    os << std::setiosflags(std::ios_base::fixed);
    os << std::setw(3) << global_prn_index
       << " " << std::setw(2) << prn_per_symbol_index
       << " " << std::setw(7) << std::setprecision(1) << magnitude / magnitude_factor
       << " " << std::setw(7) << std::setprecision(1) << magnitude_avg / magnitude_factor;

    if (magnitude_avg == 0.0)
    {
    	os << " - N/A -";
    }
    else
    {
        os << " " << std::setw(7) << std::setprecision(3) << magnitude / magnitude_avg;
    }

    os << std::endl;
}


void AutocorrelationRecord::dump_banner(std::ostringstream& os)
{
    os << "PNi Ni Mag.... Avg.... Mag/Avg" << std::endl;
}
