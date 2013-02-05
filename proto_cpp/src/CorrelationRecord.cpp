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
#include "CorrelationRecord.h"
#include <iostream>
#include <iomanip>

CorrelationRecord::CorrelationRecord()
{
	reset();
}


CorrelationRecord::~CorrelationRecord()
{}


void CorrelationRecord::reset()
{
    global_prn_index = 0;     // PRN index in the global sequence of received PRNs
    prn_per_symbol_index = 0; // PRN index in the arbitrary PRN per symbol cycle
    prn_index_max = 0;        // Index of the PRN in GC sequences having maximum correlation
    magnitude_max = 0.0;      // Magnitude of maximum correlation
    magnitude_avg = 0.0;      // Average of correlation magnitudes
    shift_index_max = 0;      // Time shift of correlation peak in PRN sequence
    phase_at_max = 0.0;       // Phase of correlation at maximum
    noise_max = 0.0;          // Noise PRN maximum correlation magnitude
    noise_avg = 0.0;          // Average of noise PRN samples magnitude
    f_rx = 0.0;               // Receiving frequency relative to input samples frequency
    frequency_locked = false; // Indicates if receiving frequency was in locked status
    pilot_shift = 0;          // Pilot sequence shift corresponding to this PRN (piloted operation only)
    selected = false;         // PRN sample has been selected as valid for message correlation by the pilot correlation analyser
}
        
void CorrelationRecord::dump(unsigned int magnitude_factor, std::ostringstream& os) const
{
    os << std::setiosflags(std::ios_base::fixed);
    //"-- %2d:%2d>" % ((bn-1)/self.N,Ni), "Fdi:", fd_i, "Ni:", Ni, "PRN:", "%03d"% prn_index_max, "Max:", "%07.1f"% value_max, "Shift:", "%04d" % bin_index_max, "Phase @ max: %5.2f" % phase_max, "Fd: %.2f" % fd, "Locked:", frequency_locked
    os << "PNi: " << std::setw(3) << global_prn_index
       << " Ni: " << std::setw(2) << prn_per_symbol_index
       << " PRN: " << std::setw(3) << std::setfill('0') << prn_index_max
       << " Max: " << std::setw(7) << std::setprecision(1) << magnitude_max / magnitude_factor;

	if (magnitude_avg != 0.0)
	{
	   os << " Max/Avg: " << std::setw(7) << std::setprecision(3) << magnitude_max / magnitude_avg;
	}

    os << " Shift: " << std::setw(4) << shift_index_max
       << " Pilot shift: " << std::setw(4) << pilot_shift
       << " Noise max: " << std::setw(7) << std::setprecision(1) << noise_max / magnitude_factor
       << " Noise avg: " << std::setw(7) << std::setprecision(1) << noise_avg / magnitude_factor
       << " Phase @ max: " << std::setw(5) << std::setprecision(2) << std::setfill(' ') << phase_at_max
       << " Frx: " << std::setw(6) << std::setprecision(2) << std::setfill(' ') << f_rx
       << " Locked: " << (frequency_locked ? "yes" : " no")
       << " Selected: " << (selected ? "yes" : " no");

    if (noise_avg != 0.0)
    {
    	os << " S/N: " << std::setw(5) << std::setprecision(2) << magnitude_max / noise_avg;
    }
    else
    {
    	os << " S/N: N/A";
    }
}


void CorrelationRecord::dump_line(wsgc_float magnitude_factor, std::ostringstream& os) const
{
    os << std::setiosflags(std::ios_base::fixed);
    os << std::setw(3) << global_prn_index
       << " " << std::setw(2) << prn_per_symbol_index
       << " " << std::setw(3) << std::setfill('0') << prn_index_max
       << " " << std::setw(7) << std::setprecision(1) << magnitude_max / magnitude_factor;

    if (magnitude_avg == 0.0)
    {
    	os << " - N/A -";
    }
    else
    {
       os << " " << std::setw(7) << std::setprecision(3) << magnitude_max / magnitude_avg;
    }

    os << " " << std::setw(4) << shift_index_max
       << " " << std::setw(4) << pilot_shift
       << " " << std::setw(7) << std::setprecision(1) << noise_max / magnitude_factor
       << " " << std::setw(7) << std::setprecision(1) << noise_avg / magnitude_factor
       << " " << std::setw(5) << std::setprecision(2) << std::setfill(' ') << phase_at_max
       << " " << std::setw(6) << std::setprecision(2) << std::setfill(' ') << f_rx
       << " " << (frequency_locked ? " Y" : " N");

    if (noise_avg != 0.0)
    {
    	os << " " << std::setw(5) << std::setprecision(2) << magnitude_max / noise_avg;
    }
    else
    {
    	os << "   N/A";
    }
    
    os << " " << (selected ? "Y" : "N") << std::endl;
}


void CorrelationRecord::dump_banner(std::ostringstream& os)
{
    os << "PNi Ni PN# Mag.Max Max/Avg Ti.. PlTi Nse.Max Nse.Avg Ph@Mx Frx... Fl S/N.. S" << std::endl;
}


unsigned int CorrelationRecord::get_time_shift() const
{
    return shift_index_max;
}


wsgc_float CorrelationRecord::get_correlation_peak() const
{
    return magnitude_max;
}


void CorrelationRecord::set_selected(bool _selected)
{
    selected = _selected;
}