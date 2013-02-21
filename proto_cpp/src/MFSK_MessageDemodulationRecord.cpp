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

     MFSK_MessageDemodulationRecord

     Class that holds messsage demodulation data

*/
#include "MFSK_MessageDemodulationRecord.h"
#include <iomanip>
#include <iostream>

MFSK_MessageDemodulationRecord::MFSK_MessageDemodulationRecord()
{
	reset();
}


MFSK_MessageDemodulationRecord::~MFSK_MessageDemodulationRecord()
{}


void MFSK_MessageDemodulationRecord::reset()
{
	_symbol_index = 0;
	_symbol_ordinal = 0;
	_max_magnitude = 0.0;
	_avg_magnitude = 0.0;
	_noise_magnitude = 0.0;
}


void MFSK_MessageDemodulationRecord::dump(std::ostringstream& os, wsgc_float magnitude_factor) const
{
    os << std::setiosflags(std::ios_base::fixed);
    os << "Si: " << std::setw(3) << _symbol_index
       << " So: " << std::setw(3) << std::setfill('0') << _symbol_ordinal
       << " Max: " << std::setw(7) << std::setprecision(1) << _max_magnitude / magnitude_factor
       << " Max/Avg: " << std::setw(7) << std::setprecision(1) << _max_magnitude / _avg_magnitude
       << " Nse: " << std::setw(7) << std::setprecision(1) << _noise_magnitude / magnitude_factor << std::endl;
}


void MFSK_MessageDemodulationRecord::dump_line(std::ostringstream& os, wsgc_float magnitude_factor) const
{
    os << std::setiosflags(std::ios_base::fixed);
    os << std::setfill('0') << std::setw(3) << _symbol_index
       << " " << std::setw(5) << _fft_index
       << " " << std::setw(3) << _symbol_ordinal
       << " " << std::setw(7) << std::setprecision(1) << _max_magnitude / magnitude_factor;

    if (_avg_magnitude == 0.0)
    {
    	os << " - N/A -";
    }
    else
    {
       os << " " << std::setw(7) << std::setprecision(3) << _max_magnitude / _avg_magnitude;
    }

    os << " " << std::setw(7) << std::setprecision(1) <<  _noise_magnitude / magnitude_factor << std::endl;
}


void MFSK_MessageDemodulationRecord::dump_banner(std::ostringstream& os)
{
    os << "Si. FFTi. So. Mag.Max Max/Avg Nse.Max" << std::endl;
}

