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
      
     Pilot correlation record

     Stores one pilot correlation round data 
*/
#include "PilotCorrelationRecord.h"
#include <iostream>
#include <iomanip>

PilotCorrelationRecord::PilotCorrelationRecord()
{
    reset();
}

PilotCorrelationRecord::~PilotCorrelationRecord()
{}
        
void PilotCorrelationRecord::dump(std::ostringstream& os) const
{
    os << std::setiosflags(std::ios_base::fixed);
    os << "Bi ...........: " << std::setw(2) << block_count << std::endl
       << " PRNi ........: " << std::setw(3) << prn_index << std::endl
       << " PRN .........: " << std::setw(3) << std::setfill('0') << pilot_index << std::endl
       << " Max .........: " << std::setw(7) << std::setprecision(1) << magnitude_max << std::endl
       << " Phase @ max .: " << std::setw(5) << std::setprecision(2) << std::setfill(' ') << phase_at_max << std::endl
       << " Ti ..........: " << std::setw(4) << t_index_max << std::endl
       << " Fi ..........: " << std::setw(4) << f_index_max << std::endl
       << " dF ..........: " << std::setw(7) << std::setprecision(2) << delta_f
       << " Selected ....: " << std::setw(2) << (selected ? "Y" : "N") << std::endl;
}

void PilotCorrelationRecord::dump_oneline(std::ostringstream& os) const
{
    os << std::setiosflags(std::ios_base::fixed);
    os << std::setw(2) << block_count  
       << " " << std::setw(3) << prn_index 
       << " " << std::setw(3) << std::setfill('0') << pilot_index
       << " " << std::setw(7) << std::setprecision(1) << magnitude_max
       << " " << std::setw(5) << std::setprecision(2) << std::setfill(' ') << phase_at_max
       << " " << std::setw(4) << t_index_max
       << " " << std::setw(4) << f_index_max
       << " " << std::setw(7) << std::setprecision(2) << delta_f
       << " " << std::setw(1) << (selected ? "Y" : "N") << std::endl;
}

void PilotCorrelationRecord::dump_oneline_banner(std::ostringstream& os)
{
    os << "Ai PNi PN# Mag.Max Ph@Mx Ti.. Fi.. delta.F S" << std::endl;
}



void PilotCorrelationRecord::reset()
{
    prn_index = 0;            // Global index of PRN
    block_count = 0;          // Number of averaging blocks processed
    pilot_index = 0;          // Pilot PRN index in GC codes
    magnitude_max = 0.0;      // Maximum magnitude of correlation
    phase_at_max = 0.0;       // Phase at the correlation maximum
    t_index_max = 0;          // Time index (i.e. delay) of the correlation maximum
    f_index_max = 0;          // Frequency index of the correlation maximum
    delta_f = 0.0;            // Corresponding frequency displacement of the correlation maximum relative to zero IF
    selected = true;          // PRN sample has been selected as valid for message correlation by the pilot correlation analyser
}
