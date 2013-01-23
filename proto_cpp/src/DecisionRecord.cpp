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
      
     Decision record

     Stores one decision round data 
*/
#include "DecisionRecord.h"
#include <iostream>
#include <iomanip>

DecisionRecord::DecisionRecord()
{
	reset();
}


DecisionRecord::~DecisionRecord()
{}


void DecisionRecord::reset()
{
	validated = false;
    symbol_index = 0;         // Symbol index in the message for which decison applies
    global_prn_index = 0;     // PRN index in the global sequence of received PRNs at which decision is taken
    prn_per_symbol_index = 0; // PRN index in the arbitrary PRN per symbol cycle at which decision is taken
    prn_index_max = 0;        // Index of the PRN in GC sequences having maximum correlation
    select_count = 0;
    magnitude_max = 0.0;      // Magnitude of maximum correlation for PRN at which decision is taken
    magnitude_avg = 0.0;      // Average of correlation magnitudes for PRN at which decision is taken
    noise_avg = 0.0;          // Average of noise PRN samples magnitude for PRN at which decision is taken
    shift_index_max = 0;      // Time shift of correlation peak in PRN sequence for PRN at which decision is taken
    f_rx = 0.0;               // Receiving frequency relative to input samples frequency for PRN at which decision is taken
}
        
void DecisionRecord::dump(std::ostringstream& os, wsgc_float magnitude_factor) const
{
    os << std::setiosflags(std::ios_base::fixed);
    os << "Si: " << std::setw(2) << symbol_index
       << " PNi: " << std::setw(3) << global_prn_index
       << " Ni: " << std::setw(2) << prn_per_symbol_index
       << " PRN: " << std::setw(3) << std::setfill('0') << prn_index_max
       << " Max: " << std::setw(7) << std::setprecision(1) << magnitude_max / magnitude_factor;

	if (magnitude_avg != 0.0)
	{
	   os << " Max/Avg: " << std::setw(7) << std::setprecision(3) << magnitude_max / magnitude_avg;
	}
    else
    {
	   os << " Max/Avg: N/A";
    }

    os << " Noise avg: " << std::setw(7) << std::setprecision(1) << noise_avg / magnitude_factor;
    
    if (noise_avg != 0.0)
    {
    	os << " S/N: " << std::setw(5) << std::setprecision(2) << magnitude_max / noise_avg;
    }
    else
    {
    	os << " S/N: N/A";
    }
    
    os << " Shift: " << std::setw(4) << shift_index_max
       << " Frx: " << std::setw(6) << std::setprecision(2) << std::setfill(' ') << f_rx << " >";
       
    dump_decision_type(os);
    os << std::endl;
}


void DecisionRecord::dump_line(std::ostringstream& os, wsgc_float magnitude_factor) const
{
    os << std::setiosflags(std::ios_base::fixed);
    os << std::setw(3) << symbol_index
       << " " << std::setw(3) << global_prn_index
       << " " << std::setw(2) << prn_per_symbol_index
       << " " << std::setw(3) << prn_index_max
       << " " << std::setw(2) << select_count
       << " " << std::setw(7) << std::setfill('0') << std::setprecision(1) << magnitude_max / magnitude_factor;

    if (magnitude_avg == 0.0)
    {
    	os << " " << std::setw(7) << std::setprecision(3) << 0.0;
    }
    else
    {
    	os << " " << std::setw(7) << std::setprecision(3) << magnitude_max / magnitude_avg;
    }

    os << " " << std::setw(7) << std::setprecision(1) << noise_avg / magnitude_factor;   
    
    if (noise_avg != 0.0)
    {
    	os << " " << std::setw(5) << std::setprecision(2) << magnitude_max / noise_avg;
    }
    else
    {
    	os << " " << std::setw(5) << std::setprecision(2) << 0.0;
    }

    os << " " << std::setw(4) << shift_index_max
       << " " << std::setw(6) << std::setprecision(2) << std::setfill(' ') << f_rx << " ";
       
    dump_decision_type(os);
    os << std::endl;
}


void DecisionRecord::dump_banner(std::ostringstream& os)
{
    os << "PNi  Ni Pi PN# S# Mag.Max Max/Avg Nse.Avg S/N.. Ti.. Frx... Decision" << std::endl;
}


void DecisionRecord::dump_decision_type(std::ostringstream& os) const
{
    switch(decision_type)
    {
        case decision_ok_strong:
            os << "OK_++";
            break;
        case decision_ok_not_enough_rec:
            os << "OK_NE";
            break;
        case decision_ko_too_weak:
            os << "KO_TW";
            break;
        case decision_ko_not_enough_rec:
            os << "KO_NE";
            break;
        case decision_ko_no_valid_rec:
            os << "KO_NR";
            break;
        default:
            os << "UN";
            break;
    }
}
