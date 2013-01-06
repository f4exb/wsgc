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
#ifndef __CORRELATION_RECORD_H__
#define __CORRELATION_RECORD_H__

#include "WsgcTypes.h"
#include <sstream>

class CorrelationRecord
{
    public:
        CorrelationRecord();
        ~CorrelationRecord();
        
        void reset();
        void dump(unsigned int magnitude_factor, std::ostringstream& os) const;
        void dump_line(wsgc_float magnitude_factor, std::ostringstream& os) const;
        static void dump_banner(std::ostringstream& os);
        
        // Ni, prn_index_max, value_max, bin_index_max, phase_max, fd, frequency_locked
        unsigned int global_prn_index;     //!< PRN index in the global sequence of received PRNs
        unsigned int prn_per_symbol_index; //!< PRN index in the arbitrary PRN per symbol cycle
        unsigned int prn_index_max;        //!< Index of the PRN in GC sequences having maximum correlation
        unsigned int prn_index_max_i;      //!< Instantaneous (no averaging sum) index of the PRN in GC sequences having maximum correlation
        wsgc_float magnitude_max;          //!< Magnitude of maximum correlation
        wsgc_float magnitude_max_i;        //!< Instantaneous (no averaging sum) magnitude of maximum correlation
        wsgc_float magnitude_avg;          //!< Average of correlation magnitudes
        unsigned int shift_index_max;      //!< Time shift of correlation peak in PRN sequence
        wsgc_float phase_at_max;           //!< Phase of correlation at maximum
        wsgc_float noise_max;              //!< Noise PRN maximum correlation magnitude
        wsgc_float noise_avg;              //!< Average of noise PRN samples magnitude
        wsgc_float f_rx;                   //!< Receiving frequency relative to input samples frequency
        bool frequency_locked;             //!< Indicates if receiving frequency was in locked status
        unsigned int pilot_shift;          //!< Pilot sequence shift corresponding to this PRN (piloted operation only)
        bool selected;                     //!< PRN sample has been selected as valid for message correlation by the pilot correlation analyser
};

#endif // __CORRELATION_RECORD_H__
