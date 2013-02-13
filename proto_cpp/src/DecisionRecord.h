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
#ifndef __DECISION_RECORD_H__
#define __DECISION_RECORD_H__

#include "WsgcTypes.h"
#include <sstream>

class DecisionRecord
{
public:
    typedef enum
    {
        decision_ok_strong, //!< Symbol validated because values are strong enough
        decision_ok_not_enough_rec, //!< Symbol validated despite not enough records
        decision_ko_too_weak, //!< Symbol not validated because values are too weak for selected PRN
        decision_ko_not_enough_rec, //!< Symbol not validated because there are not enough valid PRN records in this cycle
        decision_ko_no_valid_rec //!< Symbol not validated because there are no valid records in this cycle
    } decision_type_t;

    DecisionRecord();
    ~DecisionRecord();
    
    void reset();
    void dump(std::ostringstream& os, wsgc_float magnitude_factor = 1.0) const;
    void dump_line(std::ostringstream& os, std::string tag, wsgc_float magnitude_factor = 1.0) const;
    static void dump_banner(std::ostringstream& os);
    
    // Ni, prn_index_max, value_max, bin_index_max, phase_max, fd, frequency_locked
    bool validated;                    //!< Symbol was validated
    decision_type_t decision_type;     //!< Type of decision that is taken
    unsigned int symbol_index;         //!< Symbol index in the message for which decison applies
    unsigned int global_prn_index;     //!< PRN index in the global sequence of received PRNs at which decision is taken
    unsigned int prn_per_symbol_index; //!< PRN index in the arbitrary PRN per symbol cycle at which decision is taken
    unsigned int prn_index_max;        //!< Index of the PRN in GC sequences having maximum correlation
    unsigned int select_count;         //!< Count of selected PRNs in the symbol cycle
    wsgc_float magnitude_max;          //!< Magnitude of maximum correlation for PRN at which decision is taken
    wsgc_float magnitude_avg;          //!< Average of correlation magnitudes for PRN at which decision is taken
    wsgc_float noise_avg;              //!< Average of noise PRN samples magnitude for PRN at which decision is taken
    unsigned int shift_index_max;      //!< Time shift of correlation peak in PRN sequence for PRN at which decision is taken
    wsgc_float f_rx;                   //!< Receiving frequency relative to input samples frequency for PRN at which decision is taken
    
protected:
    void dump_decision_type(std::ostringstream& os) const;
};

#endif // __CORRELATION_RECORD_H__
