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
      
     
     Decision box thresholds
     
*/
#ifndef __DECISION_BOX_THRESHOLDS_H__
#define __DECISION_BOX_THRESHOLDS_H__

#include "WsgcTypes.h"
#include <string>
#include <sstream>

class DecisionBox_Thresholds
{
public:
	DecisionBox_Thresholds();
    ~DecisionBox_Thresholds();

    void set_cuda_defaults();
    bool parse_options(std::string thresholds_params);
	void print_options(std::ostringstream& os);
	void get_help(std::ostringstream& os);
        
	wsgc_float _peak_margin_threshold;             //!< Number of selected peak occurences / Number of correlated PRNs
    wsgc_float _max_to_avg_ok_threshold;           //!< Peak to average ratio unconditional selection lower bound threshold
    wsgc_float _max_to_avg_cdt_threshold;          //!< Peak to average ratio conditional selection lower bound threshold
    wsgc_float _signal_to_noise_avg_ko_threshold;  //!< Peak to noise average ratio unconditional rejection upper bound threshold
    wsgc_float _signal_to_noise_avg_cdt_threshold; //!< Peak to noise average ratio lower bound threshold for conditional selection on peak to average
    wsgc_float _max_to_avg_ok_threshold_ner_0_5;   //!< Bypass selected PRNs in symbol < 0.75 rejection if peak to average is greater than this and selected PRNs in symbol < 0.5
    wsgc_float _max_to_avg_ok_threshold_ner_0_75;  //!< Bypass selected PRNs in symbol < 0.75 rejection if peak to average is greater than this and selected PRNs in symbol in [0.5,0.75(
};

#endif // __DECISION_BOX_THRESHOLDS_H__
