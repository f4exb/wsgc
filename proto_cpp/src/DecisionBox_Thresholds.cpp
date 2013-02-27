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
#include "DecisionBox_Thresholds.h"
#include "WsgcUtils.h"
#include <vector>
#include <iomanip>

//=================================================================================================
DecisionBox_Thresholds::DecisionBox_Thresholds() :
    _peak_margin_threshold(0.13),
    _max_to_avg_ok_threshold(1.85),
    _max_to_avg_cdt_threshold(1.6),
    _signal_to_noise_avg_ko_threshold(1.0),
    _signal_to_noise_avg_cdt_threshold(2.0),
    _max_to_avg_ok_threshold_ner_0_5(4.0),
    _max_to_avg_ok_threshold_ner_0_75(3.0)
{}

    
//=================================================================================================
DecisionBox_Thresholds::~DecisionBox_Thresholds()
{}


//=================================================================================================
void DecisionBox_Thresholds::set_cuda_defaults()
{
    _max_to_avg_ok_threshold = 3.4; // was 3.0
    _max_to_avg_cdt_threshold = 2.6;
    _signal_to_noise_avg_ko_threshold = 3.0;
    _signal_to_noise_avg_cdt_threshold = 3.5;
    _max_to_avg_ok_threshold_ner_0_5 = 7.0;
    _max_to_avg_ok_threshold_ner_0_75 = 5.0;
}


//=================================================================================================
bool DecisionBox_Thresholds::parse_options(std::string thresholds_params)
{
	bool status;
	std::vector<wsgc_float> raw_params;
	status = extract_vector<wsgc_float>(raw_params, thresholds_params);

	if (status)
    {
        if (raw_params.size() > 0)
        {
        	if (!(raw_params[0] < 0.0))
        	{
        		_peak_margin_threshold = raw_params[0];
        	}
        }
        if (raw_params.size() > 1)
        {
        	if (!(raw_params[1] < 0.0))
        	{
        		_max_to_avg_ok_threshold = raw_params[1];
        	}
        }
        if (raw_params.size() > 2)
        {
        	if (!(raw_params[2] < 0.0))
        	{
        		_max_to_avg_cdt_threshold = raw_params[2];
        	}
        }
        if (raw_params.size() > 3)
        {
        	if (!(raw_params[3] < 0.0))
        	{
        		_signal_to_noise_avg_ko_threshold = raw_params[3];
        	}
        }
        if (raw_params.size() > 4)
        {
        	if (!(raw_params[4] < 0.0))
        	{
        		_signal_to_noise_avg_cdt_threshold = raw_params[4];
        	}
        }
        if (raw_params.size() > 5)
        {
        	if (!(raw_params[5] < 0.0))
        	{
        		_max_to_avg_ok_threshold_ner_0_5 = raw_params[5];
        	}
        }
        if (raw_params.size() > 6)
        {
        	if (!(raw_params[6] < 0.0))
        	{
        		_max_to_avg_ok_threshold_ner_0_75 = raw_params[6];
        	}
        }
    }

	return status;
}


//=================================================================================================
void DecisionBox_Thresholds::print_options(std::ostringstream& os)
{
    os << std::setiosflags(std::ios_base::fixed);
    os << "Decision thresholds (not all may be used):" << std::endl;
    os << "Peak margin ....................: " << std::setw(10) << std::setprecision(3) << std::right << _peak_margin_threshold << std::endl;
    os << "Max to Avg OK ..................: " << std::setw(10) << std::setprecision(3) << std::right << _max_to_avg_ok_threshold << std::endl;
    os << "Max to Avg Conditional .........: " << std::setw(10) << std::setprecision(3) << std::right << _max_to_avg_cdt_threshold << std::endl;
    os << "SNR Avg KO .....................: " << std::setw(10) << std::setprecision(3) << std::right << _signal_to_noise_avg_ko_threshold << std::endl;
    os << "SNR Avg Conditional ............: " << std::setw(10) << std::setprecision(3) << std::right << _signal_to_noise_avg_cdt_threshold << std::endl;
    os << "Max to Avg NER bypass (r<0.5) ..: " << std::setw(10) << std::setprecision(3) << std::right << _max_to_avg_ok_threshold_ner_0_5 << std::endl;
    os << "Max to Avg NER bypass (r<0.75) .: " << std::setw(10) << std::setprecision(3) << std::right << _max_to_avg_ok_threshold_ner_0_75 << std::endl;
}


//=================================================================================================
void DecisionBox_Thresholds::get_help(std::ostringstream& os)
{
    os << "Decision box thresholds..." << std::endl;
    os << "   Comma separated list of parameters in this order (negative means default):" << std::endl;
    os << "   - peak_margin_threshold (default 0.13) all correlation schemes (not MFSK)" << std::endl;
    os << "       . Number of selected peak occurences / Number of correlated PRNs" << std::endl;
    os << "   - max_to_avg_ok_threshold (default 1.85 / CUDA: 3.4) all schemes" << std::endl;
    os << "       . Peak to average ratio unconditional selection lower bound threshold" << std::endl;
    os << "   - max_to_avg_cdt_threshold (default 1.6 / CUDA: 2.6) BPSK if noise avg available" << std::endl;
    os << "       . Peak to average ratio conditional selection lower bound threshold (see signal_to_noise_avg_cdt_threshold)" << std::endl;
    os << "   - signal_to_noise_avg_ko_threshold (default 1.0 / CUDA: 3.0) BPSK if noise avg available" << std::endl;
    os << "       . Peak to noise average ratio unconditional rejection upper bound threshold" << std::endl;
    os << "   - signal_to_noise_avg_cdt_threshold (default 2.0 / CUDA: 3.5) BPSK if noise avg available" << std::endl;
    os << "       . Peak to noise average ratio lower bound threshold for conditional selection on peak to average" << std::endl;
    os << "   - max_to_avg_ok_threshold_ner_0_5 (default 4.0 / CUDA: 7.0) BPSK" << std::endl;
    os << "       . Bypass selected PRNs in symbol < 0.75 rejection lower bound value if selected PRNs in symbol < 0.5" << std::endl;
    os << "   - max_to_avg_ok_threshold_ner_0_75 (default 3.0 / CUDA: 5.0) BPSK" << std::endl;
    os << "       . Bypass selected PRNs in symbol < 0.75 rejection lower bound value if selected PRNs in symbol in [0.5,0.75(" << std::endl;
}
