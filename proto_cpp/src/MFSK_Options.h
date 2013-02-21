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
      
     Options specific to MFSK

     Options parsing and holding
     
*/
#ifndef __MFSK_OPTIONS_H__
#define __MFSK_OPTIONS_H__

#include "WsgcTypes.h"
#include <string>
#include <sstream>

class MFSK_Options
{
public:
	MFSK_Options(wsgc_float f_sampling);
    ~MFSK_Options();

    bool parse_options(std::string& mfsk_params);
	void print_options(std::ostringstream& os);
        
	wsgc_float _f_sampling;
	unsigned int _symbol_bandwidth_log2;
	int _symbol_time_log2;
	int _zero_fft_slot;
	wsgc_float _symbol_bandwidth;
	wsgc_float _symbol_time;
	wsgc_float _zero_frequency;
	unsigned int _fft_N;
	unsigned int _nb_fft_per_symbol;
};

#endif // __MFSK_OPTIONS_H__
