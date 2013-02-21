/*
     Copyright 2012-2013 Edouard Griffiths <f4exb at free dot fr>
 
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
#include "MFSK_Options.h"
#include "WsgcUtils.h"


MFSK_Options::MFSK_Options(wsgc_float f_sampling) :
	_f_sampling(f_sampling),
	_symbol_bandwidth_log2(0),
	_symbol_time_log2(0),
	_zero_fft_slot(0),
	_fft_N(int(f_sampling)),
	_nb_fft_per_symbol(1),
	_symbol_time(1.0),
	_symbol_bandwidth(1),
	_zero_frequency(0.0)
{}


MFSK_Options::~MFSK_Options()
{}


bool MFSK_Options::parse_options(std::string& mfsk_params)
{
	bool status;
	std::vector<int> raw_mfsk_params;
	status = extract_vector<int>(raw_mfsk_params, mfsk_params);

	if (!status)
	{
		return false;
	}
	else if (raw_mfsk_params.size() < 3)
	{
		return false;
	}
	else
	{
		int symbol_bandwidth_log2 = raw_mfsk_params[0];
		_symbol_time_log2 = raw_mfsk_params[1];
		_zero_fft_slot = raw_mfsk_params[2];

		if (symbol_bandwidth_log2 < 0)
		{
			std::cout << "MFSK symbol bandwidth cannot be smaller than 1" << std::endl;
			return false;
		}
		else if (_symbol_time_log2+symbol_bandwidth_log2 < 0)
		{
			std::cout << "MFSK symbol time cannot be shorter than the inverse of its bandwidth" << std::endl;
			return false;
		}
		else
		{
			_symbol_bandwidth_log2 = symbol_bandwidth_log2;
			_symbol_bandwidth = (wsgc_float) (1 << _symbol_bandwidth_log2);
			_fft_N = int(_f_sampling / _symbol_bandwidth);

			if (abs(_zero_fft_slot) > _fft_N/2)
			{
				std::cout << "MFSK zero frequency slot cannot be larger than half the FFT size" << std::endl;
				return false;
			}

			if (_symbol_time_log2 < 0)
			{
				_symbol_time = 1.0 / (1 << -_symbol_time_log2);
			}
			else
			{
				_symbol_time = (wsgc_float) (1 << _symbol_time_log2);
			}

			_nb_fft_per_symbol = 1 << (_symbol_time_log2+symbol_bandwidth_log2);
			_zero_frequency = _zero_fft_slot * _symbol_bandwidth;

			return true;
		}
	}
}


void MFSK_Options::print_options(std::ostringstream& os)
{
    os << std::setiosflags(std::ios_base::fixed);
    os << "Sampling frequency ........: " << std::setw(8) << std::setprecision(1) << std::right << _f_sampling << std::endl;
    os << "Symbol bandwidth ..........: " << std::setw(8) << std::setprecision(1) << std::right << _symbol_bandwidth << std::endl;
    os << "Zero frequency ............: " << std::setw(8) << std::setprecision(1) << std::right << _zero_frequency << std::endl;
    os << "Symbol time ...............: " << std::setw(10) << std::setprecision(3) << std::right << _symbol_time << std::endl;
    os << "FFT size ..................: " << std::setw(6) << std::right << _fft_N << std::endl;
    os << "Zero FFT slot .............: " << std::setw(6) << std::right << _zero_fft_slot << std::endl;
    os << "Nb FFT per symbol .........: " << std::setw(6) << std::right << _nb_fft_per_symbol << std::endl;
}
