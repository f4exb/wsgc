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

     MessageCorrelationMatrices

	 This class encapsulates message correlation sparse matrices used with piloted resolution
	 The matrices are:
	 - Correlation results: 3 dimensional Delta-t x PRNi x Pi where Pi is the PRN sequence within symbol
	 	 This matrix has only one non null column (PRNi) per PRNi x Delta-t plane
	 	 It is stored as one PRNi x Pi matrix of correlation values and one PRNi vector of Delta-t values
	 - Averaging sum results: 2 dimensional PRNi x Delta-t
	 	 This matrix has only one non null column (PRNi) per PRNi x Delta-t plane
	 	 It is stored as one PRNi x Pi matrix of correlation values and one PRNi vector of Delta-t values
*/

#include "MessageCorrelationMatrices.h"
#include "WsgcUtils.h"
#include <cstring>
#include <iostream>

MessageCorrelationMatrices::MessageCorrelationMatrices(unsigned int nb_message_prns, unsigned int prn_per_symbol, unsigned int fft_N) :
	_nb_message_prns(nb_message_prns),
	_prn_per_symbol(prn_per_symbol),
	_fft_N(fft_N),
	_nb_averaging_sums(0),
	_nb_corr_vectors(0),
	_corr_item_index(0),
	_corr_pi_at_zero_index(0)
{
	_corr_values = new wsgc_complex[_nb_message_prns * _prn_per_symbol];
	_corr_delta_t = new unsigned int[_prn_per_symbol];
	_avg_values = new wsgc_complex[_nb_message_prns * _prn_per_symbol];
	_avg_delta_t = new unsigned int[_prn_per_symbol];
}


MessageCorrelationMatrices::~MessageCorrelationMatrices()
{
	delete[] _avg_delta_t;
	delete[] _avg_values;
	delete[] _corr_delta_t;
	delete[] _corr_values;
}


void MessageCorrelationMatrices::add_correlation_item(unsigned int prni, unsigned int pi, unsigned int delta_t, wsgc_complex value)
{
	_corr_values[pi*_nb_message_prns + prni] = value;
	_corr_delta_t[pi] = delta_t;
}


void MessageCorrelationMatrices::validate_prni_vector(unsigned int pi)
{
	if (_nb_corr_vectors == 0)
	{
		_corr_pi_at_zero_index = pi;
	}

	if (_nb_corr_vectors < _prn_per_symbol)
	{
		_nb_corr_vectors++;
	}
}


void MessageCorrelationMatrices::process_averaging()
{
	_avg_dict.clear();
	std::map<unsigned int, std::vector<wsgc_complex> >::iterator avg_dict_it;

	for (unsigned int pi=0; pi < _nb_corr_vectors; pi++)
	{
		avg_dict_it = _avg_dict.find(_corr_delta_t[pi]);

		if (avg_dict_it == _avg_dict.end())
		{
			for (unsigned int prni = 0; prni < _nb_message_prns; prni++)
			{
				_avg_dict[_corr_delta_t[pi]].push_back(_corr_values[pi*_nb_message_prns + prni]);
			}
		}
		else
		{
			for (unsigned int prni = 0; prni < _nb_message_prns; prni++)
			{
				_avg_dict[_corr_delta_t[pi]][prni] += _corr_values[pi*_nb_message_prns + prni];
			}
		}
	}
}


void MessageCorrelationMatrices::get_mag_max(CorrelationTuple_t& correlation_tuple, wsgc_float &max_mag, wsgc_float &avg_mag)
{
	max_mag = 0.0;
    wsgc_float sum_mag[_prn_per_symbol];
    unsigned int max_ai;
	wsgc_float mag;
	unsigned int max_prni;
	unsigned int max_delta_t;
	wsgc_complex max_value;

	std::map<unsigned int, std::vector<wsgc_complex> >::iterator avg_dict_it = _avg_dict.begin();
	const std::map<unsigned int, std::vector<wsgc_complex> >::const_iterator avg_dict_end = _avg_dict.end();

	for (unsigned int ai=0; avg_dict_it != avg_dict_end; ++avg_dict_it, ai++)
	{
        sum_mag[ai] = 0.0;
    
		for (unsigned int prni = 0; prni < _nb_message_prns-1; prni++) // exclude noise PRN
		{
			WsgcUtils::magnitude_estimation(&avg_dict_it->second[prni], &mag);
            sum_mag[ai] += mag;

			if (mag > max_mag)
			{
				max_mag = mag;
				max_delta_t = avg_dict_it->first;
				max_prni = prni;
				max_value = avg_dict_it->second[prni];
                max_ai = ai;
			}
		}
	}

    avg_mag = sum_mag[max_ai]/(_nb_message_prns-1);
	correlation_tuple.delta_t = max_delta_t;
	correlation_tuple.prni = max_prni;
	correlation_tuple.value = max_value;
}


void MessageCorrelationMatrices::get_noise_mag_max(CorrelationTuple_t& correlation_tuple, wsgc_float &max_mag)
{
	max_mag = 0.0;
	wsgc_float mag;
	unsigned int max_prni;
	unsigned int max_delta_t;
	wsgc_complex max_value;

	std::map<unsigned int, std::vector<wsgc_complex> >::iterator avg_dict_it = _avg_dict.begin();
	const std::map<unsigned int, std::vector<wsgc_complex> >::const_iterator avg_dict_end = _avg_dict.end();

	for (; avg_dict_it != avg_dict_end; ++avg_dict_it)
	{
		WsgcUtils::magnitude_estimation(&avg_dict_it->second[_nb_message_prns-1], &mag);

		if (mag > max_mag)
		{
			max_mag = mag;
			max_delta_t = avg_dict_it->first;
			max_prni = _nb_message_prns-1;
			max_value = avg_dict_it->second[_nb_message_prns-1];
		}
	}

	correlation_tuple.delta_t = max_delta_t;
	correlation_tuple.prni = max_prni;
	correlation_tuple.value = max_value;
}


wsgc_float MessageCorrelationMatrices::get_noise_avg()
{
	wsgc_float mag;
	wsgc_float mag_sum = 0.0;

	std::map<unsigned int, std::vector<wsgc_complex> >::iterator avg_dict_it = _avg_dict.begin();
	const std::map<unsigned int, std::vector<wsgc_complex> >::const_iterator avg_dict_end = _avg_dict.end();

	for (; avg_dict_it != avg_dict_end; ++avg_dict_it)
	{
		WsgcUtils::magnitude_estimation(&avg_dict_it->second[_nb_message_prns-1], &mag);
		mag_sum += mag;
	}

    return mag_sum;
}
