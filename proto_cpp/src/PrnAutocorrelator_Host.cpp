/*
 * PrnAutocorrelator_Host.cpp
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

     PrnAutocorrelator

     Do an auto-correlation of one PRN to the next. Host version.

 */

#include "PrnAutocorrelator_Host.h"
#include "WsgcUtils.h"
#include <cstring>

PrnAutocorrelator_Host::PrnAutocorrelator_Host(
			unsigned int fft_N,
			unsigned int prn_per_symbol) :
    PrnAutocorrelator(fft_N, prn_per_symbol)
{
	_prn_samples = new wsgc_complex[2*fft_N];
	_acorrs = new wsgc_complex[prn_per_symbol];
}


PrnAutocorrelator_Host::~PrnAutocorrelator_Host()
{
	delete[] _acorrs;
	delete[] _prn_samples;
}


void PrnAutocorrelator_Host::set_source_block(const wsgc_complex *source_block, unsigned int global_prn_index)
{
	_global_prn_index = global_prn_index;
	unsigned int prn_parity =  global_prn_index % 2;
	memcpy(&_prn_samples[_fft_N*prn_parity], source_block, _fft_N);
}
    

void PrnAutocorrelator_Host::make_correlation(std::vector<AutocorrelationRecord>& autocorrelation_records)
{
	static const AutocorrelationRecord tmp_acorr_record;
	wsgc_complex acorr(0.0, 0.0);
	wsgc_complex acorr_avgsum = _acorrs[0];

	for (unsigned int psi = 0; psi < _fft_N; psi++)
	{
		acorr += _prn_samples[psi] * _prn_samples[psi + _fft_N];
	}

	_acorrs[_global_prn_index%_prn_per_symbol] = acorr;

	for (unsigned int ai=1; ai<_prn_per_symbol; ai++)
	{
		acorr_avgsum += _acorrs[ai];
	}

	autocorrelation_records.push_back(tmp_acorr_record);
	autocorrelation_records.back().global_prn_index = _global_prn_index;
	autocorrelation_records.back().prn_per_symbol_index = _global_prn_index % _prn_per_symbol;
	WsgcUtils::magnitude_estimation(&acorr, &autocorrelation_records.back().magnitude);
	WsgcUtils::magnitude_estimation(&acorr_avgsum, &autocorrelation_records.back().magnitude_avg);
	autocorrelation_records.back().acorr = acorr;
}
