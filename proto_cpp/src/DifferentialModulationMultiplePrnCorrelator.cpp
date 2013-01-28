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

     DifferentialModulationMultiplePrnCorrelator

     This flavour of correlator deals with PRNs encoded with differential modulation

*/

#include "DifferentialModulationMultiplePrnCorrelator.h"
#include <assert.h>

DifferentialModulationMultiplePrnCorrelator::DifferentialModulationMultiplePrnCorrelator(
		wsgc_float f_sampling,
		wsgc_float f_chip,
		unsigned int prn_length,
		const std::vector<unsigned int>& prn_list,
		unsigned int prn_window_size) :
		_f_sampling(f_sampling),
		_f_chip(f_chip),
		_prn_length(prn_length),
		_fft_N(int((f_sampling*prn_length)/f_chip)),
		_fractional_chip_per_sample(prn_length/_fft_N),
		_fractional_samples_per_chip(_fft_N/prn_length),
		_int_samples_per_chip(int(_fractional_samples_per_chip)),
		_prn_list(prn_list),
		_samples_length(0),
		_prn_window_size(prn_window_size)
{
    assert(_fractional_samples_per_chip > 4.0); // Need at least 4 samples per chip
}
  
  
DifferentialModulationMultiplePrnCorrelator::~DifferentialModulationMultiplePrnCorrelator()
{}
