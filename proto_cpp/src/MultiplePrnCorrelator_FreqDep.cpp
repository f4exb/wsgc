/*
     Copyright 2012 Edouard Griffiths <f4exb at free dot fr>
 
     This file is part of WSGC. A Weak Signal transmission mode using Gold Codes.
 
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
      
     MultiplePrnCorrelator_FreqDep
       
     Does a IFFT batch correlation and averaging given a list of PRNs to look for over a range of frequencies.
     
*/
#include "MultiplePrnCorrelator_FreqDep.h"
#include "WsgcTypes.h"
#include "WsgcUtils.h"
#include "LocalCodesFFT.h"
#include <iostream>

const wsgc_complex MultiplePrnCorrelator_FreqDep::cone(1.0, 0.0);

MultiplePrnCorrelator_FreqDep::MultiplePrnCorrelator_FreqDep(
			unsigned int fft_N,
			wsgc_float f_sampling,
			unsigned int nb_message_symbols,
			unsigned int nb_f_bins,
			unsigned int prn_per_symbol,
			unsigned int nb_batch_prns,
			unsigned int freq_step_division) :
    _fft_N(fft_N),
    _f_sampling(f_sampling),
    _nb_message_symbols(nb_message_symbols),
    _nb_f_bins(nb_f_bins),
    _freq_step_division(freq_step_division),
    _prn_per_symbol(prn_per_symbol),
    _nb_batch_prns(nb_batch_prns),
    _storage_depth(2*nb_batch_prns),
    _batch_index(0),
    _nb_codes(0)
{}


MultiplePrnCorrelator_FreqDep::~MultiplePrnCorrelator_FreqDep()
{}
        
