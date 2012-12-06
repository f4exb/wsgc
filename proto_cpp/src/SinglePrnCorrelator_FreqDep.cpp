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
      
     SinglePrnCorrelator_FreqDep
       
     Does a IFFT batch correlation and averaging given one single PRN to look for over a range of frequencies.
     
*/
#include "SinglePrnCorrelator_FreqDep.h"
#include "WsgcTypes.h"
#include "WsgcUtils.h"
#include "LocalCodesFFT.h"
#include <iostream>

const wsgc_complex SinglePrnCorrelator_FreqDep::cone(1.0, 0.0);

SinglePrnCorrelator_FreqDep::SinglePrnCorrelator_FreqDep(
			unsigned int fft_N,
			unsigned int nb_f_bins,
			unsigned int prn_per_block,
			unsigned int nb_batch_prns,
			unsigned int freq_step_division) :
    _fft_N(fft_N),
    _nb_f_bins(nb_f_bins),
    _freq_step_division(freq_step_division),
    _prn_per_block(prn_per_block),
    _nb_batch_prns(nb_batch_prns),
    _storage_depth(2*nb_batch_prns),
    _batch_index(0),
    _batch_max_magnitudes(nb_batch_prns,0.0),
    _batch_max_composite_indexes(nb_batch_prns,0),
    _batch_complex_values_max(nb_batch_prns, cone)
{
}


SinglePrnCorrelator_FreqDep::~SinglePrnCorrelator_FreqDep()
{
}
        
