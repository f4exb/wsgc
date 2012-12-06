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
      
     Source FFT
      
     Multiplies input samples by each of the sub-frequency LOs
     Do the FFT of the result
     
*/
#include "SourceFFT.h"

SourceFFT::SourceFFT(wsgc_float f_sampling,
		  wsgc_float f_chip,
		  unsigned int fft_N,
		  unsigned int freq_step_division) :
    _f_sampling(f_sampling),
    _f_chip(f_chip),
    _fft_N(fft_N),
    _freq_step_division(freq_step_division)
{}

SourceFFT::~SourceFFT()
{}
