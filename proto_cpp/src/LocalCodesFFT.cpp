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
      
     LocalCodes
      
     Creates a local copy of the complex conjugate of the FFT transform of all symbols gold codes
     This pre-calculates the code to be inserted in the final IFFT.
     Used with BPSK complex signals.
     
*/
#include "LocalCodesFFT.h"
#include "GoldCodeGenerator.h"
#include "CodeModulator.h"
#include <string.h>
#include <assert.h>

LocalCodesFFT::LocalCodesFFT(
			CodeModulator& code_modulator,
			GoldCodeGenerator& gc_generator,
			wsgc_float f_sampling,
			wsgc_float f_chip,
			std::vector<unsigned int>& symbols) :
	_code_modulator(code_modulator),
	_gc_generator(gc_generator),
	_f_sampling(f_sampling),
	_f_chip(f_chip),
	_nb_code_samples(gc_generator.get_nb_code_samples(f_sampling,f_chip)),
	_symbols(symbols)
{
}

LocalCodesFFT::~LocalCodesFFT()
{
}
