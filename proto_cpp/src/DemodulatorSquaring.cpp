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
      
     DemodulatorSquaring
      
*/
#include "DemodulatorSquaring.h"

DemodulatorSquaring::DemodulatorSquaring() {}

DemodulatorSquaring::~DemodulatorSquaring() {}

void DemodulatorSquaring::demodulate_in_place(wsgc_complex *samples_in, unsigned int samples_length)
{
    for (unsigned int i = 0; i < samples_length; i++)
    {
        samples_in[i].real() = std::norm(samples_in[i]);
        samples_in[i].imag() = 0.0;
    }
}

void DemodulatorSquaring::demodulate_out_of_place(wsgc_complex *samples_in, wsgc_complex *samples_out, unsigned int samples_length)
{
    for (unsigned int i = 0; i < samples_length; i++)
    {
        samples_out[i].real() = std::norm(samples_in[i]);
        samples_out[i].imag() = 0.0;
    }
}
