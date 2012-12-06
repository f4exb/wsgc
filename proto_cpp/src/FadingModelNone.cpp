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
      
     Fading Model
      
     Applies a fading model to the test signal. Actual Fading model is instantiated in derivated classes.
     This parent class supports the AWGN addition which is common to all.
*/

#include "FadingModelNone.h"
#include <iostream>


FadingModelNone::FadingModelNone(wsgc_float f_sampling) :
    FadingModel::FadingModel(f_sampling, false)
{
}


FadingModelNone::~FadingModelNone() 
{}


void FadingModelNone::apply_fading(const wsgc_complex *samples_in, wsgc_complex *samples_out, unsigned int nb_samples)
{
    for (unsigned int sample_i=0; sample_i<nb_samples; sample_i++)
    {
        samples_out[sample_i] = samples_in[sample_i];
    }    
}


void FadingModelNone::print_fading_data(std::ostringstream& os) const
{
    os << "None";
}
