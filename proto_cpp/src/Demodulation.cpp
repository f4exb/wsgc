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
      
     Demodulation
      
     Common part of all demodulation classes.
     
*/

#include "Demodulation.h" 
#include "DemodulatorDifferential.h"
#include "DemodulatorSquaring.h"
#include "Options.h"
#include "GoldCodeGenerator.h"

#include <iostream>

//=================================================================================================
void Demodulation::demodulate_generic(Options& options, const GoldCodeGenerator& gc_generator, wsgc_complex *faded_source_samples, unsigned int nb_faded_source_samples)
{
    unsigned int fft_N = gc_generator.get_nb_code_samples(options.f_sampling, options.f_chip);

    if (options.modulation.getScheme() == Modulation::Modulation_OOK)
    {
        std::cout << "Simulate AM power detection" << std::endl;
        DemodulatorSquaring demodulator;
        demodulator.demodulate_in_place(faded_source_samples, nb_faded_source_samples);
    }
    else if (options.modulation.isDifferential())
    {
        unsigned int int_samples_per_chip = ((wsgc_float) fft_N) /  gc_generator.get_code_length();
        static const wsgc_complex c_zero(0.0, 0.0);

        DemodulatorDifferential demodulator(int_samples_per_chip);
        demodulator.set_value_at_origin(c_zero);
        demodulator.demodulate_in_place(faded_source_samples, nb_faded_source_samples);
    }
}
