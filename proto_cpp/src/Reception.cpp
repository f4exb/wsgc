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
      
     Reception
      
     Common part of all recepion classes.
     
*/

#include "Reception.h"
#include "DemodulatorDifferential.h"
#include "DemodulatorSquaring.h"
#include <iostream>

//=================================================================================================
void Reception::demodulate_before_correlate(wsgc_complex *faded_source_samples, unsigned int nb_faded_source_samples)
{
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

//=================================================================================================
void Reception::generate_training_prn_list(std::vector<unsigned int>& prn_list)
{
	for (unsigned int prni=0; prni < gc_generator.get_nb_message_codes(); prni++) 
	{
		prn_list.push_back(prni);
	}
}

//=================================================================================================
void Reception::generate_training_prn_list_unpiloted(std::vector<unsigned int>& prn_list)
{
	unsigned int start_code = gc_generator.get_nb_message_codes() + gc_generator.get_nb_service_codes();

	for (unsigned int prni=start_code; prni < start_code + gc_generator.get_nb_training_codes(); prni++)
	{
		prn_list.push_back(prni);
	}
}

//=================================================================================================
void Reception::generate_message_prn_list(std::vector<unsigned int>& prn_list)
{
	for (unsigned int prni=0; prni <= gc_generator.get_nb_message_codes(); prni++) // incudes the noise PRN that is the last hence the <= comparison
	{
		prn_list.push_back(prni);
	}
}


//=================================================================================================
void Reception::generate_pilot_prn_list(std::vector<unsigned int>& prn_list, unsigned int pilot_prni)
{
    prn_list.push_back(gc_generator.get_nb_message_codes()+pilot_prni);
}

