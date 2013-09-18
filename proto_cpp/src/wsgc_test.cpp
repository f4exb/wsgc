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

     Main test program

*/
#include "WsgcTypes.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <time.h>
#include <cmath>
#include "Options.h"
#include "GoldCodeGenerator.h"
#include "WsgcUtils.h"
#include "Transmission.h"
#include "Reception_WSGC.h"
#include "Reception_WSGCE.h"
#include "Reception_WSGCO.h"
#include "Reception_WSGCD.h"
#include "Reception_MFSK.h"

#ifdef _CCSOFT
#include "CC_Encoding_base.h"
#endif

//=================================================================================================
int main(int argc, char *argv[])
{
	std::string binary_name(argv[0]);
    Options options(binary_name);

    srand(WsgcUtils::timenow_usec_hour());

    if (options.get_options(argc,argv))
    {
        std::ostringstream os;
        options.print_options(os);
        std::cout << os.str() << std::endl;

        unsigned int nb_message_symbols;

#ifdef _CCSOFT
        if (options.fec_scheme == Options::OptionFEC_CCSoft)
        {
            unsigned int tmp_n, tmp_k, tmp_m;
            ccsoft::get_cc_parameters<unsigned int>(options.cc_k_constraints, 
                options.cc_generator_polys,
                tmp_n,
                tmp_k,
                tmp_m);
            nb_message_symbols = 1<<tmp_n;
        }
        else
        {
            nb_message_symbols = options.nb_message_symbols;
        }
#else
        nb_message_symbols = options.nb_message_symbols;
#endif

        GoldCodeGenerator gc_generator(options.gc_nb_stages, nb_message_symbols, options.nb_service_symbols, options.nb_training_symbols, options.g1_poly_powers, options.g2_poly_powers);
        Transmission transmission(options, gc_generator);
        transmission.generate_samples(options.prns);
        wsgc_complex *faded_source_samples = transmission.get_samples();
        unsigned int nb_faded_source_samples = transmission.get_nb_samples();

        if (faded_source_samples)
        {
            if (options.transmission_scheme == Options::OptionTrans_MFSK)
            {
                Reception_MFSK reception_MFSK(options, gc_generator);
                reception_MFSK.message_processing(faded_source_samples, nb_faded_source_samples);
            }
            else if (options.transmission_scheme == Options::OptionTrans_WSGC)
            {
                Reception_WSGC reception_WSGC(options, gc_generator);

                if (options.simulate_training)
                {
                    reception_WSGC.training_processing(faded_source_samples, nb_faded_source_samples);
                }
                else
                {
                    reception_WSGC.message_processing(faded_source_samples, nb_faded_source_samples);
                }
            }
            else if (options.transmission_scheme == Options::OptionTrans_WSGCE)
            {
                Reception_WSGCE reception_WSGCE(options, gc_generator);

                if (options.simulate_training)
                {
                    std::cerr << "Simulating training sequence is not implemented" << std::endl;
                }
                else
                {
                    reception_WSGCE.message_processing(faded_source_samples, nb_faded_source_samples);
                }
            }
            else if (options.transmission_scheme == Options::OptionTrans_WSGCO)
            {
                Reception_WSGCO reception_WSGCO(options, gc_generator);
                
                if (options.simulate_training)
                {
                    reception_WSGCO.training_processing(faded_source_samples, nb_faded_source_samples);
                }
                else
                {
                    reception_WSGCO.message_processing(faded_source_samples, nb_faded_source_samples);
                }
            }
            else if (options.transmission_scheme == Options::OptionTrans_WSGCD)
            {
                Reception_WSGCD reception_WSGCD(options, gc_generator);
                
                if (options.simulate_training)
                {
                    reception_WSGCD.training_processing(faded_source_samples, nb_faded_source_samples);
                }
                else
                {
                    reception_WSGCD.message_processing(faded_source_samples, nb_faded_source_samples);
                }
            }
            else
            {
                std::cout << "Unknown or undefined transmission scheme" << std::endl;
            }
        }
        else
        {
            std::cout << "No samples" << std::endl;
        }

        return 0;
    }
    else
    {
        std::cout << "Incorrect options. Please correct and re-submit" << std::endl;
        return -1;
    }
}
