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
#ifndef __RECEPTION_H__
#define __RECEPTION_H__

#include "WsgcTypes.h"
#include "Options.h"
#include "GoldCodeGenerator.h"

class Reception
{
public:
    Reception(Options& _options, const GoldCodeGenerator& _gc_generator) :
        options(_options),
        gc_generator(_gc_generator),
        fft_N(gc_generator.get_nb_code_samples(options.f_sampling, options.f_chip))
    {}
    
    ~Reception()
    {}

    void demodulate_before_correlate(wsgc_complex *faded_source_samples, unsigned int nb_faded_source_samples);
    
protected:
    void generate_training_prn_list(std::vector<unsigned int>& prn_list);
    void generate_training_prn_list_unpiloted(std::vector<unsigned int>& prn_list);
    void generate_message_prn_list(std::vector<unsigned int>& prn_list);
    void generate_pilot_prn_list(std::vector<unsigned int>& prn_list, unsigned int pilot_prni);
    
    Options& options;
    const GoldCodeGenerator& gc_generator;
    unsigned int fft_N;
};

#endif // __RECEPTION_H__
