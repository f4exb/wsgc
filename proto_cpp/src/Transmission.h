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
      
     Transmission

     Handles transmission operation including channel simulation
     
*/

#ifndef __TRANSMISSION_H__
#define __TRANSMISSION_H__

#include "WsgcTypes.h"
#include <vector>

class Options;
class GoldCodeGenerator;

class Transmission
{
public:
    Transmission(Options& _options,
        const GoldCodeGenerator& _gc_generator);
    
    ~Transmission();
    
    wsgc_complex *get_samples()
    {
        return signal_samples;
    }
    
    unsigned int get_nb_samples() const
    {
        return nb_signal_samples;
    }
    
    void generate_samples();
    void normalize();
    
protected:
    void apply_fir(wsgc_complex *inout, unsigned int& nb_samples, const std::vector<wsgc_float>& fir_coef);
    void apply_channel();
    void apply_fec();

    Options& options; //!< Reference to the options
    const GoldCodeGenerator& gc_generator; //!< Reference to the Gold Code generator
    wsgc_complex *signal_samples; //!< Signal samples
    unsigned int nb_signal_samples; //!< Number of signal samples generated
};

#endif
