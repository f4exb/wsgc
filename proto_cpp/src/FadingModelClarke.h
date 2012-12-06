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
      
     Applies a fading model to the test signal. It follows the model known as the 
     Clarke model and implemented using Matlab by Athuranathan Viswanathan at:
     http://gaussianwaves.blogspot.com/2011/05/simulation-of-rayleigh-fading-clarkes.html
      
*/
#ifndef __FADING_MODEL_CLARKE_H__
#define __FADING_MODEL_CLARKE_H__

#include "WsgcTypes.h"
#include "FadingModel.h"
#include <tr1/random> // c++0x extension (needs -std=c++0x to compile)
#include <vector>


class FadingModelClarke : public FadingModel
{
    public:
        FadingModelClarke(unsigned int nb_paths, wsgc_float spread_frequency, wsgc_float f_sampling);
        virtual ~FadingModelClarke();
        virtual void apply_fading(const wsgc_complex *samples_in, wsgc_complex *samples_out, unsigned int samples_length);
        virtual void print_fading_data(std::ostringstream& os) const;
        
        virtual unsigned int get_output_size(unsigned int input_size) const
        {
            return input_size;
        }
        
    private:      
        unsigned int _nb_paths; // M
        wsgc_float _spread_frequency; // fd
        std::tr1::uniform_real<wsgc_float> _unif;
        std::vector<wsgc_float> alpha_factors;
        std::vector<wsgc_float> beta_factors;
        std::vector<wsgc_float> theta_factors;
        unsigned int _global_sample_index;
};

#endif // __FADING_MODEL_CLARKE_H__
