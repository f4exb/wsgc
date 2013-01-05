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
#include "FadingModelClarke.h"
#include <time.h>
#include <assert.h>
#include <math.h>


FadingModelClarke::FadingModelClarke(unsigned int nb_paths, wsgc_float spread_frequency, wsgc_float f_sampling) :
    FadingModel::FadingModel(f_sampling, true),
    _nb_paths(nb_paths),
    _spread_frequency(spread_frequency),
    _unif(0.0, 2*M_PI),
    _global_sample_index(0)
{
    assert(nb_paths > 0);
    
    for (unsigned int i=0; i<nb_paths; i++)
    {
        alpha_factors.push_back(_unif(_random_engine));
        beta_factors.push_back(_unif(_random_engine));
        theta_factors.push_back(_unif(_random_engine));        
    }
}


FadingModelClarke::~FadingModelClarke()
{
}


void FadingModelClarke::apply_fading(const wsgc_complex *samples_in, wsgc_complex *samples_out, unsigned int samples_length)
{
    wsgc_complex h;
    static const wsgc_float inv_sqrt_2 = 1.0/sqrt(2.0);
    
    for (unsigned int sample_i=0; sample_i < samples_length; sample_i++, _global_sample_index++)
    {
        h = (0.0, 0.0); // fading term at one sample
        
        for (unsigned int path_i=0; path_i < _nb_paths; path_i++)
        {
            wsgc_float x = cos(((2.0*path_i)*M_PI+theta_factors[path_i])/(4*_nb_paths));
            h += (cos(2.0*M_PI*_spread_frequency*x*(_global_sample_index/_f_sampling)+alpha_factors[path_i]),
                  sin(2.0*M_PI*_spread_frequency*x*(_global_sample_index/_f_sampling)+beta_factors[path_i]));
        }
        
        //h /= 1.0/sqrt(_nb_paths);
        h /= (_nb_paths*inv_sqrt_2);
        samples_out[sample_i] = h*samples_in[sample_i]; // multiply by fading term
    }
}
     
     
void FadingModelClarke::print_fading_data(std::ostringstream& os) const
{
    os << "Clarke [nb_paths=" << _nb_paths << ",FDoppler=" << _spread_frequency << "]";
}
     
