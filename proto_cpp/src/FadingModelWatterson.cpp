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
      
     Applies a fading model to the test signal. It follows the model implemented
     in Moe Wheatley's (AE4JY) PathSim simulation program. It is known as the 
     Watterson's model. See: http://www.moetronix.com/ae4jy/pathsim.htm
      
     Input signal is duplicated in several paths. On each path a delay and attenuation
     factor is applied. Each path is then passed through a Doppler spread simulation
     consisting in the multiplication by a low pass filtered random AWG noise. The bandwidth
     of this lowpass filter determines the amount of Doppler spread.
*/
#include "FadingModelWatterson.h"
#include <iostream>
#include <iomanip>
#include <assert.h>

FadingModelWatterson::FadingModelWatterson(wsgc_float f_sampling) :
    FadingModel::FadingModel(f_sampling, true)
{
	/*
    // Presently works for sample rates multiple of 5^2, 5^3 or 5^4
    unsigned int f_sampling_int = int(f_sampling);
    assert ((f_sampling_int % (5*5*5*5) == 0) || (f_sampling_int % (5*5*5) == 0) || (f_sampling_int % (5*5) == 0));
    */
}

FadingModelWatterson::~FadingModelWatterson()
{
}


void FadingModelWatterson::add_path_description(wsgc_float delay, wsgc_float amplitude_factor, wsgc_float spread_frequency, wsgc_float offset_frequency)
{
    static const FadingModelPath_t tmp_path_description = {0.0, 0.0, 0.0};
    
    _fading_model_description.push_back(tmp_path_description);
    _fading_model_description.back().delay = delay;
    _fading_model_description.back().amplitude_factor = amplitude_factor;
    _fading_model_description.back().spread_frequency = spread_frequency;
    _fading_model_description.back().offset_frequency = offset_frequency;
    _nb_paths = _fading_model_description.size();
}


void FadingModelWatterson::apply_fading(const wsgc_complex *samples_in, wsgc_complex *samples_out, unsigned int samples_length)
{
    wsgc_complex path_sample;
    
    for (unsigned int sample_i = 0; sample_i < (samples_length+_max_delay_samples); sample_i++)
    {
        samples_out[sample_i] = (0.0, 0.0);
        
        for (unsigned int path_i = 0; path_i < _fading_model_description.size(); path_i++)
        {
            if (sample_i-_delay_samples[path_i] < samples_length)
            {
                if (sample_i < _delay_samples[path_i]) // signal has not arrived yet
                {
                    path_sample = (0.0, 0.0);
                }
                else
                {
                    _paths[path_i].CalcPathSample(&samples_in[sample_i-_delay_samples[path_i]], &path_sample);
                    path_sample *= _fading_model_description[path_i].amplitude_factor;
                }
            }
            else // signal is finished but last has not arrived yet
            {
                path_sample = (0.0, 0.0);
            }
            
            samples_out[sample_i] += path_sample;
        }     
    }
}

void FadingModelWatterson::print_fading_data(std::ostringstream& os) const
{
    os << std::setw(6) << std::setprecision(4) << "Watterson [";
    std::vector<FadingModelPath_t>::const_iterator it = _fading_model_description.begin();
    const std::vector<FadingModelPath_t>::const_iterator itEnd = _fading_model_description.end();
    
    for (;it != itEnd; ++it)
    {
        if (it == _fading_model_description.begin())
        {
            os << "(";
        }
        else
        {
            os << ",(";
        }
        
        os << it->delay << "," << it->amplitude_factor << "," << it->spread_frequency << "," << it->offset_frequency << ")";
    }    
    
    os << "]";
}

void FadingModelWatterson::calculate_paths_description()
{
    std::vector<FadingModelPath_t>::iterator it = _fading_model_description.begin();
    const std::vector<FadingModelPath_t>::iterator itEnd = _fading_model_description.end();
    wsgc_float max_delay = 0.0;
    
    for (;it != itEnd; ++it)
    {
        static const Path tmp_path(_f_sampling);
        _paths.push_back(tmp_path);
        _paths.back().InitPath(it->spread_frequency, it->offset_frequency, 1, _fading_model_description.size(), true);
        
        if (it == _fading_model_description.begin())
        {
            _delay_samples.push_back(0); // first path is direct path
        }
        else
        {
            _delay_samples.push_back(int(it->delay*_f_sampling)); // path delay in number of samples
        }
        
        if (it->delay > max_delay)
        {
            max_delay = it->delay;
        }        
    }
    
    _max_delay_samples = int(max_delay * _f_sampling);
}


unsigned int FadingModelWatterson::get_output_size(unsigned int input_size) const
{
    return input_size + _max_delay_samples;
}
        
