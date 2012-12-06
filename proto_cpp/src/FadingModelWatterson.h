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
#ifndef __FADING_MODEL_WATTERSON_H__
#define __FADING_MODEL_WATTERSON_H__

#include "WsgcTypes.h"
#include "FadingModel.h"
#include "Path.h"
#include <vector>

class FadingModelWatterson : public FadingModel
{
    public:
        typedef struct FadingModelPath
        {
            wsgc_float delay;
            wsgc_float amplitude_factor;
            wsgc_float spread_frequency;
        } FadingModelPath_t;
    
        FadingModelWatterson(wsgc_float f_sampling);
        virtual ~FadingModelWatterson();
        virtual void apply_fading(const wsgc_complex *samples_in, wsgc_complex *samples_out, unsigned int samples_length);
        virtual void print_fading_data(std::ostringstream& os) const;
        virtual unsigned int get_output_size(unsigned int input_size) const;
        void add_path_description(wsgc_float delay, wsgc_float amplitude_factor, wsgc_float spread_frequency);
        void calculate_paths_description();
        
    private:
        std::vector<FadingModelPath_t> _fading_model_description;
        std::vector<Path> _paths;
        std::vector<unsigned int> _delay_samples;
        unsigned int _nb_paths;
        unsigned int _max_delay_samples;        
};

#endif // __FADING_MODEL_WATTERSON_H__
