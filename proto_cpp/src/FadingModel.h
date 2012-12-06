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
      
     Fading model class. 
     This parent class supports the AWGN addition which is common to all.
     An actual active Fading model is instantiated in derivated classes.
*/
#ifndef __FADING_MODEL_H__
#define __FADING_MODEL_H__

#include "WsgcTypes.h"
#include <sstream>
#include <tr1/random> // c++0x extension (needs -std=c++0x to compile)


class FadingModel
{
    public:
        FadingModel(wsgc_float f_sampling, bool active);
        
        virtual ~FadingModel();
        virtual void apply_fading(const wsgc_complex *samples_in, wsgc_complex *samples_out, unsigned int nb_samples) = 0;
        virtual void print_fading_data(std::ostringstream& os) const = 0;        
        virtual unsigned int get_output_size(unsigned int input_size) const = 0;
        
        bool is_fading_active()
        {
            return _active;
        }
        
        void set_verbose(bool verbose)
        {
            _verbose = verbose;
        }
        
        void set_random_seed(unsigned int seed)
        {
            _random_engine.seed(seed);
        }
        
        static wsgc_float get_mean_signal_power(wsgc_complex *samples, unsigned int nb_samples);
        void apply_awgn(wsgc_complex *samples, unsigned int nb_samples, unsigned int signal_shift, wsgc_float snr_db);
        
    protected:
        typedef std::tr1::ranlux64_base_01 RandomEngine; 
        typedef std::tr1::normal_distribution<wsgc_float> NormalDistribution; 

        RandomEngine _random_engine;
        NormalDistribution _AWGN_distribution_unit;

        wsgc_float _f_sampling;
        bool _active;
        bool _verbose;
};

#endif // __FADING_MODEL_WATTERSON_H__
