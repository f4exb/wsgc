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
      
     Simulated source

     Generates samples for a simulated noisy WSGC source
*/
#ifndef __SIMULATED_SOURCE_H__
#define __SIMULATED_SOURCE_H__

#include "WsgcTypes.h"
#include "ContinuousPhaseCarrier.h"
#include "CodeModulator_BPSK.h"
#include <vector>
#include <tr1/random> // c++0x extension (needs -std=c++0x to compile)

class GoldCodeGenerator;

class SimulatedSource
{
    public:
        SimulatedSource(GoldCodeGenerator& gc_generator, std::vector<unsigned int>& prn_list, wsgc_float f_sampling, wsgc_float f_chip, wsgc_float f_tx, unsigned int code_shift,
                        unsigned int prns_per_symbol=4, wsgc_float start_phase=0.0);
        ~SimulatedSource();
        void create_samples();

        void set_noisy(bool make_noise)
        {
            _make_noise = make_noise;
        }
        
        bool get_next_code_samples(wsgc_complex **samples);
        
        wsgc_complex *get_samples()
        {
            return _samples;
        }
        
        unsigned int get_nb_samples() // facility, the caller should normally know it
        {
            return _nb_samples;
        }
        
        void set_code_modulator(CodeModulator *code_modulator)
        {
            _code_modulator = code_modulator;
        }
        
        void set_prn_list(std::vector<unsigned int>& prn_list)
        {
        	_prn_list = prn_list;
        }

    private:
        typedef std::tr1::ranlux64_base_01 RandomEngine; 
        typedef std::tr1::normal_distribution<wsgc_float> NormalDistribution; 
        
        RandomEngine _randomEngine;
        NormalDistribution _noiseDistributionUnit;
        ContinuousPhaseCarrier _localOscillator;
        
        GoldCodeGenerator& _gc_generator;
        CodeModulator *_code_modulator;
        std::vector<unsigned int>& _prn_list;
        wsgc_float _f_sampling;
        wsgc_float _f_chip;
        wsgc_float _f_tx;
        unsigned int _code_shift;
        unsigned int _prns_per_symbol;
        wsgc_float _snr_db;
        bool _make_noise;
        wsgc_float _start_phase; 
        unsigned int _nb_code_samples;
        
        wsgc_complex *_samples;
        unsigned int _nb_samples;
        unsigned int _serviced_samples_index;        
};

#endif // __SIMULATED_SOURCE_H__
