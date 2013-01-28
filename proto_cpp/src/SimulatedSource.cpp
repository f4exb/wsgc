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
#include "SimulatedSource.h"
#include "GoldCodeGenerator.h"
#include <iostream>
#include <time.h>
#include <math.h>
#include <assert.h>

SimulatedSource::SimulatedSource(GoldCodeGenerator& gc_generator, std::vector<unsigned int>& prn_list, wsgc_float f_sampling, wsgc_float f_chip, wsgc_float f_tx, unsigned int code_shift,
                                 unsigned int prns_per_symbol, wsgc_float start_phase) :
    _noiseDistributionUnit(0.0,1.0),
    _localOscillator(f_sampling, gc_generator.get_nb_code_samples(f_sampling, f_chip), start_phase),
    _gc_generator(gc_generator),
    _code_modulator(0),
    _prn_list(prn_list),
    _f_sampling(f_sampling),
    _f_chip(f_chip),
    _f_tx(f_tx),
    _code_shift(code_shift % gc_generator.get_nb_code_samples(f_sampling, f_chip)),
    _prns_per_symbol(prns_per_symbol),
    _snr_db(0),
    _make_noise(false),
    _start_phase(start_phase),
    _nb_code_samples(gc_generator.get_nb_code_samples(f_sampling, f_chip)),
    _nb_samples(gc_generator.get_nb_code_samples(f_sampling, f_chip)*prns_per_symbol*prn_list.size() + code_shift), // symbol samples + shift preamble zeroes
    _serviced_samples_index(0)
{
    _randomEngine.seed(time(0));
    _samples = (wsgc_complex *) WSGC_FFTW_MALLOC(_nb_samples*sizeof(wsgc_fftw_complex)); 
}    
     
SimulatedSource::~SimulatedSource()
{
    WSGC_FFTW_FREE(_samples);
}

void SimulatedSource::create_samples()
{
    wsgc_fftw_complex *samples_fftw = reinterpret_cast<wsgc_fftw_complex *>(_samples);
    const wsgc_complex *lo_samples;
    unsigned int sample_i = 0;
    
    assert(_code_modulator != 0);
    
    // preamble zeros for shift
    for (; sample_i < _code_shift; sample_i++)
    {
        samples_fftw[sample_i][0] = 0.0; // real (I)
        samples_fftw[sample_i][1] = 0.0; // imaginary (Q)
    }
    
    // symbols samples
    std::vector<char> code;
    std::vector<unsigned int>::const_iterator prn_it = _prn_list.begin();
    const std::vector<unsigned int>::const_iterator prns_end = _prn_list.end();
    
    for (; prn_it < prns_end; ++prn_it)
    {
        _gc_generator.make_code(code, *prn_it);
    
        for (unsigned int n=0; n < _prns_per_symbol; n++)
        {
            _localOscillator.make_next_samples(_f_tx);
            lo_samples = _localOscillator.get_samples();
            const wsgc_fftw_complex *lo_samples_fftw = reinterpret_cast<const wsgc_fftw_complex *>(lo_samples);
            
            _code_modulator->modulate(lo_samples_fftw, &samples_fftw[sample_i], code, _f_sampling, _f_chip);
            
            sample_i += _nb_code_samples;
        }
    }
}


bool SimulatedSource::get_next_code_samples(wsgc_complex **samples) 
{
    if (_serviced_samples_index + _nb_code_samples < _nb_samples)
    {
        *samples = &_samples[_serviced_samples_index];
        _serviced_samples_index += _nb_code_samples;
        return true;
    }
    else
    {
        return false;
    }
}
