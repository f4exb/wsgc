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
      
     Applies a fading model to the test signal. Actual Fading model is instantiated in derivated classes.
     This parent class supports the AWGN addition which is common to all.
*/

#include "FadingModel.h"
#include <iostream>
#include <cmath>
#include <assert.h>


FadingModel::FadingModel(wsgc_float f_sampling, bool active) :
    _AWGN_distribution_unit(0.0,1.0),
    _f_sampling(f_sampling),
    _active(active),
    _verbose(false)
{
    _random_engine.seed(time(0));
}


FadingModel::~FadingModel() 
{}


void apply_fading(const wsgc_complex *samples_in, wsgc_complex *samples_out, unsigned int nb_samples)
{
    for (unsigned int sample_i=0; sample_i<nb_samples; sample_i++)
    {
        samples_out[sample_i] = samples_in[sample_i];
    }    
}


void FadingModel::print_fading_data(std::ostringstream& os) const
{
    os << "None";
}

        
wsgc_float FadingModel::get_mean_signal_power(wsgc_complex *samples, unsigned int nb_samples)
{
    wsgc_float signal_power = 0;
    
    for (unsigned int i=0; i<nb_samples; i++)
    {
        signal_power += std::norm(samples[i]);
    }
    
    return signal_power / nb_samples;
}


void FadingModel::apply_awgn(wsgc_complex *samples, unsigned int nb_samples, unsigned int signal_shift, wsgc_float snr_db)
{
    assert(signal_shift < nb_samples);
    
    wsgc_float mean_signal_power = get_mean_signal_power(&samples[signal_shift], nb_samples - signal_shift);
    wsgc_fftw_complex *samples_fftw = reinterpret_cast<wsgc_fftw_complex *>(samples);
    wsgc_float noise_amplitude = sqrt(mean_signal_power / pow(10.0, (snr_db/10.0)));
    
    std::cout << "Signal power: " << mean_signal_power << " Noise amplitude: " << noise_amplitude << std::endl;
    
    noise_amplitude /= sqrt(2.0); // because of complex (analytical) signal
    
    if (!_verbose)
    {
        for (unsigned int sample_i=0; sample_i<nb_samples; sample_i++)
        {
            //samples_fftw[sample_i][0] += noise_amplitude * _AWGN_distribution_unit(_random_engine);
            //samples_fftw[sample_i][1] += noise_amplitude * _AWGN_distribution_unit(_random_engine);
        	// Added normalization
            samples_fftw[sample_i][0] = (samples_fftw[sample_i][0] + noise_amplitude * _AWGN_distribution_unit(_random_engine)) / (noise_amplitude * sqrt(2.0));
            samples_fftw[sample_i][1] = (samples_fftw[sample_i][1] + noise_amplitude * _AWGN_distribution_unit(_random_engine)) / (noise_amplitude * sqrt(2.0));
        }
    }
    else
    {
        
        wsgc_complex noise_samples[nb_samples];
        wsgc_fftw_complex *noise_samples_fftw = reinterpret_cast<wsgc_fftw_complex *>(noise_samples);
        
        for (unsigned int sample_i=0; sample_i<nb_samples; sample_i++)
        {
            wsgc_float noise_i, noise_q;
            
            noise_i = noise_amplitude * _AWGN_distribution_unit(_random_engine);
            noise_q = noise_amplitude * _AWGN_distribution_unit(_random_engine);
            
            noise_samples_fftw[sample_i][0] = noise_i;
            noise_samples_fftw[sample_i][1] = noise_q;
            
            //samples_fftw[sample_i][0] += noise_i;
            //samples_fftw[sample_i][1] += noise_q;
        	// Added normalization
            samples_fftw[sample_i][0] = (samples_fftw[sample_i][0]+noise_i)/(noise_amplitude * sqrt(2.0));
            samples_fftw[sample_i][1] = (samples_fftw[sample_i][1]+noise_q)/(noise_amplitude * sqrt(2.0));
        }
        
        wsgc_float noise_power = get_mean_signal_power(noise_samples, nb_samples);
        
        std::cout << "Noise power: " << noise_power << " S/N (dB): " << 10.0*log10(mean_signal_power/noise_power) << std::endl;
        std::cout << "New signal power: " << get_mean_signal_power(&samples[signal_shift], nb_samples - signal_shift) << std::endl;
    }
}
        

