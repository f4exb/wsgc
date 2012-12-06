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
      
     MultipleChannelCorrelator
      
     Takes samples for the length of a PRN and processes it in order to:
     - process correlation over all PRNs symbols in the alphabet for one code sequence
     - do the frequency tracking
     
     Used with BPSK complex signals
     
*/

#include "MultiplePrnCorrelator_FreqDep.h"
#include "GoldCodeGenerator.h"
#include "WsgcUtils.h"
#include <string.h>

#ifdef WSGC_VERBOSE_CORRELATOR
#include <iostream>
#include <iomanip>
#endif

// frequency tracking constants
const wsgc_float MultiplePrnCorrelator_FreqDep::lock_free_delta_phase_threshold = 0.1;     // no attempt to correct if phase shift from one round to the other is lower
const wsgc_float MultiplePrnCorrelator_FreqDep::unlock_delta_phase_threshold = M_PI / 8.0; // declare frequency locked if phase shift from one round to the other is lower
const wsgc_float MultiplePrnCorrelator_FreqDep::ftrack_noise_margin_limit = 1.4;           // reject phase shift sensing if noise margin too small (below this ratio)

MultiplePrnCorrelator_FreqDep::MultiplePrnCorrelator_FreqDep(GoldCodeGenerator& gc_generator, LocalCodesFFT_Host& local_codes, wsgc_float f_sampling, wsgc_float f_chip, wsgc_float f_rx_init, unsigned int prn_per_symbol,
                                                           unsigned int prn_per_symbol_i_init, unsigned int frequency_tracking_averaging_limit, bool file_debugging) :
    _code_modulator(),
    MultiplePrnCorrelator::MultiplePrnCorrelator(local_codes, gc_generator, f_sampling, f_chip, prn_per_symbol, prn_per_symbol_i_init, file_debugging),
    _f_rx(f_rx_init),
    _local_oscillator(f_sampling, gc_generator.get_nb_code_samples(f_sampling,f_chip)),
    _frequency_tracking_averaging_limit(frequency_tracking_averaging_limit)
{
}


MultiplePrnCorrelator_FreqDep::~MultiplePrnCorrelator_FreqDep()
{
}


void MultiplePrnCorrelator_FreqDep::set_source_block(wsgc_fftw_complex *fft_source_block)
{
	MultiplePrnCorrelator::set_source_block(fft_source_block);
    
    // mix with LO to obtain zero IF
    _local_oscillator.make_next_samples(-_f_rx);
    const wsgc_complex *lo_source_block = _local_oscillator.get_samples();
    wsgc_complex *fft_source_block_x = reinterpret_cast<wsgc_complex *>(_fft_source_block);

    for (unsigned int i=0; i<_fft_N; i++)
    {
        fft_source_block_x[i] *= lo_source_block[i];
    }
}
   

void MultiplePrnCorrelator_FreqDep::peak_and_track()
{
    wsgc_float phase_at_max;
    
    // calculate phase at maximum
    phase_at_max = std::arg(_last_prns_correlation[_prn_index_max*_fft_N + _shift_index_max]);
#ifdef WSGC_VERBOSE_CORRELATOR   
    wsgc_float magnitude_at_max;
    WsgcUtils::magnitude_estimation(&_last_prns_correlation[_prn_index_max*_fft_N + _shift_index_max], &magnitude_at_max);
#endif
    
    // frequency tracking
    if (_tracking_count > 0)
    {
        wsgc_float delta_phase = phase_at_max - _last_phase_at_max;
        
#ifdef WSGC_VERBOSE_CORRELATOR   
        std::cout << _last_symbol_count << ":" << _last_prn_per_symbol_i 
                  << "> Phase at max [" << _prn_index_max << "," << _shift_index_max << "] = " << phase_at_max 
                  << " Delta: " << delta_phase 
                  << " Mag: " << magnitude_at_max  
                  << " Mag Max: " << _magnitude_max  
                  << " Noise Max: " << _noise_mag_max << std::endl;
#endif            
        if (_magnitude_max / _noise_mag_max < ftrack_noise_margin_limit)
        {
#ifdef WSGC_VERBOSE_CORRELATOR   
            std::cout << "-    rejected ftrack: 6. noise margin too small" << std::endl;
#endif            
        }
        else if (!(_last_prn_index_max == _prn_index_max))
        {
#ifdef WSGC_VERBOSE_CORRELATOR   
            std::cout << "-    rejected ftrack: 2. selected prn change" << std::endl;
#endif            
        }
        else if (!((_last_magnitude_max < _magnitude_max)))
        {
#ifdef WSGC_VERBOSE_CORRELATOR   
            std::cout << "-    rejected ftrack: 3. not increasing correlation peak value" << std::endl;
#endif            
        }
        else if (_shift_index_max != _max_shift_hit)
        {
#ifdef WSGC_VERBOSE_CORRELATOR   
            std::cout << "-    rejected ftrack: 4. current correlation peak shift is not the most hit " << _shift_index_max << ":" << _max_shift_hit << std::endl;
#endif            
        }
        else if (_last_shift_index_max != _shift_index_max)
        {
#ifdef WSGC_VERBOSE_CORRELATOR   
            std::cout << "-    rejected ftrack: 5. correlation peak shift has changed" << std::endl;
#endif            
        }
        else // valid candidate values for frequency tracking (but we're not there yet!)
        {
#ifdef WSGC_VERBOSE_CORRELATOR   
            wsgc_float delta_frequency = (delta_phase/(2*M_PI)) / _code_period;
#endif            
            if (_frequency_lock_done) // frequency lock was done last valid round
            {
#ifdef WSGC_VERBOSE_CORRELATOR   
                std::cout << "+    skipped ftrack: " << std::setw(6) << std::setprecision(2) << delta_phase << " (" << std::setw(7) << std::setprecision(3) << delta_frequency
                    << "): lock just done skip phase averaging" << std::endl;
#endif            
                _frequency_lock_done = false;
            }
            else if ((delta_phase > -M_PI) and (delta_phase < M_PI))
            {
                _delta_phase_nb += 1;
                _delta_phase_sum += delta_phase;
                wsgc_float phase_average = _delta_phase_sum / _delta_phase_nb;
                _frequency_locked = (fabs(phase_average) < unlock_delta_phase_threshold);
#ifdef WSGC_VERBOSE_CORRELATOR   
                std::cout << "+    VALID ftrack: dPh:" << delta_phase << " Avg:" << phase_average << " dF:" << delta_frequency << " Flock: " << (_frequency_locked ? "yes" : "no") << std::endl;
#endif          
                if (_delta_phase_nb % _frequency_tracking_averaging_limit == 0)
                {
                    adjust_frequency(phase_average);
                }
            }
            else
            {
#ifdef WSGC_VERBOSE_CORRELATOR   
                std::cout << "+    skipped ftrack: delta_phase out of [-pi,pi]:" << delta_phase << std::endl;
#endif          
            }
        }
    }
    else
    {
#ifdef WSGC_VERBOSE_CORRELATOR 
        std::cout << _last_symbol_count << ":" << _last_prn_per_symbol_i 
                  << "> Phase at max [" << _prn_index_max << "," << _shift_index_max << "] = " << phase_at_max 
                  << " Mag: " << magnitude_at_max << std::endl;
        std::cout << "-    rejected ftrack: 1. first time round" << std::endl;
#endif
    }
    
    // record for tracking
    MultiplePrnCorrelator::peak_and_track();

    // BPSK specific (frequency tracking)
    _last_phase_at_max = phase_at_max;
    _last_frequency_locked = _frequency_locked;
    _last_f_rx = _f_rx;
}


void MultiplePrnCorrelator_FreqDep::adjust_frequency(wsgc_float delta_phase)
{
    if (fabs(delta_phase) > lock_free_delta_phase_threshold)
    {
        wsgc_float delta_frequency = (delta_phase/(2*M_PI)) / _code_period;
#ifdef WSGC_VERBOSE_CORRELATOR   
        std::cout << "+    -> Frequency adjust, delta f:" << delta_frequency << " New Frx:" <<  _f_rx+delta_frequency << std::endl;
#endif
        _f_rx += delta_frequency;
        // reset averaging
        _delta_phase_nb = 0; 
        _delta_phase_sum = 0.0;
        //_frequency_locked = true; // should be decided by phase measurement
        _frequency_lock_done = true;
    }
    else
    {
#ifdef WSGC_VERBOSE_CORRELATOR   
        std::cout << "+    -> Delta phase within tolerance, no frequency adjust" <<  std::endl;
#endif
    }
}
