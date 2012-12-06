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
     - do peak estimation and possible tracking
     
*/

#include "MultiplePrnCorrelator.h"
#include "GoldCodeGenerator.h"
#include "WsgcUtils.h"
#include <string.h>

#ifdef WSGC_VERBOSE_CORRELATOR
#include <iostream>
#include <iomanip>
#endif

MultiplePrnCorrelator::MultiplePrnCorrelator(LocalCodesFFT_Host& local_codes, GoldCodeGenerator& gc_generator, wsgc_float f_sampling, wsgc_float f_chip, unsigned int prn_per_symbol,
                                                   unsigned int prn_per_symbol_i_init, unsigned int frequency_tracking_averaging_limit, bool file_debugging) :
    _f_sampling(f_sampling),
    _f_chip(f_chip),
    _prn_per_symbol(prn_per_symbol),
    _prn_per_symbol_i(prn_per_symbol_i_init % prn_per_symbol),
    _symbol_count(0),
    _fft_N(gc_generator.get_nb_code_samples(f_sampling,f_chip)),
    _local_codes(local_codes),
    _nb_prns(gc_generator.get_nb_codes()),
    _gc_generator(gc_generator),
    _code_period(gc_generator.get_code_length()/f_chip),
    _fft_source_block(0),
    _tracking_count(0),
    _file_debugging(file_debugging)
{
    _fft_sample_in = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_fftw_complex));
    _fft_sample_out = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_fftw_complex));
    _fft_sample_plan = WSGC_FFTW_PLAN(_fft_N, reinterpret_cast<wsgc_fftw_complex *>(_fft_sample_in), reinterpret_cast<wsgc_fftw_complex *>(_fft_sample_out), FFTW_FORWARD, FFTW_ESTIMATE);
    
    _ifft_code_in = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_fftw_complex));
    _ifft_code_out = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_fftw_complex));
    _ifft_code_plan = WSGC_FFTW_PLAN(_fft_N, reinterpret_cast<wsgc_fftw_complex *>(_ifft_code_in), reinterpret_cast<wsgc_fftw_complex *>(_ifft_code_out), FFTW_BACKWARD, FFTW_ESTIMATE);

    _last_prns_correlation = (wsgc_complex *) WSGC_FFTW_MALLOC(_nb_prns*_fft_N*sizeof(wsgc_fftw_complex));
    _prn_modules = new wsgc_float[_nb_prns*_fft_N*_prn_per_symbol];
    
    for (unsigned int i=0; i<_nb_prns*_fft_N*_prn_per_symbol; i++)
    {
        _prn_modules[i] = 0.0;
    }
    
    reset_peak_tracking();
    
#ifdef WSGC_VERBOSE_CORRELATOR       
    if (_file_debugging)
    {
        unsigned int code_length = (1<<_gc_generator.get_nb_stages())-1;
        unsigned int nb_samples_per_code = (unsigned int) ((_f_sampling / _f_chip) * ((1<<_gc_generator.get_nb_stages())-1));

        _debug_file.open("wsgc_corr_debug.log", std::ios::out | std::ios::app);
        _debug_file << std::endl << "******** New run ********" << std::endl;
        _debug_file << std::endl << "Run with options:" << std::endl;
        _debug_file << std::endl << "----------------" << std::endl;
        _debug_file << std::setiosflags(std::ios_base::fixed);
        _debug_file << "Sampling frequency ........: " << std::setw(8) << std::setprecision(1) << std::right << _f_sampling << std::endl;
        _debug_file << "Chip frequency ............: " << std::setw(8) << std::setprecision(1) << std::right << _f_chip << std::endl;
        _debug_file << "Code length ...............: " << std::setw(6) << std::right << code_length << std::endl;
        _debug_file << "Samples/code = FFT size ...: " << std::setw(6) << std::right << nb_samples_per_code << std::endl;
        _debug_file << "PRNs per symbol ...........: " << std::setw(6) << std::right << _prn_per_symbol << std::endl;
        _debug_file << "Nb message symbols ........: " << std::setw(6) << std::right << _gc_generator.get_nb_message_codes() << std::endl;
        _debug_file << "Nb service symbols ........: " << std::setw(6) << std::right << _gc_generator.get_nb_service_codes() << std::endl; 
        _debug_file << std::endl;
        _debug_file << "Generator polynomials:" << std::endl;
        std::ostringstream os;
        os << "G1 = ";
        WsgcUtils::print_polynomial(_gc_generator.get_nb_stages(), _gc_generator.get_g1_powers(), os); os << std::endl;
        os << "G2 = ";
        WsgcUtils::print_polynomial(_gc_generator.get_nb_stages(), _gc_generator.get_g2_powers(), os); os << std::endl;
        _debug_file << os.str() << std::endl << std::endl;
        std::cout << "wsgc_corr_debug.log opened" << std::endl;
    }
#endif
}

MultiplePrnCorrelator::~MultiplePrnCorrelator()
{
    WSGC_FFTW_DESTROY_PLAN(_fft_sample_plan);
    WSGC_FFTW_FREE(_fft_sample_in);
    WSGC_FFTW_FREE(_fft_sample_out);

    WSGC_FFTW_DESTROY_PLAN(_ifft_code_plan);
    WSGC_FFTW_FREE(_ifft_code_in);
    WSGC_FFTW_FREE(_ifft_code_out);
    
#ifdef WSGC_VERBOSE_CORRELATOR       
    if (_file_debugging)
    {
        _debug_file.close();
        std::cout << "wsgc_corr_debug.log closed" << std::endl;
    }
#endif
}


void MultiplePrnCorrelator::set_source_block(wsgc_fftw_complex *fft_source_block)
{
    _fft_source_block = fft_source_block;
    wsgc_complex *fft_source_block_x = reinterpret_cast<wsgc_complex *>(_fft_source_block);
    
#ifdef WSGC_VERBOSE_CORRELATOR       
    if (_file_debugging)
    {
        std::ostringstream os;
        os << "=== (" << _symbol_count << "," << _prn_per_symbol_i << ") Source data ===" << std::endl;
        print_array_x( _fft_N, fft_source_block_x, os, 8);
        _debug_file << os.str() << std::endl;
    }
#endif
}
    
    
void MultiplePrnCorrelator::make_correlation()
{
    wsgc_complex *fft_source_block_x = reinterpret_cast<wsgc_complex *>(_fft_source_block);
    wsgc_float magnitude_sum;
    std::map<unsigned int, unsigned int, std::greater<unsigned int> >::iterator shift_occurences_it;
    const wsgc_complex *local_code_block;
    
    assert(_fft_source_block != 0);

    // reset correlation run variables
    _prn_index_max = 0;
    _shift_index_max = 0;
    _magnitude_max = 0.0;
    _noise_mag_sum = 0.0;
    _noise_mag_max = 0.0;
    
#ifdef WSGC_VERBOSE_CORRELATOR       
    if (_file_debugging)
    {
        std::ostringstream os;
        os << "=== (" << _symbol_count << "," << _prn_per_symbol_i << ") FFT input data ===" << std::endl;
        print_array_x( _fft_N, fft_source_block_x, os, 8);
        _debug_file << os.str() << std::endl;
    }
#endif

    // do the FFT of the zero IF input signal
    memcpy(_fft_sample_in, _fft_source_block, _fft_N*sizeof(wsgc_fftw_complex));
    WSGC_FFTW_EXECUTE(_fft_sample_plan);
    
#ifdef WSGC_VERBOSE_CORRELATOR       
    if (_file_debugging)
    {
        std::ostringstream os;
        os << "=== (" << _symbol_count << "," << _prn_per_symbol_i << ") FFT output data ===" << std::endl;
        print_array_x( _fft_N, _fft_sample_out, os, 8);
        _debug_file << os.str() << std::endl;
    }
#endif

    // do the IFFT of the product of the FFT of the zero IF input signal with the complex conjugate of the FFT of the local code for each local code
    for (unsigned int prn_i=0; prn_i<_nb_prns; prn_i++)
    {
        local_code_block = _local_codes.get_local_code(prn_i); // get the complex conjugate of the FFT of the local code

#ifdef WSGC_VERBOSE_CORRELATOR       
        if ((_file_debugging) && (prn_i == 0))
        {
            std::ostringstream os;
            os << "=== (" << _symbol_count << "," << _prn_per_symbol_i << ") LC block PRN 0 ===" << std::endl;
            print_array_x( _fft_N, local_code_block, os, 8);
            _debug_file << os.str() << std::endl;
        }
#endif
        
        // do the product
        for (unsigned int fft_i=0; fft_i<_fft_N; fft_i++)
        {
            _ifft_code_in[fft_i] = _fft_sample_out[fft_i] * local_code_block[fft_i]; 
        }
        
#ifdef WSGC_VERBOSE_CORRELATOR       
        if ((_file_debugging) && (prn_i == 0))
        {
            std::ostringstream os;
            os << "=== (" << _symbol_count << "," << _prn_per_symbol_i << ") IFFT input data PRN 0 ===" << std::endl;
            print_array_x( _fft_N, _ifft_code_in, os, 8);
            _debug_file << os.str() << std::endl;
        }
#endif

        WSGC_FFTW_EXECUTE(_ifft_code_plan);
        
#ifdef WSGC_VERBOSE_CORRELATOR       
        if ((_file_debugging) && (prn_i == 0))
        {
            std::ostringstream os;
            os << "=== (" << _symbol_count << "," << _prn_per_symbol_i << ") IFFT output data PRN 0 ===" << std::endl;
            print_array_x( _fft_N, _ifft_code_out, os, 8);
            _debug_file << os.str() << std::endl;
        }
#endif

        // store output
        for (unsigned int ifft_i=0; ifft_i < _fft_N; ifft_i++)
        {
            _last_prns_correlation[prn_i*_fft_N + ifft_i] = _ifft_code_out[ifft_i]; // re + im
            WsgcUtils::magnitude_estimation(&_ifft_code_out[ifft_i], &_prn_modules[prn_i*_fft_N*_prn_per_symbol+ifft_i*_prn_per_symbol+_prn_per_symbol_i]); // module

            // find magnitude grand maximum location and value
            
            /* Possible optimization with hardcoding (template worker class) instead of loop on [1,_prn_per_symbol[
            magnitude_sum = _prn_modules[prn_i*_fft_N*_prn_per_symbol+ifft_i*_prn_per_symbol]; 
            magnitude_sum += _prn_modules[prn_i*_fft_N*_prn_per_symbol+ifft_i*_prn_per_symbol+1];
            magnitude_sum += _prn_modules[prn_i*_fft_N*_prn_per_symbol+ifft_i*_prn_per_symbol+2];
            magnitude_sum += _prn_modules[prn_i*_fft_N*_prn_per_symbol+ifft_i*_prn_per_symbol+3];
            */
            
            magnitude_sum = _prn_modules[prn_i*_fft_N*_prn_per_symbol+ifft_i*_prn_per_symbol]; 
            
            for (unsigned int i = 1; i < _prn_per_symbol; i++)
            {
                magnitude_sum += _prn_modules[prn_i*_fft_N*_prn_per_symbol+ifft_i*_prn_per_symbol+i];
            }
            
            if (prn_i == _gc_generator.get_nb_message_codes()) // noise PRN is the first service PRN hence right after last message PRN
            {
                if (magnitude_sum > _noise_mag_max)
                {
                    _noise_mag_max = magnitude_sum;
                }
                _noise_mag_sum += magnitude_sum;
            }
            
            if (magnitude_sum > _magnitude_max)
            {
                _magnitude_max = magnitude_sum;
                _prn_index_max = prn_i;
                _shift_index_max = ifft_i;
            }
        }
    }
    
    // store current shift maximum in shifts dictionnary
    shift_occurences_it = _shift_occurences.find(_shift_index_max);
    
    if (shift_occurences_it == _shift_occurences.end())
    {
        _shift_occurences[_shift_index_max] = 1;
    }
    else
    {
        _shift_occurences[_shift_index_max] += 1;
    }
    
    // store most occurent shift
    if (_shift_occurences[_shift_index_max] > _max_shift_count)
    {
        _max_shift_count = _shift_occurences[_shift_index_max];
        _max_shift_hit = _shift_index_max;
    }    
    
    // updates for next run
    _last_prn_per_symbol_i = _prn_per_symbol_i;
    _last_symbol_count = _symbol_count;
    
    if (_prn_per_symbol_i < _prn_per_symbol-1)
    {
        _prn_per_symbol_i += 1;
    }
    else
    {
        _symbol_count++;
        _prn_per_symbol_i = 0;
    }
    
    _fft_source_block = 0; // indicate input block has been consumed
}


void MultiplePrnCorrelator::peak_and_track()
{
    // record for tracking
    _last_prn_index_max = _prn_index_max;
    _last_magnitude_max = _magnitude_max;
    _last_shift_index_max = _shift_index_max;
    _last_noise_max = _noise_mag_max;
    _last_noise_avg = _noise_mag_sum / _fft_N;
    _tracking_count++;
}


void MultiplePrnCorrelator::print_prn_correlation_module(unsigned int prn_i, unsigned int prn_shift_i, unsigned int fft_length, unsigned int prn_per_symbol, wsgc_float *module, std::ostringstream& os)
{
    print_array_xzy<wsgc_float>(_nb_prns, _fft_N, _prn_per_symbol, module, os, 8);
}


void MultiplePrnCorrelator::print_prn_complex_correlation(unsigned int prn_i, unsigned int fft_length, wsgc_complex *correlation, std::ostringstream& os)
{
    wsgc_float magnitude;
    
    for (unsigned int ifft_i=0; ifft_i<fft_length; ifft_i++)
    {
        WsgcUtils::magnitude_estimation(&correlation[prn_i*fft_length + ifft_i], &magnitude);
        os << "[" << ifft_i << "] = " << std::fixed << std::setw(11) << std::setprecision(8) << correlation[prn_i*fft_length + ifft_i] << " -> (" << magnitude << "," << std::arg(correlation[prn_i*fft_length + ifft_i]) << ")" << std::endl;
    }
}

