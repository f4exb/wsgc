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
#ifndef __MULTIPLE_CHANNEL_CORRELATOR_H__
#define __MULTIPLE_CHANNEL_CORRELATOR_H__

#include "WsgcTypes.h"
#include "ContinuousPhaseCarrier.h"
#include "LocalCodesFFT_Host.h"
#include "CorrelationRecord.h"
#include <map>
#include <fftw3.h>
#include <fstream>

class GoldCodeGenerator;

class MultiplePrnCorrelator
{
    public:
        MultiplePrnCorrelator(LocalCodesFFT_Host& local_codes, GoldCodeGenerator& gc_generator, wsgc_float f_sampling, wsgc_float f_chip, unsigned int prn_per_symbol=4,
                                  unsigned int prn_per_symbol_i_init=0, unsigned int frequency_tracking_averaging_limit=4, bool file_debugging=false);
        virtual ~MultiplePrnCorrelator();
        
        virtual void set_source_block(wsgc_fftw_complex *fft_source_block);
        void make_correlation();
        virtual void peak_and_track();
        
        unsigned int get_tracking_count() const
        {
            return _tracking_count;
        }
        
        const std::map<unsigned int, unsigned int>& get_shift_occurences() const
        {
            return _shift_occurences;
        }
        
        const wsgc_complex *get_last_prns_correlation() const
        {
            return _last_prns_correlation;
        }
        
        virtual void reset_peak_tracking()
        {
            _shift_occurences.clear();
            _max_shift_count = 0;
            _tracking_count = 0;
        }
        
        virtual void get_correlation_record(CorrelationRecord& record) const
        {
            record.prn_per_symbol_index = _last_prn_per_symbol_i;
            record.prn_index_max = _last_prn_index_max;
            record.magnitude_max = _last_magnitude_max;
            record.shift_index_max = _last_shift_index_max;
            record.noise_max = _last_noise_max;
            record.noise_avg = _last_noise_avg;
            // Overridden in case of BPSK:
            record.phase_at_max = 0.0;
            record.f_rx = 0.0;
            record.frequency_locked = true;            
        }
        
        unsigned int get_max_shift_hit()
        {
            return _max_shift_hit;
        }
        
    protected:
        wsgc_float _f_sampling;
        wsgc_float _f_chip;
        unsigned int _prn_per_symbol;
        unsigned int _prn_per_symbol_i;
        unsigned int _symbol_count;
        unsigned int _fft_N;
        unsigned int _nb_prns;
        GoldCodeGenerator& _gc_generator;
        wsgc_float _code_period;
        
        wsgc_fftw_complex *_fft_source_block;
        wsgc_fftw_plan _fft_sample_plan;
        wsgc_complex *_fft_sample_in, *_fft_sample_out;
        wsgc_fftw_plan _ifft_code_plan;
        wsgc_complex *_ifft_code_in, *_ifft_code_out;
        
        LocalCodesFFT_Host& _local_codes;
        
        wsgc_complex *_last_prns_correlation;
        wsgc_float *_prn_modules;
        
        unsigned int _max_shift_count;
        unsigned int _max_shift_hit;
        unsigned int _prn_index_max;
        unsigned int _shift_index_max;
        wsgc_float   _magnitude_max;
        unsigned int _tracking_count;
        
        unsigned int _last_prn_per_symbol_i;
        unsigned int _last_symbol_count;
        unsigned int _last_prn_index_max;
        wsgc_float   _last_magnitude_max;
        wsgc_float   _noise_mag_sum;
        wsgc_float   _noise_mag_max;
        unsigned int _last_shift_index_max;
        wsgc_float   _last_noise_max;
        wsgc_float   _last_noise_avg;

        std::map<unsigned int, unsigned int> _shift_occurences;
        
        bool _file_debugging;
        std::ofstream _debug_file;

        void adjust_frequency(wsgc_float phase_average);
        void print_prn_correlation_module(unsigned int prn_i, unsigned int prn_shift_i, unsigned int fft_length, unsigned int prn_per_symbol, wsgc_float *module, std::ostringstream& os);
        static void print_prn_complex_correlation(unsigned int prn_i, unsigned int fft_length, wsgc_complex *correlation, std::ostringstream& os);
};

#endif // __MULTIPLE_CHANNEL_CORRELATOR_BPSK_H__
