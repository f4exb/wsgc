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
#ifndef __MULTIPLE_PRN_CORRELATOR_BPSK_H__
#define __MULTIPLE_PRN_CORRELATOR_BPSK_H__

#include "WsgcTypes.h"
#include "ContinuousPhaseCarrier.h"
#include "CodeModulator_BPSK.h"
#include "CorrelationRecord.h"
#include "MultiplePrnCorrelator.h"

class GoldCodeGenerator;
class CodeModulator_BPSK;

/**
 * \brief MultiplePrnCorrelator_FreqDep engine to acquire and track many possible PRN(s)
 *
 * This class of Multiple PRN Correlator is used when the modulation implies a frequency dependance. The correlation can be done best when the
 * receiving frequency makes it match exactly at zero IF. As oscillators are never infinitely stable there should be a frequency tracking
 * mechanism.
 *
*/
class MultiplePrnCorrelator_FreqDep : public MultiplePrnCorrelator
{
    public:
        MultiplePrnCorrelator_FreqDep(GoldCodeGenerator& gc_generator, LocalCodesFFT_Host& local_codes, wsgc_float f_sampling, wsgc_float f_chip, wsgc_float f_rx_init, unsigned int prn_per_symbol=4,
                                      unsigned int prn_per_symbol_i_init=0, unsigned int frequency_tracking_averaging_limit=4, bool file_debugging=false);
        virtual ~MultiplePrnCorrelator_FreqDep();
        
        virtual void set_source_block(wsgc_fftw_complex *fft_source_block);
        virtual void peak_and_track();
        
        virtual void reset_peak_tracking()
        {
        	MultiplePrnCorrelator::reset_peak_tracking();
            // BPSK:
            _delta_phase_sum = 0.0;
            _delta_phase_nb = 0;
            _frequency_locked = false;
            _frequency_lock_done = false;
        }
        
        virtual void get_correlation_record(CorrelationRecord& record) const
        {
        	MultiplePrnCorrelator::get_correlation_record(record);
            // BPSK:
            record.phase_at_max = _last_phase_at_max;
            record.f_rx = _last_f_rx;
            record.frequency_locked = _last_frequency_locked;
        }
        
               
    private:
        wsgc_float _f_rx;
        ContinuousPhaseCarrier _local_oscillator;
        CodeModulator_BPSK _code_modulator;
        
        wsgc_float   _last_phase_at_max;
        wsgc_float   _last_f_rx;
        bool         _last_frequency_locked;

        unsigned int _frequency_tracking_averaging_limit;
        bool _frequency_lock_done;
        wsgc_float _delta_phase_sum;
        unsigned int _delta_phase_nb;
        bool _frequency_locked;
        
        static const wsgc_float unlock_delta_phase_threshold;
        static const wsgc_float lock_free_delta_phase_threshold;
        static const wsgc_float ftrack_noise_margin_limit;
        
        void adjust_frequency(wsgc_float phase_average);
};

#endif // __MULTIPLE_PRN_CORRELATOR_BPSK_H__
