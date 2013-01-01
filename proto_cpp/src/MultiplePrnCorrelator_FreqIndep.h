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
     - process correlation over all PRNs symbols in the alphabet for one code sequence (superclass)
     - tracking (superclass)
     
     Used with OOK complex signals
*/
#ifndef __MULTIPLE_CHANNEL_CORRELATOR_OOK_H__
#define __MULTIPLE_CHANNEL_CORRELATOR_OOK_H__

#include "WsgcTypes.h"
#include "ContinuousPhaseCarrier.h"
#include "CodeModulator_OOK.h"
#include "CorrelationRecord.h"
#include "MultiplePrnCorrelator.h"

class GoldCodeGenerator;

/**
 * \brief MultiplePrnCorrelator_FreqIndep engine to acquire and track many possible PRN(s)
 *
 * This class of Multiple PRN Correlator is used when the modulation has no or out of scope dependance on the receiving frequency. The correlation
 * can be done reliably even if the receiving frequency shifts by a certain tolerable amount. It is assumed that it stays within this amount throughout
 * the decoding process.
 *
 * \deprecated Use UnpilotedMessageCorrelator classes for frequency independant modulations
 *
*/
class MultiplePrnCorrelator_FreqIndep : public MultiplePrnCorrelator
{
    public:
        MultiplePrnCorrelator_FreqIndep(GoldCodeGenerator& gc_generator, LocalCodesFFT_Host& local_codes, wsgc_float f_sampling, wsgc_float f_chip, unsigned int prn_per_symbol=4,
                                      unsigned int prn_per_symbol_i_init=0, bool file_debugging=false);
        virtual ~MultiplePrnCorrelator_FreqIndep();
        
        virtual void set_source_block(wsgc_fftw_complex *fft_source_block)
        {
        	MultiplePrnCorrelator::set_source_block(fft_source_block);
        }

        virtual void peak_and_track()
        {
        	MultiplePrnCorrelator::peak_and_track();
        }
        
        virtual void reset_peak_tracking()
        {
        	MultiplePrnCorrelator::reset_peak_tracking();
        }
        
        virtual void get_correlation_record(CorrelationRecord& record) const
        {
        	MultiplePrnCorrelator::get_correlation_record(record);
        }
               
    private:
        CodeModulator_OOK _code_modulator;
};

#endif // __MULTIPLE_CHANNEL_CORRELATOR_BPSK_H__
