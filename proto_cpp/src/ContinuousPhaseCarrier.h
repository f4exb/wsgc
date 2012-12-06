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
      
     Continuous phase carrier
      
     Makes pure sine carrier samples of specified length. The next batch of samples
     is in continuous phase with the previous
*/
#ifndef __CONTINUOUS_PHASE_CARRIER_H__
#define __CONTINUOUS_PHASE_CARRIER_H__

#include "WsgcTypes.h"

/**
 * \brief Makes samples batches with continuous phase between batches
 */
class ContinuousPhaseCarrier
{
    public:
        /**
        * \brief Constructor
        * Makes a new ContinuousPhaseCarrier object
        * \param f_sampling sampling frequency
        * \param length number of samples per batch
        * \param init_phase initial phase
        */
        ContinuousPhaseCarrier(wsgc_float f_sampling, unsigned int length, wsgc_float init_phase = 0.0);
        virtual ~ContinuousPhaseCarrier();
        /*
         * Get a new batch of samples
         * \return A pointer to the first element of the array
         */
        const wsgc_complex *get_samples() const;
        /*
         * Prepare the next batch of samples
         * \param f Frequency for the next batch
         */
        void make_next_samples(wsgc_float f);
        /*
         * Set current phase
         * \param phase New phase
         */
        void set_phase(wsgc_float phase);
        
    private:
        wsgc_float _f_sampling;
        unsigned int _length;
        wsgc_complex *_samples;
        wsgc_float _phase_mod_2pi;
        wsgc_float _phase_intpart;
};

#endif // __CONTINUOUS_PHASE_CARRIER_H__
