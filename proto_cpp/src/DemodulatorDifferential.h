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
      
     DemodulatorDifferential
      
*/
#ifndef __DEMODULATOR_DIFFERENTIAL_H__
#define __DEMODULATOR_DIFFERENTIAL_H__

#include "Demodulator.h"
#include "WsgcTypes.h"
#include <vector>

/**
 * \brief Differentialy demodulate input samples. This is typically used for DBPSK
 */
class DemodulatorDifferential : public Demodulator
{
    public:
        /**
         * Constructs a differential demodulator
         * \param differential_length Number of samples in the differential delay
         */
        DemodulatorDifferential(unsigned int differential_length);
        
        virtual ~DemodulatorDifferential();

        /**
         * Demodulate samples in place
         * \param samples_inout Samples in input and output
         * \param samples_length Number of samples
         */
        virtual void demodulate_in_place(wsgc_complex *samples_in, unsigned int samples_length);

        /**
         * Demodulate samples out of place
         * \param samples_in Samples in input and output
         * \param samples_out Samples in input and output
         * \param samples_length Number of samples
         */
        virtual void demodulate_out_of_place(wsgc_complex *samples_in, wsgc_complex *samples_out, unsigned int samples_length);
        
        /**
         * Sets the origin samples at a value
         * \param samples_value Samples value
         */
        void set_value_at_origin(wsgc_complex value);

        /**
         * Copy given samples to the origin samples. There must be at least the _differential_length of samples in input
         * \param samples_value Samples value
         */
        void set_samples_at_origin(wsgc_complex *samples);

    protected:
        unsigned int _differential_length; //!< Number of samples in the differential delay
        wsgc_complex *_differential_samples_at_origin; //!< Keep in memory the differential samples not processed in case of successive calls
};

#endif // __DEMODULATOR_DIFFERENTIAL_H__
