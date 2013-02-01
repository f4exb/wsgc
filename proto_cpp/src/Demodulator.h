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
      
     Demodulator
      
*/
#ifndef __DEMODULATOR_H__
#define __DEMODULATOR_H__

#include "WsgcTypes.h"

/**
 * \brief Demodulator abstract class
 */
class Demodulator
{
    public:
        Demodulator() {}
        virtual ~Demodulator() {}

        /**
         * Demodulate samples in place
         * \param samples_inout Samples in input and output
         * \param samples_length Number of samples
         */
        virtual void demodulate_in_place(wsgc_complex *samples_in, unsigned int samples_length) = 0;

        /**
         * Demodulate samples out of place
         * \param samples_in Samples in input and output
         * \param samples_out Samples in input and output
         * \param samples_length Number of samples
         */
        virtual void demodulate_out_of_place(wsgc_complex *samples_in, wsgc_complex *samples_out, unsigned int samples_length) = 0;
};

#endif // __DEMODULATOR_H__
