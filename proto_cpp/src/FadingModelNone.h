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
      
     Applies a no fading model to the test signal. 
     Useful when only AWGN addition from parent class is necessary.
     An actual active Fading model is instantiated in derivated classes.
*/
#ifndef __FADING_MODEL_NONE_H__
#define __FADING_MODEL_NONE_H__

#include "WsgcTypes.h"
#include "FadingModel.h"
#include <sstream>


class FadingModelNone : public FadingModel
{
    public:
        FadingModelNone(wsgc_float f_sampling);
        virtual ~FadingModelNone();
        virtual void apply_fading(const wsgc_complex *samples_in, wsgc_complex *samples_out, unsigned int nb_samples);
        virtual void print_fading_data(std::ostringstream& os) const;
        
        virtual unsigned int get_output_size(unsigned int input_size) const
        {
            return input_size;
        }
};

#endif // __FADING_MODEL_NONE_H__
