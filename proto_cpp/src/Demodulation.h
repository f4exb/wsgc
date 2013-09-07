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
      
     Demodulation
      
     Common part of all demodulation classes.
     
*/
#ifndef __DEMODULATION_H__
#define __DEMODULATION_H__

#include "WsgcTypes.h"

class Options;
class GoldCodeGenerator;

class Demodulation
{
public:
    static void demodulate_generic(Options& options, const GoldCodeGenerator& gc_generator, wsgc_complex *faded_source_samples, unsigned int nb_faded_source_samples);    
};

#endif // __DEMODULATION_H__
