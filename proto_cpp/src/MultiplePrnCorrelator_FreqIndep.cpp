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

#include "MultiplePrnCorrelator_FreqIndep.h"
#include "GoldCodeGenerator.h"
#include "WsgcUtils.h"
#include <string.h>

MultiplePrnCorrelator_FreqIndep::MultiplePrnCorrelator_FreqIndep(GoldCodeGenerator& gc_generator, LocalCodesFFT_Host& local_codes, wsgc_float f_sampling, wsgc_float f_chip, unsigned int prn_per_symbol,
                                                            unsigned int prn_per_symbol_i_init, bool file_debugging) :
    _code_modulator(),
    MultiplePrnCorrelator::MultiplePrnCorrelator(local_codes, gc_generator, f_sampling, f_chip, prn_per_symbol, prn_per_symbol_i_init, file_debugging)
{
}


MultiplePrnCorrelator_FreqIndep::~MultiplePrnCorrelator_FreqIndep()
{
}
