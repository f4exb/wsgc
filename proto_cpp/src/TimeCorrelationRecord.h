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
      
     Time Correlation record

     Interface for time shift related methods on correlation records
*/
#ifndef __TIME_CORRELATION_RECORD_H__
#define __TIME_CORRELATION_RECORD_H__

#include "WsgcTypes.h"
#include <sstream>

class TimeCorrelationRecord
{
    public:
        TimeCorrelationRecord() {};
        virtual ~TimeCorrelationRecord() {};
        
        virtual unsigned int get_time_shift() const = 0;
        virtual wsgc_float get_correlation_peak() const = 0;
        virtual void set_selected(bool selected) = 0;
};

#endif // __TIME_CORRELATION_RECORD_H__
