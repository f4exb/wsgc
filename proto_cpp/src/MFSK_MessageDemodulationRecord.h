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

     MFSK_MessageDemodulationRecord

     Class that holds messsage demodulation data

*/
#ifndef __MFSK_MESSAGE_DEMODULATION_RECORD__
#define __MFSK_MESSAGE_DEMODULATION_RECORD__

#include "WsgcTypes.h"
#include <sstream>

/**
 * \brief Class that holds messsage demodulation data
 */
class MFSK_MessageDemodulationRecord
{
public:
	MFSK_MessageDemodulationRecord();
	~MFSK_MessageDemodulationRecord();

    void reset();
    void dump(std::ostringstream& os, wsgc_float magnitude_factor = 1.0) const;
    void dump_line(std::ostringstream& os, wsgc_float magnitude_factor = 1.0) const;
    static void dump_banner(std::ostringstream& os);

	unsigned int _symbol_index;
	unsigned int _symbol_ordinal;
	wsgc_float   _max_magnitude;
	wsgc_float   _avg_magnitude;
	wsgc_float   _noise_magnitude;
};

#endif // __MFSK_MESSAGE_DEMODULATION_RECORD__
