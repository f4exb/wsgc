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

     MFSK_MessageDemodulator

     Class to incoherently demodulate MFSK message. This is not using a correlation scheme
     and is there only for comparison to incoherent MFSK

     Abstract class to support Host or CUDA implementation

*/

#include "MFSK_MessageDemodulator.h"
#include "MFSK_MessageDemodulationRecord.h"


MFSK_MessageDemodulator::MFSK_MessageDemodulator(
		unsigned int fft_N,
		unsigned int nb_fft_per_symbol) :
	_fft_N(fft_N),
	_nb_fft_per_symbol(nb_fft_per_symbol)
{}


MFSK_MessageDemodulator::~MFSK_MessageDemodulator()
{}


void MFSK_MessageDemodulator::dump_demodulation_records(std::ostringstream& os, wsgc_float magnitude_factor) const
{
	std::vector<MFSK_MessageDemodulationRecord>::const_iterator it = _demodulation_records.begin();
	const std::vector<MFSK_MessageDemodulationRecord>::const_iterator it_end = _demodulation_records.end();

	MFSK_MessageDemodulationRecord::dump_banner(os);

	for (; it != it_end; ++it)
	{
		it->dump_line(os, magnitude_factor);
	}
}
