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

     FIRCoefGenerator

     Generic FIR (real) coefficients generator

*/

#include <FIRCoefGenerator.h>
#include <cmath>
#include <iostream>
#include <iomanip>


//=================================================================================================
FIRCoefGenerator::FIRCoefGenerator(
			wsgc_float f_sampling,
			wsgc_float f_cutoff,
			unsigned int nb_taps) :
	_f_sampling(f_sampling),
	_f_cutoff(f_cutoff),
	_nb_taps(nb_taps),
	_half_taps((nb_taps-1) / 2)
{
	// verify input

	if ((_f_cutoff <= 0.0) || (_f_cutoff > _f_sampling/2.0))
	{
		_f_cutoff = _f_sampling/2.0;
	}

	if (_nb_taps < 2)
	{
		_nb_taps = 2;
	}

	if (_nb_taps % 2 == 0)
	{
		_nb_taps += 1;
	}
}


//=================================================================================================
FIRCoefGenerator::~FIRCoefGenerator()
{
}


//=================================================================================================
void FIRCoefGenerator::dump_common(std::ostringstream& os) const
{
	os << "Sampling frequency ........: " << std::setw(8) << std::setprecision(1) << std::right << _f_sampling << std::endl;
    os << "Cutoff frequency ..........: " << std::setw(8) << std::setprecision(1) << std::right << _f_cutoff << std::endl;
    os << "Nb taps ...................: " << std::setw(6) << std::right << _nb_taps << std::endl;
}


//=================================================================================================
void FIRCoefGenerator::dump_coefs(std::ostringstream& os) const
{
	std::vector<wsgc_float>::const_iterator it = _tap_coefs.begin();
	const std::vector<wsgc_float>::const_iterator it_end = _tap_coefs.end();

    os << "Tap coefficients:" << std::endl;

	int i = 0;

	for (; it != it_end; ++it, i++)
	{
		os << std::setw(2) << i << " (" << std::setw(3) << i-(int)_half_taps << "): " << std::setw(8) << std::setprecision(6) << *it;
		os << std::endl;
	}
}
