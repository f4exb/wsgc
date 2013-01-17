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

     FIR_RCoef

     Complex FIR filter with real coefficients

*/

#include "FIR_RCoef.h"
#include <assert.h>

const wsgc_complex FIR_RCoef::_c_zero = (0.0, 0.0);

FIR_RCoef::FIR_RCoef(std::vector<wsgc_float>& tap_coefs) :
	_tap_coefs(tap_coefs),
	_last_index(0)
{
	assert((tap_coefs.size() % 2 == 1) && (tap_coefs.size() > 1));
	_nb_taps = tap_coefs.size();
	_delay_line = new wsgc_complex[_nb_taps];

	for (unsigned int i=0; i < _nb_taps; i++)
	{
		_delay_line[i] = _c_zero;
	}
}


FIR_RCoef::~FIR_RCoef()
{
	delete[] _delay_line;
}


wsgc_complex FIR_RCoef::calc(wsgc_complex in)
{
	_delay_line[_last_index] = in;

	wsgc_complex acc = _c_zero;

	for (unsigned int i = 0; i < _tap_coefs.size(); ++i)
	{
		acc += _delay_line[(_last_index+i)%_tap_coefs.size()] * _tap_coefs[i];
	}

	_last_index = (_last_index + 1) % _tap_coefs.size();

	return acc;
}

