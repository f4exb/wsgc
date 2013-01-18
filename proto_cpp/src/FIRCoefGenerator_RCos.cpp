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

     FIRCoefGenerator_RCos

     FIR (real) coefficients generator for raised cosine filters

*/

#include "FIRCoefGenerator_RCos.h"
#include <cmath>
#include <deque>
#include <iostream>
#include <iomanip>

//=================================================================================================
FIRCoefGenerator_RCos::FIRCoefGenerator_RCos(
			wsgc_float f_sampling,
			wsgc_float f_cutoff,
			wsgc_float alpha,
			unsigned int nb_taps) :
	FIRCoefGenerator(f_sampling, f_cutoff, nb_taps),
	_alpha(alpha)
{
	// verify input

	if ((_alpha < 0.0) || (_alpha > 1.0))
	{
		_alpha = 1.0;
	}

	// make filter coefficients
	std::deque<wsgc_float> tmp_coefs;

	for (unsigned int tap_i = 0; tap_i <= _half_taps; tap_i++)
	{
		wsgc_float coef = rcos((tap_i*_f_cutoff)/_f_sampling);

		if (tap_i == 0)
		{
			tmp_coefs.push_back(coef);
		}
		else
		{
			tmp_coefs.push_back(coef);
			tmp_coefs.push_front(coef);
		}
	}

	_tap_coefs.assign(tmp_coefs.begin(), tmp_coefs.end());
}


//=================================================================================================
FIRCoefGenerator_RCos::~FIRCoefGenerator_RCos()
{
}


//=================================================================================================
wsgc_float FIRCoefGenerator_RCos::sinc(wsgc_float x)
{
	if (x == 0.0)
	{
		return 1.0;
	}
	else
	{
		return sin(x)/x;
	}
}


//=================================================================================================
wsgc_float FIRCoefGenerator_RCos::rcos(wsgc_float t_tau)
{
	wsgc_float two_alpha_t_tau = 2.0*_alpha*t_tau;

	if (two_alpha_t_tau == 1.0)
	{
		return (M_PI_4)*sinc(t_tau);
	}
	else
	{
		return (sinc(t_tau)*cos(_alpha*M_PI*t_tau)) / (1.0-(two_alpha_t_tau*two_alpha_t_tau));
	}
}


//=================================================================================================
void FIRCoefGenerator_RCos::dump(std::ostringstream& os) const
{
	os << "Type ......................: Raised cosine" << std::endl;
	dump_common(os);
    os << "Rolloff factor ............: " << std::setw(9) << std::setprecision(2) << std::right << _alpha << std::endl;
    dump_coefs(os);
}
