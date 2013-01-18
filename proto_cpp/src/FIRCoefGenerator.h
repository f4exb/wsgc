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

#ifndef __FIR_COEF_GENERATOR__
#define __FIR_COEF_GENERATOR__

#include "WsgcTypes.h"
#include <sstream>
#include <vector>

/**
 * \brief FIR (real) coefficients generator
 */
class FIRCoefGenerator
{
public:
	/**
	 * Make a new set of tap coefficients
	 * \param f_sampling Sampling frequency
	 * \param f_cutoff Cutoff frequency. Must be twice the chip frequency for complex signals
	 * \param nb_taps Number of taps. Must be an odd number > 1. Defaults to next odd number if even.
	 */
	FIRCoefGenerator(
			wsgc_float f_sampling,
			wsgc_float f_cutoff,
			unsigned int nb_taps);

	virtual ~FIRCoefGenerator();

	/**
	 * Get the filter coefficients
	 * \return reference to the coefficients vector
	 */
	const std::vector<wsgc_float>& get_coefs() const
	{
		return _tap_coefs;
	}

	/**
	 * Dump filter characteristics to output stream
	 * \param os Output stream
	 */
	virtual void dump(std::ostringstream& os) const = 0;


protected:
	wsgc_float _f_sampling;  //!< Sampling frequency
	wsgc_float _f_cutoff;    //!< Cutoff frequency. Twice the chip frequency for complex signals.
	unsigned int _nb_taps;   //!< Number of taps. Must be an odd number > 1.
	unsigned int _half_taps; //!< Integer half of the number of taps ex: 11 -> 5
	std::vector<wsgc_float> _tap_coefs; //!< Tap coefficients in filter order

	/**
	 * Dump filter common characteristics to output stream
	 * \param os Output stream
	 */
	void dump_common(std::ostringstream& os) const;

	/**
	 * Dump filter coefficients to output stream
	 * \param os Output stream
	 */
	void dump_coefs(std::ostringstream& os) const;
};


#endif // __FIR_COEF_GENERATOR_RCOS__
