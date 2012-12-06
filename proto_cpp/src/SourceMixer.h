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

     Source mixer

     Mixes two simulated source outputs
*/

#ifndef __SOURCE_MIXER_H__
#define __SOURCE_MIXER_H__

#include "WsgcTypes.h"

class SimulatedSource;

class SourceMixer
{
	public:
		SourceMixer(SimulatedSource *source_A, SimulatedSource *source_B, wsgc_float b_gain);
		~SourceMixer();

		wsgc_complex *get_samples()
		{
			return _samples;
		}

		unsigned int get_nb_samples()
		{
			return _nb_samples;
		}


	protected:
		SimulatedSource *_source_A;
		SimulatedSource *_source_B;
		wsgc_float _B_gain;
		wsgc_complex *_samples;
		unsigned int _nb_samples;
};

#endif /* __SOURCE_MIXER_H__ */
