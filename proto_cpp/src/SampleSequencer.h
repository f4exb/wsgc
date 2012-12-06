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

     Sample sequencer

     Manages the sequence of given samples so that samples corresponding to one PRN are given at a time
*/

#ifndef __SAMPLE_SEQUENCER_H__
#define __SAMPLE_SEQUENCER_H__

#include "WsgcTypes.h"

/**
 * \brief SampleSequencer Class to handle PRN samples sequencing
 *
 * This class manages the whole sequence of samples created by the source(s) simulator(s) so that samples for
 * one PRN are given at a time
 *
*/
class SampleSequencer
{
	public:
		/**
		 * Constructor
		 * \param samples Pointer to the array of samples
		 * \param nb_samples Number of samples
		 * \param nb_code_samples Number of samples in one PRN length
		 */
		SampleSequencer(wsgc_complex *samples, unsigned int nb_samples, unsigned int nb_code_samples);

		/**
		 * Gives a next bunch of PRN samples
		 * \param samples Pointer to the beginning of the PRN samples array
		 * \return True if the pointer to the PRN samples is valid else false
		 */
		bool get_next_code_samples(wsgc_complex **samples);


	protected:
		wsgc_complex *_samples; //!< Source samples
        unsigned int _nb_samples; //!< Number of source samples
        unsigned int _nb_code_samples; //!< Number of samples in one PRN
        unsigned int _serviced_samples_index; //!< Index of last sample serviced
};

#endif /* __SAMPLE_SEQUENCER_H__ */
