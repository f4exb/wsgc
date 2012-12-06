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

     Correlator

     Takes samples for the length of a PRN and processes it in order to:
     - process correlation over all PRNs symbols in the alphabet for one code sequence
     - do the frequency tracking
*/
#ifndef __CORRELATOR_H__
#define __CORRELATOR_H__

class PilotCorrelator;

/*
 * Takes samples for the length of a PRN and processes it in order to:
 * - process correlation over all PRNs symbols in the alphabet for one code sequence
 * - do the frequency tracking
 */
class Correlator
{
public:
    /**
    * Runs on a source samples batch (one PRN length)
    * \param gc_generator Gold Code generator
    * \param local_codes PRN signals local copy
    * \param f_sampling Sampling frequency
    * \param f_chip Chip rate (frequency)
    * \param pilot1 Pilot #1 PRN index
    * \param pilot2 Pilot #2 PRN index
    * \param prn_per_block Number of PRNs per (symbol) block
    * \param nb_pilot_f_bins Number of frequency bins explored for pilot acquisition and tracking
    * \param correlation_blocking_factor Number of PRN cycles processed by the correlator in one batch
    */
    Correlator(GoldCodeGenerator& gc_generator, LocalCodesFFT& local_codes, wsgc_float f_sampling, wsgc_float f_chip,
                    unsigned int pilot1, unsigned int pilot2, unsigned int prn_per_block=4, unsigned int nb_pilot_f_bins=3,
                    unsigned int correlation_blocking_factor=1);
    virtual ~Correlator();
    /**
    * Process correlation over samples for the length of a PRN. This is implementation (CPU/GPU) dependent
    * \param source_block Pointer to first element in source samples array
    */
    virtual void execute(wsgc_fftw_complex *source_block) = 0;

};

#endif // __CORRELATOR_H__

