/*
 * PrnAutocorrelator.h
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

     PrnAutocorrelator

     Do an auto-correlation of one PRN to the next. Host version.

 */

#ifndef __PRN_AUTOCORRELATOR_HOST_H__
#define __PRN_AUTOCORRELATOR_HOST_H__

#include "PrnAutocorrelator.h"

/**
 * \brief Autocorrelator engine that makes the correlation of one PRN to the next. Averages over one symbol length.
 */
class PrnAutocorrelator_Host : public PrnAutocorrelator
{
public:
    
    /**
    * Correlator engine constructor
    * \param fft_N FFT size (number of PRN samples)
    * \param prn_per_symbol Number of PRNs per symbols (= nb of averaging PRNs)
    */
	PrnAutocorrelator_Host(
			unsigned int fft_N,
			unsigned int prn_per_symbol);

    virtual ~PrnAutocorrelator_Host();
    /**
     * Set the PRN source block pointer. It is assumed the samples are one PRN length.
     * \param source_block source samples
     */
    virtual void set_source_block(const wsgc_complex *source_block, unsigned int global_prn_index);
    
    /**
     * Execute one autocorrelation operation with the current PRN samples.
     * \param autcorrelation_records Reference to the (message) autocorrelation records
     */
    virtual void make_correlation(std::vector<AutocorrelationRecord>& autocorrelation_records);
    

protected:
    wsgc_complex *_prn_samples; //!< PRN samples double buffer
    wsgc_complex *_acorrs;      //!< Autocorrelation values over one symbol length
};

#endif // __PRN_AUTOCORRELATOR_HOST_H__
