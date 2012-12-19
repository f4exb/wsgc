/*
     Copyright 2012 Edouard Griffiths <f4exb at free dot fr>
 
     This file is part of WSGC. A Weak Signal transmission mode using Gold Codes.
 
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
      
     SinglePrnCorrelator_FreqDep
      
     Does a IFFT batch correlation and averaging given one single PRN to look for over a range of frequencies.
     
*/
#ifndef __SINGLE_PRN_CORRELATOR_FREQ_DEP_CUDA_H__
#define __SINGLE_PRN_CORRELATOR_FREQ_DEP_CUDA_H__

#include "SinglePrnCorrelator_FreqDep.h"
#include "WsgcTypes.h"
#include "LocalCodesFFT_Cuda.h"

#include <cufft.h>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <cublas_v2.h>

#include <vector>

class GoldCodeGenerator;
class LocalCodesFFT_Cuda;
class CudaManager;
class PilotCorrelationAnalyzer;

/**
 * \brief Correlator engine to look for one PRN(s) in a range of shift frequencies
 *
 * Original one PRN length samples have been multiplied by a range of sub-frequencies in a previous step
 * Takes the succession of sub-frequency mixed PRN samples and processes it in order to:
 * - process correlation using multiplication in the frequency domain and IFFT method
 * - do the complex averaging 
 *
 * Does correlation in the frequency domain using multiplication by the local code (times the number of frequency shift steps) followed by IFFT.
 * Frequency shift is performed by rotating the local code. This gives discrete frequencies which step corresponds to the FFT bin size that is
 * also the recurrence frequency of the PRN code. This leaves gaps where correlation cannot be made as it has to be at most 10% of the FFT bin size
 * away. Therefore at each step the bunch of pre-mixed PRN copies is processed (and not only just one PRN at zero IF). The pre-mixing is done at
 * an earlier stage before this correlator is invoked.
 * In addition it takes a block of inputs to be processed in one batch (see prn_per_block parameter of constructor).
 *
 * This is a GPU implementation using CUDA.
 *
 */
class SinglePrnCorrelator_FreqDep_Cuda : public SinglePrnCorrelator_FreqDep
{
public:
	/**
	* Runs on a source samples batch (one PRN length)
	*
	* IFFT organization
	* For IFFT elements place elements of successive batches next to each other so as to optimize average (sum) calculation later => as much as possible in the cache
	* <-- fpi ---------------><- fpi+1 --- ...   Elements for fpi then elements for fpi+1 etc... (nb_pilot_f_bins frequencies)
	* [......|......| ...    ][......|...
	*    ||      /\            |\        IFFT element j for batch #0 then IFFT element j for batch #1 etc... (_nb_batch_prns=3 => _storage_depth = 6)
	*    ||     /  \          |ti+1,0,0|...
	*    /\    /    \
	*   /  \   |ti,1,0|ti,1,1|...
	*  /    \
	* [ti,0,0|ti,0,1|ti,0,2||ti,0,3|ti,0,4|ti,0,5|
	* <- even -------------><- odd --------------><- even ...
	*
	* ti,j,k with i = frequency index; j = FFT sample index; k = batch index
	*
	* In the next documentation concerning this class:
	* - B is the batch size
	* - T is the FFT/IFFT size
	* - Isf is the number of sub-frequency steps or sub-frequencies
	* - Ihf is the number of frequency steps or harmonic frequencies
	*
	* \param gc_generator Gold Code generator used to build the codes
	* \param code_modulator Modulator used to build the codes
	* \param f_sampling Sampling frequency
	* \param f_chip Chip rate (frequency)
    * \param _pilot_symbols Reference to the list of pilot symbol PRNs
	* \param nb_f_bins Number of frequency bins explored around IF=0
	* \param prn_per_block Number of PRNs per (symbol) block
	* \param nb_batch_prns Number of PRNs processed in one batch ("PRN batch factor")
	* \param frequency_step_division Frequency step division
	* \param cuda_device CUDA device number to use
	*/
	SinglePrnCorrelator_FreqDep_Cuda(
			GoldCodeGenerator& gc_generator,
			CodeModulator& code_modulator,
			wsgc_float f_sampling,
			wsgc_float f_chip,
			std::vector<unsigned int>& pilot_symbols,
			unsigned int nb_f_bins,
			unsigned int prn_per_block=4,
			unsigned int nb_batch_prns=3,
			unsigned int frequency_step_division=1,
			unsigned int cuda_device=0);

	virtual ~SinglePrnCorrelator_FreqDep_Cuda();

	/**
	 * Runs multiplication and IFFT on a source samples one PRN length.
	 * \param source_block Pointer to vector of FFT of source samples array
	 * \param pilot_prn_index Index of pilot PRN in the local codes matrix
	 * \param PRN position in a two batch (storage depth) cycle
	 */
	void multiply_and_ifft(const thrust::device_vector<cuComplex>& source_block, unsigned int pilot_prn_index, unsigned int prn_position);

	/**
	 * Do the averaging sum of one half of the IFFT frames
	 * \param first_half true if the first half is to be processed else false
	 */
	virtual void execute_averaging(bool first_half);

	void set_pilot_correlation_analyzer(PilotCorrelationAnalyzer *pilot_correlation_analyzer)
	{
		_pilot_correlation_analyzer = pilot_correlation_analyzer;
	}

protected:
	LocalCodesFFT_Cuda _local_codes;                 //!< Local copy of pilot PRNs conjugate FFT codes
	cufftHandle _ifft_plan;                          //!< CUFFT transform plan for IFFT
	int _n[1];                                       //!< CUFFT Plan FFT size parameter
	int _inembed[1];                                 //!< CUFFT Plan parameter
	int _onembed[1];                                 //!< CUFFT Plan parameter
	thrust::device_vector<cuComplex> _d_ifft_in;     //!< Input area for IFFT
	thrust::device_vector<cuComplex> _d_ifft_out;    //!< Output area for IFFT
	thrust::device_vector<int> _d_avg_keys;          //!< Result keys for averaging
	unsigned int _nb_pilot_prns;                     //!< Number of possible pilot PRNs
	unsigned int _cuda_device;                       //!< CUDA device number on which to run
	PilotCorrelationAnalyzer *_pilot_correlation_analyzer;
	cublasHandle_t _cublas_handle;
};

#endif // __SINGLE_PRN_CORRELATOR_FREQ_DEP_CUDA_H__
