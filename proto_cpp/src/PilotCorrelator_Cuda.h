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
      
     PilotCorrelator - CUDA version
      
     Takes samples for the length of a PRN and processes it in order to:
     - process correlation over all frequency bins for one or two pilot code sequences
     - do peak estimation therefore delay and frequency tracking
     
*/
#ifndef __PILOT_CORRELATOR_CUDA_H__
#define __PILOT_CORRELATOR_CUDA_H__

#include "WsgcTypes.h"
#include "PilotCorrelationRecord.h"
#include "PilotCorrelator.h"
#include "SourceFFT_Cuda.h"
#include "LocalCodesFFT_Cuda.h"
#include "SinglePrnCorrelator_FreqDep_Cuda.h"

class GoldCodeGenerator;
class PilotCorrelationAnalyzer;

/**
 * \brief Correlator engine to acquire and track one or two pilot PRN(s) - Host version
 *
 * It assumes the search for correlation is done in two dimensions:
 * - One is the time shift of the origin of the PRN sequence in the time referential of the receiving system (of course)
 * - The other is the frequency shift from the zero IF on which the receiving system is set. This also assumes that the 
 *   modulation used has some frequency dependence and the exact frequency of carrier must be acquired and tracked.
 *
 * Takes samples for the length of a PRN and processes it in order to:
 * - process correlation over all frequency bins for one or two pilot code sequences
 * - do peak estimation therefore delay and frequency tracking
 *
 * Pilot PRN(s) is (are) taken among service PRNs normally the second and third service PRNs respectively.
 *
 * Does correlation in the frequency domain using FFT+IFFT on one or two pilot PRNs (two is for autonomous message synchronization).
 * It looks in a number of frequency intervals or "bins" around the center receiving frequency. The length of this frequency interval
 * is half the chip frequency divided by the code length. For example with a code length of 1023 bits and a chip frequency of 1023 Hz
 * this is 0.5 Hz. See SinglePrnCorrelator class for the IFFT and averaging operations using a pilot PRN.
 *
 * Batch processing allows for some parallelization (good for GPU implementation). The result at each PRN code length
 * is obtained by a moving average over the number of PRNs per batch (PRN batch factor) last results. We keep two times this number
 * of [fi,ti] frames in memory. This way the moving average can be always be calculated with a delay vs IFFT. Example with 4 PRNs
 * per block and 3 PRNs per batch. This yields 3 correlation cycles per batch. This is the minimum number of correlation cycles
 * as there should always be at least the number of PRNs in an averaging cycle minus one (see documentation for details):
 *
 * |  |  |  |  |  |  |
 *  --------              IFFT 0
 *           --------     IFFT 1
 *  -----------           AVG 0
 *  --------              IFFT 2
 *           --------     |AVG1
 *  --                    |
 *
 *  Averaging (sum) is calculated as follows:
 *
 *  |  |  |  |  |  |  |
 *   ^  :  :  :    Sum 4 next cells and put result in first cell
 *   +--+--+--+             |
 *      ^  :  :  :          |
 *      +--+--+--+          V
 *         ^  :  :  :       Repeat 3 times (Averaging length - 1)
 *         +--+--+--+
 *
 * Each time the number of "PRN per block" (or averaging length) minus one results are calculated
 *
 * This fits CUDA very well and is a bit CUDA oriented. Nevertheless it allows performance improvement for CPU based calculation also
 *
 */

class PilotCorrelator_Cuda : public PilotCorrelator
{
public:
	/**
	* Runs on a source samples batch (one PRN length)
	* \param gc_generator Gold Code generator used to build the codes
	* \param code_modulator Modulator used to build the codes
	* \param f_sampling Sampling frequency
	* \param f_chip Chip rate (frequency)
    * \param _pilot_symbols Reference to the list of pilot symbol PRNs
    * \param cuda_device CUDA GPU# on which to run
	* \param prn_per_symbol Number of PRNs per symbol thus per averaging block
	* \param nb_pilot_f_bins Number of frequency bins explored for pilot acquisition and tracking
	* \param nb_batch_prns Number of PRNs processed in one batch ("PRN batch factor")
	* \param freq_step_division Frequency step division
	*/
	PilotCorrelator_Cuda(
        		const GoldCodeGenerator& gc_generator,
        		CodeModulator& code_modulator,
        		wsgc_float f_sampling,
        		wsgc_float f_chip,
        		std::vector<unsigned int>& _pilot_symbols,
        		unsigned int cuda_device,
                unsigned int prn_per_symbol=4,
                unsigned int nb_pilot_f_bins=3,
                unsigned int nb_batch_prns=3,
                unsigned int freq_step_division=1);

    virtual ~PilotCorrelator_Cuda();
        
	/**
	 * Runs on a source samples of one PRN length. Depending on the batch processing settings
	 * this may or may not trigger the actual calculation immediately.
	 * The pilot correlation analyzer manages the source samples and pilot correlation records
	 * \param pilot_correlation_analyzer Reference of the pilot correlation analyzer
	 * \param pilot_prn_index Index of pilot PRN code
	 */
	virtual void execute(PilotCorrelationAnalyzer& pilot_correlation_analyzer, unsigned int pilot_prn_index);

    protected:
        SourceFFT_Cuda _source_fft; //!< Source samples FFT manager
        SinglePrnCorrelator_FreqDep_Cuda _ifft_correlator_pilot; //!< IFFT correlator for pilot
};

#endif // __PILOT_CORRELATOR_HOST_H__
