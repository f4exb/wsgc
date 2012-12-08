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
      
     PilotCorrelator
      
     Takes samples for the length of a PRN and processes it in order to:
     - process correlation over all frequency bins for one or two pilot code sequences
     - do peak estimation therefore delay and frequency tracking
     
*/
#ifndef __PILOT_CORRELATOR_H__
#define __PILOT_CORRELATOR_H__

#include "WsgcTypes.h"
#include "PilotCorrelationRecord.h"
#include <vector>

class GoldCodeGenerator;
class CodeModulator;
class LocalCodesFFT;
class SinglePrnCorrelator_FreqDep;

/**
 * \brief Correlator engine to acquire and track one or two pilot PRN(s)
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
class PilotCorrelationAnalyzer;

class PilotCorrelator
{
public:
	/**
	* Runs on a source samples batch (one PRN length)
	* \param f_sampling Sampling frequency
	* \param f_chip Chip rate (frequency)
	* \param fft_N FFT/IFFT size or number of samples per PRN
	* \param pilot_symbols Reference to the list of pilot symbol PRNs
	* \param prn_per_symbol Number of PRNs per symbol thus per averaging block
	* \param nb_pilot_f_bins Number of frequency bins explored for pilot acquisition and tracking
	* \param nb_batch_prns Number of PRNs processed in one batch ("PRN batch factor")
	* \param freq_step_division Frequency step division
	*/
	PilotCorrelator(
			wsgc_float f_sampling,
			wsgc_float f_chip,
			unsigned int fft_N,
			std::vector<unsigned int>& pilot_symbols,
			unsigned int prn_per_symbol=4,
			unsigned int nb_pilot_f_bins=3,
			unsigned int nb_batch_prns=3,
			unsigned int freq_step_division=1);

	virtual ~PilotCorrelator();

	/**
	 * Runs on a source samples of one PRN length. Depending on the batch processing settings
	 * this may or may not trigger the actual calculation immediately.
	 * The pilot correlation analyzer manages the source samples and pilot correlation records
	 * \param pilot_correlation_analyzer Reference of the pilot correlation analyzer
	 * \param pilot_prn_index Index of pilot PRN code
	 */
	virtual void execute(PilotCorrelationAnalyzer& pilot_correlation_analyzer, unsigned int pilot_prn_index) = 0;
	/**
	 * Get the result of one correlation. Will be successful only after the last "execute" method invocation of the second
	 * correlation batch. Results are those of the last complete correlation batch.
	 * \param index in the correlation batch
	 * \param alternate true if alternate pilot (pilot #2) is requested
	 * \return reference to Correlation Record
	 */
	const PilotCorrelationRecord& get_correlation_record(unsigned int index=0, bool alternate=false) const;
	/**
	 * Get the result of a batch correlation. Will be successful only after the last "execute" method invocation of the second
	 * correlation batch. Results are those of the last complete correlation batch.
	 * \param alternate true if alternate pilot (pilot #2) is requested
	 * \return reference to Correlation Records vector
	 */
	const std::vector<PilotCorrelationRecord>& get_correlation_records(bool alternate=false) const;

	/**
	 * Get the pipeline length i.e the number of PRNs to process in input before one correlation record is available in output.
	 * In other terms the result obtained after execution on one PRN is the one corresponding to that many PRNs back from the PRN
	 * that has just been processed. This is normally the PRN batch factor multiplied by two.
	 * \return the pipeline length in number of PRNs
	 */
	unsigned int get_pipeline_length() const;

	/**
	 * Get the batch size in number of PRNs.
	 * \return the batch size in number of PRNs.
	 */
	unsigned int get_batch_size()
	{
		return _nb_batch_prns;
	}

	 /**
	 * Test if valid correlation records are available
	 * \return True if valid correlation records are available else false
	 */
	bool are_correlation_records_available() const
	{
		return _result_available;
	}

	/**
	 * Test if a new batch has been processed after last call to execute() method
	 * \return True if a new batch has been processed else false
	 */
	bool new_batch_processed()
	{
		return _new_batch_processed;
	}


protected:
	wsgc_float _f_sampling; //!< Sampling frequency
	wsgc_float _f_chip; //!< Chip rate
	std::vector<unsigned int>& _pilot_symbols; //!< Reference to the list of pilot symbol PRNs
	unsigned int _prn_per_symbol; //!< Number of PRNs per averaging block is also the number of PRNs per symbol
	unsigned int _nb_pilot_f_bins; //!< Number of FFT frequency bins explored. Odd number greater or equal to 3.
	unsigned int _prn_per_avg_i; //!< Index of the current PRN being processed in the averaging block
	unsigned int _freq_step_division; //!< Frequency step division factor
	unsigned int _fft_N; //!< Size of FFT
	unsigned int _nb_batch_prns; //!< Number of PRNs processed in one batch (pipelining) or "PRN batch factor"
	bool _result_available; //!< Is true as soon as valid correlation records are available
	bool _new_batch_processed; //!< A new batch has been processed following the last execute() command. Thus new results are available

	unsigned int _storage_depth; //!< Storage depth factor: 2*(nb_batch_prns)
	unsigned int _batch_index; //!< Batch index. Incremented after each batch processing. Ever increasing.
	std::vector<wsgc_float> _batch_max_magnitudes; //!< maximum magnitudes for the current batch
	std::vector<unsigned int> _batch_max_composite_indexes; //!< maximum composite indexes for the current batch

	/**
	 * Append correlation records with batch result
	 * \param pilot_correlation_analyzer Reference of the pilot correlation analyzer
	 * \param pilot_prn Pilot PRN index in GC codes
	 * \param even_batch True if the batch sequence is even else false
	 * \param pilot_prn_index Global index of pilot PRN at batch processing point
	 * \param ifft_correlator_pilot Pointer to the IFFT correlator object
	 */
	void update_pilot_correlation_records(
			PilotCorrelationAnalyzer& pilot_correlation_analyzer,
			unsigned int pilot_prn,
			bool even_batch,
			unsigned int pilot_prn_index,
			SinglePrnCorrelator_FreqDep *ifft_correlator_pilot);
};

#endif // __PILOT_CORRELATOR_H__
