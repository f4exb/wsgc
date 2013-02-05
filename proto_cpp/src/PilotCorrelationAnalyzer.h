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

     Pilot correlation analyzer

     Stores pilot correlation records and corresponding source samples for further reference. 
     Analyzes pilot correlation records to drive message correlation. Used with piloted mode.
*/

#ifndef __PILOT_CORRELATION_ANALYZER_H__
#define __PILOT_CORRELATION_ANALYZER_H__

#include "WsgcTypes.h"
#include "PilotCorrelationRecord.h"
#include "CorrelationRecord.h"
#include "TrainingCorrelationRecord.h"
#include "TimeCorrelationAnalyzer.h"
#include <time.h>
#include <vector>
#include <map>
#include <iostream>

/**
 * \brief SampleAnalyzer Class to handle PRN samples sequencing
 *
 * Stores pilot correlation records and corresponding source samples for further reference. 
 * Analyzes pilot correlation records to drive message correlation. Used with piloted mode.
 * It is important that samples are stored in PRN sequence
 *
*/
class PilotCorrelationAnalyzer
{
	public:

		typedef struct timings_s
		{
			timespec _beg;
			timespec _end;
		} timings_t;

		/**
		 * Constructor
		 * \param analysis_window_size Size of analysis window in number symbols
		 * \param prn_per_symbol Number of PRNs per symbol
         * \param fft_N Size of the FFT, this is also the number of samples in one PRN
         * \param two_pilots True if the system supports message synchronization with two pilots
         * \param time_index_tolerance Tolerance +/- to validate the best time index
		 */
		PilotCorrelationAnalyzer(
            unsigned int analysis_window_size, 
            unsigned int prn_per_symbol, 
            unsigned int fft_N, 
            bool two_pilots = false, 
            unsigned int time_index_tolerance = 0);
        virtual ~PilotCorrelationAnalyzer();
        
		/**
		 * Stores the samples for one PRN. This increments the PRN global index.
		 * \param samples Pointer to the PRN samples
		 */
        void store_prn_samples(wsgc_complex *samples);
        
        /**
         * Get the samples for one PRN given its global index if available else null
         * \param prn_index PRN global index
         * \return Pointer to the first sample of PRN if available (still in buffer) else null
         */
        wsgc_complex *get_samples(unsigned int prn_index);
        
        /**
         * Get the samples for one PRN given its global index and sample shift index if available else null.
         * Used for training sequence
         * \param prn_index PRN global index
         * \param shift_index start sample shift index within PRN
         * \return Pointer to the first sample of PRN at start point of PRN if available (in buffer) else null
         */
        wsgc_complex *get_synchronized_samples(unsigned int prn_index, unsigned int shift_index);

        /**
         * Get a pointer to the last entered source samples
         * \return Pointer to the first sample of PRN of the last PRN signal source entered in the analyzer
         */
        const wsgc_complex *get_last_samples() const;

        /**
         * Append a new pilot correlation record
         * \return Reference to the pilot correlation record
         */
        PilotCorrelationRecord& new_pilot_correlation_record();

        /**
         * Append a new message correlation record
         * \param global_prn_index Global PRN index corresponding to this message PRN
         * \return Reference to the message correlation record
         */
        CorrelationRecord& new_message_correlation_record(unsigned int global_prn_index);

        /**
         * Append a new training correlation record
         * \param global_prn_index Global PRN index corresponding to this message PRN
         * \param sequence_length Length of training sequence
         * \param analysis_window_prn Length of analysis window in PRN numbers
         * \return Reference to the training correlation record
         */
        TrainingCorrelationRecord& new_training_correlation_record(
        		unsigned int global_prn_index,
        		unsigned int sequence_length,
        		unsigned int analysis_window_prn);
                
        /**
         * Get a reference to the pilot correlation records vector
         * \return Reference to the pilot correlation records vector
         */
        const std::vector<PilotCorrelationRecord>& get_pilot_correlation_records() const
		{
       		return _pilot_correlation_records;
		}
        
        /**
         * Get a reference to the message correlation records vector
         * \return Reference to the message correlation records vector
         */
        const std::vector<CorrelationRecord>& get_message_correlation_records() const
		{
        	return _message_correlation_records;
		}

        /**
         * Get a reference to the training correlation records vector
         * \return Reference to the training correlation records vector
         */
        std::vector<TrainingCorrelationRecord>& get_training_correlation_records() 
		{
        	return _training_correlation_records;
		}

        /**
         * Validate pilot correlation records for pilot 1 and pilot 2 that many places from the end. 
         * This can trigger the analysis process with these records when the analysis window is complete.
         * \param reverse_index Reverse index in the pilot correlation records vector
         * \return True if a cycle is completed and therefore the analysis and message processing using this pilot can take place
         */
        bool validate_pilot_correlation_records_back(unsigned int reverse_index);

        /**
         * Analyze the results of one analysis window length
         * \param best_time_shift_start Reference to an integer set to the best matching time shift peak start discovered
         * \param best_time_shift_length Reference to an integer set to the best matching time shift peak length discovered
         * \return True if the best matching time shift was discovered else false
         */
        bool analyze(unsigned int& best_time_shift_start, unsigned int& best_time_shift_length);
        
        /**
         * Post process noise data in message correlation records of the current analysis window
         * by averaging all valid samples and put back this value in all message correlation records of the current analysis window
         */
        void post_process_noise();

        /** 
         * Reset analysis objects and update counters and indexes to start a new analysis window
         */
        void reset_analysis();

        /**
         * Get a reference to the pilot correlation record that many places before the last element (0 for last element)
         * \param reverse_index Reverse index in the pilot correlation records vector
         * \return Reference to the pilot correlation record
         */
        PilotCorrelationRecord& get_pilot_correlation_record_back(unsigned int reverse_index);

        /**
         * Get index of pilot correlation records at start of analysis
         * \return Index of pilot correlation records at start of analysis
         */
        unsigned int get_start_analysis_pilot_correlation_records_index() const
        {
        	return _start_pilot_correlation_records_index;
        }

        /**
         * Get index of message correlation records at start of analysis
         * \return Index of message correlation records at start of analysis
         */
        unsigned int get_start_analysis_message_correlation_records_index() const
        {
        	return _start_message_correlation_records_index;
        }

        /**
         * Get analysis window size
         * \return Analysis window size
         */
        unsigned int get_analysis_window_size() const
        {
        	return _analysis_window_size;
        }

        /**
         * Get analysis window size in number of PRNs
         * \return Analysis window size
         */
        unsigned int get_analysis_window_size_in_prns() const
        {
        	return _analysis_window_size*_prn_per_symbol;
        }

        /**
         * Get PRN time shift analysis size
         * \return PRN time shift analysis size
         */
        unsigned int get_prn_time_shift_analysis_size() const
        {
        	return _fft_N;
        }

        /**
         * Get the current PRN global index that the one corresponding to the last entered sample
         * \return The PRN global index
         */
        unsigned int get_prn_index() const
        {
            return _global_prn_index_bot + _buffer_prn_count - 1;
        }

        /**
         * Get a reference to the time shifts dictionnary
         * \return Reference to the time shifts dictionnary
         */
        const std::map<unsigned int, unsigned int>& get_pilot_time_shift_occurences() const
        {
        	return _time_analyzer.get_pilot_time_shift_occurences();
            //return _pilot1_time_shift_occurences;
        }
        
        /**
         * Get a reference to the time shifts histogram vector 
         * \return Reference to the time shifts histogram vector
         */
        const std::vector<std::pair<unsigned int, unsigned int> >& get_histo_basic_time_shift_occurences() const
        {
        	return _time_analyzer.get_histo_basic_time_shift_occurences();
            //return _histo_basic_time_shift_occurences;
        }

        /**
         * Dumps the pilot correlation record vector to a string stream
         * \param os The output string stream
         */
        void dump_pilot_correlation_records(std::ostringstream& os) const;
        
        /**
         * Dumps the message correlation record vector to a string stream
         * \param os The output string stream
         */
        void dump_message_correlation_records(std::ostringstream& os) const;

        /**
         * Dumps the training correlation record vector to a string stream
         * \param os The output string stream
         */
        void dump_training_correlation_records(std::ostringstream& os) const;

        /**
         * Dumps the time shifts histogram vector 
         * \param os The output string stream
         */
        void dump_histo_time_shift_occurences(std::ostringstream& os) const;

        /**
         * Dump the timings for performance inspection
         * \param os Output string stream
         */
        void dump_timings(std::ostringstream& os);
        
        /**
         * Sets the pilot magnitude display factor. Magnitudes are divided by this factor for display.
         * \param mag_display_factor The magnitude display factor.
         */
        void set_pilot_mag_display_factor(wsgc_float mag_display_factor)
        {
            _pilot_mag_display_factor = mag_display_factor;
        }

        /**
         * Sets the message PRNs magnitude display factor. Magnitudes are divided by this factor for display.
         * \param mag_display_factor The magnitude display factor.
         */
        void set_message_mag_display_factor(wsgc_float mag_display_factor)
        {
            _message_mag_display_factor = mag_display_factor;
        }

        /**
         * Sets the training PRNs magnitude display factor. Magnitudes are divided by this factor for display.
         * \param mag_display_factor The magnitude display factor.
         */
        void set_training_mag_display_factor(wsgc_float mag_display_factor)
        {
            _training_mag_display_factor = mag_display_factor;
        }

        /**
         * Get the number of PRNs per symbol
         * \return The number of PRNs per symbol
         */
        unsigned int get_prn_per_symbol() const
        {
        	return _prn_per_symbol;
        }

        std::vector<timings_t> _pilot_mul_ifft_times;
        std::vector<timings_t> _pilot_avg_times;
        std::vector<timings_t> _pilot_cuda_avg_times;
        std::vector<timings_t> _pilot_cuda_reduc_times;
        std::vector<timings_t> _message_times;
        static const int _time_option;
        static const timings_t tmp_time;

	protected:

        // Time shift analyzer with pilot correlation records data
        TimeCorrelationAnalyzer<PilotCorrelationRecord> _time_analyzer;

        // analysis window
        unsigned int _analysis_window_size; //!< Size of analysis window in number symbols
        unsigned int _analysis_index; //!< Current analysis window index

        // PRN stuff
        unsigned int _prn_per_symbol; //!< Number of PRNs per symbol
        unsigned int _fft_N; //!< Size of the FFT, this is also the number of samples in one PRN
        
        // correlation records
        std::vector<PilotCorrelationRecord> _pilot_correlation_records; //!< Pilot correlation records corresponding to the stored samples for pilot 1
        std::vector<CorrelationRecord> _message_correlation_records; //!< Message correlation records filled by the message correlator
        std::vector<TrainingCorrelationRecord> _training_correlation_records; //!< Synchronization training correlation records filled by the message correlator
        unsigned int _start_pilot_correlation_records_index; //!< Index at start of pilot(s) correlation records analysis window
        unsigned int _start_message_correlation_records_index; //!< Index at start of message correlation records analysis window

        // process tuning
        unsigned int _time_index_tolerance; //!< Tolerance +/- to validate the best time index
        
        // display factors
        wsgc_float _pilot_mag_display_factor; //!< Pilot magnitudes are divided by this factor for display
        wsgc_float _message_mag_display_factor; //!< Message magnitudes are divided by this factor for display
        wsgc_float _training_mag_display_factor; //!< Message magnitudes are divided by this factor for display
        
        // memory buffer of samples management
		wsgc_complex *_samples_ext; //!< Copy of source samples. This is a double buffer twice the size of the analysis window plus one PRN extension at each end for synchronized access
		wsgc_complex *_samples; //!< Pointer to the core copy of source samples (without extension). It is the base for access
		unsigned int _sample_buffer_len; //!< Length of complete sample buffer in number of PRNs
        unsigned int _buffer_prn_index_bot; //!< Current bottom PRN index local to buffer (inclusive)
        unsigned int _buffer_prn_index_top; //!< Current top PRN index local to buffer (inclusive)
        unsigned int _global_prn_index_bot; //!< Current global bottom PRN index (absolute index value of _buffer_prn_index_bot PRN)
        unsigned int _buffer_prn_count; //!< Current count of PRNs in buffer

        /**
         * Append a new pilot correlation record
         * \param Reference to the pilot correlation records vector
         * \return Reference to the pilot correlation record
         */
        PilotCorrelationRecord& new_pilot_correlation_record(std::vector<PilotCorrelationRecord>& pilot_correlation_records);
        
        /**
         * Validate a pilot correlation record.
         * \param pilot_correlation_record Reference to the pilot correlation record
         */
        void validate_pilot_correlation_record(
        		PilotCorrelationRecord& pilot_correlation_record);
        
        /**
         * Get a reference to the pilot correlation record that many places before the last element (0 for last element)
         * \param pilot_correlation_records Reference to the pilot correlation records vector
         * \param reverse_index Reverse index in the pilot correlation records vector
         * \return Reference to the pilot correlation record
         */
        PilotCorrelationRecord& get_pilot_correlation_record_back(std::vector<PilotCorrelationRecord>& pilot_correlation_records, unsigned int reverse_index);

        
        /**
         * Return pointer to the start of synchronized samples given arbitrarily synchronized PRN start and
         * sample shift index of synchronization point
         * \param arbitrary_start_index Starting index as in get_samples method
         * \param shift_index Start sample shift index within PRN i.e. synchronization point
         * \return Pointer to the first sample
         */
        wsgc_complex *get_synchronized_samples_at_arbitrary_index(
        		unsigned int arbitrary_start_index,
        		unsigned int shift_index);




};

#endif /* __PILOT_CORRELATION_ANALYZER_H__ */
