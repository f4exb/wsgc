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
         * Get a pointer to the last entered source samples
         * \return Pointer to the first sample of PRN of the last PRN signal source entered in the analyzer
         */
        const wsgc_complex *get_last_samples() const;

        /*
         * Append a new pilot correlation record
         * \param alternate_pilot Reference to a boolean set to true if pilot 2 has been selected, else false
         * \return Reference to the pilot correlation record
         */
        PilotCorrelationRecord& new_pilot_correlation_record(bool alternate = false);

        /*
         * Append a new message correlation record
         * \param global_prn_index Global PRN index corresponding to this message PRN
         * \return Reference to the message correlation record
         */
        CorrelationRecord& new_message_correlation_record(unsigned int global_prn_index);

        /**
         * Validate pilot correlation records for pilot 1 and pilot 2 that many places from the end. 
         * This triggers the analysis process with these records.
         * \param reverse_index Reverse index in the pilot correlation records vector
         * \param alternate_pilot Reference to a boolean set to true if pilot 2 has been selected, else false
         * \return True if a cycle is completed and therefore the analysis using this pilot can take place
         */
        bool validate_pilot_correlation_records_back(unsigned int reverse_index);

        /**
         * Analyze the results of one analysis window length
         * \param alternate_pilot Reference to a boolean set to true if pilot 2 has been selected, else false
         * \param best_time_shift_start Reference to an integer set to the best matching time shift peak start discovered
         * \param best_time_shift_length Reference to an integer set to the best matching time shift peak length discovered
         * \return True if the best matching time shift was discovered else false
         */
        bool analyze(bool& alternate_pilot, unsigned int& best_time_shift_start, unsigned int& best_time_shift_length);
        
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
         * \param alternate_pilot Reference to a boolean set to true if pilot 2 has been selected, else false
         * \return Reference to the pilot correlation record
         */
        PilotCorrelationRecord& get_pilot_correlation_record_back(unsigned int reverse_index, bool alternate = false);

        /**
         * Get index of pilot correlation records at start of analysis
         * \return Index of pilot correlation records at start of analysis
         */
        unsigned int get_start_analysis_pilot_correlation_records_index()
        {
        	return _start_pilot_correlation_records_index;
        }

        /**
         * Get index of message correlation records at start of analysis
         * \return Index of message correlation records at start of analysis
         */
        unsigned int get_start_analysis_message_correlation_records_index()
        {
        	return _start_message_correlation_records_index;
        }

        /**
         * Get analysis window size
         * \return Analysis window size
         */
        unsigned int get_analysis_window_size()
        {
        	return _analysis_window_size;
        }

        /**
         * Get PRN time shift analysis size
         * \return PRN time shift analysis size
         */
        unsigned int get_prn_time_shift_analysis_size()
        {
        	return _fft_N;
        }

        /**
         * Get the current PRN global index that the one corresponding to the last entered sample
         * \return The PRN global index
         */
        unsigned int get_prn_index()
        {
            return _prn_index -1 + _samples_start_global_prn_index;
        }

        /**
         * Get a reference to the pilot correlation records vector
         * \param alternate True for pilot2, false for pilot1
         * \return Reference to the pilot correlation records vector
         */
        const std::vector<PilotCorrelationRecord>& get_pilot_correlation_records(bool alternate = false) const
		{
        	if (alternate)
        	{
        		return _pilot2_correlation_records;
        	}
        	else
        	{
        		return _pilot1_correlation_records;
        	}
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
         * Get a reference to the time shifts dictionnary
         * \param alternate True for pilot2, false for pilot1
         * \return Reference to the time shifts dictionnary
         */
        const std::map<unsigned int, unsigned int>& get_pilot_time_shift_occurences(bool alternate) const
        {
            if (alternate)
            {
                return _pilot2_time_shift_occurences;
            }
            else
            {
                return _pilot1_time_shift_occurences;
            }
        }
        
        /**
         * Get a reference to the time shifts histogram vector 
         * \return Reference to the time shifts histogram vector
         */
        const std::vector<std::pair<unsigned int, unsigned int> >& get_histo_basic_time_shift_occurences() const
        {
            return _histo_basic_time_shift_occurences;
        }

        /**
         * Dumps the pilot correlation record vector to a string stream
         * \param os The output string stream
         * \param alternate True for pilot2, false for pilot1
         */
        void dump_pilot_correlation_records(std::ostringstream& os, bool alternate = false) const;
        
        /**
         * Dumps the message correlation record vector to a string stream
         * \param os The output string stream
         */
        void dump_message_correlation_records(std::ostringstream& os) const;

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
         * Sets the magnitude display factor. Magnitudes are divided by this factor for display.
         * \param mag_display_factor The magnitude display factor.
         */
        void set_mag_display_factor(wsgc_float mag_display_factor)
        {
            _mag_display_factor = mag_display_factor;
        }

        std::vector<timings_t> _pilot_mul_ifft_times;
        std::vector<timings_t> _pilot_avg_times;
        std::vector<timings_t> _pilot_cuda_avg_times;
        std::vector<timings_t> _pilot_cuda_reduc_times;
        std::vector<timings_t> _message_times;
        static const int _time_option;
        static const timings_t tmp_time;

	protected:
    
        class PrnTimeShiftRecord
        {
        public:
            PrnTimeShiftRecord();
            void dump(std::ostringstream& os, unsigned int time_shifts_size = 0) const;
            static bool order_peak_mag_sum(const PrnTimeShiftRecord& left, const PrnTimeShiftRecord& right);

            unsigned int peak_start;
            unsigned int peak_length;
            unsigned int peak_mag_sum;
            unsigned int peak_shift_at_max;
            unsigned int peak_max;
        };

        unsigned int _analysis_window_size; //!< Size of analysis window in number symbols
        unsigned int _prn_per_symbol; //!< Number of PRNs per symbol
        unsigned int _fft_N; //!< Size of the FFT, this is also the number of samples in one PRN
        bool _two_pilots; //!< True if the system supports message synchronization with two pilots
        unsigned int _time_index_tolerance; //!< Tolerance +/- to validate the best time index
		wsgc_complex *_samples; //!< Copy of source samples. This is a double buffer twice the size of the analysis window to allow parrallel or pipelined processing
		unsigned int _sample_buffer_len; //!< Length of sample buffer in number of PRNs
        unsigned int _samples_start_global_prn_index; //!< Global PRN index at start of the copy of source samples
        unsigned int _prn_index; //!< Current PRN index
        unsigned int _start_pilot_correlation_records_index; //!< Index at start of pilot(s) correlation records analysis window
        wsgc_float _sum_max_pilot1; //!< Sum of the magnitude maxima of pilot1 over an analysis window
        wsgc_float _sum_max_pilot2; //!< Sum of the magnitude maxima of pilot2 over an analysis window
        std::vector<PilotCorrelationRecord> _pilot1_correlation_records; //!< Pilot correlation records corresponding to the stored samples for pilot 1
        std::vector<PilotCorrelationRecord> _pilot2_correlation_records; //!< Pilot correlation records corresponding to the stored samples for pilot 2
        std::map<unsigned int, unsigned int> _pilot1_time_shift_occurences; //!< Dictionnary (by time shift) of pilot1 PRN time shifts occurences
        std::map<unsigned int, unsigned int> _pilot2_time_shift_occurences; //!< Dictionnary (by time shift) of pilot2 PRN time shifts occurences
        std::vector<std::pair<unsigned int, unsigned int> > _histo_basic_time_shift_occurences; //!< Histogram of time shift occurences for selected pilot (ordered by occurence number)
        //std::vector<PrnTimeShiftRecord> _histo_time_shift_occurences; //!< Histogram of time shift occurences for selected pilot (ordered by occurence number)
        std::vector<std::vector<PrnTimeShiftRecord> > _histo_time_shift_occurences_vector; //!< Vector of histogram of time shift occurences for selected pilot (ordered by occurence number)
        std::vector<unsigned int> _pilot_numbers; //!< Vector of pilot numbers detected
        unsigned int _analysis_index; //!< Current analysis window index
        std::vector<CorrelationRecord> _message_correlation_records; //!< Message correlation records filled by the message correlator
        unsigned int _start_message_correlation_records_index; //!< Index at start of message correlation records analysis window
        unsigned int _validate_count; //!< Counter of number of pilot correlation records validated
        wsgc_float _mag_display_factor; //!< Magnitudes are divided by this factor for display

        static const wsgc_float _best_time_margin_threshold;  //<! selected correlation time delta difference with next / number of correlations ratio threshold

        
        /**
         * Append a new pilot correlation record
         * \param Reference to the pilot correlation records vector
         * \return Reference to the pilot correlation record
         */
        PilotCorrelationRecord& new_pilot_correlation_record(std::vector<PilotCorrelationRecord>& pilot_correlation_records);
        
        /**
         * Validate a pilot correlation record.
         * \param pilot_correlation_record Reference to the pilot correlation record
         * \param shift_occurences_dict Dictionnary of pilot PRN time shifts occurences
         * \param sum_max_pilot Sum of magnitude maxima of pilot PRN correlation
         */
        void validate_pilot_correlation_record(
        		PilotCorrelationRecord& pilot_correlation_record,
        	    std::map<unsigned int, unsigned int>& shift_occurences_dict,
        	    wsgc_float& sum_max_pilot);
        
        /**
         * Get a reference to the pilot correlation record that many places before the last element (0 for last element)
         * \param pilot_correlation_records Reference to the pilot correlation records vector
         * \param reverse_index Reverse index in the pilot correlation records vector
         * \return Reference to the pilot correlation record
         */
        PilotCorrelationRecord& get_pilot_correlation_record_back(std::vector<PilotCorrelationRecord>& pilot_correlation_records, unsigned int reverse_index);

        
        /**
        * Builds the time shift histogram by occurence number
        * \param shift_occurences_dict Reference to the time shifts occurences dictionnary
        */
        void make_time_shift_histogram(std::map<unsigned int, unsigned int>& shift_occurences_dict);

        /**
        * Analyzes time shift occurence to select best time shift if possible and flags the pilot correlation records accordingly
        * \param shift_occurences_dict Reference to the time shifts occurences dictionnary
        * \param pilot_correlation_records Reference to the pilot correlation records vector
        * \param best_time_shift_start Reference to an integer set to the best matching time shift peak start discovered
        * \param best_time_shift_length Reference to an integer set to the best matching time shift peak length discovered
        * \return True if the best time shift could be determined else false        
        */
        bool analyze_time_shifts(
            std::map<unsigned int, unsigned int>& shift_occurences_dict,
            std::vector<PilotCorrelationRecord>& pilot_correlation_records, 
            unsigned int& best_time_shift_start,
            unsigned int& best_time_shift_length);        

};

#endif /* __PILOT_CORRELATION_ANALYZER_H__ */
