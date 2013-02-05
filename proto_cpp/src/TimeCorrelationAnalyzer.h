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

     Time correlation analyzer

     Analyzes time shift data of correlation records.
*/

#ifndef __TIME_CORRELATION_ANALYZER_H__
#define __TIME_CORRELATION_ANALYZER_H__

#include "WsgcTypes.h"
#include "TimeCorrelationRecord.h"
#include <vector>
#include <map>
#include <iostream>


/**
 * \brief Time shift data of correlation records analyzer common data and methods
 *
 */
class TimeCorrelationAnalyzer_Common
{
public:
	TimeCorrelationAnalyzer_Common(unsigned int fft_N, unsigned int time_index_tolerance = 0);
    virtual ~TimeCorrelationAnalyzer_Common();

    /**
     * Get a reference to the time shifts dictionnary
     * \return Reference to the time shifts dictionnary
     */
    const std::map<unsigned int, unsigned int>& get_pilot_time_shift_occurences() const
    {
        return _time_shift_occurences;
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
     * Get the number of validated records so far
     */
    unsigned int get_validate_count() const
    {
    	return _time_correlations_count;
    }

    /**
     * Dumps the time shifts histogram vector
     * \param os The output string stream
     */
    void dump_histo_time_shift_occurences(std::ostringstream& os) const;

    /**
     * Reset analysis data
     */
    void reset_analysis();

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
        unsigned int nb_corr_records; //!< Number of correlation records involved in the statistics
    };

    unsigned int _fft_N; //!< FFT size
    unsigned int _time_index_tolerance; //!< Tolerate that amount of jitter samples on the correlation peak position
    std::vector<std::vector<PrnTimeShiftRecord> > _histo_time_shift_occurences_vector; //!< Vector of histogram of time shift occurences for selected pilot (ordered by occurence number)
    std::vector<std::pair<unsigned int, unsigned int> > _histo_basic_time_shift_occurences; //!< Histogram of time shift occurences for selected pilot (ordered by occurence number)
    std::map<unsigned int, unsigned int> _time_shift_occurences; //!< Dictionnary (by time shift) of PRN time shifts occurences
    unsigned int _time_correlations_count; //!< Count of time correlation records processed
    wsgc_float _sum_corr_max; //!< Sum of correlation maxima for the time correlation records processed

    // static
    static const wsgc_float _best_time_margin_threshold;  //<! selected correlation time delta difference with next / number of correlations ratio threshold

    /**
    * Builds the time shift histogram by occurence number
    */
    void make_time_shift_histogram();
};


/**
 * \brief Analyzes time shift data of correlation records for correlation records implementing the TimeCorrelationRecord interface.
 *
 */
template <typename T_CorrelationRecord>
class TimeCorrelationAnalyzer : public TimeCorrelationAnalyzer_Common
{
	public:

        /**
		 * Constructor
         * \param fft_N FFT size or number of possible time shifts
         * \param time_index_tolerance Tolerance +/- to validate the best time index
		 */
		TimeCorrelationAnalyzer(unsigned int fft_N, unsigned int time_index_tolerance = 0);

        virtual ~TimeCorrelationAnalyzer();
        
        /**
         * Validates a time correlation record in the statistics
         * \param time_correlation_record Reference of the time correlation record
         */
        void validate_time_correlation_record(T_CorrelationRecord& time_correlation_record);
            
        
        /**
        * Analyzes time shift occurence to select best time shift if possible and flags the pilot correlation records accordingly
        * \param time_correlation_records Reference to the time correlation records vector
        * \param start_time_correlation_records_index Index of the first correlation record corresponding to the validated records
        * \param best_time_shift_start Reference to an integer set to the best matching time shift peak start discovered
        * \param best_time_shift_length Reference to an integer set to the best matching time shift peak length discovered
        * \return True if the best time shift could be determined else false        
        */
        bool analyze_time_shifts(
            std::vector<T_CorrelationRecord>& time_correlation_records,
            unsigned int start_time_correlation_records_index,
            unsigned int& best_time_shift_start,
            unsigned int& best_time_shift_length);        

};


//=================================================================================================
template <typename T_CorrelationRecord>
TimeCorrelationAnalyzer<T_CorrelationRecord>::TimeCorrelationAnalyzer(unsigned int fft_N, unsigned int time_index_tolerance) :
	TimeCorrelationAnalyzer_Common::TimeCorrelationAnalyzer_Common(fft_N, time_index_tolerance)
{
	std::cout << "TimeCorrelationAnalyzer::TimeCorrelationAnalyzer" << std::endl;
}


//=================================================================================================
template <typename T_CorrelationRecord>
TimeCorrelationAnalyzer<T_CorrelationRecord>::~TimeCorrelationAnalyzer()
{}


//=================================================================================================
template <typename T_CorrelationRecord>
void TimeCorrelationAnalyzer<T_CorrelationRecord>::validate_time_correlation_record(T_CorrelationRecord& time_correlation_record)
{
    // store current shift maximum in shifts dictionnary
    unsigned int shift_index_max = time_correlation_record.get_time_shift();
    std::map<unsigned int, unsigned int>::iterator shift_occurences_it;
    shift_occurences_it = _time_shift_occurences.find(shift_index_max);

    if (shift_occurences_it == _time_shift_occurences.end())
    {
    	_time_shift_occurences[shift_index_max] = 1;
    }
    else
    {
    	_time_shift_occurences[shift_index_max] += 1;
    }

    _time_correlations_count++;
    // accumulate the sum of maximum magnitudes
    _sum_corr_max += time_correlation_record.get_correlation_peak();
}


//=================================================================================================
template <typename T_CorrelationRecord>
bool TimeCorrelationAnalyzer<T_CorrelationRecord>::analyze_time_shifts(
    std::vector<T_CorrelationRecord>& time_correlation_records,
    unsigned int start_time_correlation_records_index,
    unsigned int& best_time_shift_start,
    unsigned int& best_time_shift_length)
{
    make_time_shift_histogram();

    PrnTimeShiftRecord& best_time_shift_peak = _histo_time_shift_occurences_vector.back().front();
    best_time_shift_start = best_time_shift_peak.peak_start;
    best_time_shift_length = best_time_shift_peak.peak_length;
    unsigned int best_time_shift_margin;

    if (_histo_time_shift_occurences_vector.back().size() > 1)
    {
    	best_time_shift_margin = best_time_shift_peak.peak_mag_sum - (_histo_time_shift_occurences_vector.back().begin()+1)->peak_mag_sum;
    }
    else
    {
    	best_time_shift_margin = best_time_shift_peak.peak_mag_sum;
    }

    bool time_shift_margin_ko = ((wsgc_float) best_time_shift_margin/_time_correlations_count) < _best_time_margin_threshold;

    typename std::vector<T_CorrelationRecord>::iterator correlation_record_it = time_correlation_records.begin() + start_time_correlation_records_index;
    const typename std::vector<T_CorrelationRecord>::iterator correlation_record_end = time_correlation_records.end();

    for (; correlation_record_it != correlation_record_end; ++correlation_record_it)
    {
    	TimeCorrelationRecord& correlation_record = *correlation_record_it;

    	if (time_shift_margin_ko)
    	{
    		correlation_record.set_selected(false);
    	}
    	else
    	{
            // checks whether the PRN time shift lies in the best PRN time shift peak interval
            if (best_time_shift_peak.peak_start + best_time_shift_peak.peak_length < _fft_N) // no time shift window boundary in the middle
            {
                if ((correlation_record.get_time_shift() < best_time_shift_peak.peak_start)
                 || (correlation_record.get_time_shift() >= best_time_shift_peak.peak_start + best_time_shift_peak.peak_length))
                {
                	correlation_record.set_selected(false);
                }
                else
                {
                	correlation_record.set_selected(true);
                }
            }
            else
            {
                if ((correlation_record.get_time_shift() < best_time_shift_peak.peak_start)
                 && (correlation_record.get_time_shift() >= (best_time_shift_peak.peak_start + best_time_shift_peak.peak_length) % _fft_N))
                {
                	correlation_record.set_selected(false);
                }
                else
                {
                	correlation_record.set_selected(true);
                }
            }
    	}
    }

    return !time_shift_margin_ko;
}




#endif /* __TIME_CORRELATION_ANALYZER_H__ */
