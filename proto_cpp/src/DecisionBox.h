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
      
     Decision box

     Final correlation data analysis and message estimation. Derivated into two flavours, see:
        - DecisionBox_Piloted   : for operation with pilot(s) PRN(s) transmitted simultaneously with the message PRNs
        - DecisionBox_Unpiloted : for operation with only the message PRNs
*/
#ifndef __DECISION_BOX_H__
#define __DECISION_BOX_H__

#include "CorrelationRecord.h"
#include "DecisionRecord.h"
#include <vector>
#include <map>

/**
 * \brief Correlation data analysis and message estimation. Derivated into two flavours, see:
 *    - DecisionBox_Piloted   : for operation with pilot(s) PRN(s) transmitted simultaneously with the message PRNs
 *    - DecisionBox_Unpiloted : for operation with only the message PRNs
 */
class DecisionBox
{
    public:

		typedef enum
		{
			decision_status_true_accept,
			decision_status_true_reject,
			decision_status_false_accept,
			decision_status_false_reject
		} decision_status_t;

		/**
		 * Constructor
		 * \param prn_per_symbol Number of PRNs per symbol
         * \param fft_N Size of the FFT, this is also the number of samples in one PRN
		 */
        DecisionBox(unsigned int prn_per_symbol, unsigned int _fft_size);
        
        virtual ~DecisionBox();

		/**
		 * Returns a reference to the vector of decoded symbols
		 * \return Reference to the vector of decoded symbols
		 */
        const std::vector<int>& get_decoded_symbols() const
        {
            return _decoded_symbols;
        }
        
		/**
		 * Does a preparatory analysis phase - implementation specific
		 */
        virtual void analyze_records() = 0;
        
		/**
		 * Does a message symbols estimation - implementation specific
		 */
        virtual void estimate_symbols() = 0;
        
		/**
		 * Adjust magnitude value displayed
		 * \param mag_display_adj_factor Adjustement factor
		 */
        void set_mag_display_adj_factor(wsgc_float mag_display_adj_factor)
        {
        	_mag_display_adj_factor = mag_display_adj_factor;
        }

		/**
		 * Get magnitude value displayed adjustment factor
		 * \return Adjustement factor
		 */
        wsgc_float get_mag_display_adj_factor() const
        {
        	return _mag_display_adj_factor;
        }

        /**
         * Set CUDA indicator
         * \param use_cuda True if using CUDA
         */
        void set_use_cuda(bool use_cuda)
        {
        	_use_cuda = use_cuda;
        }

        bool is_prni_at_max_invalid() const
        {
        	return _prni_at_max_invalid;
        }

        /**
         * Append a new decision record
         * \return Reference to the decision record
         */
        DecisionRecord& new_decision_record();

        /**
         * Get a reference to the decision records
         * \return Reference to the vector of decision records
         */
        const std::vector<DecisionRecord>& get_decision_records() const
        {
            return _decision_records;
        }
        
        /**
         * Print decision records to output stream
         * \param os The output stream
         */
        void dump_decision_records(std::ostringstream& os) const;
        
        /**
         * Print decision errors
         * \param os The output stream
         * \param original_symbols The symbols that were sent
         * \param no_trivial Do not print trivial results (true accept status and no_record decision)
         */
        void dump_decision_status(std::ostringstream& os, std::vector<unsigned int>& original_symbols, bool no_trivial=true) const;

        
    protected:
        /**
        * \brief Specialized class to sort PRN index in symbol histogram
        */
        class GreaterPrnIndex
        {
            public:
                GreaterPrnIndex(unsigned int prn_per_symbol) : _prn_per_symbol(prn_per_symbol) {}
                
                bool operator()(const std::pair<unsigned int, unsigned int>& i, const std::pair<unsigned int, unsigned int>& j)
                {
                    if (j.second == i.second)
                    {
                        return ((i.first-j.first) % _prn_per_symbol) < ((j.first-i.first) % _prn_per_symbol);
                    }
                    else
                    {
                        return j.second < i.second;
                    }                    
                }
            
            private:
                unsigned int _prn_per_symbol;
        };
        
        unsigned int _prn_per_symbol; //!< Number of PRNs per symbol
        unsigned int _fft_size; //!< Size of the FFT, this is also the number of samples in one PRN
        unsigned int _preferred_symbol_prn_i; //!< Preferred PRN index in symbol for symbol estimation
        unsigned int _mag_display_adj_factor; //!< Adjustement factor for magnitude displays
        bool _prni_at_max_invalid; //!< Cannot estimate preferred PRN index in symbol with confidence
        bool _use_cuda; //!< Set if using CUDA. Calculations may be different so different criteria could be used.
        std::map<unsigned int, unsigned int> _symbol_prn_i_at_max; //!< Maxumum module PRN index in symbol occurences dictionnary
        std::vector<std::pair<unsigned int, unsigned int> > _histo_symbol_prn_i_at_max; //!< Maxumum module PRN index in symbol histogram
        std::vector<int> _decoded_symbols; //!< Symbols decoded, -1 for erasure
        std::vector<DecisionRecord> _decision_records; //!< Decision records giving details on decision process

        static const unsigned int preferred_symbol_prn_i_margin_threshold; //!< minimum margin between preferred PRN index in symbol and next preferred to have confidence
        static const unsigned int single_preferred_symbol_prn_i_threshold; //!< minimum number of occurences of preferred PRN index in symbol when single to have confidence
        static const wsgc_float peak_margin_threshold;  //!< selected correlation peak max difference with next / correlations ratio threshold

        /**
         * Do an estimation of the preferred PRN index in the symbol PRN repetition for estimation (estimates end of symbol)
         * \param correlation_records Reference to the correlation records vectors
         */
        void estimate_preferred_symbol_prn_i(const std::vector<CorrelationRecord>& correlation_records); //!< Estimate the preferred PRN index in symbol for symbol estimation

        /**
         * Common histogram data ordering method for sorting
         */
        static bool histo_order(const std::pair<unsigned int, unsigned int>& i, const std::pair<unsigned int, unsigned int>& j); //!< Method to sort simple histograms

        /**
         * Print decoding status
         */
        void dump_decoding_status(std::ostringstream& os, decision_status_t decision_status) const;
};

#endif // __DECISION_BOX_H__
