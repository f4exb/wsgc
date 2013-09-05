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
      
     Gold Code generator
      
     Creates a Gold Code sequences used in the symbol alphabet of the system. 
     It can be in binary (0,1) or (1,-1) format and optionnally sampled at given frequency
     
*/
#ifndef __GOLD_CODE_GENERATOR__
#define __GOLD_CODE_GENERATOR__

#include "WsgcTypes.h"
#include <vector>
#include <algorithm>
#include <iostream>

/**
 * \brief Gold Code generator
 * A Gold Code generator given the two generator polynomials
 */
class GoldCodeGenerator
{
    public:
        GoldCodeGenerator(
        		unsigned int nb_stages,
        		unsigned int nb_message_symbols,
        		unsigned int nb_service_symbols,
        		unsigned int nb_training_symbols,
        		std::vector<unsigned int>& g1_poly,
        		std::vector<unsigned int>& g2_poly);

        ~GoldCodeGenerator();
        
        /**
         * \return Total number of symbols (message and service)
         */
        unsigned int get_nb_codes() const
        {
            return _nb_message_symbols + _nb_service_symbols + _nb_training_symbols;
        }
        
        /**
         * \return Number of message symbols
         */
        unsigned int get_nb_message_codes() const
        {
            return _nb_message_symbols;
        }
        
        /**
         * \return Number of service symbols
         */
        unsigned int get_nb_service_codes() const
        {
            return _nb_service_symbols;
        }
        
        /**
         * \return Number of training symbols
         */
        unsigned int get_nb_training_codes() const
        {
            return _nb_training_symbols;
        }

        /**
         * \return Number of stages for shift register
         */
        unsigned int get_nb_stages() const
        {
            return _nb_stages;
        }
        
        /**
         * \return Length of code in bits
         */
        unsigned int get_code_length() const
        {
            return (1 << _nb_stages)-1;
        }
        
        /**
         * \return Number of samples for the code length
         */
        unsigned int get_nb_code_samples(wsgc_float f_sampling, wsgc_float f_chip) const
        {
            return int((get_code_length()/f_chip)*f_sampling);
        }
        
        /**
         * \return Successive LFSR register values for G1
         */
        const std::vector<unsigned int>& get_g1_regs() const
        {
            return _g1_regs;
        }
        
        /**
         * \return Successive LFSR register values for G2
         */
        const std::vector<unsigned int>& get_g2_regs() const
        {
            return _g2_regs;
        }
        
        /**
         * \return Powers of G1 polynomial
         */
        const std::vector<unsigned int>& get_g1_powers() const
        {
            return _g1_poly;
        }
        
        /**
         * \return Powers of G2 polynomial
         */
        const std::vector<unsigned int>& get_g2_powers() const
        {
            return _g2_poly;
        }
        
        /**
         * Makes code bits {0,1} at given sampling frequency and chip rate (frequency)
         * \param code Code bits are constructed there. Vector is cleared on input
         * \param prn_number Given code number
         * \param f_sample Sampling frequency. If zero (default) return raw code bits
         * \param f_chip Chip rate (frequency)
         */
        void make_code(std::vector<char>& code, unsigned int prn_number, wsgc_float f_sample=0.0, wsgc_float f_chip=0.0) const;
        /**
         * Print successive LFSR register values for G1
         * \param os Printing string stream
         */
        void print_g1_regs(std::ostringstream& os);
        /**
         * Print successive LFSR register values for G2
         * \param os Printing string stream
         */
        void print_g2_regs(std::ostringstream& os);
        
        
    private:
        unsigned int _nb_stages; //!< Number of LFSR stages
        unsigned int _nb_message_symbols; //!< Symbols used in the messages (first part of possible codes)
        unsigned int _nb_service_symbols; //!< Symbols used for service such as noise or pilot (middle part of possible codes)
        unsigned int _nb_training_symbols; //!< Training symbols used for unpiloted correlation (last part of possible codes)
        std::vector<unsigned int>& _g1_poly; // Generator polynomial powers for G1 (excluding max and 0)
        std::vector<unsigned int>& _g2_poly; // Generator polynomial powers for G2 (excluding max and 0)
        std::vector<unsigned int> _g1_regs;  // successive LFSR register values for G1
        std::vector<unsigned int> _g2_regs;  // successive LFSR register values for G2
        
        void lfsr_generator(std::vector<unsigned int>& g, std::vector<unsigned int>& polynomial_powers);
};

#endif // __GOLD_CODE_GENERATOR__
