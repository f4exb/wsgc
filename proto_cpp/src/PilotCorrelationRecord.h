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
      
     Pilot correlation record

     Stores one pilot correlation round data 
*/
#ifndef __PILOT_CORRELATION_RECORD_H__
#define __PILOT_CORRELATION_RECORD_H__

#include "WsgcTypes.h"
#include <sstream>

/**
 * \brief Holds pilot PRN correlation data.
 *
 * It is basically just a structure with a display method
 */
class PilotCorrelationRecord
{
    public:
        PilotCorrelationRecord();
        ~PilotCorrelationRecord();
        
        /**
         * Dumps the data to a string stream in frame form
         * \param os The output string stream
         * \param mag_display_factor Magnitudes are divided by thos factor for display
         */
        void dump(std::ostringstream& os, wsgc_float mag_display_factor = 1.0) const;

        /**
         * Dumps the data to a string stream in one line form
         * \param os The output string stream
         */
        void dump_oneline(std::ostringstream& os, wsgc_float mag_display_factor = 1.0) const;

        /**
         * Dumps a banner used when the data is printed in one line form
         * \param os The output string stream
         */
        static void dump_oneline_banner(std::ostringstream& os);
        
        /**
         * Resets values to initial values
         */
        void reset();
        
        unsigned int prn_index; //!< Global index of correlated PRN
        unsigned int block_count; //!< Number of averaging blocks processed
        unsigned int pilot_index; //!< PRN number of the pilot
        wsgc_float magnitude_max; //!< Maximum magnitude of correlation
        wsgc_float phase_at_max; //!< Phase at the correlation maximum
        unsigned int t_index_max; //!< Time index (i.e. delay) of the correlation maximum
        unsigned int f_index_max; //!< Frequency index of the correlation maximum
        wsgc_float delta_f; //!< Corresponding frequency displacement of the correlation maximum relative to zero IF
        bool selected; //!< PRN sample has been selected as valid for message correlation by the pilot correlation analyser
};

#endif // __PILOT_CORRELATION_RECORD_H__
