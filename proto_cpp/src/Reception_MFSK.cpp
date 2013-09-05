/*
     Copyright 2012-2013 Edouard Griffiths <f4exb at free dot fr>

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

     Reception MFSK
*/

#include "Reception_FEC.h"
#include "Reception_MFSK.h"
#include "MFSK_MessageDemodulator.h"
#include "MFSK_MessageDemodulator_Host.h"
#include "DecisionBox_MFSK.h"
#include "SampleSequencer.h"
#include "WsgcUtils.h"

#include <iostream>
#include <iomanip>

#ifdef _RSSOFT
#include "RS_ReliabilityMatrix.h"
#endif

#ifdef _CCSOFT
#include "CC_ReliabilityMatrix.h"
#endif


//=================================================================================================
Reception_MFSK::Reception_MFSK(Options& _options, const GoldCodeGenerator& _gc_generator) :
    Reception(_options, _gc_generator)
{}

//=================================================================================================
Reception_MFSK::~Reception_MFSK()
{}

//=================================================================================================
void Reception_MFSK::message_processing(wsgc_complex *faded_source_samples, unsigned int nb_faded_source_samples)
{
    timespec time1, time2;
    int time_option = CLOCK_REALTIME;
    MFSK_MessageDemodulator_Host mfsk_message_demodulator(
            options.mfsk_options._fft_N,
            options.mfsk_options._nb_fft_per_symbol,
            options.mfsk_options._zero_fft_slot,
            options.nb_message_symbols,
            options.nb_service_symbols);
    SampleSequencer sample_sequencer(faded_source_samples, nb_faded_source_samples, options.mfsk_options._fft_N*options.mfsk_options._nb_fft_per_symbol);

    wsgc_complex *signal_samples;

    if (options.fec_scheme == Options::OptionFEC_RSSoft)
    {
#ifdef _RSSOFT
        rssoft::RS_ReliabilityMatrix rs_relmat(options.rs_logq, (1<<options.rs_logq) - 1);

        clock_gettime(time_option, &time1);

        while (sample_sequencer.get_next_code_samples(&signal_samples)) // pseudo real time loop, one PRN length at a time
        {
            mfsk_message_demodulator.execute(signal_samples, rs_relmat);
        }

        clock_gettime(time_option, &time2);
        std::cout << "MFSK message sequence decoding time: " << std::setw(12) << std::setprecision(9) << WsgcUtils::get_time_difference(time2,time1) << " s" << std::endl << std::endl;

        Reception_FEC::run_rssoft_decoding(options, rs_relmat);
#else
        std::cout << "Program not linked with RSSoft, aborting decode" << std::endl;s
#endif
    }
    else if (options.fec_scheme == Options::OptionFEC_CCSoft)
    {
#ifdef _CCSOFT
        ccsoft::CC_ReliabilityMatrix cc_relmat(options.nb_message_symbols, options.prns.size()-1);

        clock_gettime(time_option, &time1);

        while (sample_sequencer.get_next_code_samples(&signal_samples)) // pseudo real time loop, one PRN length at a time
        {
            mfsk_message_demodulator.execute(signal_samples, cc_relmat);
        }

        clock_gettime(time_option, &time2);
        std::cout << "MFSK message sequence decoding time: " << std::setw(12) << std::setprecision(9) << WsgcUtils::get_time_difference(time2,time1) << " s" << std::endl << std::endl;

        Reception_FEC::run_ccsoft_decoding(options, cc_relmat);
#else
        std::cout << "Program not linked with CCSoft, aborting decode" << std::endl;s
#endif
    }
    else
    {
        clock_gettime(time_option, &time1);

        while (sample_sequencer.get_next_code_samples(&signal_samples)) // pseudo real time loop, one PRN length at a time
        {
            mfsk_message_demodulator.execute(signal_samples);
        }

        clock_gettime(time_option, &time2);
        std::cout << "MFSK message sequence decoding time: " << std::setw(12) << std::setprecision(9) << WsgcUtils::get_time_difference(time2,time1) << " s" << std::endl << std::endl;

        std::ostringstream demod_os;
        demod_os << "--- demodulation records:" << std::endl;
        mfsk_message_demodulator.dump_demodulation_records(demod_os);
        std::cout << demod_os.str() << std::endl;

        DecisionBox_MFSK decision_box(options.mfsk_options._fft_N, options.mfsk_options._nb_fft_per_symbol, options.decision_thresholds);
        decision_box.estimate_symbols(mfsk_message_demodulator.get_demodulation_records());

        options.prns.pop_back(); // pop padding symbol

        std::vector<unsigned int> symbol_indices;
        for (unsigned int i=0; i<options.prns.size(); i++)
        {
            symbol_indices.push_back(i);
        }

        std::ostringstream os_result;
        options.decision_thresholds.print_options(os_result);
        os_result << std::endl << "Decisions status:" << std::endl;
        decision_box.dump_decision_status(os_result, options.prns, mfsk_message_demodulator.get_demodulation_records());
        os_result << std::endl << "Index, original and decoded symbols (-1 denotes an erasure):";
        os_result << std::endl << "-----------------------------------------------------------" << std::endl;
        os_result << "_SIN "; print_vector<unsigned int, unsigned int>(symbol_indices, 4, os_result); os_result << std::endl;
        os_result << "_SOR "; print_vector<unsigned int, unsigned int>(options.prns, 4, os_result); os_result << std::endl;
        os_result << "_SRS "; print_vector<int, int>(decision_box.get_decoded_symbols(), 4, os_result); os_result << std::endl;
        std::cout << os_result.str() << std::endl;

        unsigned int erasure_counts = 0;
        unsigned int error_counts = 0;

        for (std::vector<int>::const_iterator it = decision_box.get_decoded_symbols().begin(); it != decision_box.get_decoded_symbols().end(); ++it)
        {
            if (*it < 0)
            {
                erasure_counts++;
            }
            else
            {
                if (*it != options.prns[it-decision_box.get_decoded_symbols().begin()])
                {
                    error_counts++;
                }
            }
        }

        std::cout << erasure_counts << " erasures (" << ((float) erasure_counts)/decision_box.get_decoded_symbols().size() << ")" << std::endl;
        std::cout << error_counts << " errors (" << ((float) error_counts)/decision_box.get_decoded_symbols().size() << ")" << std::endl;
        std::cout << erasure_counts+error_counts << " total (" << ((float) erasure_counts+error_counts)/decision_box.get_decoded_symbols().size() << ")" << std::endl;
        std::cout << "_SUM " << erasure_counts << "," << error_counts << "," << erasure_counts+error_counts << "," << erasure_counts+2*error_counts << std::endl;
    }
}
