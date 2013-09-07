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

     Reception

     Common part of all recepion classes.

*/

#include "Reception_FEC.h"
#include "Demodulation.h"
#include "Reception.h"
#include "LocalCodesFFT_Host.h"
#include "UnpilotedMultiplePrnCorrelator_Host.h"
#include "UnpilotedTrainingMessageCorrelator_Host.h"
#include "DecisionBox_Unpiloted_And_Synced.h"
#include "SampleSequencer.h"
#include "WsgcUtils.h"

#include <iostream>

//=================================================================================================
void Reception::unpiloted_message_correlation(CodeModulator& localCodeModulator,
    wsgc_complex *faded_source_samples,
    unsigned int nb_faded_source_samples,
    std::vector<CorrelationRecord>& correlation_records)
{
    timespec time1, time2;
    int time_option = CLOCK_REALTIME;

    Demodulation::demodulate_generic(options, gc_generator, faded_source_samples, nb_faded_source_samples);

    std::vector<unsigned int> message_prn_numbers;
    generate_message_prn_list(message_prn_numbers);

    std::cout << "Unpiloted correlation with frequency independant modulations" << std::endl;

    // only host for now
    LocalCodesFFT_Host local_codes_fft(localCodeModulator, gc_generator, options.f_sampling, options.f_chip, message_prn_numbers); // make local codes time domain
    UnpilotedMultiplePrnCorrelator_Host dm_correlator(
            options.f_sampling,
            options.f_chip,
            gc_generator.get_code_length(),
            options.nb_prns_per_symbol,
            message_prn_numbers,
            options.batch_size, // will be taken as number of symbols
            options.analysis_window_size,
            correlation_records,
            local_codes_fft);

    // Do the correlation
    std::cout << "Do the correlation..." << std::endl;
    wsgc_complex *signal_samples;
    SampleSequencer sample_sequencer(faded_source_samples, nb_faded_source_samples, fft_N);
    clock_gettime(time_option, &time1);

    while (sample_sequencer.get_next_code_samples(&signal_samples)) // pseudo real time loop, one PRN length at a time
    {
        if (dm_correlator.set_samples(signal_samples))
        {
            dm_correlator.execute();
        }
    }

    clock_gettime(time_option, &time2);
    std::cout << "Message correlation time: " << std::setw(12) << std::setprecision(9) << WsgcUtils::get_time_difference(time2,time1) << " s" << std::endl << std::endl;

    // Print some stats
    std::cout << std::endl;
    std::ostringstream corr_os;
    corr_os << "Message correlation analysis results:" << std::endl;
    dm_correlator.dump_correlation_records(corr_os, fft_N);
    corr_os << "Time shift analysis results:" << std::endl;
    dm_correlator.dump_time_analyzer_results(corr_os);
    corr_os << std::endl;
    std::cout << corr_os.str() << std::endl;
}


//=================================================================================================
void Reception::unpiloted_training_correlation(CodeModulator& localCodeModulator,
        wsgc_complex *faded_source_samples,
        unsigned int nb_faded_source_samples,
        std::vector<TrainingCorrelationRecord>& training_correlation_records)
{
    timespec time1, time2;
    int time_option = CLOCK_REALTIME;

    Demodulation::demodulate_generic(options, gc_generator, faded_source_samples, nb_faded_source_samples);

    std::vector<unsigned int> training_prn_numbers;
    generate_training_prn_list_unpiloted(training_prn_numbers);

    unsigned int prn_length = gc_generator.get_code_length();

    LocalCodesFFT_Host local_codes_fft(
            localCodeModulator,
            gc_generator,
            options.f_sampling,
            options.f_chip,
            training_prn_numbers);

    UnpilotedTrainingMessageCorrelator_Host unpiloted_training_message_correlator(
            options.f_sampling,
            options.f_chip,
            prn_length,
            options.prns.size(),
            options.analysis_window_size,
            training_prn_numbers,
            training_correlation_records,
            local_codes_fft);

    SampleSequencer sample_sequencer(faded_source_samples, nb_faded_source_samples, fft_N);
    wsgc_complex *signal_samples;

    clock_gettime(time_option, &time1);

    while (sample_sequencer.get_next_code_samples(&signal_samples)) // pseudo real time loop, one PRN length at a time
    {
        unpiloted_training_message_correlator.set_samples(signal_samples);
        unpiloted_training_message_correlator.execute();
    }

    clock_gettime(time_option, &time2);
    std::cout << "Training sequence correlation time: " << std::setw(12) << std::setprecision(9) << WsgcUtils::get_time_difference(time2,time1) << " s" << std::endl << std::endl;
}


//=================================================================================================
void Reception::unpiloted_and_synced_decision(std::vector<CorrelationRecord>& correlation_records)
{
    std::ostringstream corr_os;
    std::cout << "Unpiloted correlation" << std::endl;

    DecisionBox_Unpiloted_And_Synced decision_box(options.nb_prns_per_symbol, fft_N, options.decision_thresholds, correlation_records);
    decision_box.set_mag_display_adj_factor(fft_N);

    decision_box.analyze_records();

    if (decision_box.is_prni_at_max_invalid())
    {
        std::cout << "Symbol boundaries cannot be estimated with confidence - message is discarded" << std::endl;
    }
    else
    {
        decision_box.estimate_symbols();

        std::vector<unsigned int> symbol_indices;

        for (unsigned int i=0; i<options.prns.size(); i++)
        {
            symbol_indices.push_back(i);
        }

        std::ostringstream os_result;
        os_result << std::endl;
        options.decision_thresholds.print_options(os_result);
        os_result << std::endl << "Decision box records:" << std::endl;
        decision_box.dump_decision_records(os_result);
        os_result << std::endl << "Decisions status:" << std::endl;
        decision_box.dump_decision_status(os_result, options.prns);
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

//=================================================================================================
void Reception::generate_training_prn_list(std::vector<unsigned int>& prn_list)
{
	for (unsigned int prni=0; prni < gc_generator.get_nb_message_codes(); prni++)
	{
		prn_list.push_back(prni);
	}
}

//=================================================================================================
void Reception::generate_training_prn_list_unpiloted(std::vector<unsigned int>& prn_list)
{
	unsigned int start_code = gc_generator.get_nb_message_codes() + gc_generator.get_nb_service_codes();

	for (unsigned int prni=start_code; prni < start_code + gc_generator.get_nb_training_codes(); prni++)
	{
		prn_list.push_back(prni);
	}
}

//=================================================================================================
void Reception::generate_message_prn_list(std::vector<unsigned int>& prn_list)
{
	for (unsigned int prni=0; prni <= gc_generator.get_nb_message_codes(); prni++) // incudes the noise PRN that is the last hence the <= comparison
	{
		prn_list.push_back(prni);
	}
}


//=================================================================================================
void Reception::generate_pilot_prn_list(std::vector<unsigned int>& prn_list, unsigned int pilot_prni)
{
    prn_list.push_back(gc_generator.get_nb_message_codes()+pilot_prni);
}

