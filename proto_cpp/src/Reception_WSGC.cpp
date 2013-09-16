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
      
     Reception WSGC: Correlation, BPSK, Pilot
*/

#include "Reception_WSGC.h"
#include "Reception_FEC.h"
#include "WsgcUtils.h"
#include "LocalCodes.h"
#include "PilotCorrelator.h"
#include "PilotedMessageCorrelator.h"
#include "PilotCorrelationAnalyzer.h"
#include "PilotedMultiplePrnCorrelator.h"
#include "DecisionBox.h"
#include "DecisionBox_Piloted_And_Synced.h"
#include "SampleSequencer.h"
#include "PilotedTrainingMultiplePrnCorrelator.h"

#ifdef _CUDA
#include "SinglePrnCorrelator_FreqDep_Cuda.h"
#include "PilotedMessageCorrelator_Cuda.h"
#include "PilotedTrainingMessageCorrelator_Cuda.h"
#include "LocalCodes_Cuda.h"
#include "PilotCorrelator_Cuda.h"
#endif
#include "SinglePrnCorrelator_FreqDep_Host.h"
#include "PilotedMessageCorrelator_Host.h"
#include "PilotedTrainingMessageCorrelator_Host.h"
#include "LocalCodes_Host.h"
#include "PilotCorrelator_Host.h"

#ifdef _RSSOFT
#include "RS_ReliabilityMatrix.h"
#ifdef _CUDA
#include "RSSoft_ReliabilityMatrixCuda.h"
#endif
#endif

#include <iostream>
#include <sstream>
#include <time.h>

//=================================================================================================
Reception_WSGC::Reception_WSGC(Options& _options, const GoldCodeGenerator& _gc_generator) :
    Reception(_options, _gc_generator)
#ifdef _CUDA
    ,cuda_manager(
        gc_generator.get_nb_message_codes(),
        options.nb_pilot_prns,
        gc_generator.get_nb_code_samples(options.f_sampling, options.f_chip),
        options.batch_size,
        options.df_steps,
        options.nb_prns_per_symbol,
        options.f_step_division)
#endif
{
#ifdef _CUDA
    if (options.gpu_affinity_specified)
    {
        cuda_manager.set_gpu_affinity(options.gpu_affinity);
    }
    bool use_cuda_ok = cuda_manager.diagnose();
    options.use_cuda = options.use_cuda && use_cuda_ok;
    if (options.use_cuda)
    {
        std::ostringstream cuda_os;
        cuda_manager.dump(cuda_os);
        std::cout << cuda_os.str() << std::endl << std::endl;
    }
#endif
}

//=================================================================================================
Reception_WSGC::~Reception_WSGC()
{}

//=================================================================================================
void Reception_WSGC::message_processing(wsgc_complex *faded_source_samples, unsigned int nb_faded_source_samples)
{
    timespec time1, time2;
    int time_option = CLOCK_REALTIME;
    LocalCodes *local_codes = 0;
    const std::map<unsigned int, unsigned int> *prn_shift_occurences = 0;
    PilotCorrelator *pilot_correlator = 0;
    PilotedMessageCorrelator *message_correlator = 0;
    std::vector<CorrelationRecord> correlation_records;
    
#ifdef _RSSOFT
#ifdef _CUDA
    RSSoft_ReliabilityMatrixCuda *rssoft_reliability_matrix_cuda = 0;
#endif
    rssoft::RS_ReliabilityMatrix rs_relmat(options.rs_logq, (1<<options.rs_logq) - 1);
#endif

    std::vector<unsigned int> message_prn_numbers;
    std::vector<unsigned int> pilot_prn_numbers;
    
    generate_message_prn_list(message_prn_numbers);
    generate_pilot_prn_list(pilot_prn_numbers, 1);

    // Modulation should be frequency dependant. It doesn't make much sense if there is no frequency tracking necessary.
    // This pilot correlation works exclusively in the two dimensions of (frequency shift, PRN time shift).
    // this means it works only with frequency dependant modulations
    // For now this concerns BPSK modulation with pilot use
    
    std::cout << "Creating piloted multiple PRN correlator" << std::endl;
    PilotCorrelationAnalyzer pilot_correlation_analyzer(options.analysis_window_size, options.nb_prns_per_symbol, options.nb_samples_per_code);
#ifdef _CUDA
    if (options.use_cuda)
    {
        std::cout << "!!! USING CUDA !!!" << std::endl;
        unsigned int cuda_device = cuda_manager.get_pilot_device();
        local_codes = new LocalCodes_Cuda(localCodeModulator, gc_generator, options.f_sampling, options.f_chip, message_prn_numbers, cuda_device); // make local codes time domain
        pilot_correlator = new PilotCorrelator_Cuda(gc_generator, localCodeModulator, options.f_sampling, options.f_chip, pilot_prn_numbers, cuda_device, options.nb_prns_per_symbol, options.df_steps, options.batch_size, options.f_step_division);
        message_correlator = new PilotedMessageCorrelator_Cuda(*((LocalCodes_Cuda *) local_codes), options.f_sampling, options.f_chip, options.nb_prns_per_symbol, cuda_device);
#ifdef _RSSOFT
        // If Reed-Solomon is active, allocate and set RSSoft cuda reliability matrix 
        if (options.fec_scheme == Options::OptionFEC_RSSoft)
        {
            rssoft_reliability_matrix_cuda = new RSSoft_ReliabilityMatrixCuda(options.rs_logq, (1<<options.rs_logq)-1);
            ((PilotedMessageCorrelator_Cuda *)message_correlator)->set_reliability_matrix_cuda(rssoft_reliability_matrix_cuda);
        }
#endif
    }
    else
    {
        local_codes = new LocalCodes_Host(localCodeModulator, gc_generator, options.f_sampling, options.f_chip, message_prn_numbers); // make local codes time domain
        pilot_correlator = new PilotCorrelator_Host(gc_generator, localCodeModulator, options.f_sampling, options.f_chip, pilot_prn_numbers, options.nb_prns_per_symbol, options.df_steps, options.batch_size, options.f_step_division);
        message_correlator = new PilotedMessageCorrelator_Host(*((LocalCodes_Host *) local_codes), options.f_sampling, options.f_chip, options.nb_prns_per_symbol);
    }
#else
    local_codes = new LocalCodes_Host(localCodeModulator, gc_generator, options.f_sampling, options.f_chip, message_prn_numbers); // make local codes time domain
    pilot_correlator = new PilotCorrelator_Host(gc_generator, localCodeModulator, options.f_sampling, options.f_chip, pilot_prn_numbers, options.nb_prns_per_symbol, options.df_steps, options.batch_size, options.f_step_division);
    message_correlator = new PilotedMessageCorrelator_Host(*((LocalCodes_Host *) local_codes), options.f_sampling, options.f_chip, options.nb_prns_per_symbol);
#endif
    PilotedMultiplePrnCorrelator piloted_mprn_correlator(pilot_correlation_analyzer, correlation_records, *pilot_correlator, *message_correlator);

    // Do the correlation
    std::cout << "Do the correlation..." << std::endl;

    wsgc_complex *signal_samples;

    clock_gettime(time_option, &time1);
    
    // Support pipeline processing if necessary (with pilot assisted correlation using CUDA)
    bool input_samples_available = false;           // input samples (length of one PRN) are available for processing  
    bool process_next_output = false;
    
    SampleSequencer sample_sequencer(faded_source_samples, nb_faded_source_samples, fft_N);
    unsigned int prn_i = 0;

    while ((input_samples_available = sample_sequencer.get_next_code_samples(&signal_samples)) || process_next_output) // pseudo real time loop, one PRN length at a time
    {
        if (input_samples_available) // process input
        {
            // Use PilotedMultiplePrnCorrelator specialized class that encompasses pilot and message correlators and eventually produces correlation records usable directly in the Decision Box
            //   - set source block
            //   - execute the correlation / averaging process
            piloted_mprn_correlator.set_source_block(reinterpret_cast<wsgc_fftw_complex *>(signal_samples));
            piloted_mprn_correlator.make_correlation(0); // use the only pilot (code index 0)
            
        }
    }
    
    clock_gettime(time_option, &time2);
    std::cout << "Message correlation time: " << std::setw(12) << std::setprecision(9) << WsgcUtils::get_time_difference(time2,time1) << " s" << std::endl << std::endl;
    
    // Do the decoding with the decision box
    std::cout << "Do the decoding with the decision box..." << std::endl;

    std::ostringstream corr_os;

    corr_os << "Piloted correlation" << std::endl << std::endl;
    corr_os << "Pilot correlation analysis results:" << std::endl;
    pilot_correlation_analyzer.dump_histo_time_shift_occurences(corr_os);
    corr_os << std::endl;

#ifdef _CUDA
    if (options.use_cuda)
    {
        pilot_correlation_analyzer.set_pilot_mag_display_factor((fft_N/2)*(fft_N/2));
        pilot_correlation_analyzer.set_message_mag_display_factor(fft_N/2);
    }
#endif

    corr_os << "--- pilot correlation records:" << std::endl;
    pilot_correlation_analyzer.dump_pilot_correlation_records(corr_os);
    corr_os << std::endl << "--- correlation records:" << std::endl;
    pilot_correlation_analyzer.dump_message_correlation_records(corr_os);
    std::cout << corr_os.str() << std::endl;
    
#ifdef _RSSOFT
    if (options.fec_scheme == Options::OptionFEC_RSSoft)
    {
#ifdef _CUDA
        if (rssoft_reliability_matrix_cuda)
        {
            rssoft_reliability_matrix_cuda->copy_to_host(rs_relmat); // move reliability data from device to host
        }
#endif         
        Reception_FEC::run_rssoft_decoding(options, rs_relmat); // Do RSSoft decoding with RSSoft decision box
    }
#endif

    DecisionBox_Piloted_And_Synced decision_box(options.nb_prns_per_symbol, fft_N, options.decision_thresholds, pilot_correlation_analyzer);

#ifdef _CUDA
    if (options.use_cuda)
    {
        decision_box.set_mag_display_adj_factor(fft_N/2);
        decision_box.set_use_cuda(true);
        if (!options.decision_thresholds_specified)
        {
            options.decision_thresholds.set_cuda_defaults();
        }
    }
#endif
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
    
#if _RSSOFT
#if _CUDA
    if (rssoft_reliability_matrix_cuda)
    {
        delete rssoft_reliability_matrix_cuda;
    }
#endif
#endif    
    if (message_correlator)
    {
        delete message_correlator;
    }

    if (pilot_correlator)
    {
        delete pilot_correlator;
    }

    if (local_codes)
    {
        delete local_codes;
    }
}

//=================================================================================================
void Reception_WSGC::training_processing(wsgc_complex *faded_source_samples, unsigned int nb_faded_source_samples)
{
    timespec time1, time2;
    int time_option = CLOCK_REALTIME;
    LocalCodes *local_codes = 0;
    PilotCorrelator *pilot_correlator = 0;
    PilotedTrainingMessageCorrelator *tr_message_correlator = 0;
    
    // - Modulation should support code division
    // - Modulation processing should be frequency dependant
    // - There should be at least 2 pilot PRN(s) - using second
    if (options.nb_pilot_prns > 1)
    {
        std::vector<unsigned int> pilot_prn_numbers;
        generate_pilot_prn_list(pilot_prn_numbers, 2); // Pilot PRN #2 is used

        std::vector<unsigned int> training_prn_numbers;
        generate_training_prn_list(training_prn_numbers);

        std::cout << "Creating managing objects..." << std::endl;
        PilotCorrelationAnalyzer pilot_correlation_analyzer(options.analysis_window_size, options.nb_prns_per_symbol, options.nb_samples_per_code);    
#ifdef _CUDA
        if (options.use_cuda)
        {
            pilot_correlation_analyzer.set_pilot_mag_display_factor((fft_N/2)*(fft_N/2));
            pilot_correlation_analyzer.set_training_mag_display_factor(fft_N/8);
            unsigned int cuda_device = cuda_manager.get_pilot_device();
            local_codes = new LocalCodes_Cuda(localCodeModulator, gc_generator, options.f_sampling, options.f_chip, training_prn_numbers, cuda_device); // make local codes time domain
            pilot_correlator = new PilotCorrelator_Cuda(gc_generator, localCodeModulator, options.f_sampling, options.f_chip, pilot_prn_numbers, cuda_device, options.nb_prns_per_symbol, options.df_steps, options.batch_size, options.f_step_division);
        	tr_message_correlator = new PilotedTrainingMessageCorrelator_Cuda(*((LocalCodes_Cuda *) local_codes), options.f_sampling, options.f_chip, pilot_correlation_analyzer.get_analysis_window_size_in_prns(), options.nb_random_prns, cuda_device);
        }
        else
        {
            local_codes = new LocalCodes_Host(localCodeModulator, gc_generator, options.f_sampling, options.f_chip, training_prn_numbers); // make local codes time domain
            pilot_correlator = new PilotCorrelator_Host(gc_generator, localCodeModulator, options.f_sampling, options.f_chip, pilot_prn_numbers, options.nb_prns_per_symbol, options.df_steps, options.batch_size, options.f_step_division);
        	tr_message_correlator = new PilotedTrainingMessageCorrelator_Host(*((LocalCodes_Host *) local_codes), options.f_sampling, options.f_chip, options.nb_random_prns);
        }
#else
        local_codes = new LocalCodes_Host(localCodeModulator, gc_generator, options.f_sampling, options.f_chip, training_prn_numbers); // make local codes time domain
        pilot_correlator = new PilotCorrelator_Host(gc_generator, localCodeModulator, options.f_sampling, options.f_chip, pilot_prn_numbers, options.nb_prns_per_symbol, options.df_steps, options.batch_size, options.f_step_division);
        tr_message_correlator = new PilotedTrainingMessageCorrelator_Host(*((LocalCodes_Host *) local_codes), options.f_sampling, options.f_chip, options.nb_random_prns);
#endif
        PilotedTrainingMultiplePrnCorrelator piloted_tr_mprn_correlator(pilot_correlation_analyzer, *pilot_correlator, *tr_message_correlator);
        
        // Do the correlation
        std::cout << "Do the correlation..." << std::endl;
        wsgc_complex *signal_samples;
        SampleSequencer sample_sequencer(faded_source_samples, nb_faded_source_samples, fft_N);
        unsigned int prn_i = 0;

        clock_gettime(time_option, &time1);
        
        while (sample_sequencer.get_next_code_samples(&signal_samples)) // pseudo real time loop, one PRN length at a time
        {
            piloted_tr_mprn_correlator.set_source_block(reinterpret_cast<wsgc_fftw_complex *>(signal_samples));
            piloted_tr_mprn_correlator.make_correlation(0); // use the only pilot (code index 0)
        }
        
        clock_gettime(time_option, &time2);
        std::cout << "Training sequence correlation time: " << std::setw(12) << std::setprecision(9) << WsgcUtils::get_time_difference(time2,time1) << " s" << std::endl << std::endl;

        std::ostringstream corr_os;

        corr_os << "Pilot correlation analysis results:" << std::endl;
        pilot_correlation_analyzer.dump_histo_time_shift_occurences(corr_os);
        corr_os << std::endl;

        corr_os << "--- pilot correlation records:" << std::endl;
        pilot_correlation_analyzer.dump_pilot_correlation_records(corr_os);

        corr_os << "--- training correlation records:" << std::endl;
        pilot_correlation_analyzer.dump_training_correlation_records(corr_os);

        std::cout << corr_os.str() << std::endl;
        
        // TODO: conclude on the message epoch
    }    
    else
    {
        std::cout << "Need 2 pilot PRNs as using second for training identification" << std::endl;
    }
    
    if (tr_message_correlator)
    {
        delete tr_message_correlator;
    }
    
    if (pilot_correlator)
    {
        delete pilot_correlator;
    }
    
    if (local_codes)
    {
        delete local_codes;
    }
}

