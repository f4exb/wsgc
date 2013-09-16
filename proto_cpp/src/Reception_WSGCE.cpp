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

#include "Reception_WSGCE.h"
#include "Reception_FEC.h"
#include "WsgcUtils.h"
#include "LocalCodes.h"
#include "SampleSequencer.h"
#include "ExhaustivePrnCorrelator.h"

#ifdef _CUDA
#include "SinglePrnCorrelator_FreqDep_Cuda.h"
#include "LocalCodes_Cuda.h"
#endif
#include "SinglePrnCorrelator_FreqDep_Host.h"
#include "LocalCodes_Host.h"

#ifdef _CCSOFT
#include "CC_ReliabilityMatrix.h"
#endif
#include "PilotCorrelationRecord.h"

#include <iostream>
#include <sstream>
#include <time.h>

//=================================================================================================
Reception_WSGCE::Reception_WSGCE(Options& _options, const GoldCodeGenerator& _gc_generator) :
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
Reception_WSGCE::~Reception_WSGCE()
{}

//=================================================================================================
void Reception_WSGCE::message_processing(wsgc_complex *faded_source_samples, unsigned int nb_faded_source_samples)
{
    timespec time1, time2;
    int time_option = CLOCK_REALTIME;
    LocalCodes *local_codes = 0;
    SinglePrnCorrelator_FreqDep *ifft_correlator = 0;
    
#ifdef _CCSOFT
    ccsoft::CC_ReliabilityMatrix cc_relmat(options.nb_message_symbols, options.prns.size()-1);
#endif
    std::vector<PilotCorrelationRecord> correlation_records;

    std::vector<unsigned int> message_prn_numbers;
    generate_message_prn_list(message_prn_numbers);

#ifdef _CUDA
    if (options.use_cuda)
    {
        std::cout << "!!! USING CUDA !!!" << std::endl;
        unsigned int cuda_device = cuda_manager.get_pilot_device();
        local_codes = new LocalCodes_Cuda(localCodeModulator,
        		gc_generator,
        		options.f_sampling,
        		options.f_chip,
        		message_prn_numbers,
        		cuda_device); // make local codes time domain
        ifft_correlator = new SinglePrnCorrelator_FreqDep_Cuda(gc_generator,
        		localCodeModulator,
        		options.f_sampling,
        		options.f_chip,
        		message_prn_numbers,
        		options.df_steps,
        		cuda_device,
        		options.nb_prns_per_symbol,
        		options.batch_size,
        		options.f_step_division);
    }
    else
    {
        local_codes = new LocalCodes_Host(localCodeModulator,
        		gc_generator,
        		options.f_sampling,
        		options.f_chip,
        		message_prn_numbers); // make local codes time domain
        ifft_correlator = new SinglePrnCorrelator_FreqDep_Host(gc_generator,
        		localCodeModulator,
        		options.f_sampling,
        		options.f_chip,
        		message_prn_numbers,
        		options.df_steps,
        		options.nb_prns_per_symbol,
        		options.batch_size,
        		options.f_step_division);
    }
#else
    local_codes = new LocalCodes_Host(localCodeModulator,
    		gc_generator,
    		options.f_sampling,
    		options.f_chip,
    		message_prn_numbers); // make local codes time domain
    ifft_correlator = new SinglePrnCorrelator_FreqDep_Host(gc_generator,
    		localCodeModulator,
    		options.f_sampling,
    		options.f_chip,
    		message_prn_numbers,
    		options.df_steps,
    		options.nb_prns_per_symbol,
    		options.batch_size,
    		options.f_step_division);
#endif
    // Do the correlation
    ExhaustivePrnCorrelator wsgce_correlator(local_codes, ifft_correlator);

    wsgc_complex *signal_samples;

    std::cout << "Do the correlation..." << std::endl;
    clock_gettime(time_option, &time1);

    bool use_ccsoft = false;

#ifdef _CCSOFT
    if (options.fec_scheme == Options::OptionFEC_CCSoft)
    {
        use_ccsoft = true;
    }
#endif
    
    SampleSequencer sample_sequencer(faded_source_samples, nb_faded_source_samples, fft_N);
    unsigned int prn_i = 0;

    while (sample_sequencer.get_next_code_samples(&signal_samples)) // pseudo real time loop, one PRN length at a time
    {
#ifdef _CCSOFT
        if (use_ccsoft)
        {
            wsgce_correlator.make_correlation(signal_samples, cc_relmat);
        }
        else
        {
            wsgce_correlator.make_correlation(signal_samples, correlation_records);
        }
#else
        wsgce_correlator.make_correlation(signal_samples, correlation_records);
#endif
    }
    
    clock_gettime(time_option, &time2);
    std::cout << "Message correlation time: " << std::setw(12) << std::setprecision(9) << WsgcUtils::get_time_difference(time2,time1) << " s" << std::endl << std::endl;
    
    // Do the decoding with the decision box
    std::cout << "Do the decoding with the decision box..." << std::endl;

    std::ostringstream corr_os;

#ifdef _CCSOFT
    if (use_ccsoft)
    {
        Reception_FEC::run_ccsoft_decoding(options, cc_relmat); // Do CCSoft decoding with CCSoft decision box
    }
    else
    {
        print_correlation_records(corr_os, correlation_records);
    }
#else
    print_correlation_records(corr_os, correlation_records);
#endif

    if (ifft_correlator)
    {
        delete ifft_correlator;
    }

    if (local_codes)
    {
        delete local_codes;
    }
}


//=================================================================================================
void Reception_WSGCE::print_correlation_records(std::ostream& os, const std::vector<PilotCorrelationRecord>& correlation_records)
{

}

