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

#include "MultiplePrnCorrelator_FreqDep_Host.h"
#include "LocalCodes_Host.h"

#ifdef _CCSOFT
#include "CC_ReliabilityMatrix.h"
#include "CC_Encoding_base.h"
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
    MultiplePrnCorrelator_FreqDep *ifft_correlator = 0;
    
#ifdef _CCSOFT
    unsigned int tmp_n, tmp_k, tmp_m;
    ccsoft::get_cc_parameters<unsigned int>(options.cc_k_constraints, 
        options.cc_generator_polys,
        tmp_n,
        tmp_k,
        tmp_m);
    ccsoft::CC_ReliabilityMatrix cc_relmat(tmp_n, options.prns.size());
    std::ostringstream os;
    print_vector<unsigned int, unsigned int>(options.prns, 4, os); os << std::endl;
    std::cout << os.str() << std::endl;
#endif
    std::vector<PilotCorrelationRecord> correlation_records;

    std::vector<unsigned int> message_prn_numbers;
    generate_message_prn_list(message_prn_numbers);

    local_codes = new LocalCodes_Host(localCodeModulator,
    		gc_generator,
    		options.f_sampling,
    		options.f_chip,
    		message_prn_numbers); // make local codes time domain
    ifft_correlator = new MultiplePrnCorrelator_FreqDep_Host(gc_generator,
    		localCodeModulator,
    		options.f_sampling,
    		options.f_chip,
    		gc_generator.get_nb_message_codes(),
    		message_prn_numbers,
    		options.df_steps,
    		options.nb_prns_per_symbol,
    		options.batch_size,
    		options.f_step_division);
    // Do the correlation
    ExhaustivePrnCorrelator wsgce_correlator(local_codes, ifft_correlator);

    wsgc_complex *signal_samples;

    bool use_ccsoft = false;

#ifdef _CCSOFT
    if (options.fec_scheme == Options::OptionFEC_CCSoft)
    {
        use_ccsoft = true;
    }
#endif
    
    SampleSequencer sample_sequencer(faded_source_samples, nb_faded_source_samples, fft_N);
    unsigned int prn_i = 0;

    std::cout << "Do the correlation..." << std::endl;
    clock_gettime(time_option, &time1);

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
    
    // Do the decoding
    std::cout << "Do the decoding..." << std::endl;

#ifdef _CCSOFT
    if (use_ccsoft)
    {
        Reception_FEC::run_ccsoft_decoding(options, cc_relmat); // Do CCSoft decoding with CCSoft decision box
    }
    else
    {
        print_correlation_records(std::cout, correlation_records);
    }
#else
    print_correlation_records(std::cout, correlation_records);
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
    std::ostringstream oss;
    PilotCorrelationRecord::dump_oneline_banner(oss);
    std::vector<PilotCorrelationRecord>::const_iterator r_it = correlation_records.begin();

    for (; r_it != correlation_records.end(); ++r_it)
    {
        r_it->dump_oneline(oss);
    }

    os << oss.str();
}

