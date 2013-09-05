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
      
     Transmission

     Handles transmission operation including channel simulation
     
*/

#include "Transmission.h"
#include "Options.h"
#include "GoldCodeGenerator.h"
#include "SimulatedSource.h"
#include "CodeModulator_BPSK.h"
#include "CodeModulator_DBPSK.h"
#include "CodeModulator_OOK.h"
#include "CodeModulator_MFSK.h"
#include "SourceMixer.h"
#include "FIR_RCoef.h"
#include "FIRCoefGenerator.h"
#include "FadingModel.h"

#ifdef _RSSOFT
#include "RSSoft_Engine.h"
#endif

#ifdef _CCSOFT
#include "CCSoft_Engine.h"
#endif


//=================================================================================================
Transmission::Transmission(Options& _options, const GoldCodeGenerator& _gc_generator) :
    options(_options),
    gc_generator(_gc_generator),
    signal_samples(0),
    nb_signal_samples(0)
{}
    
//=================================================================================================
Transmission::~Transmission()
{
    if (signal_samples)
    {
        WSGC_FFTW_FREE(signal_samples);
    }
}
    
//=================================================================================================
void Transmission::generate_samples()
{
    apply_fec();

    if (options.transmission_scheme == Options::OptionTrans_WSGC)
    {
        CodeModulator_BPSK code_modulator;
        unsigned int nb_prns_per_symbol = options.nb_prns_per_symbol;
        unsigned int pilot_id = options.pilot1;
        unsigned int nb_pilot_prns = options.prns.size() + (options.batch_size/options.nb_prns_per_symbol);
        
        if (options.simulate_training) // only one "PRN per symbol" for training sequence
        {
            nb_prns_per_symbol = 1;
            pilot_id = options.pilot2;
            nb_pilot_prns = options.prns.size() + options.batch_size;
        }
        
        SimulatedSource message_source(gc_generator, options.prns, options.f_sampling, options.f_chip,
                                       options.f_tx, options.code_shift, nb_prns_per_symbol, 0.0);
        message_source.set_code_modulator(&code_modulator);  
        
        std::vector<unsigned int> pilot_prns;
        pilot_prns.assign(nb_pilot_prns, pilot_id);
        SimulatedSource pilot_source(gc_generator, pilot_prns, options.f_sampling, options.f_chip,
                                     options.f_tx, options.code_shift, nb_prns_per_symbol, 0.0);
        pilot_source.set_code_modulator(&code_modulator); 
        
        wsgc_float pilot_gain = pow(10.0, (options.pilot_gain_db / 10.0));
        SourceMixer source_mixer(message_source, pilot_source, pilot_gain);
        source_mixer.get_samples(&signal_samples);
        nb_signal_samples = source_mixer.get_nb_samples();
    }
    else if (options.transmission_scheme == Options::OptionTrans_WSGCD)
    {
        CodeModulator_DBPSK code_modulator;
        
        SimulatedSource message_source(gc_generator, options.prns, options.f_sampling, options.f_chip,
                                       options.f_tx, options.code_shift, options.nb_prns_per_symbol, 0.0);
        message_source.set_code_modulator(&code_modulator);  
        message_source.create_samples(&signal_samples);
        nb_signal_samples = message_source.get_nb_samples();
    }
    else if (options.transmission_scheme == Options::OptionTrans_WSGCO)
    {
        CodeModulator_OOK code_modulator;
        
        SimulatedSource message_source(gc_generator, options.prns, options.f_sampling, options.f_chip,
                                       options.f_tx, options.code_shift, options.nb_prns_per_symbol, 0.0);
        message_source.set_code_modulator(&code_modulator);  
        message_source.create_samples(&signal_samples);
        nb_signal_samples = message_source.get_nb_samples();
    }
    else if (options.transmission_scheme == Options::OptionTrans_MFSK)
    {
        CodeModulator_MFSK codeModulator_MFSK(
            options.mfsk_options._f_sampling,
            options.mfsk_options._zero_frequency + options.f_tx,
            options.mfsk_options._symbol_bandwidth,
            options.mfsk_options._symbol_time);
            
        options.prns.push_back(0); // pad with one symbol
        nb_signal_samples = codeModulator_MFSK.get_nb_symbol_samples()*options.prns.size();
        signal_samples = (wsgc_complex *) WSGC_FFTW_MALLOC(nb_signal_samples*sizeof(wsgc_fftw_complex));
        codeModulator_MFSK.modulate(reinterpret_cast<wsgc_fftw_complex*>(signal_samples), options.prns);
    }
    
    if (signal_samples)
    {
        // Apply lowpass filter if any
        if (options._fir_coef_generator != 0)
        {
            apply_fir(signal_samples, nb_signal_samples, options._fir_coef_generator->get_coefs());
        }

        apply_channel();
    }
}

//=================================================================================================
void Transmission::normalize()
{
    if (signal_samples)
    {
        FadingModel *fading = options._fading_model;
        fading->normalize(signal_samples, nb_signal_samples);
    }
}

//=================================================================================================
void Transmission::apply_fir(wsgc_complex *inout, unsigned int& nb_samples, const std::vector<wsgc_float>& fir_coef)
{
	std::cout << "Apply lowpass FIR filter" << std::endl;

	FIR_RCoef fir_filter(fir_coef);
	static const wsgc_complex c_zero = (0.0, 0.0);

	for (unsigned int i = 0; i<nb_samples; i++)
	{
		inout[i] = fir_filter.calc(inout[i]);
	}
}
    
//=================================================================================================
void Transmission::apply_channel()
{
    // get fading model
    FadingModel *fading = options._fading_model;
    wsgc_complex *faded_signal_samples;
    unsigned int nb_faded_signal_samples; // fading may add delay hence there could be more samples after the fading process
    
    // apply fading
    if (fading->is_fading_active())
    {
        std::cout << "Apply fading" << std::endl;
        nb_faded_signal_samples = fading->get_output_size(nb_signal_samples);
        faded_signal_samples = (wsgc_complex *) WSGC_FFTW_MALLOC(nb_faded_signal_samples*sizeof(wsgc_fftw_complex));
        fading->apply_fading(signal_samples, faded_signal_samples, nb_signal_samples);
        WSGC_FFTW_FREE(signal_samples);
        signal_samples = faded_signal_samples;
        nb_signal_samples = nb_faded_signal_samples;
    }
    
    // apply AWGN
    if (options.make_noise)
    {
        std::cout << "Apply AWGN" << std::endl;
        fading->apply_awgn(signal_samples, nb_signal_samples, options.code_shift, options.snr);
    }
}

//=================================================================================================
void Transmission::apply_fec()
{
    if (options.fec_scheme == Options::OptionFEC_RSSoft)
    {
#ifdef _RSSOFT
        if (1<<options.rs_logq != options.nb_message_symbols)
        {
            std::cout << "Invalid size of symbols alphabet (" << options.nb_message_symbols << ") for the Reed-Solomon parameters RS(" << (1<<options.rs_logq) - 1 << "," << options.rs_k << ")" << std::endl;
            std::cout << "FEC ignored" << std::endl;
            return;
        }

        if (options.prns.size() != options.rs_k)
        {
            std::cout << "Invalid number of message symbols for the Reed-Solomon parameters" << std::endl;
            std::cout << "FEC ignored" << std::endl;
            return;
        }

        if (options.rs_k > options.nb_message_symbols - 2)
        {
            std::cout << "Invalid k parameter for Reed-Solomon" << std::endl;
            std::cout << "FEC ignored" << std::endl;
            return;
        }

        RSSoft_Engine rssoft_engine(options.rs_logq, options.rs_k);
        rssoft_engine.set_initial_global_multiplicity(options.rs_init_M);
        rssoft_engine.set_nb_retries(options.rs_r);
        rssoft_engine.set_retry_base_increment(options.rs_inc);
        rssoft_engine.set_retry_mm_strategy(options.rs_inc_strategy);

        options.source_prns = options.prns;
        options.prns.clear();

        rssoft_engine.encode(options.source_prns, options.prns);
#else
        std::cout << "Program not linked with RSSoft library, ignoring FEC with RSSoft" << std::endl;
#endif
    }
    else if (options.fec_scheme == Options::OptionFEC_CCSoft)
    {
#ifdef _CCSOFT
        CCSoft_Engine ccsoft_engine(options.cc_k_constraints, options.cc_generator_polys);

        switch (options.cc_algorithm_type)
        {
        case CCSoft_Engine_defs::Algorithm_Stack:
            ccsoft_engine.init_decoding_stack(options.cc_edge_bias);
            break;
        case CCSoft_Engine_defs::Algorithm_Fano:
            ccsoft_engine.init_decoding_fano(options.cc_edge_bias, options.cc_fano_init_metric, options.cc_fano_delta_metric, options.cc_fano_tree_cache_size, options.cc_fano_delta_init_threshold);
            break;
        default:
            std::cout << "Unrecognised convolutional coding decoding algorithm" << std::endl;
            std::cout << "FEC ignored" << std::endl;
            return;
        }

        if (options.cc_use_node_limit)
        {
            ccsoft_engine.set_node_limit(options.cc_node_limit);
        }

        if (options.cc_use_metric_limit)
        {
            ccsoft_engine.set_metric_limit(options.cc_metric_limit);
        }

        if (1<<(ccsoft_engine.get_k()) != options.nb_message_symbols)
        {
            std::cout << "Invalid size of symbols alphabet (" << options.nb_message_symbols << ") for a (n,k,m) CC code with k = " << ccsoft_engine.get_k() << std::endl;
            std::cout << "FEC ignored" << std::endl;
            return;
        }

        options.source_prns = options.prns;
        ccsoft_engine.zero_pad(options.source_prns);
        options.prns.clear();
        ccsoft_engine.encode(options.source_prns, options.prns);
#else
        std::cout << "Program not linked with CCSoft library, ignoring FEC with CCSoft" << std::endl;
#endif
    }
}
