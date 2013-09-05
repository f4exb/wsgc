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
      
     Reception class specialized in FEC processing
*/

#include "Reception_FEC.h"
#include "Options.h"

#ifdef _RSSOFT
#include "RSSoft_Engine.h"
#include "RSSoft_DecisionBox.h"
#include "RS_ReliabilityMatrix.h"
#endif

#ifdef _CCSOFT
#include "CCSoft_Engine.h"
#include "CCSoft_DecisionBox.h"
#include "CC_ReliabilityMatrix.h"
#endif

#ifdef _RSSOFT
//=================================================================================================
void Reception_FEC::run_rssoft_decoding(Options& options, rssoft::RS_ReliabilityMatrix& relmat)
{
    RSSoft_Engine rssoft_engine(options.rs_logq, options.rs_k);
    rssoft_engine.set_initial_global_multiplicity(options.rs_init_M);
    rssoft_engine.set_nb_retries(options.rs_r);
    rssoft_engine.set_retry_base_increment(options.rs_inc);
    rssoft_engine.set_retry_mm_strategy(options.rs_inc_strategy);
    
    RSSoft_DecisionBox rssoft_decision_box(rssoft_engine, options._source_codec, relmat);

    switch (options.rs_decoding_mode)
    {
        case RSSoft_Engine_defs::RSSoft_decoding_all:
        case RSSoft_Engine_defs::RSSoft_decoding_full:
        case RSSoft_Engine_defs::RSSoft_decoding_best:
        case RSSoft_Engine_defs::RSSoft_decoding_first:
            rssoft_decision_box.run(options.rs_decoding_mode);
            break;
        case RSSoft_Engine_defs::RSSoft_decoding_regex:
            rssoft_decision_box.run_regex(options.rs_decoding_match_str);
            break;
        case RSSoft_Engine_defs::RSSoft_decoding_match:
            rssoft_decision_box.run_match(options.rs_decoding_match_str);
            break;
        case RSSoft_Engine_defs::RSSoft_decoding_binmatch:
            rssoft_decision_box.run_match(options.source_prns);
            break;
        case RSSoft_Engine_defs::RSSoft_decoding_relthr:
            rssoft_decision_box.run_reliability_threshold(options.rs_reliability_threshold);
            break;
        default:
            std::cout << "Unknown RSSoft decoding options" << std::endl;
    }

    rssoft_decision_box.print_stats(options.source_prns, options.prns, std::cout);
}
#endif

#ifdef _CCSOFT
//=================================================================================================
void Reception_FEC::run_ccsoft_decoding(Options& options, ccsoft::CC_ReliabilityMatrix& relmat)
{
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
        return;
    }
    
	//ccsoft_engine.reset();
    
    CCSoft_DecisionBox ccsoft_decision_box(ccsoft_engine, options._source_codec);
    ccsoft_decision_box.run(relmat);
    ccsoft_decision_box.print_retrieved_message(std::cout);
    ccsoft_decision_box.print_stats(options.source_prns, options.prns, std::cout);
}
#endif
