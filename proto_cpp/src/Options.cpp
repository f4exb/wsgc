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
      
     Options

     Options parsing and holding
     
*/
#include "Options.h"
#include "WsgcUtils.h"
#include "FadingModel.h"
#include "FadingModelNone.h"
#include "FadingModelClarke.h"
#include "FadingModelWatterson.h"
#include "FIRCoefGenerator.h"
#include "FIRCoefGenerator_RCos.h"
#include "SourceCodec.h"
#include "SourceCodec_JT65.h"
#include <stdlib.h>
#include <getopt.h>
#include <cstring>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <boost/lexical_cast.hpp>


//=================================================================================================
// template to extract information from getopt more easily
template<typename TOpt, typename TField> bool extract_option(TField& field, char short_option)
{
    TOpt option_value;
    
    try
    {
        option_value = boost::lexical_cast<TOpt>(optarg);
        field = option_value;
        return true;
    }
    catch (boost::bad_lexical_cast &)
    {
        std::cout << "wrong argument for -" << short_option << ": " << optarg << " leave default (" << field << ")";
        std::cout << std::endl;
        return false;
    }
}


//=================================================================================================
Options::Options(std::string& _binary_name, Options_Executable _options_executable) :
	options_executable(_options_executable),
	transmission_scheme(OptionTrans_WSGC),
	fec_scheme(OptionFEC_None),
    binary_path(_binary_name),
    f_sampling(4096.0),
    f_chip(1023.0),
    snr(0.0),
    make_noise(false),
    f_tx(100.0),
    code_shift(1020),
    noise_test(false),
    nb_prns_per_symbol(4),
    prn_shift(0),
    f_init_rx(100.1),
    nb_random_prns(4),
    tracking_phase_average_cycles(4),
    gc_nb_stages(10),
    nb_message_symbols(64),
    nb_service_symbols(3),
    nb_training_symbols(0),
    nb_samples_per_code(0),
    file_debugging(false),
    _indicator_int(0),
    _fading_model(0),
    modulation(Modulation::Modulation_BPSK),
    nb_pilot_prns(1),
    pilot_gain_db(0.0),
    pilot1(0),
    pilot2(0),
    df_steps(31),
    f_step_division(8),
    batch_size(0),
    noise_prn(0),
    use_cuda(false),
    analysis_window_size(4),
    simulate_training(false),
    simulate_demod(false),
    gpu_affinity(0),
    gpu_affinity_specified(false),
	_fir_coef_generator(0),
	mfsk_options(4096.0),
    decision_thresholds_specified(false),
    _source_codec(0),
    rs_logq(0),
    rs_k(0)
#ifdef _CCSOFT
    ,cc_algorithm_type(CCSoft_Engine_defs::Algorithm_Fano),
    cc_decoding_mode(CCSoft_Engine_defs::Decoding_normal),
    cc_edge_bias(0.0),
    cc_node_limit(0),
    cc_use_node_limit(false),
    cc_metric_limit(0.0),
    cc_use_metric_limit(false),
    cc_interleave(true),
    cc_fano_init_metric(0.0),
    cc_fano_delta_metric(1.0),
    cc_fano_tree_cache_size(0),
    cc_fano_delta_init_threshold(0.0),
    cc_nb_retries(1),
    cc_edge_bias_decrement(0.0)
#endif
{
    srand(time(0));

    std::vector<std::string> path_nodes;
    WsgcUtils::extract_string_vector(path_nodes, binary_path);
    binary_name = path_nodes.back();

    std::cout << binary_name << std::endl;
}


//=================================================================================================
Options::~Options()
{
    if (_fading_model != 0)
    {
        delete _fading_model;
    }
    if (_source_codec != 0)
    {
    	delete _source_codec;
    }
}


//=================================================================================================
bool Options::get_options(int argc, char *argv[])
{
    int c;
    bool status = true;
    bool help = false;
    unsigned int random_seed;
    bool has_random_seed = false;
    
    while (true)
    {
        static struct option long_options[] =
        {
            // these options set a flag
            {"noise-test", no_argument, &_indicator_int, 1},
            {"file-debugging", no_argument, &_indicator_int, 1},
            {"help", no_argument, &_indicator_int, 1},
            {"cuda", no_argument, &_indicator_int, 1},
            {"simulate-trn", no_argument, &_indicator_int, 1},
            {"simulate-demod", no_argument, &_indicator_int, 1},
            // these options do not set a flag
            {"f-sampling", required_argument, 0, 's'},
            {"f-chip", required_argument, 0, 'c'},
            {"snr", required_argument, 0, 'n'},
            {"f-tx", required_argument, 0, 't'},
            {"code-shift", required_argument, 0, 'C'},
            {"prns", required_argument, 0, 'p'},
            {"nb-prn-per-symbol", required_argument, 0, 'N'},
            {"prn-shift", required_argument, 0, 'I'},   
            {"f-rx", required_argument, 0, 'r'},      
            {"random-prns", required_argument, 0, 'R'},
            {"tracking-phase-average", required_argument, 0, 'T'},
            {"nb-stages", required_argument, 0, 'm'},
            {"nb-message-symbols", required_argument, 0, 'M'},
            {"nb-service-symbols", required_argument, 0, 'S'},
            {"nb-training-symbols", required_argument, 0, 'Y'},
            {"g1-poly", required_argument, 0, 'g'},
            {"g2-poly", required_argument, 0, 'G'},
            {"fading-model", required_argument, 0, 'f'},
            {"modulation-scheme", required_argument, 0, 'd'},
            {"random-seed", required_argument, 0, 'a'},
            {"pilot-prns", required_argument, 0, 'P'},
            {"pilot-gain-db", required_argument, 0, 'A'},
            {"df-steps", required_argument, 0, 'F'},
            {"df-sub-steps", required_argument, 0, 'U'},
            {"batch-size", required_argument, 0, 'B'},
            {"analysis-window-size", required_argument, 0, 'z'},
            {"samples-output-file", required_argument, 0, 'o'},
            {"fir-filter-model", required_argument, 0, 'L'},
            {"gpu-affinity", required_argument, 0, 'y'},
            {"decision-thresholds", required_argument, 0, 'H'},
            {"source-coding", required_argument, 0, 'j'},
            {"fec-option", required_argument, 0, 'O'}
        };
        
        int option_index = 0;
        
        c = getopt_long (argc, argv, "s:c:n:t:C:p:N:I:r:R:T:m:M:S:g:G:f:d:a:P:A:F:B:z:U:o:L:Y:y:H:j:O:", long_options, &option_index);
        
        if (c == -1) // end of options
        {
            break;
        }
            
        switch(c)
        {
            case 0: // set flag
                if (strcmp("noise-test", long_options[option_index].name) == 0)
                {
                    noise_test = true;
                }
                else if (strcmp("file-debugging", long_options[option_index].name) == 0)
                {
                    file_debugging = true;
                }
                else if (strcmp("help", long_options[option_index].name) == 0)
                {
                    std::cout << "Showing help..." << std::endl;
                    help = true;
                }
                else if (strcmp("cuda", long_options[option_index].name) == 0)
                {
                    use_cuda = true;
                }
                else if (strcmp("simulate-trn", long_options[option_index].name) == 0)
                {
                    simulate_training = true;
                }
                else if (strcmp("simulate-demod", long_options[option_index].name) == 0)
                {
                    simulate_demod = true;
                }
                _indicator_int = 0;
                break;
            case 's':
                status = extract_option<double, wsgc_float>(f_sampling, 's');
                break;
            case 'c':
                status = extract_option<double, wsgc_float>(f_chip, 'c');
                break;
            case 'n':
                make_noise = true;
                status = extract_option<double, wsgc_float>(snr, 'n');
                break;
            case 't':
                status = extract_option<double, wsgc_float>(f_tx, 't');
                break;
            case 'C':
                status = extract_option<int, unsigned int>(code_shift, 'C');
                break;
            case 'p':
                status = extract_vector<unsigned int>(prns, std::string(optarg));
                break;
            case 'N':
                status = extract_option<int, unsigned int>(nb_prns_per_symbol, 'N');
                break;
            case 'I':
                status = extract_option<int, unsigned int>(prn_shift, 'I');
                break;
            case 'r':
                status = extract_option<double, wsgc_float>(f_init_rx, 'r');
                break;
            case 'R':
                status = extract_option<int, unsigned int>(nb_random_prns, 'R');
                break;
            case 'T':
                status = extract_option<int, unsigned int>(tracking_phase_average_cycles, 'T');
                break;
            case 'm':
                status = extract_option<int, unsigned int>(gc_nb_stages, 'm');
                break;
            case 'M':
                status = extract_option<int, unsigned int>(nb_message_symbols, 'M');
                break;
            case 'S':
                status = extract_option<int, unsigned int>(nb_service_symbols, 'S');
                break;
            case 'Y':
                status = extract_option<int, unsigned int>(nb_training_symbols, 'Y');
                break;
            case 'g':
                status = extract_vector<unsigned int>(g1_poly_powers, std::string(optarg));
                break;
            case 'G':
                status = extract_vector<unsigned int>(g2_poly_powers, std::string(optarg));
                break;
            case 'f':
                status = parse_fading_model_data(std::string(optarg));
                break;
            case 'd':
                status = parse_modulation_data(std::string(optarg));
                break;
            case 'a':
                status = extract_option<int, unsigned int>(random_seed, 'a');
                has_random_seed = true;
                break;
            case 'P':
                status = extract_option<int, unsigned int>(nb_pilot_prns, 'P');
                status = parse_pilot_prns_data(std::string(optarg));
                break;
            case 'A':
                status = extract_option<double, wsgc_float>(pilot_gain_db, 'A');
                break;
            case 'F':
                status = extract_option<int, unsigned int>(df_steps, 'F');
                break;
            case 'B':
                status = extract_option<int, unsigned int>(batch_size, 'B');
                break;
            case 'z':
                status = extract_option<int, unsigned int>(analysis_window_size, 'z');
                break;
            case 'U':
                status = extract_option<int, unsigned int>(f_step_division, 'U');
                break;
            case 'o':
            	samples_output_file = std::string(optarg);
            	status = true;
                break;
            case 'L':
                status = parse_fir_filter_model_data(std::string(optarg));
                break;
            case 'y':
                status = extract_option<int, unsigned int>(gpu_affinity, 'y');
                gpu_affinity_specified = true;
                break;
            case 'H':
                status = decision_thresholds.parse_options(std::string(optarg));
                decision_thresholds_specified = true;
                break;
            case 'j':
                status = parse_source_coding_data(std::string(optarg));
                break;
            case 'O':
                status = parse_fec_option(std::string(optarg));
                break;
            case '?':
                std::ostringstream os;
                get_help(os);
                std::cout << os.str() << std::endl;
                status = false;
                break;
        }
    }
    
    if (help)
    {
        std::ostringstream os;
        get_help(os);
        std::cout << os.str() << std::endl;       
        return false;
    }
    else
    {
        if (status)
        {
            // post process options           
            if (f_sampling < f_chip)
            {
                std::cout << "Sampling frequency must be higher than chip frequency (-s and -c options)" << std::endl;
                return false;
            }
            
            if (_fading_model == 0)
            {
                _fading_model = new FadingModelNone(f_sampling);
            }

            _fading_model->set_verbose(true);
            
            if (has_random_seed)
            {
                _fading_model->set_random_seed(random_seed);
            }
            
            if (nb_prns_per_symbol < 2)
            {
                std::cout << "Need at least 2 PRN codes per symbol (-N option)" << std::endl;
                return false;
            }
            
            if (tracking_phase_average_cycles < 2)
            {
                std::cout << "Frequency tracking phase average must be done at least on two valid phase measurements (-T option)" << std::endl;
                return false;
            }
            
            if (gc_nb_stages < 5)
            {
                std::cout << "Need generator LFSR of at least 5 stages producing codes of at least 31 chips (-m option)" << std::endl;
                return false;
            }
            
            if ((nb_message_symbols < 1) || ((nb_message_symbols + nb_service_symbols) > (1<<gc_nb_stages)-1))
            {
                std::cout << "Number of message symbols must be between 1 and " << (1<<gc_nb_stages)-1-nb_service_symbols << " (-M and -S options)" << std::endl;
                return false;
            }
            
            if ((nb_pilot_prns < 0) || (nb_pilot_prns > 2))
            {
            	std::cout << "Number of pilot PRNs is between 0 and 2 inclusive" << std::endl;
            	return false;
            }

            if (df_steps == 0)
            {
            	df_steps = 1;
            }

            if (f_step_division == 0)
            {
            	f_step_division = 1;
            }

            if ((nb_service_symbols < 1 + nb_pilot_prns) || ((nb_message_symbols + nb_service_symbols) > (1<<gc_nb_stages)-1))
            {
                std::cout << "Number of service symbols must be between " << 1 + nb_pilot_prns << " and " << (1<<gc_nb_stages)-1-nb_message_symbols << " (-M and -S options)" << std::endl;
                return false;
            }
            else
            {
                noise_prn = nb_message_symbols;

                if (nb_pilot_prns == 1)
                {
                    pilot1 = nb_message_symbols + 1;
                    pilot2 = pilot1;
                }
                else
                {
                    pilot1 = nb_message_symbols + 1;
                    pilot2 = nb_message_symbols + 2;
                }
            }
            
            nb_samples_per_code = (unsigned int) ((f_sampling / f_chip) * ((1<<gc_nb_stages)-1));
            
            if ((code_shift < 0) || (code_shift > nb_samples_per_code))
            {
                std::cout << "PRN code shift in number of samples must be between 0 and " << nb_samples_per_code << " (-C option)" << std::endl;
                return false;
            }
            
            // Build PRN list

            if (simulate_training) // Build training sequence and extra training simulation checks
            {
            	if (modulation.isCodeDivisionCapable() && (nb_pilot_prns > 0))
            	{
            		nb_training_symbols = 0; // Force zero

					if (nb_pilot_prns < 2)
					{
						std::cout << "The number of pilot PRNs should be 2 for synchronization training sequence simulation" << std::endl;
						return false;
					}

					if (prns.size() > 0)
					{
						prns.clear();
					}

					if (prn_shift+nb_random_prns > nb_message_symbols)
					{
						std::cout << "Cannot create training sequence: exceeding number of message symbols" << std::endl;
						return false;
					}
					else if (nb_random_prns < 4)
					{
						std::cout << "Need at least 4 PRNs in the training sequence" << std::endl;
						return false;
					}
					else
					{
						for (unsigned int prni = prn_shift; prni < prn_shift+nb_random_prns; prni++)
						{
							prns.push_back(prni);
						}
					}
            	}
            	else if (modulation.demodulateBeforeCorrelate())
            	{
					if (prn_shift+nb_random_prns > nb_training_symbols)
					{
						std::cout << "Cannot create training sequence: exceeding number of training symbols" << std::endl;
						return false;
					}
					else if (nb_random_prns < 4)
					{
						std::cout << "Need at least 4 PRNs in the training sequence" << std::endl;
						return false;
					}
					else
					{
						if (prns.size() > 0)
						{
							prns.clear();
						}
						for (unsigned int prni = prn_shift+nb_message_symbols+nb_service_symbols; prni < prn_shift+nb_message_symbols+nb_service_symbols+nb_random_prns; prni++)
						{
							prns.push_back(prni);
						}
					}
            	}
            }
            else // Build PRN list for message simulation
            {
            	nb_training_symbols = 0; // Force zero

				if (prn_shift > nb_prns_per_symbol - 1)
				{
					std::cout << "Index of the PRN in symbol where the simulation starts must be between 0 and " << nb_prns_per_symbol << " (-I option)" << std::endl;
					return false;
				}

				// check given symbol PRNs
				if (prns.size() == 0)
				{
					if (nb_random_prns < 4)
					{
						std::cout << "Need at least 4 random symbols (-R option) when no symbols given (-p option)" << std::endl;
						return false;
					}

					for (int i=0; i < nb_random_prns; i++)
					{
						prns.push_back(rand() % nb_message_symbols);
					}
				}
				else
				{
					for (std::vector<unsigned int>::const_iterator it = prns.begin(); it != prns.end(); ++it)
					{
						if ((*it < 0) || (*it > nb_message_symbols + nb_service_symbols - 1))
						{
							std::cout << "Wrong PRN number " << *it << " must be between 0 and " << nb_message_symbols + nb_service_symbols - 1 << " (options -p -M and -S)" << std::endl;
						}
					}
				}
            }

            // control analysis window size vs number of PRNs in sequence
			if (analysis_window_size > prns.size())
			{
				analysis_window_size = prns.size();
			}

			// source coding (if any) adjustments
			if (_source_codec)
			{
				if (!adjust_parameters_for_source_coding())
				{
					return false;
				}
				if (!source_codec_create_message_prns())
				{
					return false;
				}
			}

			// last validation: generator polynomials
            if ((g1_poly_powers.size() > 0) && (g2_poly_powers.size() > 0)) // basic checks on given generator polynomials
            {
                std::sort(g1_poly_powers.begin(), g1_poly_powers.end(), std::greater<unsigned int>());
                std::sort(g2_poly_powers.begin(), g2_poly_powers.end(), std::greater<unsigned int>());
                
                if ((g1_poly_powers.front() > gc_nb_stages - 1) || (g1_poly_powers.back() < 1))
                {
                    std::cout << "Specified G1 polynomial powers must be between 0 and " << gc_nb_stages - 1 << " (options -g and -m)" << std::endl;
                    return false;
                }
                
                if ((g1_poly_powers.front() > gc_nb_stages - 1) || (g1_poly_powers.back() < 1))
                {
                    std::cout << "Specified G2 polynomial powers must be between 0 and " << gc_nb_stages - 1 << " (options -G and -m)" << std::endl;
                    return false;
                }
                
                return true;
            }
            else // create default generator polynomials for orders 5 to 11 with preferred pairs of m-sequences
            {
                switch(gc_nb_stages)
                {
                    case 5:
                        g1_poly_powers.push_back(3); // G1 = X^5 + X^3 + 1
                        g2_poly_powers.push_back(4); // G2 = X^5 + X^4 + X^3 + X^2 + 1
                        g2_poly_powers.push_back(3);
                        g2_poly_powers.push_back(2);
                        return true;
                    case 6:
                        g1_poly_powers.push_back(1); // G1 = X^6 + X^1 + 1
                        g2_poly_powers.push_back(5); // G2 = X^6 + X^5 + X^2 + X^1 + 1
                        g2_poly_powers.push_back(2);
                        g2_poly_powers.push_back(1);
                        return true;
                    case 7:
                        g1_poly_powers.push_back(3); // G1 = X^7 + X^3 + 1
                        g2_poly_powers.push_back(3); // G2 = X^7 + X^3 + X^2 + X^1 + 1
                        g2_poly_powers.push_back(2);
                        g2_poly_powers.push_back(1);
                        return true;
                    case 8:
                        g1_poly_powers.push_back(7); // G1 = X^8 + X^7 + X^6 + X^1 + 1
                        g1_poly_powers.push_back(6);
                        g1_poly_powers.push_back(1);
                        g2_poly_powers.push_back(7); // G2 = X^8 + X^7 + X^6 + X^5 + X^4 + X^2 + 1
                        g2_poly_powers.push_back(6);
                        g2_poly_powers.push_back(5);
                        g2_poly_powers.push_back(4);
                        g2_poly_powers.push_back(2);
                        return true;
                    case 9:
                        g1_poly_powers.push_back(4); // G1 = X^9 + X^4 + 1
                        g2_poly_powers.push_back(6); // G2 = X^9 + X^6 + X^4 + X^3 + 1
                        g2_poly_powers.push_back(4);
                        g2_poly_powers.push_back(3);
                        return true;
                    case 10:
                        g1_poly_powers.push_back(8); // was G1 = X^10 + X^3 + 1
                        g1_poly_powers.push_back(5); // G1 = X^10 + X^8 + X^5 + X^1 + 1
                        g1_poly_powers.push_back(1);
                        g2_poly_powers.push_back(9); // was G2 = X^10 + X^8 + X^3 + X^2 + 1
                        g2_poly_powers.push_back(7); // G2 = X^10 + X^9 + X^7 + X^6 + X^4 + X^1 + 1
                        g2_poly_powers.push_back(6);
                        g2_poly_powers.push_back(4);
                        g2_poly_powers.push_back(1);
                        return true;
                    case 11:
                        g1_poly_powers.push_back(2); // G1 = X^11 + X^2 + 1
                        g2_poly_powers.push_back(8); // G2 = X^11 + X^8 + X^5 + X^2 + 1
                        g2_poly_powers.push_back(5);
                        g2_poly_powers.push_back(2);
                        return true;
                    default:
                        std::cout << "No default generator polynomials for this number of stages in LFSR: " << gc_nb_stages << " (options -m -g and -G)" << std::endl;
                        return false;
                }
            }
        }
        else
        {
            return false;
        }
    }
}

//=================================================================================================
void Options::get_help(std::ostringstream& os)
{
    struct help_struct
    {
        const char *short_option;
        const char *long_option;
        const char *help_text;
        const char *type;
        const char *default_value;
        const char *binary_name;
    };
    
    static struct help_struct help_lines[] =
    {       
        {"", "--help", "Print this help and exit", "flag", "false", "any"},
        {"", "--noise-test", "Trigger noise test", "flag", "false", ""},    
        {"", "--file-debugging", "Trigger debugging on files", "flag", "false", ""},    
        {"", "--cuda", "Trigger CUDA implementation", "flag", "false", "wsgc_test"},
        {"", "--simulate-trn", "Trigger synchronization training sequence simulation", "flag", "false", "any"},
        {"", "--simulate-demod", "Trigger external symbol synchronization simulation", "flag", "false", "wsgc_generator"},
        {"-s", "--f-sampling", "Sampling frequency", "float", "4096.0", "any"},    
        {"-c", "--f-chip", "Chip frequency", "float", "1023.0", "any"},    
        {"-C", "--code-shift", "PRN code shift in number of samples", "int", "1020", "any"},    
        {"-N", "--nb-prn-per-symbol", "Number of PRN code sequences per symbol", "int", "4", "any"},
        {"-I", "--prn-shift", "Start simulation at PRN shift in symbol", "int", "0", "any"},
        {"-t", "--f-tx", "Transmission frequency", "float", "100.0", "any"},
        {"-r", "--f-rx", "Initial receiving frequency", "float", "100.1", "wsgc_test"},
        {"-n", "--snr", "Add Gaussian noise: signal to noise ratio in signal bandwidth in dB", "float", "-24.0", "any"},    
        {"-t", "--tracking-phase-average", "Number of phase averaging cycles for frequency tracking", "int", "4", ""},
        {"-m", "--nb-stages", "Number of stages of the PRN generator LFSRs. This gives 2**N-1 chips", "int", "10", "any"},
        {"-g", "--g1-poly", "Generator 1 comma separated list of polynomial powers except N and 0", "string", "\"3\"", "any"},
        {"-G", "--g2-poly", "Generator 2 comma separated list of polynomial powers except N and 0", "string", "\"8,3,2\"", "any"},
        {"-M", "--nb-message-symbols", "Number of message symbols", "int", "64", "any"},
        {"-S", "--nb-service-symbols", "Number of service symbols", "int", "3", "any"},
        {"-Y", "--nb-training-symbols", "Number of training symbols, used with --simulate-trn and unpiloted schemes", "int", "0", "any"},
        {"-P", "--pilot-prns", "Number of pilot PRNs (see documentation)", "int", "1", "any"},
        {"-A", "--pilot-gain-db", "Gain of the pilot PRN(s) vs Message PRNs", "float", "0.0", "any"},
        {"-F", "--df-steps", "Number of frequency steps to explore (step sized by FFT size and sampling frequency)", "int", "31", "wsgc_test"},
        {"-U", "--df-sub-steps", "Number of frequency sub-steps of the above frequency step to explore", "int", "8", "wsgc_test"},
        {"-B", "--batch-size", "Batch size in number of PRNs for pilot processing", "int", "3", "any"},
        {"-p", "--prns", "Symbol PRN numbers in a comma separated list (overrides -R)", "string", "null (uses random option)", "any"},
        {"-R", "--random-prns", "Number of random symbols to generate", "int", "4", "any"},
        {"-f", "--fading-model", "Fading model data (see short help below)", "string", "null (no fading)", "any"},
        {"-d", "--modulation-scheme", "Modulation scheme (see short help below)", "string", "BPSK", "any"},
        {"-z", "--analysis-window-size", "Pilot analysis window size in number of symbols", "int", "4", "wsgc_test"},
        {"-o", "--samples-output-file", "Output file for generated samples", "string", "", "wsgc_generator"},
        {"-L", "--fir-filter-model", "Lowpass FIR filter model (see short help below)", "string", "", "any"},
        {"-y", "--gpu-affinity", "Force CUDA execution on specified GPU ID", "int", "(none)", "wsgc_test"},
        {"-H", "--decision-thresholds", "Specify decision box thresholds (see decision thresholds below)", "string", "(none)", "wsgc_test"},
        {"-f", "--source-codec", "Source codec data (see short help below)", "string", "null (no source codec)", "any"},
        {"-O", "--fec-option", "Forward Error Correction using external libraries (see short help below)", "string", "null (no FEC)", "any"},
        {0,0,0,0,0,0}
    };
    
    // get maximum width of each column
    unsigned int i = 0;
    unsigned int short_option_len = 0;
    unsigned int long_option_len = 0;
    unsigned int help_text_len = 0;
    unsigned int type_len = 0;
    unsigned int default_value_len = 0;
    
    do
    {
        if (strlen(help_lines[i].short_option) > short_option_len)
            short_option_len = strlen(help_lines[i].short_option);
        if (strlen(help_lines[i].long_option) > long_option_len)
            long_option_len = strlen(help_lines[i].long_option);
        if (strlen(help_lines[i].help_text) > help_text_len)
            help_text_len = strlen(help_lines[i].help_text);
        if (strlen(help_lines[i].type) > type_len)
            type_len = strlen(help_lines[i].type);
        if (strlen(help_lines[i].default_value) > default_value_len)
            default_value_len = strlen(help_lines[i].default_value);
        i++;
    } while (help_lines[i].short_option);
    
    // display option help
    os << "Valid options (short, long, comment, type, default)" << std::endl;
    i = 0;    

    do
    {
        char anystr[] = "any";
        
        if ((strcmp(help_lines[i].binary_name, binary_name.c_str()) == 0) || (strcmp(help_lines[i].binary_name, anystr) == 0))
        {
            os << std::setw(short_option_len) << help_lines[i].short_option << " " 
               << std::setw(long_option_len) << std::left << std::setfill(' ') << help_lines[i].long_option << " "
               << std::setw(help_text_len) << std::left << std::setfill(' ') << help_lines[i].help_text << " "
               << std::setw(type_len) << std::left << std::setfill(' ') << help_lines[i].type << " "
               << std::setw(default_value_len) << std::left << std::setfill(' ') << help_lines[i].default_value << std::endl;
        }
        else if (help_lines[i].binary_name[0] == '\0')
        {
            os << std::setw(short_option_len) << help_lines[i].short_option << " " 
               << std::setw(long_option_len) << std::left << std::setfill(' ') << help_lines[i].long_option << " "
               << std::setw(help_text_len) << std::left << std::setfill(' ') << help_lines[i].help_text << " "
               << "-- deprecated -- " << std::endl;
        }
        i++;
    } while (help_lines[i].short_option);
    
    os << std::endl;
    os << "Modulation schemes:" << std::endl;
    os << " - BPSK : BPSK" << std::endl;
    os << " - DBPSK: DBPSK (incomplete and experimental for now, work is in progress...)" << std::endl;
    os << " - OOK  : OOK"  << std::endl;
    os << " - CW   : CW (no message: only to test fading models with wsgc_generator)" << std::endl;
    os << " - MFSK : MFSK (no correlation, used for benchmarking, see MFSK options)" << std::endl;
    os << std::endl;
    os << "Fading models:" << std::endl;
    os << " - Clarke:    -f \"C:<nb paths>,<frequency spread (Hz)>\"" << std::endl;
    os << " - Watterson: -f \"W:<path params>/<path params>/...\"" << std::endl;
    os << "     path params: <delay (s)>,<amplitude>,<f Doppler (Hz)>,<f offset (Hz)>" << std::endl;
    os << std::endl;
    os << "Training sequence simulation:" << std::endl;
    os << "Some options are re-used in a different way:" << std::endl;
    os << " - Random PRNs option becomes the number of PRNs in the training sequence (>= 16 recommended)" << std::endl;
    os << " - PRN shift in symbol start becomes the PRN number at which the training sequence starts" << std::endl;
    os << " => These two options should be arranged so as not to yield a PRN number higher than the number of message PRNs" << std::endl;
    os << "The number if pilot PRNs should be 2" << std::endl;
    os << std::endl;
    os << "Lowpass FIR filter models:" << std::endl;
    os << " - Raised cosine: \"RCOS:<nb_taps>,<rolloff>,<cutoff frequency>\"" << std::endl;
    os << "   Rolloff factor must be between 0.0 and 1.0. Invalid values yield 1.0" << std::endl;
    os << " Cutoff frequency must be lower or equal to half the sampling frequency. Invalid values yield fs/2" << std::endl;
    os << std::endl;
    os << "Source codecs:" << std::endl;
    os << " - JT65: classical. Uses RS(63,12): -j \"JT65::<source text message>\"" << std::endl;
    os << " - JT257: packing the 72 bytes in 9 8-byte symbols. Otherwise classical. Uses RS(255,9): -j \"JT257::<source text message>\"" << std::endl;
    os << " - JTCC: packing the 72 bytes in 72 1-byte symbols. Otherwise classical. Uses any 1/n rate CC code: -j \"JTCC::<source text message>\"" << std::endl;
    os << std::endl;
    os << "Forward Error Correction using external libraries -O <FEC scheme>/<FEC options...>:" << std::endl;
    os << std::endl;
#ifdef _RSSOFT
    os << "  Reed Solomon encoding and soft-decision decoding using RSSoft library:" << std::endl;
    os << "  -O RS/q,k,M,r,i:<multiplicity increment strategy>,<decoding mode>(,<regular expression>)" << std::endl;
    os << "    For a RS(n,k) code with n = 2^q-1 and 1<k<n," << std::endl;
    os << "    M is the initial multiplicity matrix global multiplicity," << std::endl;
    os << "    r is the maximum number of retries," << std::endl;
    os << "    i is the base increment for next retry multiplicity calculation:" << std::endl;
    os << "    Multiplicity increment strategy:" << std::endl;
    os << "     - arith: arithmetic increment for multiplicity" << std::endl;
    os << "     - iarith: arithmetic increment for the additive multiplicity" << std::endl;
    os << "     - geom: geometric increment for multiplicity" << std::endl;
    os << "     - igeom: geometric increment for the multiplicative multiplicity" << std::endl;
    os << "    Decoding modes:" << std::endl;
    os << "     - full: gives all results hence executes r retries" << std::endl;
    os << "     - best: only return the result with best probability or none" << std::endl;
    os << "     - first: returns the first result" << std::endl;
    os << "     - regex: retries until it finds a resulting textual message matching the given regular expression" << std::endl;
    os << "     - match: retries until it finds a resulting textual message matching exactly the given string" << std::endl;
    os << "     - binmatch: retries until it finds a resulting message matching exactly the sent message (test context)" << std::endl;
    os << "     - relthr: retries until it finds a candidate message with a reliability figure above the given threshold" << std::endl;
    os << std::endl;
#endif
#ifdef _CCSOFT
    os << "  Convolutional Code encoding and soft-decision decoding using CCSoft library:" << std::endl;
    os << "  -O CC/<constraint lengths>/<generator polynomials>/<decoding algorithm>" << std::endl;
    os << "    Constraint lengths: comma separated list of constraint (actually register) lengths" << std::endl;
    os << "      - there is one constraint per input symbol bit (k)" << std::endl;
    os << "    Generator polynomials: colon separated list of lists of binary representation of generator polynomials" << std::endl;
    os << "      - there is one list per input symbol bit (k)" << std::endl;
    os << "      - each list has one polynomial per output symbol bit (n)" << std::endl;
    os << "    Decoding algorithm: <algorithm code>:<algorithm arguments comma separated>" << std::endl;
    os << "      - Fano algorithm (code = fano):" << std::endl;
    os << "         Arguments: <bias>,<init threshold>,<delta threshold>,<cache size>,<delta init threshold>,<node limit>(,<threshold limit>)" << std::endl;
    os << "      - Stack algorithm (code = stack):" << std::endl;
    os << "         Arguments: <bias>,<node limit>(,<metric limit>)" << std::endl;
    os << std::endl;
#endif
    mfsk_options.get_help(os);
    os << std::endl;
    decision_thresholds.get_help(os);
    os << std::endl;

    if (binary_name == "wsgc_test")
    {
        os << "Example: wsgc_test -s 4096 -c 1023 -C 1020 -t 100.0 -r 100.1 -n -24 -N 4 -I 0 -R 6" << std::endl;
    }
    else if (binary_name == "wsgc_generator")
    {
        os << "Example: wsgc_generator -s 4096 -c 1023 -C 1020 -t 100.0 -n -24 -N 4 -I 0 -R 6 -o \"my_samples_iq_file.raw\"" << std::endl;
    }
}


void Options::print_options(std::ostringstream& os)
{
	if (modulation.getScheme() == Modulation::Modulation_MFSK)
	{
		print_mfsk_options(os);
	}
	else
	{
		print_standard_options(os);
	}
}

        
//=================================================================================================
void Options::print_standard_options(std::ostringstream& os)
{
    // additional computations
    unsigned int code_length = (1<<gc_nb_stages)-1;
    wsgc_float code_period = code_length / f_chip;
    wsgc_float symbol_period = nb_prns_per_symbol * code_period;
    wsgc_float message_time = prns.size() * symbol_period;
    
    os << "Using options:" << std::endl;
    os << "------------- " << std::endl;
    os << std::endl;
    os << std::setiosflags(std::ios_base::fixed);
    os << "Sampling frequency ........: " << std::setw(8) << std::setprecision(1) << std::right << f_sampling << std::endl;
    os << "Chip frequency ............: " << std::setw(8) << std::setprecision(1) << std::right << f_chip << std::endl;
    os << "Samples per chip ..........: " << std::setw(10) << std::setprecision(3) << std::right << f_sampling/f_chip << std::endl;
    os << "Code length ...............: " << std::setw(6) << std::right << code_length << std::endl;
    os << "Code period ...............: " << std::setw(9) << std::setprecision(2) << std::right << code_period << std::endl;
    os << "Samples/code = FFT size ...: " << std::setw(6) << std::right << nb_samples_per_code << std::endl;
    os << "Code shift ................: " << std::setw(6) << std::right << code_shift << std::endl;
    os << "Nb of generated symbols ...: " << std::setw(6) << std::right << prns.size() << std::endl;

    if (simulate_training)
    {
        os << "PRNs per pilot averaging ..: " << std::setw(6) << std::right << nb_prns_per_symbol << std::endl;
        os << "Pilot averaging period ....: " << std::setw(9) << std::setprecision(2) << std::right << symbol_period << std::endl;
    	os << "Start PRN number ..........: " << std::setw(6) << std::right << prn_shift << std::endl;
    }
    else
    {
        os << "PRNs per symbol ...........: " << std::setw(6) << std::right << nb_prns_per_symbol << std::endl;
        os << "Symbol period .............: " << std::setw(9) << std::setprecision(2) << std::right << symbol_period << std::endl;
    	os << "Start PRN shift in symbol .: " << std::setw(6) << std::right << prn_shift << std::endl;
    }

	os << "Tx frequency ..............: " << std::setw(9) << std::setprecision(2) << std::right << f_tx << std::endl;
	os << "Initial Rx frequency ......: " << std::setw(9) << std::setprecision(2) << std::right << f_init_rx << std::endl;

    os << "SNR(dB) ...................: ";
    
    if (make_noise)
    {
        os << std::setw(8) << std::setprecision(1) << std::right << snr << std::endl; 
    }
    else
    {
        os << "No noise" << std::endl;
    }
    
    if (_fading_model != 0)
    {
        os << "Fading Model ..............: "; print_fading_model_data(os); os << std::endl;
    }
    
    if (_source_codec != 0)
    {
    	os << "Source codec ..............: "; print_source_codec_data(os); os << std::endl;
    }

#ifdef _RSSOFT
    if (rs_k != 0)
    {
        os << std::endl;
    	os << "Reed Solomon (RSSoft lib)" << std::endl; 
        print_reed_solomon_data(os); 
        os << std::endl;
    }
#endif

#ifdef _CCSOFT
    if (cc_k_constraints.size() > 0)
    {
        os << std::endl;
        os << "Conv. Code (CCSoft lib)" << std::endl; 
        print_convolutional_code_data(os); 
        os << std::endl;
    }
#endif

    os << "Modulation ................: "; modulation.print_modulation_data(os); os << std::endl;
    os << "Transmission scheme .......: "; print_transmission_scheme(os); os << std::endl;
    //os << "Nb phase averaging cycles .: " << std::setw(6) << std::right << tracking_phase_average_cycles << std::endl;
    os << "Nb message symbols ........: " << std::setw(6) << std::right << nb_message_symbols << std::endl;
    os << "Nb service symbols ........: " << std::setw(6) << std::right << nb_service_symbols << std::endl;
    os << "Nb training symbols .......: " << std::setw(6) << std::right << nb_training_symbols << std::endl;
    os << "Noise PRN .................: " << std::setw(6) << std::right << noise_prn << " (" << noise_prn - nb_message_symbols << ")" << std::endl;

    if (modulation.isCodeDivisionCapable() && nb_pilot_prns > 0)
    {
    	wsgc_float freq_step_size = f_sampling/nb_samples_per_code;
    	int f_step_low  = -(df_steps/2);
    	int f_step_high = df_steps-(df_steps/2)-1;

        os << "Pilot PRN gain ............: " << std::setw(8) << std::setprecision(1) << std::right << pilot_gain_db << " dB" << std::endl;
        os << "Nb frequency steps ........: " << std::setw(6) << std::right << df_steps << std::endl;
        os << "Nb frequency sub-steps ....: " << std::setw(6) << std::right << f_step_division << std::endl;
        os << "Frequency range ...........: " << std::setprecision(1) << "[" << f_step_low*freq_step_size << ":" << f_step_high*freq_step_size << "]" << std::endl;
        os << "Minor freq step size ......: " << std::setw(9) << std::setprecision(3) << std::right << (freq_step_size/f_step_division) << std::endl;

        if (pilot1 == pilot2)
        {
            os << "Pilot PRN .................: " << std::setw(6) << std::right << pilot1 << " (" << pilot1 - nb_message_symbols << ")" << std::endl;
        }
        else
        {
            os << "Pilot PRN 1................: " << std::setw(6) << std::right << pilot1 << " (" << pilot1 - nb_message_symbols << ")" << std::endl;
            os << "Pilot PRN 2................: " << std::setw(6) << std::right << pilot2 << " (" << pilot2 - nb_message_symbols << ")" << std::endl;
        }
    }

    if (modulation.isCodeDivisionCapable() && (nb_pilot_prns > 0))
    {
        os << "Batch size ................: " << std::setw(6) << std::right << batch_size << " PRNs" << std::endl;
    }
    else
    {
        os << "Batch size ................: " << std::setw(6) << std::right << batch_size << " Symbols" << std::endl;
    }

    if (simulate_training)
    {
    	if (modulation.isCodeDivisionCapable() && (nb_pilot_prns > 0))
    	{
    		os << "Analysis window size ......: " << std::setw(6) << std::right << analysis_window_size*nb_prns_per_symbol << " PRNs" << std::endl;
    	}
    	else
    	{
    		os << "Analysis window size ......: " << std::setw(6) << std::right << analysis_window_size << " PRNs" << std::endl;
    	}
    }
    else
    {
    	os << "Analysis window size ......: " << std::setw(6) << std::right << analysis_window_size << " Symbols" << std::endl;
    }

	os << "Analysis window time ......: " << std::setw(9) << std::setprecision(2) << std::right << symbol_period << std::endl;
	os << "Message time ..............: " << std::setw(9) << std::setprecision(2) << std::right << message_time << std::endl;

	if (binary_name == "wsgc_generator")
	{
		os << "Samples output file .......: " << samples_output_file << std::endl;
		os << "Demodulation ..............: " << (simulate_demod ? "Yes" : "No") << std::endl;
	}

    os << std::endl;
    
    if (_fir_coef_generator != 0)
    {
    	os << "Lowpass FIR filter model:" << std::endl;
    	_fir_coef_generator->dump(os);
    	os << std::endl;
    }

    os << "Processing options:" << std::endl;
    
    if (use_cuda)
    {
        os << " - Using CUDA implementation" << std::endl;
        
        if (gpu_affinity_specified)
        {
            os << " - Running on GPU #" << gpu_affinity << " requested" << std::endl;
        }
    }
    else
    {
        std::cout << " - Using Host implementation" << std::endl;
    }
    
    if (simulate_training)
    {
        os << " - " << "Simulate synchronization training sequence" << std::endl;
    }
    else
    {
        os << " - " << "Simulate message correlation" << std::endl;
    }
    os << std::endl;

    os << "Generator polynomials:" << std::endl;
    os << "G1 = ";
    WsgcUtils::print_polynomial(gc_nb_stages, g1_poly_powers, os); os << std::endl;
    os << "G2 = ";
    WsgcUtils::print_polynomial(gc_nb_stages, g2_poly_powers, os); os << std::endl;
    
    if (noise_test)
    {
        os << " - Performing noise measurements" << std::endl;
    }
    
    os << std::endl;
    os << "Processing " << prns.size() << " symbols: ";
    print_vector<unsigned int, unsigned int>(prns, 4, os); os << std::endl;
}

//=================================================================================================
void Options::print_transmission_scheme(std::ostringstream& os)
{
    switch (transmission_scheme)
    {
    case OptionTrans_None:
        os << "None";
        break;
    case OptionTrans_WSGC:
        os << "WSGC";
        break;
    case OptionTrans_WSGCE:
        os << "WSGCE";
        break;
    case OptionTrans_WSGCD:
        os << "WSGCD";
        break;
    case OptionTrans_WSGCO:
        os << "WSGCO";
        break;
    case OptionTrans_MFSK:
        os << "MFSK";
        break;
    default:
        os << "Unknown";
        break;
    }
}

//=================================================================================================
void Options::print_mfsk_options(std::ostringstream& os)
{
    os << "Using options:" << std::endl;
    os << "------------- " << std::endl;
    os << std::endl;
    os << "Modulation ................: "; modulation.print_modulation_data(os); os << std::endl;
    mfsk_options.print_options(os);
    os << "Nb message symbols ........: " << std::setw(6) << std::right << nb_message_symbols << std::endl;
    os << "Nb service symbols ........: " << std::setw(6) << std::right << nb_service_symbols << std::endl;
    os << "Nb training symbols .......: " << std::setw(6) << std::right << nb_training_symbols << std::endl;
	os << "Max symbol frequency.......: " << std::setw(9) << std::setprecision(2) << std::right << (nb_message_symbols+nb_service_symbols)*mfsk_options._symbol_bandwidth + mfsk_options._zero_frequency << std::endl;
    os << "Nb of generated symbols ...: " << std::setw(6) << std::right << prns.size() << std::endl;
	os << "Message time ..............: " << std::setw(9) << std::setprecision(2) << std::right <<  prns.size() * mfsk_options._symbol_time << std::endl;
	os << "Tx shift frequency ........: " << std::setw(9) << std::setprecision(2) << std::right << f_tx << std::endl;
    os << "SNR(dB) ...................: ";

    if (make_noise)
    {
        os << std::setw(8) << std::setprecision(1) << std::right << snr << std::endl;
    }
    else
    {
        os << "No noise" << std::endl;
    }

    if (_fading_model != 0)
    {
        os << "Fading Model ..............: "; print_fading_model_data(os); os << std::endl;
    }

    if (_source_codec != 0)
    {
    	os << "Source codec ..............: "; print_source_codec_data(os); os << std::endl;
    }

#ifdef _RSSOFT
    if (rs_k != 0)
    {
        os << std::endl;
    	os << "Reed Solomon (RSSoft lib)" << std::endl; 
        print_reed_solomon_data(os); 
        os << std::endl;
    }
#endif

#ifdef _CCSOFT
    if (cc_k_constraints.size() > 0)
    {
        os << std::endl;
        os << "Conv. Code (CCSoft lib)" << std::endl; 
        print_convolutional_code_data(os); 
        os << std::endl;
    }
#endif

	if (binary_name == "wsgc_generator")
	{
		os << "Samples output file .......: " << samples_output_file << std::endl;
		os << "Demodulation ..............: " << (simulate_demod ? "Yes" : "No") << std::endl;
	}

	os << "Processing options:" << std::endl;

    if (use_cuda)
    {
        os << " - Using CUDA implementation" << std::endl;

        if (gpu_affinity_specified)
        {
            os << " - Running on GPU #" << gpu_affinity << " requested" << std::endl;
        }
    }
    else
    {
        std::cout << " - Using Host implementation" << std::endl;
    }

    os << std::endl;
    os << "Processing " << prns.size() << " symbols: ";
    print_vector<unsigned int, unsigned int>(prns, 4, os); os << std::endl;
}


//=================================================================================================
void Options::print_fading_model_data(std::ostringstream& os)
{
    if (_fading_model)
    {
        _fading_model->print_fading_data(os);
    }
    else
    {
        os << "None";
    }
}


//=================================================================================================
void Options::print_source_codec_data(std::ostringstream& os)
{
    if (_source_codec)
    {
        _source_codec->print_source_codec_data(os);
        os << "; Message: \"" << _source_message_str << "\"";
    }
    else
    {
        os << "None";
    }
}

#ifdef _RSSOFT
//=================================================================================================
void Options::print_reed_solomon_data(std::ostringstream& os)
{
    if (rs_k)
    {
    	if (options_executable == Options_wsgc_generator)
    	{
    		os << "RS(" << (1<<rs_logq)-1 << "," << rs_k << ")";
    	}
    	else if (options_executable == Options_wsgc_test)
    	{
			os << "RS(" << (1<<rs_logq)-1 << "," << rs_k << "), M = " << rs_init_M << ", r = " << rs_r << ", i = " << rs_inc;

			os << ", increment strategy = ";
			switch (rs_inc_strategy)
			{
			case RSSoft_Engine_defs::MMatrix_retry_arithmetic:
				os << "arithmetic";
				break;
			case RSSoft_Engine_defs::MMatrix_retry_arithmetic_increment:
				os << "arithmetic increment";
				break;
			case RSSoft_Engine_defs::MMatrix_retry_geometric:
				os << "geometric";
				break;
			case RSSoft_Engine_defs::MMatrix_retry_geometric_increment:
				os << "geometric increment";
				break;
			default:
				os << "none";
				break;
			}

			os << ", decoding mode = ";
			switch (rs_decoding_mode)
			{
			case RSSoft_Engine_defs::RSSoft_decoding_all:
				os << "all";
				break;
			case RSSoft_Engine_defs::RSSoft_decoding_full:
				os << "full";
				break;
			case RSSoft_Engine_defs::RSSoft_decoding_best:
				os << "best";
				break;
			case RSSoft_Engine_defs::RSSoft_decoding_first:
				os << "first";
				break;
			case RSSoft_Engine_defs::RSSoft_decoding_regex:
				os << "regex: \"" << rs_decoding_match_str << "\"";
				break;
			case RSSoft_Engine_defs::RSSoft_decoding_match:
				os << "match: \"" << rs_decoding_match_str << "\"";
				break;
			case RSSoft_Engine_defs::RSSoft_decoding_binmatch:
				os << "binmatch";
				break;
			case RSSoft_Engine_defs::RSSoft_decoding_relthr:
				os << "relthr: " << rs_reliability_threshold << " dB/Symbol";
				break;
			default:
				os << "none";
				break;
			}
        }
    }
    else
    {
        os << "None";
    }
}
#endif

#ifdef _CCSOFT
//=================================================================================================
void Options::print_convolutional_code_data(std::ostringstream& os)
{
    os << "Edge metric bias ..................: " << std::setw(5) << std::setprecision(2) << cc_edge_bias << std::endl;
    os << "Symbol interleaving ...............: " << (cc_interleave ? "Yes" : "No") << std::endl;

    if (cc_use_node_limit)
    {
        os << "Node number limit .................: " << cc_node_limit << std::endl;
    }

    if (cc_use_metric_limit)
    {
        os << "Path metric number limit ..........: " << cc_metric_limit << std::endl;
    }

    switch (cc_algorithm_type)
    {
    case CCSoft_Engine_defs::Algorithm_Stack:
        os << "Using stack algorithm" << std::endl;
        break;
    case CCSoft_Engine_defs::Algorithm_Fano:
        os << "Using Fano algorithm" << std::endl;
        os << "  Initial metric threshold ........: " << cc_fano_init_metric << std::endl;
        os << "  Metric threshold delta ..........: " << cc_fano_delta_metric << std::endl;
        os << "  Nb nodes cached .................: " << cc_fano_tree_cache_size << std::endl;
        os << "  Initial metric threshold delta ..: " << cc_fano_delta_init_threshold << std::endl;
        break;
    default:
        os << "Unknown algorithm" << std::endl;
        break;
    }

    switch (cc_decoding_mode)
    {
    case CCSoft_Engine_defs::Decoding_normal:
        os << "Using normal decoding" << std::endl;
        break;
    case CCSoft_Engine_defs::Decoding_regex:
        os << "Using regexp decoding" << std::endl;
        os << "  Nb retries ......................: " << cc_nb_retries << std::endl;
        os << "  Edge bias decrement .............: " << cc_edge_bias_decrement << std::endl;
        os << "  Regexp ..........................: " << cc_match_str << std::endl;
        break;
    case CCSoft_Engine_defs::Decoding_match_str:
        os << "Using match exact string decoding" << std::endl;
        os << "  Nb retries ......................: " << cc_nb_retries << std::endl;
        os << "  Edge bias decrement .............: " << cc_edge_bias_decrement << std::endl;
        os << "  Matching string .................: " << cc_match_str << std::endl;
        break;
    case CCSoft_Engine_defs::Decoding_match_msg:
        os << "Using match source message decoding" << std::endl;
        os << "  Nb retries ......................: " << cc_nb_retries << std::endl;
        os << "  Edge bias decrement .............: " << cc_edge_bias_decrement << std::endl;
        break;
    default:
        os << "Unknown decoding mode" << std::endl;
        break;
    }
}
#endif

//=================================================================================================
bool Options::parse_fir_filter_model_data(std::string fir_data_str)
{
    bool status = false;

    size_t colon_pos = fir_data_str.find(":");

    if ((fir_data_str.length() > 3) && (colon_pos != std::string::npos) && (colon_pos > 0))
    {
        std::vector<wsgc_float> raw_float_parameters;
        std::vector<std::string> raw_string_parameters;
        std::string parameter_str = fir_data_str.substr(colon_pos+1);
        std::string filter_type(fir_data_str.substr(0,colon_pos));
        std::transform(filter_type.begin(), filter_type.end(), filter_type.begin(), ::toupper);

        if (filter_type == "RCOS")
        {
        	status = extract_vector<wsgc_float>(raw_float_parameters, parameter_str);

        	if (raw_float_parameters.size() > 2)
        	{
        		unsigned int nb_taps = int(raw_float_parameters[0]);
        		_fir_coef_generator = new FIRCoefGenerator_RCos(f_sampling, raw_float_parameters[2], raw_float_parameters[1], nb_taps);
        	}

            status = true;
        }
    }

    if (!status)
    {
        std::cout << "Invalid FIR filter model specification" << std::endl;
    }

    return status;
}

//=================================================================================================
bool Options::parse_fading_model_data(std::string fading_data_str)
{
    bool status = false;

    size_t colon_pos = fading_data_str.find(":");

    if ((fading_data_str.length() > 3) && (colon_pos != std::string::npos) && (colon_pos > 0))
    {
        std::vector<wsgc_float> raw_float_parameters;
        std::vector<std::string> raw_string_parameters;
        std::string parameter_str = fading_data_str.substr(colon_pos+1);
        status = true;
        
        switch(fading_data_str[0])
        {
            case 'C': // Clarke
            case 'c':
                status = extract_vector<wsgc_float>(raw_float_parameters, parameter_str);
                
                if (raw_float_parameters.size() > 1)
                {
                    unsigned int nb_paths = int(raw_float_parameters[0]);
                    _fading_model = new FadingModelClarke(nb_paths, raw_float_parameters[1], f_sampling);
                    status = true;
                }
                
                break;
            
            case 'W': // Watterson
            case 'w':
                status = WsgcUtils::extract_string_vector(raw_string_parameters, parameter_str);
                
                if (status)
                {
                    _fading_model = new FadingModelWatterson(f_sampling);
                
                    for(std::vector<std::string>::iterator it = raw_string_parameters.begin(); it !=  raw_string_parameters.end(); ++it)
                    {
                        status = extract_vector<wsgc_float>(raw_float_parameters, *it);
                        
                        if ((status) && (raw_float_parameters.size() > 2))
                        {
                            if (raw_float_parameters[2] > f_sampling * (5*5) * 10) // Doppler spread frequency > 10 times maximum Doppler sampling frequency (which is 1/25th of the sampling frequency)
                            {
                                status = false;
                                break;
                            }
                            else
                            {
                                static_cast<FadingModelWatterson*>(_fading_model)->add_path_description(raw_float_parameters[0], raw_float_parameters[1], raw_float_parameters[2], raw_float_parameters[3]);
                            }
                        }
                        else
                        {
                            status = false;
                            break;
                        }    

                        raw_float_parameters.clear();
                    }
                    
                    if (status)
                    {
                        static_cast<FadingModelWatterson*>(_fading_model)->calculate_paths_description();
                    }
                }
                                
                break;

            default:
                break;
        }
    }
    
    if (!status)
    {
        std::cout << "Invalid fading model specification" << std::endl;
    }
    
    return status;
}


//=================================================================================================
bool Options::parse_modulation_data(std::string modulation_data_str)
{
	std::transform(modulation_data_str.begin(), modulation_data_str.end(), modulation_data_str.begin(), ::toupper);

	if (modulation_data_str == "BPSK")
	{
		modulation.setScheme(Modulation::Modulation_BPSK);
        transmission_scheme = OptionTrans_WSGC;
		return true;
	}
    if (modulation_data_str == "WSGCE")
    {
        modulation.setScheme(Modulation::Modulation_BPSK);
        transmission_scheme = OptionTrans_WSGCE;
        return true;
    }
	if (modulation_data_str == "DBPSK")
	{
		modulation.setScheme(Modulation::Modulation_DBPSK);
        transmission_scheme = OptionTrans_WSGCD;
		return true;
	}
	else if (modulation_data_str == "OOK")
	{
        modulation.setScheme(Modulation::Modulation_OOK);
        transmission_scheme = OptionTrans_WSGCO;
        return true;
	}
	else if (modulation_data_str == "CW")
	{
        modulation.setScheme(Modulation::Modulation_CW);
        transmission_scheme = OptionTrans_None;
        return true;
	}
	else if (modulation_data_str.compare(0,4,"MFSK") == 0)
	{
        transmission_scheme = OptionTrans_MFSK;
        
		if (modulation_data_str.size() < 6)
		{
			std::cout << "MFSK modulation specification is incorrect" << std::endl;
			return false;
		}
		else
		{
			size_t colon_pos = modulation_data_str.find(":");

			if (colon_pos == std::string::npos)
			{
				std::cout << "MFSK modulation specification is incorrect" << std::endl;
				return false;
			}
			else
			{
				modulation.setScheme(Modulation::Modulation_MFSK);
				std::string parameter_str = modulation_data_str.substr(colon_pos+1);
                mfsk_options._f_sampling = f_sampling;
				bool status = mfsk_options.parse_options(parameter_str);
				return status;
			}
		}
	}
	else
	{
        return false;
	}
}



//=================================================================================================
bool Options::parse_pilot_prns_data(std::string pilot_data_str)
{
    size_t comma_pos = pilot_data_str.find(",");

    try
    {
        if (comma_pos == std::string::npos) // one pilot
        {
            pilot1 = boost::lexical_cast<unsigned int>(pilot_data_str);
        }
        else if (comma_pos > 0)// two pilots
        {
        	pilot1 = boost::lexical_cast<unsigned int>(pilot_data_str.substr(0,comma_pos));
        	pilot2 = boost::lexical_cast<unsigned int>(pilot_data_str.substr(comma_pos+1));
        }
        else
        {
            std::cout << "Invalid pilot PRN number(s) specification" << std::endl;
            return false;
        }

        return true;
    }
    catch (boost::bad_lexical_cast &)
    {
        std::cout << "Invalid pilot PRN number(s) specification" << std::endl;
    }

    return false;
}

//=================================================================================================
bool Options::parse_source_coding_data(std::string source_coding_data_str)
{
    bool status = false;

    size_t colon1_pos = source_coding_data_str.find(":");
    size_t colon2_pos = source_coding_data_str.find(":", colon1_pos+1);

    if ((source_coding_data_str.length() > 3) && (colon1_pos != std::string::npos) && (colon1_pos > 0))
    {
    	if ((colon2_pos != std::string::npos) && (colon2_pos > colon1_pos))
    	{
    		_source_codec_type_str = source_coding_data_str.substr(0,colon1_pos);
    		std::transform(_source_codec_type_str.begin(), _source_codec_type_str.end(), _source_codec_type_str.begin(), ::toupper);
    		std::string source_coding_parameters_str = source_coding_data_str.substr(colon1_pos+1, colon2_pos-colon1_pos-1);
    		_source_message_str = source_coding_data_str.substr(colon2_pos+1);

    		if (_source_codec_type_str == "JT65")
    		{
    			// source_coding_data_str unused
    			_source_codec = new SourceCodec_JT65();
    			status = true;
    		}
            else if (_source_codec_type_str == "JT257")
            {
    			_source_codec = new SourceCodec_JT65(SourceCodec_JT65::JT257);
    			status = true;
            }
            else if (_source_codec_type_str == "JTCC")
            {
                _source_codec = new SourceCodec_JT65(SourceCodec_JT65::JTCC);
                status = true;
            }
    	}
    }

    if (!status)
    {
        std::cout << "Invalid source coding specification : " << source_coding_data_str << std::endl;
    }

    return status;
}

#ifdef _RSSOFT
//=================================================================================================
bool Options::parse_reed_solomon_data_generator(std::string parameter_str)
{
	std::vector<unsigned int> rs_num_parameters;
	bool status = extract_vector<unsigned int>(rs_num_parameters, parameter_str);

	if (rs_num_parameters.size() < 2)
	{
		std::cout << "Reed Solomon parameters take 2 values: q, k" << std::endl;
		return false;
	}

    rs_logq = rs_num_parameters[0];
    rs_k = rs_num_parameters[1];

    return true;
}

//=================================================================================================
bool Options::parse_reed_solomon_data(std::string parameter_str)
{
    bool status = false;
    std::vector<unsigned int> rs_num_parameters;
    std::vector<std::string> rs_str_parameters;
    size_t colon_pos = parameter_str.find(":");

    if (colon_pos == std::string::npos)
    {
    	std::cout << "Invalid Reed-Solomon specification" << std::endl;
    	return false;
    }

    std::string num_parameter_str = parameter_str.substr(0, colon_pos);
    std::string str_parameter_str = parameter_str.substr(colon_pos+1);

    status = extract_vector<unsigned int>(rs_num_parameters, num_parameter_str);

    if (rs_num_parameters.size() < 5)
    {
    	std::cout << "Reed Solomon parameters take 5 values: q, k, M, r, i" << std::endl;
    	return false;
    }

    unsigned int rs_n;
    rs_n = rs_num_parameters[0];
    rs_k = rs_num_parameters[1];
    rs_init_M = rs_num_parameters[2];
    rs_r = rs_num_parameters[3];
    rs_inc = rs_num_parameters[4];
    
    if (!is_power_of_2(rs_n+1, rs_logq))
    {
        std::cout << "Invalid Reed-Solomon size " << rs_n << std::endl;
        return false;
    }
    
    if (gc_nb_stages < rs_logq+1)
    {
        gc_nb_stages = rs_logq+1;
    }

    status = extract_vector<std::string>(rs_str_parameters, str_parameter_str);

    if (rs_str_parameters.size() < 2)
    {
    	std::cout << "Invalid Reed-Solomon specification" << std::endl;
    	return false;
    }

    if (rs_str_parameters[0] == "arith")
    {
    	rs_inc_strategy = RSSoft_Engine_defs::MMatrix_retry_arithmetic;
    }
    else if (rs_str_parameters[0] == "iarith")
    {
    	rs_inc_strategy = RSSoft_Engine_defs::MMatrix_retry_arithmetic_increment;
    }
    else if (rs_str_parameters[0] == "geom")
    {
    	rs_inc_strategy = RSSoft_Engine_defs::MMatrix_retry_geometric;
    }
    else if (rs_str_parameters[0] == "igeom")
    {
    	rs_inc_strategy = RSSoft_Engine_defs::MMatrix_retry_geometric_increment;
    }
    else
    {
    	std::cout << "Unrecognized Multiplicity Matrix increment strategy" << std::endl;
    	return false;
    }

    if (rs_str_parameters[1] == "all")
    {
    	rs_decoding_mode = RSSoft_Engine_defs::RSSoft_decoding_all;
    }
    else if (rs_str_parameters[1] == "full")
    {
    	rs_decoding_mode = RSSoft_Engine_defs::RSSoft_decoding_full;
    }
    else if (rs_str_parameters[1] == "best")
    {
    	rs_decoding_mode = RSSoft_Engine_defs::RSSoft_decoding_best;
    }
    else if (rs_str_parameters[1] == "first")
    {
    	rs_decoding_mode = RSSoft_Engine_defs::RSSoft_decoding_first;
    }
    else if (rs_str_parameters[1] == "regex")
    {
    	rs_decoding_mode = RSSoft_Engine_defs::RSSoft_decoding_regex;
    	if (rs_str_parameters.size() < 3)
    	{
    		std::cout << "Regular expression expected for RSSoft decoding mode with regular expression" << std::endl;
    		return false;
    	}
    	rs_decoding_match_str = rs_str_parameters[2];
    }
    else if (rs_str_parameters[1] == "match")
    {
    	rs_decoding_mode = RSSoft_Engine_defs::RSSoft_decoding_match;
    	if (rs_str_parameters.size() < 3)
    	{
    		std::cout << "Exact match string expected for RSSoft decoding mode with exact match" << std::endl;
    		return false;
    	}
    	rs_decoding_match_str = rs_str_parameters[2];
    }
    else if (rs_str_parameters[1] == "binmatch")
    {
    	rs_decoding_mode = RSSoft_Engine_defs::RSSoft_decoding_binmatch;
    }
    else if (rs_str_parameters[1] == "relthr")
    {
    	rs_decoding_mode = RSSoft_Engine_defs::RSSoft_decoding_relthr;
    	if (rs_str_parameters.size() < 3)
    	{
    		std::cout << "Regular expression expected for RSSoft decoding mode with regular expression" << std::endl;
    		return false;
    	}
        else
        {
            try
            {
                rs_reliability_threshold = boost::lexical_cast<float>(rs_str_parameters[2]);
            }
            catch (boost::bad_lexical_cast &)
            {
                std::cout << "Invalid format for reliability threshold (decimal number expected)" << std::endl;
                return false;
            }
        }
    }
    else
    {
    	std::cout << "Unrecognised RSSoft decoding mode" << std::endl;
    	return false;
    }

    return true;
}
#endif

#ifdef _CCSOFT
//=================================================================================================
bool Options::parse_convolutional_code_data_generator(std::vector<std::string> coding_data_strings)
{
    bool status = false;
    
    if (coding_data_strings.size() < 3)
    {
        std::cout << "Invalid specification for convolutional code" << std::endl;
        status = false;
    }
    else
    {
        status = parse_convolutional_code_constraints(coding_data_strings[1]);
        
        if (status)
        {
            status = parse_convolutional_code_genpolys(coding_data_strings[2]);
        }
    }
    
    return status;
}

//=================================================================================================
bool Options::parse_convolutional_code_data(std::vector<std::string> coding_data_strings)
{
    bool status = false;
    
    if (coding_data_strings.size() < 4)
    {
        std::cout << "Invalid specification for convolutional code" << std::endl;
        return false;
    }
    else
    {
        status = parse_convolutional_code_constraints(coding_data_strings[1]);
        
        if (status)
        {
            status = parse_convolutional_code_genpolys(coding_data_strings[2]);
        }
        else
        {
            std::cout << "Invalid convolutional code constraints specification" << std::endl;
            return false;
        }
        
        if (status)
        {
            std::vector<std::string> cc_specs;
            
            if (extract_vector(cc_specs, ":", coding_data_strings[3]))
            {
                if (cc_specs.size() < 2)
                {
                    std::cout << "Invalid convolutional code specification" << std::endl;
                    return false;
                }
                else
                {
                    std::vector<std::string> algo_specs;

                    if (extract_vector(algo_specs, ",", cc_specs[0], true))
                    {
                        if (algo_specs.size() < 1)
                        {
                            std::cout << "Invalid convolutional code decoding algorithm specification" << std::endl;
                            return false;
                        }
                        else
                        {
                            std::transform(algo_specs[0].begin(), algo_specs[0].end(), algo_specs[0].begin(), ::toupper);

                            if (algo_specs[0] == "STACK")
                            {
                                cc_algorithm_type = CCSoft_Engine_defs::Algorithm_Stack;
                            }
                            else if (algo_specs[0] == "FANO")
                            {
                                cc_algorithm_type = CCSoft_Engine_defs::Algorithm_Fano;
                            }

                            if (algo_specs.size() > 1)
                            {
                                std::transform(algo_specs[1].begin(), algo_specs[1].end(), algo_specs[1].begin(), ::toupper);

                                if (algo_specs[1] == "NI")
                                {
                                    cc_interleave = false;
                                }
                            }

                            if (algo_specs.size() > 2)
                            {
                                std::transform(algo_specs[2].begin(), algo_specs[2].end(), algo_specs[2].begin(), ::toupper);

                                if (algo_specs[2] == "REGEX")
                                {
                                    cc_decoding_mode = CCSoft_Engine_defs::Decoding_regex;
                                }
                                else if (algo_specs[2] == "MATCHSTR")
                                {
                                    cc_decoding_mode = CCSoft_Engine_defs::Decoding_match_str;
                                }
                                else if (algo_specs[2] == "MATCHMSG")
                                {
                                    cc_decoding_mode = CCSoft_Engine_defs::Decoding_match_msg;
                                }
                                else
                                {
                                    cc_decoding_mode = CCSoft_Engine_defs::Decoding_normal;
                                }
                            }

                            if (algo_specs.size() > 3)
                            {
                                cc_match_str = algo_specs[3];
                            }
                        }
                    }



                    std::vector<float> algo_parms;
                    
                    if (!extract_vector<float>(algo_parms, ",", cc_specs[1]))
                    {
                        std::cout << "Invalid convolutional code decoding algorithm parameters specification" << std::endl;
                        return false;
                    }
                    else
                    {
                        if (cc_algorithm_type == CCSoft_Engine_defs::Algorithm_Stack)
                        {
                            if (algo_parms.size() > 0)
                            {
                                cc_edge_bias = algo_parms[0];
                            }
                            if (algo_parms.size() > 1)
                            {
                                cc_use_node_limit = true;
                                cc_node_limit = algo_parms[1];
                            }
                            if (algo_parms.size() > 2)
                            {
                                cc_metric_limit = algo_parms[2];
                                cc_use_metric_limit = true;
                            }
                            if (algo_parms.size() > 3)
                            {
                                cc_nb_retries = algo_parms[3];
                            }
                            if (algo_parms.size() > 4)
                            {
                                cc_edge_bias_decrement = algo_parms[4];
                            }
                        }
                        else if (cc_algorithm_type == CCSoft_Engine_defs::Algorithm_Fano)
                        {
                            if (algo_parms.size() > 0)
                            {
                                cc_edge_bias = algo_parms[0];
                            }
                            if (algo_parms.size() > 1)
                            {
                                cc_fano_init_metric = algo_parms[1];
                            }
                            if (algo_parms.size() > 2)
                            {
                                cc_fano_delta_metric = algo_parms[2];
                            }
                            if (algo_parms.size() > 3)
                            {
                                cc_fano_tree_cache_size = algo_parms[3];
                            }
                            if (algo_parms.size() > 4)
                            {
                                cc_fano_delta_init_threshold = algo_parms[4];
                            }
                            if (algo_parms.size() > 5)
                            {
                                cc_nb_retries = algo_parms[5];
                            }
                            if (algo_parms.size() > 6)
                            {
                                cc_edge_bias_decrement = algo_parms[6];
                            }
                            if (algo_parms.size() > 7)
                            {
                                cc_use_node_limit = true;
                                cc_node_limit = algo_parms[7];
                            }
                            if (algo_parms.size() > 8)
                            {
                                cc_use_metric_limit = true;
                                cc_metric_limit = algo_parms[8];
                            }
                        }
                        else
                        {
                            std::cout << "Unrecognised convolutional code decoding algorithm code " << algo_specs[0] << std::endl;
                            return false;
                        }
                    }
                }
            }
            else
            {
                std::cout << "Invalid convolutional code decoding algorithm specification" << std::endl;
                return false;
            }
        }
        else
        {
            std::cout << "Invalid convolutional code generator polynomials specification" << std::endl;
            return false;
        }

        return true;
    }
}

//=================================================================================================
bool Options::parse_convolutional_code_constraints(std::string constraints_str)
{
    bool status = extract_vector<unsigned int>(cc_k_constraints, ",", constraints_str);
    return status;
}

//=================================================================================================
bool Options::parse_convolutional_code_genpolys(std::string genpolys_str)
{
    std::vector<std::string> g_strings;

    if (!extract_vector(g_strings, ":", genpolys_str))
    {
    	std::cout << "Invalid generator polynomials specification" << std::endl;
        return false;
    }

    std::vector<std::string>::const_iterator gs_it = g_strings.begin();

    for (; gs_it != g_strings.end(); ++gs_it)
    {
        std::vector<unsigned int> g;

        if (extract_vector<unsigned int>(g, ",", *gs_it))
        {
            cc_generator_polys.push_back(g);
        }
        else
        {
            std::cout << "Invalid generator polynomial specification" << std::endl;
            return false;
        }
    }

    return true;
}

#endif

//=================================================================================================
bool Options::adjust_parameters_for_source_coding()
{
	bool status = false;

	if (_source_codec_type_str == "JT65")
	{
		nb_message_symbols = 64;

		if (gc_nb_stages < 7)
		{
			gc_nb_stages = 7;
		}

		if (rs_k != 0)
		{
			rs_k = 12;
			rs_logq = 6;
		}

		status = true;
	}
	else if (_source_codec_type_str == "JT257")
	{
		nb_message_symbols = 256;

		if (gc_nb_stages < 9)
		{
			gc_nb_stages = 9;
		}

		if (rs_k != 0)
		{
			rs_k = 9;
			rs_logq = 8;
		}

		status = true;
    }
    else if (_source_codec_type_str == "JTCC")
    {
        nb_message_symbols = 2;
        status = true;
    }

	return status;
}

//=================================================================================================
bool Options::source_codec_create_message_prns()
{
	prns.clear(); // clear anything constructed before
	bool status = _source_codec->encode(_source_message_str, prns);

	if (!status)
	{
		std::cout << "Cannot encode source message" << std::endl;
	}

	return status;
}

//=================================================================================================
bool Options::is_power_of_2(unsigned int x, unsigned int& log2)
{
    unsigned int i = 0;

    while (((x & 1) == 0) && x > 1) /* While x is even and > 1 */
    {
        x >>= 1;
        i++;
    }
    
    log2 = i;
    return (x == 1);
}

//=================================================================================================
bool Options::parse_fec_option(std::string fec_option_str)
{
    std::vector<std::string> fec_elements;
    bool status = extract_vector(fec_elements, "/", fec_option_str);

    if (status)
    {
        if (fec_elements.size() > 1)
        {
            if (fec_elements[0] == "RS")
            {
#ifdef _RSSOFT
                fec_scheme = OptionFEC_RSSoft;

                if (options_executable == Options_wsgc_test)
                {
                    status = parse_reed_solomon_data(fec_elements[1]);
                }
                else if (options_executable == Options_wsgc_generator)
                {
                    status = parse_reed_solomon_data_generator(fec_elements[1]);
                }
                else
                {
                    std::cout << "Unrecognised options executable" << std::endl;
                    status = false;
                }
#else
                std::cout << "RSSoft library not supported" << std::endl;
                status = false;
#endif
            }
            else if (fec_elements[0] == "CC")
            {
#ifdef _CCSOFT
                fec_scheme = OptionFEC_CCSoft;

                if (options_executable == Options_wsgc_test)
                {
                    status = parse_convolutional_code_data(fec_elements);
                }
                else if (options_executable == Options_wsgc_generator)
                {
                    status = parse_convolutional_code_data_generator(fec_elements);
                }
                else
                {
                    std::cout << "Unrecognised options executable" << std::endl;
                    status = false;
                }
#else
                std::cout << "CCSoft library not supported" << std::endl;
                status = false;
#endif            
            }
            else
            {
                std::cout << "unrecognised FEC scheme" << std::endl;
                status = false;
            }
        }
        else
        {
            std::cout << "invalid FEC option" << std::endl;
            status = false;
        }
    }

    return status;
}
