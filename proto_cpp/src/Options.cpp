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
#include "FadingModelNone.h"
#include "FadingModelClarke.h"
#include "FadingModelWatterson.h"
#include "FIRCoefGenerator_RCos.h"
#include <stdlib.h>
#include <getopt.h>
#include <cstring>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <boost/lexical_cast.hpp>


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


Options::Options(std::string& _binary_name) :
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
	_fir_coef_generator(0)
{
    srand(time(0));

    std::vector<std::string> path_nodes;
    WsgcUtils::extract_string_vector(path_nodes, binary_path);
    binary_name = path_nodes.back();

    std::cout << binary_name << std::endl;
}


Options::~Options()
{
    if (_fading_model == 0)
    {
        delete _fading_model;
    }
}


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
        };
        
        int option_index = 0;
        
        c = getopt_long (argc, argv, "s:c:n:t:C:p:N:I:r:R:T:m:M:S:g:G:f:d:a:P:A:F:B:z:U:o:L:", long_options, &option_index);
        
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
                    std::cout << "Using CUDA implementation" << std::endl;
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
            
            if (nb_prns_per_symbol < 4)
            {
                std::cout << "Need at least 4 PRN codes per symbol (-N option)" << std::endl;
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
            else // Build PRN list for message simulation
            {
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
    
    os << "Modulation ................: "; modulation.print_modulation_data(os); os << std::endl;
    //os << "Nb phase averaging cycles .: " << std::setw(6) << std::right << tracking_phase_average_cycles << std::endl;
    os << "Nb message symbols ........: " << std::setw(6) << std::right << nb_message_symbols << std::endl;
    os << "Nb service symbols ........: " << std::setw(6) << std::right << nb_service_symbols << std::endl;
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
    	os << "Analysis window size ......: " << std::setw(6) << std::right << analysis_window_size*nb_prns_per_symbol << " PRNs" << std::endl;
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
    os << " - " << (use_cuda ? "Using CUDA implementation" : "Using Host implementation") << std::endl;
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
        os << "Performing noise measurements" << std::endl;
    }
    
    os << std::endl;
    os << "Processing " << prns.size() << " symbols: ";
    print_vector<unsigned int, unsigned int>(prns, 4, os); os << std::endl;
}


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


bool Options::parse_modulation_data(std::string modulation_data_str)
{
	std::transform(modulation_data_str.begin(), modulation_data_str.end(), modulation_data_str.begin(), ::toupper);

	if (modulation_data_str == "BPSK")
	{
		modulation.setScheme(Modulation::Modulation_BPSK);
		return true;
	}
	if (modulation_data_str == "DBPSK")
	{
		modulation.setScheme(Modulation::Modulation_DBPSK);
		return true;
	}
	else if (modulation_data_str == "OOK")
	{
        modulation.setScheme(Modulation::Modulation_OOK);
        return true;
	}
	else if (modulation_data_str == "CW")
	{
        modulation.setScheme(Modulation::Modulation_CW);
        return true;
	}
	else
	{
        return false;
	}

	/*
    switch(modulation_data_str[0])
    {
        case 'B':
        case 'b':
            modulation.setScheme(Modulation::Modulation_BPSK);
            return true;
        case 'O':
        case 'o':
            modulation.setScheme(Modulation::Modulation_OOK);
            return true;
        default:
            return false;
    }
    */
}



bool Options::parse_pilot_prns_data(std::string pilot_data_str)
{
    bool status = false;

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
        return false;
    }

}
