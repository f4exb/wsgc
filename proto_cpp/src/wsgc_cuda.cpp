#include "test_cuda.h"
#include "ContinuousPhaseCarrier.h"
#include "GoldCodeGenerator.h"
#include "SimulatedSource.h"
#include "CodeModulator_BPSK.h"

#include <iostream>
#include <sstream>

void fill_options(options_t& options);
bool fill_poly_powers(unsigned int gc_nb_stages, std::vector<unsigned int>& g1_poly_powers, std::vector<unsigned int>& g2_poly_powers);

int main(int argc, char *argv[])
{
	options_t options;
	fill_options(options);

	test_cuda test_cuda_instance(options);

	std::cout << "-- test1 : Source FFT --" << std::endl;
	test_cuda_instance.test1();

    std::cout << "-- test2 : Frequency domain correlation --" << std::endl;

    GoldCodeGenerator gc_generator(options.gc_nb_stages, options.nb_message_symbols, options.nb_service_symbols, options.g1_poly_powers, options.g2_poly_powers);
    SimulatedSource message_source(gc_generator, options.prn_list, options.f_sampling, options.f_chip, options.f_tx, options.code_shift);
    CodeModulator_BPSK code_modulator;
    message_source.set_code_modulator(&code_modulator);
    message_source.create_samples();

    test_cuda_instance.test2(message_source.get_samples(), gc_generator, code_modulator);

    std::cout << "-- test3 : mul_ifft algo --" << std::endl;

    test_cuda_instance.test3(message_source.get_samples(), gc_generator, code_modulator);

	return 0;
}


void fill_options(options_t& options)
{
	options.nb_message_symbols = 16;
	options.nb_service_symbols = 2;
	options.nb_pilot_prns = 1;
	options.f_sampling = 128;
	options.f_chip = 31;
	options.fft_N = 128;
	options.freq_step_division = 8;
	options.gc_nb_stages = 5;

	fill_poly_powers(options.gc_nb_stages, options.g1_poly_powers, options.g2_poly_powers);

	options.prn_list.push_back(0);
	options.f_tx = 0.0;
	options.code_shift = 10;
	options.nb_batch_prns = 3;
	options.nb_f_bins = 2;
}


bool fill_poly_powers(unsigned int gc_nb_stages, std::vector<unsigned int>& g1_poly_powers, std::vector<unsigned int>& g2_poly_powers)
{
    switch(gc_nb_stages)
    {
        case 5:
            g1_poly_powers.push_back(2); // G1 = X^5 + X^2 + 1
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
            g1_poly_powers.push_back(3); // G1 = X^10 + X^3 + 1
            g2_poly_powers.push_back(8); // G2 = X^10 + X^8 + X^3 + X^2 + 1
            g2_poly_powers.push_back(3);
            g2_poly_powers.push_back(2);
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
