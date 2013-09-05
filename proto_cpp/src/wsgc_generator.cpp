#include "WsgcTypes.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <time.h>
#include <cmath>
#include "Options.h"
#include "GoldCodeGenerator.h"
#include "Transmission.h"
#include "WsgcUtils.h"
#include "Demodulator.h"
#include "DemodulatorDifferential.h"
#include "DemodulatorSquaring.h"

void apply_fir(wsgc_complex *inout, unsigned int& nb_samples, const std::vector<wsgc_float>& fir_coef);

int main(int argc, char *argv[])
{
    timespec time1, time2;
    int time_option = CLOCK_REALTIME;
    std::string binary_name(argv[0]);
    Options options(binary_name, Options::Options_wsgc_generator);

    if (options.get_options(argc,argv))
    {
        std::ostringstream os;
        options.print_options(os);
        std::cout << os.str() << std::endl;

        GoldCodeGenerator gc_generator(options.gc_nb_stages, options.nb_message_symbols, options.nb_service_symbols, options.nb_training_symbols, options.g1_poly_powers, options.g2_poly_powers);

        unsigned int fft_N = gc_generator.get_nb_code_samples(options.f_sampling, options.f_chip);

        Transmission transmission(options, gc_generator);
        transmission.generate_samples();
        wsgc_complex *faded_source_samples = transmission.get_samples();
        unsigned int nb_faded_source_samples = transmission.get_nb_samples();

        if (faded_source_samples)
        {
            std::cout << "Normalize" << std::endl;
            transmission.normalize();

            // demodulate OOK
            if (options.simulate_demod)
            {
            	Demodulator *demodulator;

            	if (options.modulation.getScheme() == Modulation::Modulation_OOK)
				{
					demodulator = new DemodulatorSquaring();
					std::cout << "Simulate AM power detection" << std::endl;
				}
				else if (options.modulation.isDifferential())
				{
					unsigned int int_samples_per_chip = ((wsgc_float) fft_N) /  gc_generator.get_code_length();
					static const wsgc_complex c_zero(0.0, 0.0);

					demodulator = new DemodulatorDifferential(int_samples_per_chip);
				    ((DemodulatorDifferential *) demodulator)->set_value_at_origin(c_zero);
				}

			    demodulator->demodulate_in_place(faded_source_samples, nb_faded_source_samples);

			    delete demodulator;
            }
        
            // Write out samples
            std::ofstream ofs;
            ofs.open(options.samples_output_file.c_str(), std::ios::out | std::ios::binary);

            for (unsigned int si = 0; si < nb_faded_source_samples; si++)
            {
                float sample_r = faded_source_samples[si].real();
                float sample_i = faded_source_samples[si].imag();
            	ofs.write((const char *) &sample_r, sizeof(wsgc_complex));
            	ofs.write((const char *) &sample_i, sizeof(wsgc_complex));
            }

            // close file, that's it!
            ofs.close();
        }
        return 0;
    }
    else
    {
        std::cout << "Incorrect options. Please correct and re-submit" << std::endl;
        return -1;
    }
}
