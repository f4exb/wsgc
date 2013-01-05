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
#include "SimulatedSource.h"
#include "CodeModulator_BPSK.h"
#include "CodeModulator_OOK.h"
#include "CodeModulator_CW_Test.h"
#include "LocalCodesFFT.h"
#include "LocalCodes.h"
#include "WsgcUtils.h"
#include "SampleSequencer.h"
#include "SourceMixer.h"

int main(int argc, char *argv[])
{
    timespec time1, time2;
    int time_option = CLOCK_REALTIME;
    Options options;

    if (options.get_options(argc,argv))
    {
        std::ostringstream os;
        options.print_options(os);
        std::cout << os.str() << std::endl;

        GoldCodeGenerator gc_generator(options.gc_nb_stages, options.nb_message_symbols, options.nb_service_symbols, options.g1_poly_powers, options.g2_poly_powers);
        CodeModulator *codeModulator = 0;
        SimulatedSource *message_source = 0;
        SourceMixer *source_mixer = 0;
        std::vector<unsigned int> pilot_prns;

        // Allocate appropriate objects depending on modulation options

        // create code modulator
        if (options.modulation.getScheme() == Modulation::Modulation_BPSK)
        {
            codeModulator = new CodeModulator_BPSK();
        }
        else if (options.modulation.getScheme() == Modulation::Modulation_OOK)
        {
            codeModulator = new CodeModulator_OOK();
        }
        else if (options.modulation.getScheme() == Modulation::Modulation_CW)
        {
            codeModulator = new CodeModulator_CW_Test();
        }

        if (codeModulator) // if a code modulator is available then the actual signal processing can take place
        {
            // Produce signal samples
            std::cout << "Produce signal samples..." << std::endl;

            LocalCodesFFT *local_codes_fft = 0;
            LocalCodes *local_codes = 0;
            wsgc_complex *source_samples = 0;
            unsigned int nb_source_samples = 0;

            message_source = new SimulatedSource(gc_generator, options.prns, options.f_sampling, options.f_chip,
                                                 options.f_tx, options.code_shift, options.nb_prns_per_symbol, 0.0);
            message_source->set_code_modulator(codeModulator);
            message_source->create_samples();

            if (options.modulation.isCodeDivisionCapable() && options.nb_pilot_prns > 0)
            {
            	pilot_prns.assign(options.prns.size() + (options.batch_size/options.nb_prns_per_symbol), options.pilot1); // use only pilot 1 for the simulation
            	SimulatedSource pilot_source(gc_generator, pilot_prns, options.f_sampling, options.f_chip,
                                             options.f_tx, options.code_shift, options.nb_prns_per_symbol, 0.0);
            	pilot_source.set_code_modulator(codeModulator);
            	pilot_source.create_samples();

            	wsgc_float pilot_gain = pow(10.0, (options.pilot_gain_db / 10.0));
            	source_mixer = new SourceMixer(message_source, &pilot_source, pilot_gain);

            	source_samples = source_mixer->get_samples();
            	nb_source_samples = source_mixer->get_nb_samples();
            }
            else
            {
            	source_samples = message_source->get_samples();
            	nb_source_samples = message_source->get_nb_samples();
            }

            // get fading model
            FadingModel *fading = options._fading_model;
            wsgc_complex *faded_source_samples;
            unsigned int nb_faded_source_samples; // fading may add delay hence there could be more samples after the fading process

            // apply fading
            if (fading->is_fading_active())
            {
            	std::cout << "Apply fading" << std::endl;
                nb_faded_source_samples = fading->get_output_size(nb_source_samples);
                faded_source_samples = (wsgc_complex *) WSGC_FFTW_MALLOC(nb_faded_source_samples*sizeof(wsgc_fftw_complex));
                fading->apply_fading(source_samples, faded_source_samples, nb_source_samples);
            }
            else
            {
            	std::cout << "Do not apply fading" << std::endl;
                faded_source_samples = source_samples;
                nb_faded_source_samples = nb_source_samples;
            }

            // apply AWGN
            if (options.make_noise)
            {
            	std::cout << "Apply AWGN" << std::endl;
                fading->apply_awgn(faded_source_samples, nb_faded_source_samples, options.code_shift, options.snr);
            }

            // apply power detection for OOK
            if (options.modulation.getScheme() == Modulation::Modulation_OOK)
            { // trim imaginary part
            	std::cout << "Simulate AM power detection" << std::endl;
                for (unsigned int i = 0; i<nb_faded_source_samples; i++)
                {
                    faded_source_samples[i].real() = std::norm(faded_source_samples[i]);
                    faded_source_samples[i].imag() = 0.0;
                }
            }

            std::cout << "Normalize" << std::endl;
            fading->normalize(faded_source_samples, nb_faded_source_samples);

            // Generate samples
            std::ofstream ofs;
            ofs.open(options.samples_output_file.c_str(), std::ios::out | std::ios::binary);

            //ContinuousPhaseCarrier test_lo(options.f_sampling, options.nb_samples_per_code);

            for (unsigned int si = 0; si < nb_faded_source_samples; si++)
            {
            	/*
            	if (si % options.nb_samples_per_code ==0)
            	{
            		if (si != 0)
            		{
            			const wsgc_complex *samples = test_lo.get_samples();

            			for (unsigned int i=0; i<options.nb_samples_per_code; i++)
            			{
            				ofs.write((const char *) &samples[i].real(), sizeof(samples[i].real()));
            				ofs.write((const char *) &samples[i].imag(), sizeof(samples[i].imag()));
            			}
            		}

            		test_lo.make_next_samples(options.f_tx);
            	}
            	*/
            	ofs.write((const char *) &faded_source_samples[si].real(), sizeof(faded_source_samples[si].real()));
            	ofs.write((const char *) &faded_source_samples[si].imag(), sizeof(faded_source_samples[si].imag()));
            }

            // close file, that's it!
            ofs.close();

            if (local_codes_fft)
            {
                delete local_codes_fft;
            }

            if (local_codes)
            {
                delete local_codes;
            }

            if (codeModulator)
            {
            	delete codeModulator;
            }

        }
        else
        {
            std::cout << "Code modulator not implemented for this modulation scheme" << std::endl;
        }

        return 0;
    }
    else
    {
        std::cout << "Incorrect options. Please correct and re-submit" << std::endl;
        return -1;
    }
}

