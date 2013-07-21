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
#include "CodeModulator_DBPSK.h"
#include "CodeModulator_OOK.h"
#include "CodeModulator_CW_Test.h"
#include "LocalCodesFFT.h"
#include "LocalCodes.h"
#include "WsgcUtils.h"
#include "SampleSequencer.h"
#include "SourceMixer.h"
#include "FadingModel.h"
#include "FIR_RCoef.h"
#include "FIRCoefGenerator.h"
#include "Demodulator.h"
#include "DemodulatorDifferential.h"
#include "DemodulatorSquaring.h"
#include "CodeModulator_MFSK.h"

void apply_fir(wsgc_complex *inout, unsigned int& nb_samples, const std::vector<wsgc_float>& fir_coef);

int main(int argc, char *argv[])
{
    timespec time1, time2;
    int time_option = CLOCK_REALTIME;
    std::string binary_name(argv[0]);
    Options options(binary_name);

    if (options.get_options(argc,argv))
    {
        std::ostringstream os;
        options.print_options(os);
        std::cout << os.str() << std::endl;

        GoldCodeGenerator gc_generator(options.gc_nb_stages, options.nb_message_symbols, options.nb_service_symbols, options.nb_training_symbols, options.g1_poly_powers, options.g2_poly_powers);
        CodeModulator *codeModulator = 0;
        CodeModulator_MFSK *codeModulator_MFSK = 0;
        wsgc_complex *source_samples = 0;
        unsigned int nb_source_samples = 0;
        SimulatedSource *message_source = 0;
        SourceMixer *source_mixer = 0;
        std::vector<unsigned int> pilot_prns;
        unsigned int fft_N = gc_generator.get_nb_code_samples(options.f_sampling, options.f_chip);

        // Allocate appropriate objects depending on modulation options

        // create code modulator
        if (options.modulation.getScheme() == Modulation::Modulation_BPSK)
        {
            codeModulator = new CodeModulator_BPSK();
        }
        if (options.modulation.getScheme() == Modulation::Modulation_DBPSK)
        {
            codeModulator = new CodeModulator_DBPSK();
        }
        else if (options.modulation.getScheme() == Modulation::Modulation_OOK)
        {
            codeModulator = new CodeModulator_OOK();
        }
        else if (options.modulation.getScheme() == Modulation::Modulation_CW)
        {
            codeModulator = new CodeModulator_CW_Test();
        }
        else if (options.modulation.getScheme() == Modulation::Modulation_MFSK)
        {
            codeModulator_MFSK = new CodeModulator_MFSK(
            		options.mfsk_options._f_sampling,
            		options.mfsk_options._zero_frequency + options.f_tx,
            		options.mfsk_options._symbol_bandwidth,
            		options.mfsk_options._symbol_time);
        }

        if (codeModulator) // if a code modulator is available then the actual signal processing can take place
        {
            // Produce signal samples
            std::cout << "Produce signal samples..." << std::endl;

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
        }
        else if (codeModulator_MFSK != 0)
        {
        	nb_source_samples = codeModulator_MFSK->get_nb_symbol_samples()*options.prns.size();
        	source_samples = (wsgc_complex *) WSGC_FFTW_MALLOC(nb_source_samples*sizeof(wsgc_fftw_complex));
        	codeModulator_MFSK->modulate(reinterpret_cast<wsgc_fftw_complex*>(source_samples), options.prns);
        }
        else
        {
            std::cout << "Code modulator not implemented for this modulation scheme" << std::endl;
        }

        if (source_samples)
        {
            LocalCodesFFT *local_codes_fft = 0;
            LocalCodes *local_codes = 0;

            // Apply lowpass filter if any
            if (options._fir_coef_generator != 0)
            {
            	apply_fir(source_samples, nb_source_samples, options._fir_coef_generator->get_coefs());
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
            /*
            if (options.modulation.getScheme() == Modulation::Modulation_OOK)
            { // trim imaginary part
            	std::cout << "Simulate AM power detection" << std::endl;
                for (unsigned int i = 0; i<nb_faded_source_samples; i++)
                {
                    faded_source_samples[i].real() = std::norm(faded_source_samples[i]);
                    faded_source_samples[i].imag() = 0.0;
                }
            }
            */

            std::cout << "Normalize" << std::endl;
            fading->normalize(faded_source_samples, nb_faded_source_samples);

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


            // Generate samples
            std::ofstream ofs;
            ofs.open(options.samples_output_file.c_str(), std::ios::out | std::ios::binary);

            for (unsigned int si = 0; si < nb_faded_source_samples; si++)
            {
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
            std::cout << "No source samples were generated" << std::endl;
        }

        return 0;
    }
    else
    {
        std::cout << "Incorrect options. Please correct and re-submit" << std::endl;
        return -1;
    }
}


//=================================================================================================
void apply_fir(wsgc_complex *inout, unsigned int& nb_samples, const std::vector<wsgc_float>& fir_coef)
{
	std::cout << "Apply lowpass FIR filter" << std::endl;

	FIR_RCoef fir_filter(fir_coef);
	static const wsgc_complex c_zero = (0.0, 0.0);

	for (unsigned int i = 0; i<nb_samples; i++)
	{
		inout[i] = fir_filter.calc(inout[i]);
	}

	/* tail
	for (unsigned int i = 0; i<nb_taps; i++)
	{
		inout[i] = fir_filter.calc(c_zero);
	}

	nb_samples += nb_taps;
	*/
}
