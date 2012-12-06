#include "WsgcTypes.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <time.h>
#include "Options.h"
#include "ContinuousPhaseCarrier.h"
#include "WsgcUtils.h"
#include "FadingModel.h"

// Playback with SOX's play: play -t raw -e float -b 32 -r 8000 -L samples.raw
// assuming sampling rate is 8000 Hz

int main(int argc, char *argv[])
{
    Options options;

    if (options.get_options(argc,argv))
    {
        std::ostringstream os;
        options.print_options(os);
        std::cout << os.str() << std::endl;
        
        if (options.has_fading())
        {
            std::cout << "Processing..." << std::endl;
            std::cout << sizeof(wsgc_float) << std::endl;
            
            std::ofstream samples_file;
            samples_file.open("samples.raw", std::ios::binary);
            
            FadingModel *fading_model = options.get_fading_model();
            unsigned int samples_per_symbol = int(options.f_sampling*options.nb_prns_per_symbol);
            unsigned int output_samples_per_symbol = fading_model->get_output_size(samples_per_symbol);
            wsgc_complex *faded_signal_samples;
            const wsgc_complex *signal_samples;
            
            if (samples_per_symbol == output_samples_per_symbol) // output synchronous from input, OK to process by block
            {
                ContinuousPhaseCarrier cw(options.f_sampling, samples_per_symbol);
                faded_signal_samples = new wsgc_complex[output_samples_per_symbol]; 
                
                for (unsigned int symbol_i=0; symbol_i<options.prns.size(); symbol_i++)
                {
                    cw.make_next_samples(options.f_tx);
                    signal_samples = cw.get_samples();
                    fading_model->apply_fading(signal_samples, faded_signal_samples, samples_per_symbol);
                    
                    //write real part to raw file
                    for (unsigned int s_i=0; s_i<samples_per_symbol; s_i++)
                    {
                        static const wsgc_float test = 0.0;
                        samples_file.write((const char*) &faded_signal_samples[s_i].real(), sizeof(wsgc_float));
                        samples_file.write((const char*) &faded_signal_samples[s_i].imag(), sizeof(wsgc_float)); // stereo binaural
                    }
                }
            }
            else // need to process globally
            {
                unsigned int nb_samples = samples_per_symbol*options.prns.size();
                unsigned int nb_out_samples = fading_model->get_output_size(nb_samples);
                faded_signal_samples = new wsgc_complex[nb_out_samples]; 
                ContinuousPhaseCarrier cw(options.f_sampling, nb_samples);
                cw.make_next_samples(options.f_tx);
                signal_samples = cw.get_samples();
                fading_model->apply_fading(signal_samples, faded_signal_samples, nb_samples);
                
                //write real part to raw file
                for (unsigned int s_i=0; s_i<nb_out_samples; s_i++)
                {
                    samples_file.write((const char*) &faded_signal_samples[s_i].real(), sizeof(wsgc_float));
                }              
            }
            
            samples_file.close();
            delete[] faded_signal_samples;
        }
        else
        {
            std::cout << "No fading option, nothing to do..." << std::endl;
        }
    }    
    else
    {
        std::cout << "Incorrect options. Please correct and re-submit" << std::endl;
        return -1;
    }
}

