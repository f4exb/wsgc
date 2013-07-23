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

     Main test program

*/
#include "WsgcTypes.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <time.h>
#include <cmath>
#include "Options.h"
#include "GoldCodeGenerator.h"
#include "SimulatedSource.h"
#include "CodeModulator_BPSK.h"
#include "CodeModulator_DBPSK.h"
#include "CodeModulator_OOK.h"
#include "CodeModulator_MFSK.h"
#include "LocalCodes.h"
#include "LocalCodes_Host.h"
#include "WsgcUtils.h"
#include "PilotCorrelationRecord.h"
#include "AutocorrelationRecord.h"
#include "PilotedTrainingMultiplePrnCorrelator.h"
#include "PilotedTrainingMessageCorrelator.h"
#ifdef _CUDA
#include "SinglePrnCorrelator_FreqDep_Cuda.h"
#include "PilotedMessageCorrelator_Cuda.h"
#include "PilotedTrainingMessageCorrelator_Cuda.h"
#include "LocalCodes_Cuda.h"
#include "LocalCodesFFT_Cuda.h"
#include "PilotCorrelator_Cuda.h"
#include "SourceFFT_Cuda.h"
#endif
#include "SinglePrnCorrelator_FreqDep_Host.h"
#include "LocalCodesFFT_Host.h"
#include "PilotCorrelator_Host.h"
#include "SourceFFT_Host.h"
#include "PilotCorrelator.h"
#include "PilotedMessageCorrelator.h"
#include "PilotedMessageCorrelator_Host.h"
#include "PrnAutocorrelator.h"
#include "PrnAutocorrelator_Host.h"
#include "PilotCorrelationAnalyzer.h"
#include "PilotedMultiplePrnCorrelator.h"
#include "PilotedTrainingMessageCorrelator_Host.h"
#include "DecisionBox.h"
#include "DecisionBox_Piloted_And_Synced.h"
#include "DecisionBox_Unpiloted_And_Synced.h"
#include "SampleSequencer.h"
#include "SourceMixer.h"
#include "FadingModel.h"
#include "FIR_RCoef.h"
#include "FIRCoefGenerator.h"
#include "UnpilotedMultiplePrnCorrelator.h"
#include "UnpilotedMultiplePrnCorrelator_Host.h"
#include "UnpilotedTrainingMessageCorrelator.h"
#include "UnpilotedTrainingMessageCorrelator_Host.h"
#include "DemodulatorDifferential.h"
#include "DemodulatorSquaring.h"
#include "MFSK_MessageDemodulator.h"
#include "MFSK_MessageDemodulator_Host.h"
#include "DecisionBox_MFSK.h"
#ifdef _RSSOFT
#include "RSSoft_DecisionBox.h"
#endif

#ifdef _CUDA
#include "CudaManager.h"
#endif

void message_processing(
#ifdef _CUDA
    CudaManager& cuda_manager,
#endif    
    Options& options, 
    GoldCodeGenerator& gc_generator,
    CodeModulator *localCodeModulator,
    wsgc_complex *faded_source_samples,
    unsigned int nb_faded_source_samples
);

void training_processing(
#ifdef _CUDA
    CudaManager& cuda_manager,
#endif    
    Options& options, 
    GoldCodeGenerator& gc_generator,
    CodeModulator *localCodeModulator,
    wsgc_complex *faded_source_samples,
    unsigned int nb_faded_source_samples
);

void message_processing_MFSK(
#ifdef _CUDA
    CudaManager& cuda_manager,
#endif
    Options& options,
    wsgc_complex *faded_source_samples,
    unsigned int nb_faded_source_samples
);



void generate_training_prn_list(std::vector<unsigned int>& prn_list, GoldCodeGenerator& gc_generator);
void generate_training_prn_list_unpiloted(std::vector<unsigned int>& prn_list, GoldCodeGenerator& gc_generator);
void generate_message_prn_list(std::vector<unsigned int>& prn_list, GoldCodeGenerator& gc_generator);
void generate_pilot_prn_list(std::vector<unsigned int>& prn_list, GoldCodeGenerator& gc_generator, unsigned int pilot_prni);
void apply_fir(wsgc_complex *inout, unsigned int& nb_samples, const std::vector<wsgc_float>& fir_coef);

//=================================================================================================
int main(int argc, char *argv[])
{
	std::string binary_name(argv[0]);
    Options options(binary_name);

    srand(WsgcUtils::timenow_usec_hour());

    if (options.get_options(argc,argv))
    {
        std::ostringstream os;
        options.print_options(os);
        std::cout << os.str() << std::endl;

        GoldCodeGenerator gc_generator(options.gc_nb_stages, options.nb_message_symbols, options.nb_service_symbols, options.nb_training_symbols, options.g1_poly_powers, options.g2_poly_powers);
        // get and print CUDA mapping
#ifdef _CUDA
        CudaManager cuda_manager(
        		gc_generator.get_nb_message_codes(),
        		options.nb_pilot_prns,
        		gc_generator.get_nb_code_samples(options.f_sampling, options.f_chip),
        		options.batch_size,
        		options.df_steps,
        		options.nb_prns_per_symbol,
        		options.f_step_division
        		);
        if (options.gpu_affinity_specified)
        {
        	cuda_manager.set_gpu_affinity(options.gpu_affinity);
        }
        bool use_cuda_ok = cuda_manager.diagnose();
        options.use_cuda = options.use_cuda && use_cuda_ok;
        std::ostringstream cuda_os;
        cuda_manager.dump(cuda_os);
        std::cout << cuda_os.str() << std::endl << std::endl;
#endif

        CodeModulator *codeModulator = 0;
        CodeModulator *localCodeModulator = 0;
        CodeModulator_MFSK *codeModulator_MFSK = 0;

        unsigned int fft_N = gc_generator.get_nb_code_samples(options.f_sampling, options.f_chip);
        
        // Allocate appropriate objects depending on modulation options
                                      
        // create code modulator                                      
        if (options.modulation.getScheme() == Modulation::Modulation_BPSK) 
        {
            codeModulator = new CodeModulator_BPSK();
            localCodeModulator = new CodeModulator_BPSK();
        }
        else if (options.modulation.getScheme() == Modulation::Modulation_OOK)
        {
            codeModulator = new CodeModulator_OOK();
            //localCodeModulator = new CodeModulator_OOK_detection();
            localCodeModulator = new CodeModulator_OOK();
        }
        else if (options.modulation.getScheme() == Modulation::Modulation_DBPSK)
        {
            codeModulator = new CodeModulator_DBPSK();
            localCodeModulator = new CodeModulator_BPSK();
        }
        else if (options.modulation.getScheme() == Modulation::Modulation_MFSK)
        {
            codeModulator_MFSK = new CodeModulator_MFSK(
            		options.mfsk_options._f_sampling,
            		options.mfsk_options._zero_frequency + options.f_tx,
            		options.mfsk_options._symbol_bandwidth,
            		options.mfsk_options._symbol_time);
        }
        
        wsgc_complex *source_samples = 0;
        unsigned int nb_source_samples = 0;
        SourceMixer *source_mixer = 0;
        SimulatedSource *message_source = 0;
        SimulatedSource *pilot_source = 0;

        if (codeModulator) // if a code modulator is available then the actual signal processing can take place
        {
            // Produce signal samples
            std::cout << "Produce signal samples..." << std::endl;

            std::vector<unsigned int> pilot_prns;
            
            if (options.simulate_training)
            {
                // only one "PRN per symbol" for training sequence
                message_source = new SimulatedSource(gc_generator, options.prns, options.f_sampling, options.f_chip,
                                                     options.f_tx, options.code_shift, 1, 0.0);
            }
            else
            {
                message_source = new SimulatedSource(gc_generator, options.prns, options.f_sampling, options.f_chip,
                                                     options.f_tx, options.code_shift, options.nb_prns_per_symbol, 0.0);
            }

            message_source->set_code_modulator(codeModulator);
            message_source->create_samples();

            if (options.modulation.isCodeDivisionCapable() && options.nb_pilot_prns > 0)
            {
                if (options.simulate_training)
                {
                    pilot_prns.assign(options.prns.size() + options.batch_size, options.pilot2);
                    // only one "PRN per symbol" for training sequence
                    pilot_source = new SimulatedSource(gc_generator, pilot_prns, options.f_sampling, options.f_chip,
                                                       options.f_tx, options.code_shift, 1, 0.0);
                }
                else
                {
                    pilot_prns.assign(options.prns.size() + (options.batch_size/options.nb_prns_per_symbol), options.pilot1); 
            	    pilot_source= new SimulatedSource(gc_generator, pilot_prns, options.f_sampling, options.f_chip,
                                                      options.f_tx, options.code_shift, options.nb_prns_per_symbol, 0.0);
                }

            	pilot_source->set_code_modulator(codeModulator);
            	pilot_source->create_samples();

            	wsgc_float pilot_gain = pow(10.0, (options.pilot_gain_db / 10.0));
            	source_mixer = new SourceMixer(message_source, pilot_source, pilot_gain);

            	source_samples = source_mixer->get_samples();
            	nb_source_samples = source_mixer->get_nb_samples();
            }
            else
            {
            	source_samples = message_source->get_samples();
            	nb_source_samples = message_source->get_nb_samples();
            }
        }
        else if (codeModulator_MFSK)
        {
        	options.prns.push_back(0); // pad with one symbol
        	nb_source_samples = codeModulator_MFSK->get_nb_symbol_samples()*options.prns.size();
        	source_samples = (wsgc_complex *) WSGC_FFTW_MALLOC(nb_source_samples*sizeof(wsgc_fftw_complex));
        	codeModulator_MFSK->modulate(reinterpret_cast<wsgc_fftw_complex*>(source_samples), options.prns);
        }

        if (source_samples)
        {
            // Apply lowpass filter if any
            if (options._fir_coef_generator != 0)
            {
            	std::cout << "Apply Lowpass FIR" << std::endl;
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
                faded_source_samples = source_samples;
                nb_faded_source_samples = nb_source_samples;
            }
            
            // apply AWGN
            if (options.make_noise)
            {
            	std::cout << "Apply AWGN" << std::endl;
                fading->apply_awgn(faded_source_samples, nb_faded_source_samples, options.code_shift, options.snr);
            }

            // demodulate OOK
            if (options.modulation.demodulateBeforeCorrelate())
            {
                Demodulator *demodulator = 0;

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

	            if (demodulator)
	            {
	                delete demodulator;
	            }
            }

            if (codeModulator_MFSK)
            {
				message_processing_MFSK(
#ifdef _CUDA
				cuda_manager,
#endif
				options,
				faded_source_samples,
				nb_faded_source_samples
				);
            }
            else
            {
				// Implement correlator(s)
				if (options.simulate_training)
				{
					training_processing(
	#ifdef _CUDA
						cuda_manager,
	#endif
						options,
						gc_generator,
						localCodeModulator,
						faded_source_samples,
						nb_faded_source_samples
						);
				}
				else
				{
					message_processing(
	#ifdef _CUDA
					cuda_manager,
	#endif
					options,
					gc_generator,
					localCodeModulator,
					faded_source_samples,
					nb_faded_source_samples
					);
				}
            }
            
            if (source_mixer)
            {
                delete source_mixer;
            }

            if (pilot_source)
            {
                delete pilot_source;
            }

            if (message_source)
            {
                delete message_source;
            }
        }
        else
        {
            std::cout << "Code modulator not implemented for this modulation scheme" << std::endl;
        }
          
        if (localCodeModulator)
        {
            delete localCodeModulator;
        }
          
        if (codeModulator)
        {
            delete codeModulator;
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
void message_processing(
#ifdef _CUDA
    CudaManager& cuda_manager,
#endif    
    Options& options, 
    GoldCodeGenerator& gc_generator,
    CodeModulator *localCodeModulator,
    wsgc_complex *faded_source_samples,
    unsigned int nb_faded_source_samples
)
{
    timespec time1, time2;
    int time_option = CLOCK_REALTIME;
    LocalCodes *local_codes = 0;
    LocalCodesFFT *local_codes_fft = 0;
    //MultiplePrnCorrelator *mprn_correlator = 0;
    PilotedMultiplePrnCorrelator *piloted_mprn_correlator = 0;
    PilotCorrelationAnalyzer *pilot_correlation_analyzer = 0;
    const std::map<unsigned int, unsigned int> *prn_shift_occurences = 0;
    PilotCorrelator *pilot_correlator = 0;
    PilotedMessageCorrelator *message_correlator = 0;
    PrnAutocorrelator *prn_autocorrelator = 0;
    DecisionBox *decision_box = 0;
    std::vector<CorrelationRecord> correlation_records;
    std::vector<AutocorrelationRecord> autocorrelation_records;
    //UnpilotedMessageCorrelator *unpiloted_message_correlator = 0;
    //UnpilotedMultiplePrnCorrelator *unpiloted_mprn_correlator = 0;
    UnpilotedMultiplePrnCorrelator *dm_correlator = 0;
    unsigned int fft_N = gc_generator.get_nb_code_samples(options.f_sampling, options.f_chip);
    
    std::vector<unsigned int> message_prn_numbers;
    std::vector<unsigned int> pilot_prn_numbers;
    
    generate_message_prn_list(message_prn_numbers, gc_generator);
    generate_pilot_prn_list(pilot_prn_numbers, gc_generator, 1);

    if (options.modulation.isCodeDivisionCapable() && (options.nb_pilot_prns > 0))
    {
        // Modulation should be frequency dependant. It doesn't make much sense if there is no frequency tracking necessary.
        // This pilot correlation works exclusively in the two dimensions of (frequency shift, PRN time shift).
        // this means it works only with frequency dependant modulations
        // For now this concerns BPSK modulation with pilot use
        if (options.modulation.isFrequencyDependant())
        {
            std::cout << "Creating piloted multiple PRN correlator" << std::endl;
            pilot_correlation_analyzer = new PilotCorrelationAnalyzer(options.analysis_window_size, options.nb_prns_per_symbol, options.nb_samples_per_code);
#ifdef _CUDA
            if (options.use_cuda)
            {
                std::cout << "!!! USING CUDA !!!" << std::endl;
                unsigned int cuda_device = cuda_manager.get_pilot_device();
                local_codes = new LocalCodes_Cuda(*localCodeModulator, gc_generator, options.f_sampling, options.f_chip, message_prn_numbers, cuda_device); // make local codes time domain
                pilot_correlator = new PilotCorrelator_Cuda(gc_generator, *localCodeModulator, options.f_sampling, options.f_chip, pilot_prn_numbers, cuda_device, options.nb_prns_per_symbol, options.df_steps, options.batch_size, options.f_step_division);
                message_correlator = new PilotedMessageCorrelator_Cuda(*((LocalCodes_Cuda *) local_codes), options.f_sampling, options.f_chip, options.nb_prns_per_symbol, cuda_device);
            }
            else
            {
                local_codes = new LocalCodes_Host(*localCodeModulator, gc_generator, options.f_sampling, options.f_chip, message_prn_numbers); // make local codes time domain
                pilot_correlator = new PilotCorrelator_Host(gc_generator, *localCodeModulator, options.f_sampling, options.f_chip, pilot_prn_numbers, options.nb_prns_per_symbol, options.df_steps, options.batch_size, options.f_step_division);
                message_correlator = new PilotedMessageCorrelator_Host(*((LocalCodes_Host *) local_codes), options.f_sampling, options.f_chip, options.nb_prns_per_symbol);
            }
#else
            local_codes = new LocalCodes_Host(*localCodeModulator, gc_generator, options.f_sampling, options.f_chip, message_prn_numbers); // make local codes time domain
            pilot_correlator = new PilotCorrelator_Host(gc_generator, *localCodeModulator, options.f_sampling, options.f_chip, pilot_prn_numbers, options.nb_prns_per_symbol, options.df_steps, options.batch_size, options.f_step_division);
            message_correlator = new PilotedMessageCorrelator_Host(*((LocalCodes_Host *) local_codes), options.f_sampling, options.f_chip, options.nb_prns_per_symbol);
#endif
            prn_autocorrelator = new PrnAutocorrelator_Host(fft_N, options.nb_prns_per_symbol);
            piloted_mprn_correlator = new PilotedMultiplePrnCorrelator(*pilot_correlation_analyzer, correlation_records, *pilot_correlator, *message_correlator, *prn_autocorrelator);
        }
    }
    else
    {
        std::cout << "Unpiloted correlation with frequency independant modulations" << std::endl;

        // only host for now
        local_codes_fft = new LocalCodesFFT_Host(*localCodeModulator, gc_generator, options.f_sampling, options.f_chip, message_prn_numbers); // make local codes time domain
        dm_correlator = new UnpilotedMultiplePrnCorrelator_Host(
        		options.f_sampling,
        		options.f_chip,
        		gc_generator.get_code_length(),
        		options.nb_prns_per_symbol,
        		message_prn_numbers,
                options.batch_size, // will be taken as number of symbols
        		options.analysis_window_size, 
        		correlation_records,
        		*((LocalCodesFFT_Host *) local_codes_fft));
    }

    // Do the correlation
    std::cout << "Do the correlation..." << std::endl;

    wsgc_complex *signal_samples;

    clock_gettime(time_option, &time1);
    
    // Support pipeline processing if necessary (with pilot assisted correlation using CUDA)
    bool input_samples_available = false;           // input samples (length of one PRN) are available for processing  
    bool process_next_output = false;
    
    SampleSequencer sample_sequencer(faded_source_samples, nb_faded_source_samples, fft_N);
    unsigned int prn_i = 0;

    while ((input_samples_available = sample_sequencer.get_next_code_samples(&signal_samples)) || process_next_output) // pseudo real time loop, one PRN length at a time
    {
        if (piloted_mprn_correlator != 0) // Correlation with pilot PRN(s) - pipelined processing
        {
            //std::cout << "Correlation with pilot PRN(s) - new correlator with pilot(s), PRNi=" << piloted_mprn_correlator->get_prn_index() 
            //          << ", CORRi=" << piloted_mprn_correlator->get_correlation_index() << std::endl;
        
            if (input_samples_available) // process input
            {
                // Use PilotedMultiplePrnCorrelator specialized class that encompasses pilot and message correlators and eventually produces correlation records usable directly in the Decision Box
                //   - set source block
                //   - execute the correlation / averaging process
                piloted_mprn_correlator->set_source_block(reinterpret_cast<wsgc_fftw_complex *>(signal_samples));
                piloted_mprn_correlator->make_correlation(0); // use the only pilot (code index 0)
            }
        }
        else if (dm_correlator != 0)
        {
        	if (dm_correlator->set_samples(signal_samples))
        	{
        		dm_correlator->execute();
        	}
        }
        // else just do nothing as there are no other options
    }
    
    clock_gettime(time_option, &time2);
    std::cout << "Message correlation time: " << std::setw(12) << std::setprecision(9) << WsgcUtils::get_time_difference(time2,time1) << " s" << std::endl << std::endl;
    
    /*
    // Dump timing data
    if (pilot_correlation_analyzer != 0)
    {
        std::ostringstream corr_os;

        pilot_correlation_analyzer->dump_timings(corr_os);
        corr_os << std::endl;

        std::cout << corr_os.str() << std::endl;
    }
    */

    // Do the decoding with the decision box
    std::cout << "Do the decoding with the decision box..." << std::endl;

    if (piloted_mprn_correlator != 0) // Correlation with pilot PRN(s)
    {
        std::ostringstream corr_os;

        corr_os << "Piloted correlation" << std::endl << std::endl;
        corr_os << "Pilot correlation analysis results:" << std::endl;
        pilot_correlation_analyzer->dump_histo_time_shift_occurences(corr_os);
        corr_os << std::endl;

#ifdef _CUDA
        if (options.use_cuda)
        {
            pilot_correlation_analyzer->set_pilot_mag_display_factor((fft_N/2)*(fft_N/2));
            pilot_correlation_analyzer->set_message_mag_display_factor(fft_N/2);
        }
#endif

        corr_os << "--- pilot correlation records:" << std::endl;
        pilot_correlation_analyzer->dump_pilot_correlation_records(corr_os);
        corr_os << std::endl << "--- correlation records:" << std::endl;
        pilot_correlation_analyzer->dump_message_correlation_records(corr_os);
        std::cout << corr_os.str() << std::endl;
        
        if (options.modulation.isCodeDivisionCapable() && (options.nb_pilot_prns > 0))
        {
            decision_box = new DecisionBox_Piloted_And_Synced(options.nb_prns_per_symbol, fft_N, options.decision_thresholds, *pilot_correlation_analyzer);
        }
    }
    else if (dm_correlator != 0)
    {
        std::ostringstream corr_os;
        corr_os << "Unpiloted correlation" << std::endl << std::endl;
        corr_os << "Message correlation analysis results:" << std::endl;
        dm_correlator->dump_correlation_records(corr_os, fft_N);
    	corr_os << "Time shift analysis results:" << std::endl;
    	dm_correlator->dump_time_analyzer_results(corr_os);
        corr_os << std::endl;
        std::cout << corr_os.str() << std::endl;
        
        decision_box = new DecisionBox_Unpiloted_And_Synced(options.nb_prns_per_symbol, fft_N, options.decision_thresholds, correlation_records);
        decision_box->set_mag_display_adj_factor(fft_N);
    }

    if (decision_box)
    {
#ifdef _CUDA
        if (options.use_cuda)
        {
        	decision_box->set_mag_display_adj_factor(fft_N/2);
        	decision_box->set_use_cuda(true);
            if (!options.decision_thresholds_specified)
            {
                options.decision_thresholds.set_cuda_defaults();
            }
        }
#endif
        decision_box->analyze_records();

        if (decision_box->is_prni_at_max_invalid())
        {
            std::cout << "Symbol boundaries cannot be estimated with confidence - message is discarded" << std::endl;
        }
        else
        {
            decision_box->estimate_symbols();

            std::vector<unsigned int> symbol_indices;
            for (unsigned int i=0; i<options.prns.size(); i++)
            {
                symbol_indices.push_back(i);
            }
            
            std::ostringstream os_result;
            os_result << std::endl;
            options.decision_thresholds.print_options(os_result);
            os_result << std::endl << "Decision box records:" << std::endl;
            decision_box->dump_decision_records(os_result);
            os_result << std::endl << "Decisions status:" << std::endl;
            decision_box->dump_decision_status(os_result, options.prns);
            os_result << std::endl << "Index, original and decoded symbols (-1 denotes an erasure):";
            os_result << std::endl << "-----------------------------------------------------------" << std::endl;
            os_result << "_SIN "; print_vector<unsigned int, unsigned int>(symbol_indices, 4, os_result); os_result << std::endl;
            os_result << "_SOR "; print_vector<unsigned int, unsigned int>(options.prns, 4, os_result); os_result << std::endl;
            os_result << "_SRS "; print_vector<int, int>(decision_box->get_decoded_symbols(), 4, os_result); os_result << std::endl;
            std::cout << os_result.str() << std::endl;

            unsigned int erasure_counts = 0;
            unsigned int error_counts = 0;

            for (std::vector<int>::const_iterator it = decision_box->get_decoded_symbols().begin(); it != decision_box->get_decoded_symbols().end(); ++it)
            {
                if (*it < 0)
                {
                    erasure_counts++;
                }
                else
                {
                    if (*it != options.prns[it-decision_box->get_decoded_symbols().begin()])
                    {
                        error_counts++;
                    }
                }
            }

            std::cout << erasure_counts << " erasures (" << ((float) erasure_counts)/decision_box->get_decoded_symbols().size() << ")" << std::endl;
            std::cout << error_counts << " errors (" << ((float) error_counts)/decision_box->get_decoded_symbols().size() << ")" << std::endl;
            std::cout << erasure_counts+error_counts << " total (" << ((float) erasure_counts+error_counts)/decision_box->get_decoded_symbols().size() << ")" << std::endl;
            std::cout << "_SUM " << erasure_counts << "," << error_counts << "," << erasure_counts+error_counts << "," << erasure_counts+2*error_counts << std::endl;
        }

        // Delete dynamically allocated objects

        delete decision_box;
    }

    if (piloted_mprn_correlator)
    {
        delete piloted_mprn_correlator;
    }
    
    if (prn_autocorrelator)
    {
        delete prn_autocorrelator;
    }

    if (message_correlator)
    {
        delete message_correlator;
    }

    if (pilot_correlator)
    {
        delete pilot_correlator;
    }

    if (local_codes_fft)
    {
        delete local_codes_fft;
    }

    if (local_codes)
    {
        delete local_codes;
    }
    
    if (pilot_correlation_analyzer)
    {
        delete pilot_correlation_analyzer;
    }

}


//=================================================================================================
void training_processing(
#ifdef _CUDA
    CudaManager& cuda_manager,
#endif    
    Options& options, 
    GoldCodeGenerator& gc_generator,
    CodeModulator *localCodeModulator,
    wsgc_complex *faded_source_samples,
    unsigned int nb_faded_source_samples
)
{
    timespec time1, time2;
    int time_option = CLOCK_REALTIME;
    LocalCodes *local_codes = 0;
    PilotedTrainingMultiplePrnCorrelator *piloted_tr_mprn_correlator = 0;
    PilotCorrelationAnalyzer *pilot_correlation_analyzer = 0;
    const std::map<unsigned int, unsigned int> *prn_shift_occurences = 0;
    PilotCorrelator *pilot_correlator = 0;
    PilotedTrainingMessageCorrelator *tr_message_correlator = 0;
    unsigned int fft_N = gc_generator.get_nb_code_samples(options.f_sampling, options.f_chip);
    
    
    
    // - Modulation should support code division
    // - Modulation processing should be frequency dependant
    // - There should be at least 2 pilot PRN(s) - using second
    if (options.modulation.isCodeDivisionCapable() && options.modulation.isFrequencyDependant() && (options.nb_pilot_prns > 1))
    {
        std::vector<unsigned int> pilot_prn_numbers;
        generate_pilot_prn_list(pilot_prn_numbers, gc_generator, 2); // Pilot PRN #2 is used

        std::vector<unsigned int> training_prn_numbers;
        generate_training_prn_list(training_prn_numbers, gc_generator);

        std::cout << "Creating managing objects..." << std::endl;
        pilot_correlation_analyzer = new PilotCorrelationAnalyzer(options.analysis_window_size, options.nb_prns_per_symbol, options.nb_samples_per_code);    
#ifdef _CUDA
        if (options.use_cuda)
        {
            pilot_correlation_analyzer->set_pilot_mag_display_factor((fft_N/2)*(fft_N/2));
            pilot_correlation_analyzer->set_training_mag_display_factor(fft_N/8);
            unsigned int cuda_device = cuda_manager.get_pilot_device();
            local_codes = new LocalCodes_Cuda(*localCodeModulator, gc_generator, options.f_sampling, options.f_chip, training_prn_numbers, cuda_device); // make local codes time domain
            pilot_correlator = new PilotCorrelator_Cuda(gc_generator, *localCodeModulator, options.f_sampling, options.f_chip, pilot_prn_numbers, cuda_device, options.nb_prns_per_symbol, options.df_steps, options.batch_size, options.f_step_division);
        	tr_message_correlator = new PilotedTrainingMessageCorrelator_Cuda(*((LocalCodes_Cuda *) local_codes), options.f_sampling, options.f_chip, pilot_correlation_analyzer->get_analysis_window_size_in_prns(), options.nb_random_prns, cuda_device);
        }
        else
        {
            local_codes = new LocalCodes_Host(*localCodeModulator, gc_generator, options.f_sampling, options.f_chip, training_prn_numbers); // make local codes time domain
            pilot_correlator = new PilotCorrelator_Host(gc_generator, *localCodeModulator, options.f_sampling, options.f_chip, pilot_prn_numbers, options.nb_prns_per_symbol, options.df_steps, options.batch_size, options.f_step_division);
        	tr_message_correlator = new PilotedTrainingMessageCorrelator_Host(*((LocalCodes_Host *) local_codes), options.f_sampling, options.f_chip, options.nb_random_prns);
        }
#else
        local_codes = new LocalCodes_Host(*localCodeModulator, gc_generator, options.f_sampling, options.f_chip, training_prn_numbers); // make local codes time domain
        pilot_correlator = new PilotCorrelator_Host(gc_generator, *localCodeModulator, options.f_sampling, options.f_chip, pilot_prn_numbers, options.nb_prns_per_symbol, options.df_steps, options.batch_size, options.f_step_division);
        tr_message_correlator = new PilotedTrainingMessageCorrelator_Host(*((LocalCodes_Host *) local_codes), options.f_sampling, options.f_chip, options.nb_random_prns);
#endif
        piloted_tr_mprn_correlator = new PilotedTrainingMultiplePrnCorrelator(*pilot_correlation_analyzer, *pilot_correlator, *tr_message_correlator);
        
        // Do the correlation
        std::cout << "Do the correlation..." << std::endl;
        wsgc_complex *signal_samples;
        SampleSequencer sample_sequencer(faded_source_samples, nb_faded_source_samples, fft_N);
        unsigned int prn_i = 0;

        clock_gettime(time_option, &time1);
        
        while (sample_sequencer.get_next_code_samples(&signal_samples)) // pseudo real time loop, one PRN length at a time
        {
            if (piloted_tr_mprn_correlator != 0) // Correlation with pilot PRN(s) - pipelined processing
            {
                piloted_tr_mprn_correlator->set_source_block(reinterpret_cast<wsgc_fftw_complex *>(signal_samples));
                piloted_tr_mprn_correlator->make_correlation(0); // use the only pilot (code index 0)
            }
        }
        
        clock_gettime(time_option, &time2);
        std::cout << "Training sequence correlation time: " << std::setw(12) << std::setprecision(9) << WsgcUtils::get_time_difference(time2,time1) << " s" << std::endl << std::endl;

        std::ostringstream corr_os;

        corr_os << "Pilot correlation analysis results:" << std::endl;
        pilot_correlation_analyzer->dump_histo_time_shift_occurences(corr_os);
        corr_os << std::endl;

        corr_os << "--- pilot correlation records:" << std::endl;
        pilot_correlation_analyzer->dump_pilot_correlation_records(corr_os);

        corr_os << "--- training correlation records:" << std::endl;
        pilot_correlation_analyzer->dump_training_correlation_records(corr_os);

        std::cout << corr_os.str() << std::endl;
        
        // TODO: conclude on the message epoch
    }    
    else if (options.modulation.demodulateBeforeCorrelate())
    {
    	LocalCodesFFT *local_codes_fft = 0;
    	UnpilotedTrainingMessageCorrelator *unpiloted_training_message_correlator = 0;

        std::vector<unsigned int> training_prn_numbers;
        generate_training_prn_list_unpiloted(training_prn_numbers, gc_generator);

    	unsigned int prn_length = gc_generator.get_code_length();
    	std::vector<TrainingCorrelationRecord> training_correlation_records;

    	local_codes_fft = new LocalCodesFFT_Host(
    			*localCodeModulator,
    			gc_generator,
    			options.f_sampling,
    			options.f_chip,
    			training_prn_numbers);

    	unpiloted_training_message_correlator =	new UnpilotedTrainingMessageCorrelator_Host(
				options.f_sampling,
				options.f_chip,
				prn_length,
				options.prns.size(),
				options.analysis_window_size,
				training_prn_numbers,
				training_correlation_records,
				*((LocalCodesFFT_Host *)local_codes_fft));

    	SampleSequencer sample_sequencer(faded_source_samples, nb_faded_source_samples, fft_N);
    	wsgc_complex *signal_samples;

        clock_gettime(time_option, &time1);

        while (sample_sequencer.get_next_code_samples(&signal_samples)) // pseudo real time loop, one PRN length at a time
        {
            if (unpiloted_training_message_correlator != 0) // Correlation without pilot PRN(s) - pipelined processing
            {
            	unpiloted_training_message_correlator->set_samples(signal_samples);
            	unpiloted_training_message_correlator->execute();
            }
        }

        clock_gettime(time_option, &time2);
        std::cout << "Training sequence correlation time: " << std::setw(12) << std::setprecision(9) << WsgcUtils::get_time_difference(time2,time1) << " s" << std::endl << std::endl;

        std::ostringstream corr_os;
        corr_os << "--- training correlation records:" << std::endl;
        TrainingCorrelationRecord::dump_banner(corr_os);

        for (std::vector<TrainingCorrelationRecord>::const_iterator it = training_correlation_records.begin(); it != training_correlation_records.end(); ++it)
        {
        	it->dump_line(corr_os);
        }

        std::cout << corr_os.str() << std::endl;
    }
    else
    {
        std::cout << "Synchronization training is not implemented under these conditions" << std::endl;
    }
    
    if (piloted_tr_mprn_correlator)
    {
        delete piloted_tr_mprn_correlator;
    }
    
    if (tr_message_correlator)
    {
        delete tr_message_correlator;
    }
    
    if (pilot_correlator)
    {
        delete pilot_correlator;
    }
    
    if (local_codes)
    {
        delete local_codes;
    }
    
    if (pilot_correlation_analyzer)
    {
        delete pilot_correlation_analyzer;
    }
}


//=================================================================================================
void message_processing_MFSK(
#ifdef _CUDA
    CudaManager& cuda_manager,
#endif
    Options& options,
    wsgc_complex *faded_source_samples,
    unsigned int nb_faded_source_samples
)
{
    timespec time1, time2;
    int time_option = CLOCK_REALTIME;
    MFSK_MessageDemodulator *mfsk_message_demodulator;

    mfsk_message_demodulator = new MFSK_MessageDemodulator_Host(
    		options.mfsk_options._fft_N,
    		options.mfsk_options._nb_fft_per_symbol,
    		options.mfsk_options._zero_fft_slot,
    		options.nb_message_symbols,
    		options.nb_service_symbols);
    SampleSequencer sample_sequencer(faded_source_samples, nb_faded_source_samples, options.mfsk_options._fft_N*options.mfsk_options._nb_fft_per_symbol);

	wsgc_complex *signal_samples;

    clock_gettime(time_option, &time1);

    while (sample_sequencer.get_next_code_samples(&signal_samples)) // pseudo real time loop, one PRN length at a time
    {
    	mfsk_message_demodulator->execute(signal_samples, options._rssoft_engine);
    }

    clock_gettime(time_option, &time2);
    std::cout << "MFSK message sequence decoding time: " << std::setw(12) << std::setprecision(9) << WsgcUtils::get_time_difference(time2,time1) << " s" << std::endl << std::endl;
    
#ifdef _RSSOFT    
    if (options._rssoft_engine)
    {
        RSSoft_DecisionBox rssoft_decision_box(*options._rssoft_engine, options._source_codec);
        
        switch (options.rs_decoding_mode)
        {
            case Options::RSSoft_decoding_full:
            case Options::RSSoft_decoding_best:
            case Options::RSSoft_decoding_first:
                rssoft_decision_box.run(options.rs_decoding_mode);
                break;
            case Options::RSSoft_decoding_regex:
                rssoft_decision_box.run_regex(options.rs_decoding_regex);
                break;
            default:
                std::cout << "Unknown RSSoft decoding options" << std::endl;
        }
    }
    else
    {
#endif
        std::ostringstream demod_os;
        demod_os << "--- demodulation records:" << std::endl;
        mfsk_message_demodulator->dump_demodulation_records(demod_os);
        std::cout << demod_os.str() << std::endl;

        DecisionBox_MFSK decision_box(options.mfsk_options._fft_N, options.mfsk_options._nb_fft_per_symbol, options.decision_thresholds);
        decision_box.estimate_symbols(mfsk_message_demodulator->get_demodulation_records());

        options.prns.pop_back(); // pop padding symbol

        std::vector<unsigned int> symbol_indices;
        for (unsigned int i=0; i<options.prns.size(); i++)
        {
            symbol_indices.push_back(i);
        }

        std::ostringstream os_result;
        options.decision_thresholds.print_options(os_result);
        os_result << std::endl << "Decisions status:" << std::endl;
        decision_box.dump_decision_status(os_result, options.prns, mfsk_message_demodulator->get_demodulation_records());
        os_result << std::endl << "Index, original and decoded symbols (-1 denotes an erasure):";
        os_result << std::endl << "-----------------------------------------------------------" << std::endl;
        os_result << "_SIN "; print_vector<unsigned int, unsigned int>(symbol_indices, 4, os_result); os_result << std::endl;
        os_result << "_SOR "; print_vector<unsigned int, unsigned int>(options.prns, 4, os_result); os_result << std::endl;
        os_result << "_SRS "; print_vector<int, int>(decision_box.get_decoded_symbols(), 4, os_result); os_result << std::endl;
        std::cout << os_result.str() << std::endl;

        unsigned int erasure_counts = 0;
        unsigned int error_counts = 0;

        for (std::vector<int>::const_iterator it = decision_box.get_decoded_symbols().begin(); it != decision_box.get_decoded_symbols().end(); ++it)
        {
            if (*it < 0)
            {
                erasure_counts++;
            }
            else
            {
                if (*it != options.prns[it-decision_box.get_decoded_symbols().begin()])
                {
                    error_counts++;
                }
            }
        }

        std::cout << erasure_counts << " erasures (" << ((float) erasure_counts)/decision_box.get_decoded_symbols().size() << ")" << std::endl;
        std::cout << error_counts << " errors (" << ((float) error_counts)/decision_box.get_decoded_symbols().size() << ")" << std::endl;
        std::cout << erasure_counts+error_counts << " total (" << ((float) erasure_counts+error_counts)/decision_box.get_decoded_symbols().size() << ")" << std::endl;
        std::cout << "_SUM " << erasure_counts << "," << error_counts << "," << erasure_counts+error_counts << "," << erasure_counts+2*error_counts << std::endl;
#ifdef _RSSOFT        
    }
#endif
}


//=================================================================================================
void generate_training_prn_list(std::vector<unsigned int>& prn_list, GoldCodeGenerator& gc_generator)
{
	for (unsigned int prni=0; prni < gc_generator.get_nb_message_codes(); prni++) 
	{
		prn_list.push_back(prni);
	}
}

//=================================================================================================
void generate_training_prn_list_unpiloted(std::vector<unsigned int>& prn_list, GoldCodeGenerator& gc_generator)
{
	unsigned int start_code = gc_generator.get_nb_message_codes() + gc_generator.get_nb_service_codes();

	for (unsigned int prni=start_code; prni < start_code + gc_generator.get_nb_training_codes(); prni++)
	{
		prn_list.push_back(prni);
	}
}

//=================================================================================================
void generate_message_prn_list(std::vector<unsigned int>& prn_list, GoldCodeGenerator& gc_generator)
{
	for (unsigned int prni=0; prni <= gc_generator.get_nb_message_codes(); prni++) // incudes the noise PRN that is the last hence the <= comparison
	{
		prn_list.push_back(prni);
	}
}


//=================================================================================================
void generate_pilot_prn_list(std::vector<unsigned int>& prn_list, GoldCodeGenerator& gc_generator, unsigned int pilot_prni)
{
    prn_list.push_back(gc_generator.get_nb_message_codes()+pilot_prni);
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
