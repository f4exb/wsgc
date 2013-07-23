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
#ifndef __OPTIONS_H__
#define __OPTIONS_H__

#include "WsgcTypes.h"
#include "Modulation.h"
#include "MFSK_Options.h"
#include "DecisionBox_Thresholds.h"
#include <vector>
#include <string>
#include <sstream>

class FadingModel;
class FIRCoefGenerator;
class SourceCodec;
class RSSoft_Engine;

class Options
{
public:
	typedef enum
	{
		RSSoft_decoding_full,
		RSSoft_decoding_best,
		RSSoft_decoding_first,
		RSSoft_decoding_regex
	} RSSoft_decoding_mode;

	Options(std::string& _binary_name);
	~Options();

	bool get_options(int argc, char *argv[]);
	void print_options(std::ostringstream& os);

	bool has_fading() const
	{
		return _fading_model != 0;
	}

	FadingModel *get_fading_model() const
	{
		return _fading_model;
	}

	std::string binary_path;
	std::string binary_name;
	wsgc_float f_sampling;
	wsgc_float f_chip;
	wsgc_float snr;
	bool make_noise;
	wsgc_float f_tx;
	unsigned int code_shift;
	std::vector<unsigned int> prns;
	bool noise_test;
	unsigned int nb_prns_per_symbol;
	unsigned int prn_shift;
	wsgc_float f_init_rx;
	unsigned int nb_random_prns;
	unsigned int tracking_phase_average_cycles;
	unsigned int gc_nb_stages;
	unsigned int nb_message_symbols;
	unsigned int nb_service_symbols;
	unsigned int nb_training_symbols;
	std::vector<unsigned int> g1_poly_powers;
	std::vector<unsigned int> g2_poly_powers;
	unsigned int nb_samples_per_code;
	bool file_debugging;
	FadingModel *_fading_model;
	Modulation modulation;
	unsigned int nb_pilot_prns;
	wsgc_float pilot_gain_db;
	unsigned int pilot1;
	unsigned int pilot2;
	unsigned int df_steps;
	unsigned int f_step_division;
	unsigned int batch_size;
	unsigned int noise_prn;
	bool use_cuda;
	unsigned int analysis_window_size;
	std::string samples_output_file;
	bool simulate_training;
	bool simulate_demod;
	unsigned int gpu_affinity;
	bool gpu_affinity_specified;
	FIRCoefGenerator *_fir_coef_generator;
	MFSK_Options mfsk_options;
	DecisionBox_Thresholds decision_thresholds;
	bool decision_thresholds_specified;
	SourceCodec *_source_codec;
	std::string _source_codec_type_str;
	std::string _source_message_str;
	RSSoft_Engine *_rssoft_engine;
	unsigned int rs_logq;
	unsigned int rs_k;
	unsigned int rs_M;
	unsigned int rs_r;
	RSSoft_decoding_mode rs_decoding_mode;
	std::string rs_decoding_regex;

private:
	int _indicator_int;

	void get_help(std::ostringstream& os);
	void print_fading_model_data(std::ostringstream& os);
	void print_modulation_data(std::ostringstream& os);
	void print_source_codec_data(std::ostringstream& os);
	bool parse_fading_model_data(std::string fading_data_str);
	bool parse_modulation_data(std::string modulation_data_str);
	bool parse_pilot_prns_data(std::string pilot_data_str);
	bool parse_fir_filter_model_data(std::string fir_data_str);
	bool parse_source_coding_data(std::string coding_data_str);
	void print_standard_options(std::ostringstream& os);
	void print_mfsk_options(std::ostringstream& os);
	bool adjust_parameters_for_source_coding();
	bool source_codec_create_message_prns();
#ifdef _RSSOFT
	bool parse_reed_solomon_data(std::string coding_data_str);
	void print_reed_solomon_data(std::ostringstream& os);
	bool encode_reed_solomon();
#endif

};

#endif // __OPTIONS_H__
