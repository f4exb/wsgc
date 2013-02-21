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

     MFSK_MessageDemodulator

     Class to incoherently demodulate MFSK message. This is not using a correlation scheme
     and is there only for comparison to incoherent MFSK

     This is the Host implementation

*/
#include "MFSK_MessageDemodulator_Host.h"
#include "MFSK_MessageDemodulationRecord.h"
#include "WsgcUtils.h"
#include <cstring>

//=================================================================================================
MFSK_MessageDemodulator_Host::MFSK_MessageDemodulator_Host(
			unsigned int fft_N,
			unsigned int nb_fft_per_symbol,
			int zero_frequency_slot,
			unsigned int nb_message_symbols,
			unsigned int nb_service_symbols) :
    MFSK_MessageDemodulator::MFSK_MessageDemodulator(fft_N, nb_fft_per_symbol, zero_frequency_slot, nb_message_symbols, nb_service_symbols),
    _symbol_i(0),
    _fft_i(0)
{
    _samples = (wsgc_complex *) WSGC_FFTW_MALLOC(_nb_fft_per_symbol*_fft_N*sizeof(wsgc_fftw_complex));
    _src_fft = (wsgc_complex *) WSGC_FFTW_MALLOC(_nb_fft_per_symbol*_fft_N*sizeof(wsgc_fftw_complex));
    _magsum_fft = new wsgc_float[_fft_N];
    
    // Do the multiple FFT of input samples for the length of one symbol
    int N = _fft_N;
    int howmany = _nb_fft_per_symbol;
    _fft_plan = WSGC_FFTW_PLAN_MANY(1, &N, howmany,
    		reinterpret_cast<wsgc_fftw_complex *>(_samples),
    		0, 1, _fft_N,
    		reinterpret_cast<wsgc_fftw_complex *>(_src_fft),
    		0, 1, _fft_N,
    		FFTW_FORWARD, FFTW_ESTIMATE);
}


//=================================================================================================
MFSK_MessageDemodulator_Host::~MFSK_MessageDemodulator_Host()
{
    delete[] _magsum_fft;
    WSGC_FFTW_FREE(_src_fft);
    WSGC_FFTW_FREE(_samples);
}


//=================================================================================================
void MFSK_MessageDemodulator_Host::execute(wsgc_complex *symbol_samples)
{
    memcpy(_samples, symbol_samples, _nb_fft_per_symbol*_fft_N*sizeof(wsgc_fftw_complex));
    WSGC_FFTW_EXECUTE(_fft_plan);
    clean_magsum();
    
    for (unsigned int fft_i = 0; fft_i < _nb_fft_per_symbol; fft_i++)
    {
        cumulate_magsum(fft_i);
    }

    estimate_magpeak();
    _symbol_i++;
}


//=================================================================================================
void MFSK_MessageDemodulator_Host::clean_magsum()
{
    for (unsigned int ffti = 0; ffti < _fft_N; ffti++)
    {
        _magsum_fft[ffti] = 0.0;
    }
}


//=================================================================================================
void MFSK_MessageDemodulator_Host::cumulate_magsum(unsigned int fft_index)
{
    unsigned int mffti = 0;
    wsgc_float magnitude;
    _fft_i = 0;

    for (unsigned int ffti = fft_index*_fft_N; mffti < _fft_N; ffti++, mffti++, _fft_i++)
    {
        WsgcUtils::magnitude_estimation(&_src_fft[ffti], &magnitude);
        _magsum_fft[mffti] += magnitude;
    }
}


//=================================================================================================
void MFSK_MessageDemodulator_Host::estimate_magpeak()
{
	static const MFSK_MessageDemodulationRecord tmp_demodulation_record;
    wsgc_float max_mag = 0.0;
    wsgc_float mag_sum = 0.0;
    wsgc_float nse_mag;
    unsigned int max_ffti;

    for (unsigned int si = 0; si < _nb_message_symbols+_nb_service_symbols; si++)
    {
    	unsigned int mffti = get_fft_slot(si);
        mag_sum += _magsum_fft[mffti];

    	if (_magsum_fft[mffti] > max_mag)
        {
            max_mag = _magsum_fft[mffti];
            max_ffti = mffti;
        }

    	// noise symbol is second service symbol
    	if ((si > _nb_message_symbols) && (si - _nb_message_symbols == 1))
    	{
    		nse_mag = _magsum_fft[mffti];
    	}
    }

    _demodulation_records.push_back(tmp_demodulation_record);
    MFSK_MessageDemodulationRecord& demodulation_record = _demodulation_records.back();

    demodulation_record._symbol_index = _symbol_i;
    demodulation_record._fft_index = max_ffti;
    demodulation_record._symbol_ordinal = get_symbol_ordinal(max_ffti);
    demodulation_record._max_magnitude = max_mag;
    demodulation_record._noise_magnitude = nse_mag;
    demodulation_record._avg_magnitude = mag_sum / (_nb_message_symbols+_nb_service_symbols);
}

