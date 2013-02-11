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
      
     LocalCodes
      
     Creates a local copy of the complex conjugate of the FFT transform of all symbols gold codes
     This pre-calculates the code to be inserted in the final IFFT.
     Used with BPSK complex signals.
     
*/
#include "LocalCodesFFT_Host.h"
#include "GoldCodeGenerator.h"
#include "CodeModulator.h"
#include "WsgcException.h"
#include <string.h>
#include <assert.h>
#include <iostream>

LocalCodesFFT_Host::LocalCodesFFT_Host(
		CodeModulator& code_modulator,
		GoldCodeGenerator& gc_generator,
		wsgc_float f_sampling,
		wsgc_float f_chip,
		std::vector<unsigned int>& symbols) :
    LocalCodesFFT(code_modulator, gc_generator, f_sampling, f_chip, symbols)
{
    fill_codes_matrix();
}

LocalCodesFFT_Host::~LocalCodesFFT_Host()
{
	std::vector<wsgc_complex*>::iterator codes_it = _codes_matrix.begin();
	std::vector<wsgc_complex*>::iterator codes_end =  _codes_matrix.end();

	for (; codes_it != codes_end; ++codes_it)
	{
		delete[] *codes_it;
	}
}

void LocalCodesFFT_Host::fill_codes_matrix()
{
    _fft_code_in = (wsgc_complex *) WSGC_FFTW_MALLOC(_nb_code_samples*sizeof(wsgc_fftw_complex));
    _fft_code_out = (wsgc_complex *) WSGC_FFTW_MALLOC(_nb_code_samples*sizeof(wsgc_fftw_complex));
    _fft_code_plan = WSGC_FFTW_PLAN(_nb_code_samples, reinterpret_cast<wsgc_fftw_complex *>(_fft_code_in), reinterpret_cast<wsgc_fftw_complex *>(_fft_code_out), FFTW_FORWARD, FFTW_ESTIMATE);

    wsgc_fftw_complex *fftw_code_in = reinterpret_cast<wsgc_fftw_complex *>(_fft_code_in);
    std::vector<char> code;
    
    std::vector<unsigned int>::iterator prni_it = _symbols.begin();
    const std::vector<unsigned int>::iterator prni_end = _symbols.end();
    unsigned int symbol_index = 0;
    
    for (; prni_it != prni_end; ++prni_it, symbol_index++)
    {
        assert(*prni_it < _gc_generator.get_nb_codes());
    
        _gc_generator.make_code(code, *prni_it); // 0/1 bits
        _code_modulator.fill_code_samples(fftw_code_in, code, _f_sampling, _f_chip); // This is the modulation specific part
        
        WSGC_FFTW_EXECUTE(_fft_code_plan);
        
        wsgc_complex *code_array = new wsgc_complex[_nb_code_samples];
        
        for (unsigned int i=0; i<_nb_code_samples; i++) // copy and do the conjugate
        {
        	code_array[i].real() = _fft_code_out[i].real();
        	code_array[i].imag() = -_fft_code_out[i].imag();
        }
        
        _codes_matrix.push_back(code_array);
        //std::cout << "Symbol " << *prni_it << std::endl;
        index_symbol(symbol_index, *prni_it);
    }

	WSGC_FFTW_DESTROY_PLAN(_fft_code_plan);
	WSGC_FFTW_FREE(reinterpret_cast<wsgc_fftw_complex *>(_fft_code_in));
	WSGC_FFTW_FREE(reinterpret_cast<wsgc_fftw_complex *>(_fft_code_out));
}
        
const wsgc_complex *LocalCodesFFT_Host::get_local_code(unsigned int symbol) const
{
	std::map<unsigned int,unsigned int>::const_iterator it = _symbols_index_dictionnary.find(symbol);

	if (it == _symbols_index_dictionnary.end())
	{
		std::ostringstream eos;
		eos << "Symbol " << symbol << " not found in dictionnary" << std::endl;
		throw WsgcException(eos.str());
	}
	else
	{
		return _codes_matrix[it->second];
	}
}

const wsgc_complex *LocalCodesFFT_Host::get_local_code_by_index(unsigned int prni) const
{
	return _codes_matrix[prni];
}
