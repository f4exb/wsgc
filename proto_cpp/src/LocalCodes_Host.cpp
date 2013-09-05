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

     Creates a local copy of all symbols gold codes in the time domain
     Used for plain time domain correlation

*/

#include "LocalCodes_Host.h"
#include "GoldCodeGenerator.h"
#include "CodeModulator.h"
#include <string.h>
#include <assert.h>
#include <fftw3.h>


LocalCodes_Host::LocalCodes_Host(
            CodeModulator& code_modulator, 
            const GoldCodeGenerator& gc_generator, 
            wsgc_float f_sampling, 
            wsgc_float f_chip,
            std::vector<unsigned int>& symbols) :
    LocalCodes(code_modulator, gc_generator, f_sampling, f_chip, symbols)
{
	fill_codes_matrix();
}


LocalCodes_Host::~LocalCodes_Host()
{
    std::vector<wsgc_complex*>::iterator codes_it = _codes_matrix.begin();
    const std::vector<wsgc_complex*>::iterator codes_end =  _codes_matrix.end();

    for (; codes_it != codes_end; ++codes_it)
    {
        WSGC_FFTW_FREE(*codes_it);
    }
}


void LocalCodes_Host::fill_codes_matrix()
{
    std::vector<char> code;

    std::vector<unsigned int>::iterator prni_it = _symbols.begin();
    const std::vector<unsigned int>::iterator prni_end = _symbols.end();
    
    for (; prni_it != prni_end; ++prni_it)
    {
        wsgc_complex *code_array = (wsgc_complex *) WSGC_FFTW_MALLOC(_nb_code_samples*sizeof(wsgc_fftw_complex));
        wsgc_fftw_complex *code_in = reinterpret_cast<wsgc_fftw_complex *>(code_array);
        _gc_generator.make_code(code, *prni_it); // 0/1 bits
        _code_modulator.fill_code_samples(code_in, code, _f_sampling, _f_chip); // This is the modulation specific part
        _codes_matrix.push_back(code_array);
    }
}


const wsgc_complex *LocalCodes_Host::get_local_code(unsigned int prni) const
{
    assert(prni < _symbols.size());
    return _codes_matrix[prni];
}


