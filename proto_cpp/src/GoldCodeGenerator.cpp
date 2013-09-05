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
      
     Gold Code generator
      
     Creates a Gold Code sequences used in the symbol alphabet of the system. 
     It can be in binary (0,1) or (1,-1) format and optionally sampled at given frequency
     
*/
#include <assert.h>
#include "GoldCodeGenerator.h"
#include "WsgcUtils.h" 

#define GC_PRN_SHIFT 3

GoldCodeGenerator::GoldCodeGenerator(
		unsigned int nb_stages,
		unsigned int nb_message_symbols,
		unsigned int nb_service_symbols,
		unsigned int nb_training_symbols,
		std::vector<unsigned int>& g1_poly,
		std::vector<unsigned int>& g2_poly) :
    _nb_stages(nb_stages),
    _nb_message_symbols(nb_message_symbols),
    _nb_service_symbols(nb_service_symbols),
    _nb_training_symbols(nb_training_symbols),
    _g1_poly(g1_poly),
    _g2_poly(g2_poly)
{
	assert (nb_message_symbols+nb_service_symbols+nb_training_symbols < (1<<nb_stages));

    lfsr_generator(_g1_regs, _g1_poly);
    lfsr_generator(_g2_regs, _g2_poly);
}

GoldCodeGenerator::~GoldCodeGenerator()
{}

void GoldCodeGenerator::lfsr_generator(std::vector<unsigned int>& regs, std::vector<unsigned int>& polynomial_powers)
{
    unsigned int lfsr = (1<<_nb_stages)-1; // initial value of LFSR register (all ones)
    unsigned int feedback;                 // feedback bit (when masked by 1)

    std::vector<unsigned int>::const_iterator pow_it;
    const std::vector<unsigned int>::const_iterator pow_end = polynomial_powers.end();
    
    do
    {
        regs.push_back(lfsr);
        feedback = lfsr;

        for (pow_it = polynomial_powers.begin(); pow_it < pow_end; ++pow_it)
        {
            feedback ^= (lfsr >> (_nb_stages - *pow_it));
        }
        
        feedback &= 1; // feedback bit
        lfsr = (lfsr >> 1) | (feedback << (_nb_stages-1));
        
    } while(lfsr != (1<<_nb_stages)-1);
}

void GoldCodeGenerator::make_code(std::vector<char>& code, unsigned int prn_number, wsgc_float f_sampling, wsgc_float f_chip) const
{
    assert(prn_number+GC_PRN_SHIFT < 1<<_nb_stages);
    
    unsigned int number_of_samples;
    wsgc_float pseudo_index_increment;

    if (f_sampling > 0.0)
    {
        assert(!(f_sampling < f_chip));
        number_of_samples = get_nb_code_samples(f_sampling, f_chip);
        pseudo_index_increment = ((wsgc_float) ((1<<_nb_stages)-2)) / (number_of_samples-1); // (code_length-1) / (nb_samples-1) magic formula
    }
    else
    {
        number_of_samples = (1<<_nb_stages)-1;
        pseudo_index_increment = 1.0;
    }

    wsgc_float pseudo_index = 0.0;
    unsigned int index = 0;
    code.clear();
    
    for (int i=0; i<number_of_samples; i++)
    {
        index = int(pseudo_index);
        unsigned int g2_shift_i = (index + ((1<<_nb_stages)-1) - (prn_number+GC_PRN_SHIFT)) % _g2_regs.size();
        unsigned int bit = (_g1_regs[index] >> (_nb_stages-1)) ^ (_g2_regs[g2_shift_i] >> (_nb_stages-1));
        code.push_back(bit);
        
        pseudo_index += pseudo_index_increment;
    }
}


void GoldCodeGenerator::print_g1_regs(std::ostringstream& os)
{
    os << "G1r = ";
    print_vector<unsigned int, unsigned int>(_g1_regs, 3, os);
}


void GoldCodeGenerator::print_g2_regs(std::ostringstream& os)
{
    os << "G2r = ";
    print_vector<unsigned int, unsigned int>(_g2_regs, 3, os);
}

