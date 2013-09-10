/*
     Copyright 2012-2013 Edouard Griffiths <f4exb at free dot fr>

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

     CCSoft_Engine

     Convolutional coding soft decision engine based on ccsoft library
     https://code.google.com/p/rssoft/#Convolutional_codes_library

*/

#include "CCSoft_Engine.h"
#include "CCSoft_Exception.h"
#include "CC_StackDecoding_FA.h"
#include "CC_FanoDecoding_FA.h"

#include <cmath>

//=================================================================================================
CCSoft_Engine::CCSoft_Engine(const std::vector<unsigned int>& _k_constraints,
        const std::vector<std::vector<unsigned int> >& _generator_polys) :
            k_constraints(_k_constraints),
            generator_polys(_generator_polys),
            k(_k_constraints.size()),
            n(_generator_polys.size()),
            algorithm_type(CCSoft_Engine::Algorithm_Fano),
            cc_decoding(0),
            fano_init_metric(0.0),
            fano_delta_metric(1.0),
            fano_tree_cache_size(0),
            edge_bias(0.0),
            verbosity(0)
{}


//=================================================================================================
CCSoft_Engine::~CCSoft_Engine()
{
    if (cc_decoding)
    {
        delete cc_decoding;
    }
}


//=================================================================================================
void CCSoft_Engine::init_decoding_stack(float edge_bias)
{
    cc_decoding = new ccsoft::CC_StackDecoding_FA<unsigned int, unsigned int, 1>(k_constraints, generator_polys);
    cc_decoding->set_edge_bias(edge_bias);
}

//=================================================================================================
void CCSoft_Engine::init_decoding_fano(float edge_bias,
        float init_threshold,
        float delta_threshold,
        unsigned int tree_cache_size,
        float init_threshold_delta)
{
    cc_decoding = new ccsoft::CC_FanoDecoding_FA<unsigned int, unsigned int, 1>(k_constraints,
            generator_polys,
            init_threshold,
            delta_threshold,
            tree_cache_size,
            init_threshold_delta);
    cc_decoding->set_edge_bias(edge_bias);
}

//=================================================================================================
void CCSoft_Engine::print_convolutional_code_data(std::ostream& os)
{
    cc_decoding->get_encoding().print(os);
}

//=================================================================================================
void CCSoft_Engine::print_stats(std::ostream& os, bool decode_ok)
{
    cc_decoding->print_stats(os, decode_ok);
}

//=================================================================================================
void CCSoft_Engine::encode(const std::vector<unsigned int>& in_msg, std::vector<unsigned int>& out_codeword)
{
    std::vector<unsigned int>::const_iterator in_it = in_msg.begin();
    out_codeword.clear();
    
    for (; in_it != in_msg.end(); ++in_it)
    {
        out_codeword.push_back(0);
        cc_decoding->get_encoding().encode(*in_it, out_codeword.back());
    }
}
    
//=================================================================================================
void CCSoft_Engine::reset()
{
    if (cc_decoding)
    {
        cc_decoding->reset();
    }
}

//=================================================================================================
bool CCSoft_Engine::decode(ccsoft::CC_ReliabilityMatrix& relmat, std::vector<unsigned int>& retrieved_msg, float& score)
{
    bool success = false;
    relmat.normalize();
    
    if (cc_decoding)
    {
        success = cc_decoding->decode(relmat, retrieved_msg);
        score = cc_decoding->get_score();
        
        if (verbosity > 0)
        {
            if (success)
            {
                std::cout << "decoding successful with retrieved message score of " << score << std::endl;
            }
            else
            { 
                std::cout << "decoding unsuccessful" << std::endl;
            }
        }
    }
    
    return success;
}

//=================================================================================================
void CCSoft_Engine::interleave(std::vector<unsigned int>& symbols, bool forward)
{
	if (cc_decoding)
	{
		cc_decoding->interleave(symbols, forward);
	}
}
