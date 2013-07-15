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

     RSSoft_Engine

     Reed Solomon soft decision engine based on rssoft library
     https://code.google.com/p/rssoft

*/
#include "RSSoft_Engine.h"
#include "WsgcException.h"
#include <algorithm>

//=================================================================================================
RSSoft_PPolys::RSSoft_PPolys()
{
        rssoft::gf::GF2_Element pp_gf8[4]   = {1,1,0,1};
        rssoft::gf::GF2_Element pp_gf16[5]  = {1,0,0,1,1};
        rssoft::gf::GF2_Element pp_gf32[6]  = {1,0,0,1,0,1};
        rssoft::gf::GF2_Element pp_gf64[7]  = {1,0,0,0,0,1,1};
        rssoft::gf::GF2_Element pp_gf128[8] = {1,0,0,0,0,0,1,1};
        rssoft::gf::GF2_Element pp_gf256[9] = {1,0,0,0,1,1,1,0,1};
    
        ppolys.push_back(rssoft::gf::GF2_Polynomial(4,pp_gf8));
        ppolys.push_back(rssoft::gf::GF2_Polynomial(5,pp_gf16));
        ppolys.push_back(rssoft::gf::GF2_Polynomial(6,pp_gf32));
        ppolys.push_back(rssoft::gf::GF2_Polynomial(7,pp_gf64));
        ppolys.push_back(rssoft::gf::GF2_Polynomial(8,pp_gf128));
        ppolys.push_back(rssoft::gf::GF2_Polynomial(9,pp_gf256));    
}


//=================================================================================================
RSSoft_PPolys::~RSSoft_PPolys()
{}


//=================================================================================================
RSSoft_generic_codeword::RSSoft_generic_codeword() :
		reliability(0.0)
{}


//=================================================================================================
RSSoft_generic_codeword::~RSSoft_generic_codeword()
{}


//=================================================================================================
const rssoft::gf::GF2_Polynomial& RSSoft_PPolys::get_ppoly(unsigned int order)
{
    if ((order < 3) || (order > 8))
    {
        throw WsgcException("Galois Field size unsupported");
    }
    else
    {
        return ppolys[order-3];
    }
}


//=================================================================================================
RSSoft_Engine::RSSoft_Engine(unsigned int _m, unsigned int _k) :
	m(_m),
	n((1<<m)-1),
	q(1<<m),
	k(_k),
	nb_retries(5),
	M(1<<m),
    gf(m, ppolys.get_ppoly(m)),
    evaluation_values(gf),
    rs_encoding(gf, k, evaluation_values),
    mat_Pi(m,n),
    gskv(gf, k, evaluation_values),
    rr(gf, k),
    final_evaluation(gf, k, evaluation_values)
{
}


//=================================================================================================
RSSoft_Engine::~RSSoft_Engine()
{}


//=================================================================================================
void RSSoft_Engine::encode(const std::vector<int>& in_msg, std::vector<int>& out_codeword)
{
	rs_encoding.run(in_msg, out_codeword);
}


//=================================================================================================
void RSSoft_Engine::record_magnitudes(wsgc_float *magnitudes)
{
	 mat_Pi.enter_symbol_data(magnitudes);
}


//=================================================================================================
void RSSoft_Engine::decode(std::vector<RSSoft_generic_codeword>& candidate_messages)
{
	mat_Pi.normalize();
	unsigned int global_multiplicity = M;

    for (unsigned int ni=1; (ni<=nb_retries); ni++)
    {
    	rssoft::MultiplicityMatrix mat_M(mat_Pi, global_multiplicity);
    	const rssoft::gf::GFq_BivariatePolynomial& Q = gskv.run(mat_M);
    	if (!Q.is_in_X()) // Interpolation successful
    	{
    		std::vector<rssoft::gf::GFq_Polynomial>& res_polys = rr.run(Q);

    		if (res_polys.size() > 0) // Factorization successful
    		{
    			final_evaluation.run(res_polys, mat_Pi);
    			const std::vector<rssoft::ProbabilityCodeword>& messages = final_evaluation.get_messages();
    			std::vector<rssoft::ProbabilityCodeword>::const_iterator msg_it = messages.begin();

    			for (; msg_it != messages.end(); ++ msg_it) // Explore results
    			{
    				std::vector<RSSoft_generic_codeword>::iterator candidates_it = candidate_messages.begin();
    				bool msg_found = false;

    				for (; candidates_it != candidate_messages.end(); ++candidates_it)
    				{
    					if (candidates_it->get_symbols() == msg_it->get_codeword())
    					{
    						if (msg_it->get_probability_score() > candidates_it->get_reliability())
    						{
    							candidates_it->set_reliability(msg_it->get_probability_score());
        						msg_found = true;
    						}
    					}
    				}

    				if (!msg_found)
    				{
    					static RSSoft_generic_codeword tmp_codeword;
    					candidate_messages.push_back(tmp_codeword);
    					candidate_messages.back().set_reliability(msg_it->get_probability_score());
    					candidate_messages.back().get_symbols() = msg_it->get_codeword();
    				}
    			} // Explore results
    		} // Factorization successful
    	} // Interpolation successful
    } // Retry loop

    std::sort(candidate_messages.begin(), candidate_messages.end());
    std::reverse(candidate_messages.begin(), candidate_messages.end());
}
