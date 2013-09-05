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

     Requires -std=c++0x compiler option for regex support

*/
#include "RSSoft_Engine.h"
#include "WsgcException.h"
#include <algorithm>
#include <regex>

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
		reliability(0.0),
        retry_nb(0),
        mm_cost(0)
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
	init_M(1<<m),
    gf(m, ppolys.get_ppoly(m)),
    evaluation_values(gf),
    rs_encoding(gf, k, evaluation_values),
    gskv(gf, k, evaluation_values),
    rr(gf, k),
    final_evaluation(gf, k, evaluation_values),
    retry_base_increment(1),
    retry_increment(0),
    retry_mm_strategy(RSSoft_Engine_defs::MMatrix_retry_arithmetic)
{
}


//=================================================================================================
RSSoft_Engine::~RSSoft_Engine()
{}


//=================================================================================================
void RSSoft_Engine::encode(const std::vector<unsigned int>& in_msg, std::vector<unsigned int>& out_codeword)
{
	rs_encoding.run(in_msg, out_codeword);
}


//=================================================================================================
void RSSoft_Engine::decode(rssoft::RS_ReliabilityMatrix& mat_Pi,
		std::vector<RSSoft_generic_codeword>& candidate_messages,
		bool unique)
{
	mat_Pi.normalize();
	M = init_M;

    for (unsigned int ni=1; (ni<=nb_retries); ni++) // Retry loop
    {
    	rssoft::MultiplicityMatrix mat_M(mat_Pi, M);
        try
        {
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
                        bool msg_found = false;

                        if (unique)
                        {
                            std::vector<RSSoft_generic_codeword>::iterator candidates_it = candidate_messages.begin();

                            for (; candidates_it != candidate_messages.end(); ++candidates_it)
                            {
                                if (candidates_it->get_symbols() == msg_it->get_codeword())
                                {
                                    if (msg_it->get_probability_score() > candidates_it->get_reliability())
                                    {
                                        candidates_it->set_reliability(msg_it->get_probability_score());
                                        candidates_it->set_retry_nb(ni);
                                        candidates_it->set_mm_cost(mat_M.cost());
                                    }

                                    msg_found = true;
                                    break;
                                }
                            }
                        }

                        if (!msg_found)
                        {
                            static RSSoft_generic_codeword tmp_codeword;
                            candidate_messages.push_back(tmp_codeword);
                            candidate_messages.back().set_retry_nb(ni);
                            candidate_messages.back().set_mm_cost(mat_M.cost());
                            candidate_messages.back().set_reliability(msg_it->get_probability_score());
                            candidate_messages.back().get_symbols() = msg_it->get_codeword();
                        }
                    } // Explore results
                } // Factorization successful
            } // Interpolation successful
        }
        catch (rssoft::gf::GF_Exception& e)
        {
            std::cerr << e.what() << std::endl;
        }

        new_multiplicity();
        gskv.init();
        rr.init();
        final_evaluation.init();
    } // Retry loop

    if (unique)
    {
        std::sort(candidate_messages.begin(), candidate_messages.end());
        std::reverse(candidate_messages.begin(), candidate_messages.end());
    }
}


//=================================================================================================
bool RSSoft_Engine::decode(rssoft::RS_ReliabilityMatrix& mat_Pi,
		RSSoft_generic_codeword& retrieved_message,
		RSSoft_generic_codeword& retrieved_codeword,
		const RSSoft_generic_codeword& sent_message)
{
	mat_Pi.normalize();
    M = init_M;
    bool found = false;

    for (unsigned int ni=1; (ni<=nb_retries) && !found; ni++) // Retry loop
    {
    	rssoft::MultiplicityMatrix mat_M(mat_Pi, M);
    	const rssoft::gf::GFq_BivariatePolynomial& Q = gskv.run(mat_M);
    	if (!Q.is_in_X()) // Interpolation successful
    	{
    		std::vector<rssoft::gf::GFq_Polynomial>& res_polys = rr.run(Q);

    		if (res_polys.size() > 0) // Factorization successful
    		{
    			final_evaluation.run(res_polys, mat_Pi);
    			const std::vector<rssoft::ProbabilityCodeword>& messages = final_evaluation.get_messages();
    			const std::vector<rssoft::ProbabilityCodeword>& codewords = final_evaluation.get_codewords();
    			std::vector<rssoft::ProbabilityCodeword>::const_iterator msg_it = messages.begin();
    			std::vector<rssoft::ProbabilityCodeword>::const_iterator cdw_it = codewords.begin();

    			for (; (msg_it != messages.end()) && (cdw_it != codewords.end()); ++ msg_it, ++cdw_it) // Explore results
    			{
                    if (msg_it->get_codeword() == sent_message.get_symbols())
                    {
                        retrieved_message.get_symbols() = msg_it->get_codeword();
                        retrieved_message.set_retry_nb(ni);
                        retrieved_message.set_mm_cost(mat_M.cost());
                        retrieved_message.set_reliability(msg_it->get_probability_score());
                        retrieved_codeword.get_symbols() = cdw_it->get_codeword();
                        retrieved_codeword.set_retry_nb(ni);
                        retrieved_codeword.set_mm_cost(mat_M.cost());
                        retrieved_codeword.set_reliability(cdw_it->get_probability_score());
                        found = true;
                        break;
                    }
    			} // Explore results
    		} // Factorization successful
    	} // Interpolation successful
        new_multiplicity();
        gskv.init();
        rr.init();
        final_evaluation.init();
    } // Retry loop

    return found;
}


//=================================================================================================
bool RSSoft_Engine::decode(rssoft::RS_ReliabilityMatrix& mat_Pi,
		RSSoft_generic_codeword& first_message)
{
	mat_Pi.normalize();
    M = init_M;
    bool found = false;

    for (unsigned int ni=1; (ni<=nb_retries) && !found; ni++) // Retry loop
    {
    	rssoft::MultiplicityMatrix mat_M(mat_Pi, M);
        unsigned int mm_cost = mat_M.cost();
    	const rssoft::gf::GFq_BivariatePolynomial& Q = gskv.run(mat_M);
    	if (!Q.is_in_X()) // Interpolation successful
    	{
    		std::vector<rssoft::gf::GFq_Polynomial>& res_polys = rr.run(Q);

    		if (res_polys.size() > 0) // Factorization successful
    		{
    			final_evaluation.run(res_polys, mat_Pi);
    			const std::vector<rssoft::ProbabilityCodeword>& messages = final_evaluation.get_messages();
                first_message.get_symbols() = messages[0].get_codeword();
                first_message.set_retry_nb(ni);
                first_message.set_mm_cost(mat_M.cost());
                first_message.set_reliability(messages[0].get_probability_score());
                return true;
    		} // Factorization successful
    	} // Interpolation successful
        new_multiplicity();
        gskv.init();
        rr.init();
        final_evaluation.init();
    } // Retry loop

    return false;
}


//=================================================================================================
bool RSSoft_Engine::decode_regex(rssoft::RS_ReliabilityMatrix& mat_Pi,
		std::string& retrieved_text_msg,
        RSSoft_generic_codeword& retrieved_message,
        const SourceCodec& src_codec,
        const std::string& regexp)
{
	mat_Pi.normalize();
    M = init_M;
    bool found = false;

    for (unsigned int ni=1; (ni<=nb_retries) && !found; ni++) // Retry loop
    {
    	rssoft::MultiplicityMatrix mat_M(mat_Pi, M);
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
                    src_codec.decode(msg_it->get_codeword(), retrieved_text_msg);
                    if (regexp_match(retrieved_text_msg,regexp))
                    {
                        retrieved_message.get_symbols() = msg_it->get_codeword();
                        retrieved_message.set_retry_nb(ni);
                        retrieved_message.set_mm_cost(mat_M.cost());
                        retrieved_message.set_reliability(msg_it->get_probability_score());
                        found = true;
                        break;
                    }
    			} // Explore results
    		} // Factorization successful
    	} // Interpolation successful
        new_multiplicity();
        gskv.init();
        rr.init();
        final_evaluation.init();
    } // Retry loop

    return found;
}


//=================================================================================================
bool RSSoft_Engine::decode_match(rssoft::RS_ReliabilityMatrix& mat_Pi,
		std::string& retrieved_text_msg,
        RSSoft_generic_codeword& retrieved_message,
        const SourceCodec& src_codec,
        const std::string& match_str)
{
	mat_Pi.normalize();
    M = init_M;
    bool found = false;

    for (unsigned int ni=1; (ni<=nb_retries) && !found; ni++) // Retry loop
    {
    	rssoft::MultiplicityMatrix mat_M(mat_Pi, M);
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
                    src_codec.decode(msg_it->get_codeword(), retrieved_text_msg);
                    if (retrieved_text_msg == match_str)
                    {
                        retrieved_message.get_symbols() = msg_it->get_codeword();
                        retrieved_message.set_retry_nb(ni);
                        retrieved_message.set_mm_cost(mat_M.cost());
                        retrieved_message.set_reliability(msg_it->get_probability_score());
                        found = true;
                        break;
                    }
    			} // Explore results
    		} // Factorization successful
    	} // Interpolation successful
        new_multiplicity();
        gskv.init();
        rr.init();
        final_evaluation.init();
    } // Retry loop

    return found;
}


//=================================================================================================
bool RSSoft_Engine::decode_match(rssoft::RS_ReliabilityMatrix& mat_Pi,
		RSSoft_generic_codeword& retrieved_message,
        const std::vector<unsigned int>& match_message)
{
	mat_Pi.normalize();
    M = init_M;
    bool found = false;

    for (unsigned int ni=1; (ni<=nb_retries) && !found; ni++) // Retry loop
    {
    	rssoft::MultiplicityMatrix mat_M(mat_Pi, M);
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
                    if (msg_it->get_codeword() == match_message)
                    {
                        retrieved_message.get_symbols() = msg_it->get_codeword();
                        retrieved_message.set_retry_nb(ni);
                        retrieved_message.set_mm_cost(mat_M.cost());
                        retrieved_message.set_reliability(msg_it->get_probability_score());
                        found = true;
                        break;
                    }
    			} // Explore results
    		} // Factorization successful
    	} // Interpolation successful
        new_multiplicity();
        gskv.init();
        rr.init();
        final_evaluation.init();
    } // Retry loop

    return found;
}


//=================================================================================================
bool RSSoft_Engine::decode(rssoft::RS_ReliabilityMatrix& mat_Pi,
		RSSoft_generic_codeword& retrieved_message,
		float reliability_threshold)
{
	mat_Pi.normalize();
    M = init_M;
    bool found = false;

    for (unsigned int ni=1; (ni<=nb_retries) && !found; ni++) // Retry loop
    {
    	rssoft::MultiplicityMatrix mat_M(mat_Pi, M);
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
                    if (msg_it->get_probability_score() > reliability_threshold)
                    {
                        retrieved_message.get_symbols() = msg_it->get_codeword();
                        retrieved_message.set_retry_nb(ni);
                        retrieved_message.set_mm_cost(mat_M.cost());
                        retrieved_message.set_reliability(msg_it->get_probability_score());
                        found = true;
                        break;
                    }
    			} // Explore results
    		} // Factorization successful
    	} // Interpolation successful
        new_multiplicity();
        gskv.init();
        rr.init();
        final_evaluation.init();
    } // Retry loop

    return found;
}


//=================================================================================================
bool RSSoft_Engine::regexp_match(const std::string& value,
		const std::string& regexp) const
{
	try
	{
		return std::regex_match(value, std::regex(regexp));
	}
	catch (std::regex_error& e)
	{
		std::cout << "Regular expression error: " << e.what() << " : " << e.code() << std::endl;
	}
	return false;
}

//=================================================================================================
void RSSoft_Engine::new_multiplicity()
{
	switch (retry_mm_strategy)
	{
	case RSSoft_Engine_defs::MMatrix_retry_arithmetic:
		M += retry_base_increment;
		break;
	case RSSoft_Engine_defs::MMatrix_retry_arithmetic_increment:
		retry_increment += retry_base_increment;
		M += retry_increment;
		break;
	case RSSoft_Engine_defs::MMatrix_retry_geometric:
		M *= retry_base_increment;
		break;
	case RSSoft_Engine_defs::MMatrix_retry_geometric_increment:
		retry_increment *= retry_base_increment;
		M *= retry_increment;
		break;
	default:
		break;
	}
}

//=================================================================================================
float RSSoft_Engine::calculate_reliability(rssoft::RS_ReliabilityMatrix& mat_Pi,
		const std::vector<unsigned int>& codeword)
{
	if (codeword.size() < n)
	{
		throw WsgcException("RSSoft: incorrect size of codeword to calculate its reliability");
	}
	else
	{
		float proba_score = 0.0; // We will use log scale in dB/symbol
		unsigned int proba_count = 0; // number of individual symbol probabilities considered

		for (unsigned int i_pt = 0; i_pt < n; i_pt++)
		{
			float p_ij = mat_Pi(codeword[i_pt], i_pt);

			if (p_ij != 0) // not erased
			{
				proba_score += 10.0 * log10(p_ij); // Accumulate probability (dB)
				proba_count++;
			}
		}

		return proba_score/proba_count;
	}
}
