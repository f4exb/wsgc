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
#ifndef __RSSOFT_ENGINE_H__
#define __RSSOFT_ENGINE_H__

#include "GFq.h"
#include "GF2_Element.h"
#include "GF2_Polynomial.h"
#include "GF_Utils.h"
#include "EvaluationValues.h"
#include "ReliabilityMatrix.h"
#include "MultiplicityMatrix.h"
#include "GSKV_Interpolation.h"
#include "RR_Factorization.h"
#include "FinalEvaluation.h"
#include "RS_Encoding.h"
#include "WsgcTypes.h"
#include "SourceCodec.h"
#include <string>

/**
 * \brief Class to handle primitive polynomials for RSSoft engine
 */
class RSSoft_PPolys
{
public:
    RSSoft_PPolys();
    ~RSSoft_PPolys();
    const rssoft::gf::GF2_Polynomial& get_ppoly(unsigned int order);
private:
    std::vector<rssoft::gf::GF2_Polynomial> ppolys; //!< Predefined primitive polynomials
};

/**
 * \brief Generic codeword class, can be either a message or an encoded codeword
 */
class RSSoft_generic_codeword
{
public:
	RSSoft_generic_codeword();
	~RSSoft_generic_codeword();

	std::vector<unsigned int>& get_symbols()
	{
		return symbols;
	}

	const std::vector<unsigned int>& get_symbols() const
	{
		return symbols;
	}

	float get_reliability() const
	{
		return reliability;
	}

	void set_reliability(float _reliability)
	{
		reliability = _reliability;
	}
    
	unsigned int get_retry_nb() const
    {
        return retry_nb;
    }

    void set_retry_nb(unsigned int _retry_nb)
    {
        retry_nb = _retry_nb;
    }

	unsigned int get_mm_cost() const
    {
        return mm_cost;
    }

    void set_mm_cost(unsigned int _mm_cost)
    {
        mm_cost = _mm_cost;
    }
    
	/**
	 * reliability increasing order
	 */
	bool operator<(const RSSoft_generic_codeword& other) const
	{
		return reliability < other.reliability;
	}

protected:
	std::vector<unsigned int> symbols;
	float reliability;
    unsigned int retry_nb;
    unsigned int mm_cost;
};

/**
 * \brief Interface engine with the RSSoft Reed Solomon soft decision librarys
 */
class RSSoft_Engine
{
public:
	typedef enum
	{
		MMatrix_retry_arithmetic,            //!< Mul(Mn+1) = Mul(Mn) + inc
		MMatrix_retry_arithmetic_increment,  //!< Mul(Mn+1) = Mul(Mn) + (n+1)*inc
		MMatrix_retry_geometric,             //!< Mul(Mn+1) = Mul(Mn) * inc
		MMatrix_retry_geometric_increment    //!< Mul(Mn+1) = Mul(Mn) * (inc^(n+1))
	} MultiplicityMatrix_RetryStrategy;


	RSSoft_Engine(unsigned int _m, unsigned int _k);

	~RSSoft_Engine();

	unsigned int get_m() const
	{
		return m;
	}

	unsigned int get_n() const
	{
		return n;
	}

	unsigned int get_q() const
	{
		return q;
	}

	unsigned int get_k() const
	{
		return k;
	}
    
    void set_initial_global_multiplicity(unsigned int _init_M)
    {
        M = _init_M;
        init_M = _init_M;
    }
    
    void set_nb_retries(unsigned int _nb_retries)
    {
        nb_retries = _nb_retries;
    }
    
    void set_retry_base_increment(unsigned int _retry_base_increment)
    {
    	retry_base_increment = _retry_base_increment;
    }

    void set_retry_mm_strategy(MultiplicityMatrix_RetryStrategy _retry_mm_strategy)
    {
    	retry_mm_strategy = _retry_mm_strategy;
    }

    /**
     * Get reference to the reliability matrix for direct update
     * \return r/w reference to the reliability matrix
     */
    rssoft::ReliabilityMatrix& get_reliability_matrix()
    {
        return mat_Pi;
    }
    
    /**
     * Encode message into codeword
     * \param in_msg Message
     * \param out_codeword Codeword
     */
    void encode(const std::vector<unsigned int>& in_msg, std::vector<unsigned int>& out_codeword);
    
    /**
     * Record symbol magnitudes for one symbol position
     * \param magnitudes Pointer to magnitudes array. It is assumed to be of 2^m size.
     */
    void record_magnitudes(wsgc_float *magnitudes);
    
    /**
     * Decode one codeword based on given magnitudes for the length of one codeword. Return candidate
     * messages (only unique messages with best reliability if unique is true) along with their best value of reliability. 
     * It will systematically loop the given number of retries.
     * \param candidate_messages Vector of candidate messages filled in in decreasing value of reliability
     */
    void decode(std::vector<RSSoft_generic_codeword>& candidate_messages, bool unique=true);
    
    /**
     * Decode one codeword based on given magnitudes for the length of one codeword. Returns only the first candidate
     * message along with its reliability value. 
     * \param first_message First message found
     */
    bool decode(RSSoft_generic_codeword& first_message);
    
    /**
     * Decode one codeword based on given magnitudes for the length of one codeword. Tries to find the message that was sent
     * and returns the corresponding decoded message and codeword. Of course this is of no use in real life conditions when you 
     * don't know which message was sent. This is to be used only for prototyping.
     * \param retrieved_message Retrieved message
     * \param retrieved_codeword Retrieved codeword
     * \param sent_message Message that was sent
     * \return true if the message sent was found
     */
    bool decode(RSSoft_generic_codeword& retrieved_message, RSSoft_generic_codeword& retrieved_codeword, const RSSoft_generic_codeword& sent_message);
    
    /**
     * Decode one codeword based on given magnitudes for the length of one codeword. Tries to match found textual messages with the given
     * regular expression and returns the first match
     * \param retrieved_text_msg Retrieved textual message
     * \param retrieved_message Retrieved message symbols and associated reliability data
     * \param src_codec Source codec being used
     * \param regexp Regular expression to match with the textual message
     * \return true if a matching message is found
     */
    bool decode(std::string& retrieved_text_msg,
        RSSoft_generic_codeword& retrieved_message, 
        const SourceCodec& src_codec,
        const std::string& regexp);

    /**
     * Decode one codeword based on given magnitudes for the length of one codeword. Returns the first candidate showing a reliability figure above a given threshold
     * \param retrieved_message Retrieved message symbols and associated reliability data
     * \param src_codec Source codec being used
     * \param regexp Regular expression to match with the textual message
     * \return true if a matching message is found
     */
    bool decode(RSSoft_generic_codeword& retrieved_message, float reliability_threshold);

    /**
     * Calculate the reliability of a codeword once the reliability matrix has been normalized
     * i.e. after decode
     * \codeword codeword for which to calculate the reliability
     * \return codeword reliability with respect to current reliability matrix
     */
    float calculate_reliability(std::vector<unsigned int> codeword);
        
protected:
    bool regexp_match(const std::string& value, const std::string& regexp) const;
    void new_multiplicity();
	unsigned int m; //!< GF(2^m)
	unsigned int n; //!< 2^m-1
	unsigned int q; //!< 2^m
	unsigned int k; //!< RS(n,k), n = 2^m-1
	unsigned int nb_retries; //!< Number of soft decision retries
	unsigned int M; //!< Global multiplicity for soft decision multiplicity matrix
	unsigned int init_M; //!< Initial value for global multiplicity for soft decision multiplicity matrix at first try
	unsigned int retry_base_increment; //!< multiplicity base increment for retry
	unsigned int retry_increment; //!< current multiplicity increment for retry
	MultiplicityMatrix_RetryStrategy retry_mm_strategy; //!< strategy for multiplicity increment
    RSSoft_PPolys ppolys; //!< Collection of pre-defined primitive polynomials
    rssoft::gf::GFq gf; //!< Galois Field being used
    rssoft::EvaluationValues evaluation_values; //!< Evaluation values for RS
    rssoft::RS_Encoding rs_encoding; //!< Encoder
    rssoft::ReliabilityMatrix mat_Pi; //!< Reliability matrix
    rssoft::GSKV_Interpolation gskv; //!< Guruswami-Sudan-Koetter-Vardy interpolation engine
    rssoft::RR_Factorization rr; //!< Roth-Ruckensteil factorization engine
    rssoft::FinalEvaluation final_evaluation; //!< Evaluation engine
};



#endif // __RSSOFT_ENGINE_H__
