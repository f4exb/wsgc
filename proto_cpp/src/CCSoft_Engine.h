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
#ifndef __CCSOFT_ENGINE_H__
#define __CCSOFT_ENGINE_H__

#include "CC_ReliabilityMatrix.h"
#include "CC_SequentialDecoding.h"

/**
 * \brief Interface engine with the RSSoft Reed Solomon soft decision librarys
 */
class CCSoft_Engine
{
public:
    typedef enum
    {
        Algorithm_Stack,
        Algorithm_Fano
    } AlgoritmType;

    CCSoft_Engine(const std::vector<unsigned int>& _k_constraints,
            const std::vector<std::vector<unsigned int> >& _generator_polys);
           
    ~CCSoft_Engine();

    /**
     * Initializes decoding and encoding objects with stack algorithm
     */
    void init_decoding_stack(float edge_bias);

    /**
     * Initializes decoding and encoding objects with stack algorithm
     */
    void init_decoding_fano(float edge_bias,
        float init_threshold,
        float delta_threshold,
        unsigned int tree_cache_size,
        float init_threshold_delta);

    /**
     * Set a metric limit for the decoder
     */
    void set_metric_limit(float metric_limit)
    {
        if (cc_decoding)
        {
            cc_decoding->set_metric_limit(metric_limit);
        }
    }
    
    /**
     * Set a node limit for the decoder
     */
    void set_node_limit(unsigned int nb_of_nodes)
    {
        if (cc_decoding)
        {
            cc_decoding->set_node_limit(nb_of_nodes);
        }
    }
    
    /**
     * Get the size of input symbols in bits (k)
     */ 
    unsigned int get_k() const
    {
        return k;
    }
    
    /**
     * Get the size of outpuy symbols in bits (n)
     */ 
    unsigned int get_n() const
    {
        return n;
    }
    
    /**
     * Get the largest register length (m)
     */ 
    unsigned int get_m() const
    {
        if (cc_decoding)
        {
            return cc_decoding->get_encoding().get_m();
        }
        else
        {
            return 0;
        }
    }
    
    /**
     * Print convolutional code data to output stream
     * \param os Output stream
     */
    void print_convolutional_code_data(std::ostream& os);
    
    /**
     * Print statistics about the decoing to output stream
     * \param os Output stream
     * \param decode_ok True to show success
     */
    void print_stats(std::ostream& os, bool decode_ok);
    
    /**
     * Suffix message with m-1 zeros
     * \param msg Message
     */
    void zero_pad(std::vector<unsigned int>& msg)
    {
        for (unsigned int i=0; i < get_m() - 1; i++)
        {
            msg.push_back(0);
        }
    }

    /**
     * Encode message into codeword
     * \param in_msg Message
     * \param out_codeword Codeword
     */
    void encode(const std::vector<unsigned int>& in_msg, std::vector<unsigned int>& out_codeword);
    
    /**
     * Reset state to begin a new decoding cycle. Resets decoding object.
     */
    void reset();
    
    /**
     * Decodes message following reliability matrix data
     * \param relmat Reference to the reliability matrix
     * \param retrieved_msg Reference to the vector of retrieved input symbols
     * \param score Metric of the retrieved input message
     */
    bool decode(ccsoft::CC_ReliabilityMatrix& relmat, std::vector<unsigned int>& retrieved_msg, float& score);

protected:
    const std::vector<unsigned int>& k_constraints;
    const std::vector<std::vector<unsigned int> >& generator_polys;
    unsigned int k;
    unsigned int n;
    AlgoritmType algorithm_type; //!< Type of algorithm used
    ccsoft::CC_SequentialDecoding<unsigned int, unsigned int> *cc_decoding;
    float fano_init_metric;
    float fano_delta_metric;
    unsigned int fano_tree_cache_size;
    float edge_bias;
    unsigned int verbosity; //!< Verbosity level
};

#endif // __CCSOFT_ENGINE_H__
