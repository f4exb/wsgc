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
      
     Decision box specialized in Reed-Solomon soft decision decoding with RSSoft library

*/
#include "Options.h"
#include "RSSoft_DecisionBox.h"
#include "SourceCodec.h"
#include "WsgcUtils.h"

//=================================================================================================
RSSoft_DecisionBox::RSSoft_DecisionBox(RSSoft_Engine& _rssoft_engine, const SourceCodec *_source_codec) :
    rssoft_engine(_rssoft_engine), source_codec(_source_codec)
{}

//=================================================================================================
RSSoft_DecisionBox::~RSSoft_DecisionBox()
{}

    
//=================================================================================================
void RSSoft_DecisionBox::run(Options::RSSoft_decoding_mode rs_decoding_mode)
{
    switch (rs_decoding_mode)
    {
        case Options::RSSoft_decoding_all:
            full_scan_all();
            list_messages(std::cout);
            break;
            
        case Options::RSSoft_decoding_full:
            full_scan_unique();
            list_messages(std::cout);
            break;
            
        case Options::RSSoft_decoding_best:
            full_scan_unique();
            
            if (candidate_messages.size() > 0)
            {
                show_message(candidate_messages[0], std::cout);
            }
            else
            {
                std::cout << "No solution found" << std::endl;
            }
            break;
            
        case Options::RSSoft_decoding_first:
            if (find_first())
            {
                show_message(candidate_messages[0], std::cout);
            }
            else
            {
                std::cout << "No solution found" << std::endl;
            }
            break;
            
        default:
            std::cout << "Unrecognized RS decoding mode" << std::endl;
            break;
    }
}

//=================================================================================================
void RSSoft_DecisionBox::run_regex(const std::string& rs_decoding_regex)
{
    std::string text_message;

    if (source_codec)
    {
        if (regex_scan(text_message, rs_decoding_regex))
        {
            show_message(candidate_messages[0], std::cout);
        }
        else
        {
            std::cout << "No solution found" << std::endl;
        }
    }
    else
    {
    	std::cout << "Cannot use RS soft decision decoding with regex without source codec" << std::endl;
    }
}

//=================================================================================================
void RSSoft_DecisionBox::run_match(const std::string& rs_decoding_matching_source_message)
{
    std::string text_message;

    if (source_codec)
    {
        if (match_scan(text_message, rs_decoding_matching_source_message))
        {
            show_message(candidate_messages[0], std::cout);
        }
        else
        {
            std::cout << "No solution found" << std::endl;
        }
    }
    else
    {
    	std::cout << "Cannot use RS soft decision decoding with match string without source codec" << std::endl;
    }
}

//=================================================================================================
void RSSoft_DecisionBox::run_match(const std::vector<unsigned int>& matching_codeword)
{
    if (match_scan(matching_codeword))
    {
        show_message(candidate_messages[0], std::cout);
    }
    else
    {
        std::cout << "No solution found" << std::endl;
    }
}

//=================================================================================================
void RSSoft_DecisionBox::run_reliability_threshold(float reliability_threshold)
{
    if (find_first_above_reliability_threshold(reliability_threshold))
    {
        show_message(candidate_messages[0], std::cout);
    }
    else
    {
        std::cout << "No solution found" << std::endl;
    }
}

//=================================================================================================
void RSSoft_DecisionBox::full_scan_unique()
{
    candidate_messages.clear();
    
    rssoft_engine.decode(candidate_messages, true);
}

//=================================================================================================
void RSSoft_DecisionBox::full_scan_all()
{
    candidate_messages.clear();
    
    rssoft_engine.decode(candidate_messages, false);
}

//=================================================================================================
bool RSSoft_DecisionBox::find_first()
{
    RSSoft_generic_codeword unique_message;
    candidate_messages.clear();
    
    if (rssoft_engine.decode(unique_message))
    {
        candidate_messages.push_back(unique_message);
        return true;
    }
    else
    {
        return false;
    }
}

//=================================================================================================
bool RSSoft_DecisionBox::find_first_above_reliability_threshold(float reliability_threshold)
{
    RSSoft_generic_codeword unique_message;
    candidate_messages.clear();
    
    if (rssoft_engine.decode(unique_message, reliability_threshold))
    {
        candidate_messages.push_back(unique_message);
        return true;
    }
    else
    {
        return false;
    }
}

//=================================================================================================
bool RSSoft_DecisionBox::regex_scan(std::string& decoded_text, const std::string& rs_decoding_regex)
{
    RSSoft_generic_codeword unique_message;
    candidate_messages.clear();
    
    if (rssoft_engine.decode_regex(decoded_text, unique_message, *source_codec, rs_decoding_regex))
    {
        candidate_messages.push_back(unique_message);
        return true;
    }
    else
    {
        return false;
    }
}

//=================================================================================================
bool RSSoft_DecisionBox::match_scan(std::string& decoded_text, const std::string& rs_decoding_matching_source_message)
{
    RSSoft_generic_codeword unique_message;
    candidate_messages.clear();

    if (rssoft_engine.decode_match(decoded_text, unique_message, *source_codec, rs_decoding_matching_source_message))
    {
        candidate_messages.push_back(unique_message);
        return true;
    }
    else
    {
        return false;
    }
}

//=================================================================================================
bool RSSoft_DecisionBox::match_scan(const std::vector<unsigned int>& matching_codeword)
{
    RSSoft_generic_codeword unique_message;
    candidate_messages.clear();

    if (rssoft_engine.decode_match(unique_message, matching_codeword))
    {
        candidate_messages.push_back(unique_message);
        return true;
    }
    else
    {
        return false;
    }
}

//=================================================================================================
void RSSoft_DecisionBox::list_messages(std::ostream& os)
{
    std::vector<RSSoft_generic_codeword>::const_iterator msg_it = candidate_messages.begin();
    const std::vector<RSSoft_generic_codeword>::const_iterator msg_end = candidate_messages.end();
    unsigned int i = 0;
    std::string decoded_text;
    
    for (; msg_it != msg_end; ++msg_it, i++)
    {
        os << i << ": try #" << msg_it->get_retry_nb() << " mm_cost " << msg_it->get_mm_cost() << " rel " << msg_it->get_reliability() << " dB/Symbol ";
        print_vector<unsigned int, unsigned int, std::ostream>(msg_it->get_symbols(), 3, os);
        if (source_codec)
        {
        	if (source_codec->decode(msg_it->get_symbols(), decoded_text))
        	{
        		os << " \"" << decoded_text << "\"";
        	}
        	else
        	{
        		os << " <cannot decode source message>" << std::endl;
        	}
        }
        std::cout << std::endl;
    }
}

//=================================================================================================
void RSSoft_DecisionBox::show_message(const RSSoft_generic_codeword& message, std::ostream& os)
{
    os << "try #" << message.get_retry_nb() << " mm_cost " << message.get_mm_cost() << " rel " << message.get_reliability() << " dB/Symbol ";
    print_vector<unsigned int, unsigned int, std::ostream>(message.get_symbols(), 3u, os);
    std::string decoded_text;

    if (source_codec)
    {
        if (source_codec->decode(message.get_symbols(), decoded_text))
        {
            os << " \"" << decoded_text << "\"";
        }
        else
        {
            os << " <cannot decode source message>" << std::endl;
        }
    }

    std::cout << std::endl;
}

//=================================================================================================
void RSSoft_DecisionBox::print_stats(RSSoft_Engine& rssoft_engine, const std::vector<unsigned int>& sent_message, const std::vector<unsigned int>& sent_codeword, std::ostream& os)
{
    float source_codeword_score = rssoft_engine.calculate_reliability(sent_codeword);
    
    if (candidate_messages.size() == 0) // no result found
    {
        os << "_RES " << "0,0," << source_codeword_score << ",0,0" << std::endl;
    }
    else
    {
        std::vector<RSSoft_generic_codeword>::const_iterator msg_it = candidate_messages.begin();
        const std::vector<RSSoft_generic_codeword>::const_iterator msg_end = candidate_messages.end();
        
        for (; msg_it != msg_end; ++msg_it)
        {
            unsigned int found = (msg_it->get_symbols() == sent_message ? 1 : 0);
            float score = msg_it->get_reliability();
            unsigned int try_nb = msg_it->get_retry_nb();
            unsigned int mm_cost = msg_it->get_mm_cost();
            os << "_RES 1," << found << "," << score << "," << try_nb << "," << mm_cost << std::endl;
        }
    }
}
