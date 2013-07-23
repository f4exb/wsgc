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
#ifndef __RSSOFT_DECISION_BOX_H__
#define __RSSOFT_DECISION_BOX_H__

#include "RSSoft_Engine.h"
#include <iostream>

class SourceCodec;

class RSSoft_DecisionBox
{
public:
    RSSoft_DecisionBox(RSSoft_Engine& _rssoft_engine, const SourceCodec *_source_codec);
    ~RSSoft_DecisionBox();
    
    void run(Options::RSSoft_decoding_mode rs_decoding_mode);
    void run_regex(const std::string& rs_decoding_regex);
    
protected:
    void full_scan(std::vector<RSSoft_generic_codeword>& candidate_messages);
    bool find_first(RSSoft_generic_codeword& unique_message);
    bool regex_scan(std::string& decoded_text, RSSoft_generic_codeword& unique_message, const std::string& rs_decoding_regex);
    void list_messages(const std::vector<RSSoft_generic_codeword>& candidate_messages, std::ostream& os);
    void show_message(const RSSoft_generic_codeword& message, std::ostream& os);

    RSSoft_Engine& rssoft_engine;
    const SourceCodec *source_codec;
};


#endif // __RSSOFT_DECISION_BOX_H__
