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
      
     Decision box specialized in Convlutional Coding soft decision decoding with CCSoft library

*/
#include "CCSoft_DecisionBox.h"
#include "SourceCodec.h"
#include "CC_ReliabilityMatrix.h"

//=================================================================================================
CCSoft_DecisionBox::CCSoft_DecisionBox(CCSoft_Engine& _ccsoft_engine, const SourceCodec *_source_codec) :
    ccsoft_engine(_ccsoft_engine),
    source_codec(_source_codec),
    retrieved_message_score(0.0),
    decode_success(false)
{}

//=================================================================================================
CCSoft_DecisionBox::~CCSoft_DecisionBox()
{}
  
//=================================================================================================
void CCSoft_DecisionBox::run(ccsoft::CC_ReliabilityMatrix& relmat)
{
    decode_success = ccsoft_engine.decode(relmat, retrieved_message, retrieved_message_score);
}

//=================================================================================================
void CCSoft_DecisionBox::print_retrieved_message(std::ostream& os)
{
    if (decode_success)
    {
        os << "score = " << retrieved_message_score << std::endl;
        print_vector<unsigned int, unsigned int, std::ostream>(retrieved_message, 3u, os);
        std::string decoded_text;

        if (source_codec)
        {
            if (source_codec->decode(retrieved_message, decoded_text))
            {
                os << " \"" << decoded_text << "\"";
            }
            else
            {
                os << " <cannot decode source message>" << std::endl;
            }
        }

        os << std::endl;
    }
    else
    {
        os << "Decoding failed, no retrieved message" << std::endl;
    }
}

//=================================================================================================
void CCSoft_DecisionBox::print_stats(const std::vector<unsigned int>& sent_message, const std::vector<unsigned int>& sent_codeword, std::ostream& os)
{
    if (decode_success)
    {
        bool decode_ok = (sent_message == retrieved_message);
        ccsoft_engine.print_stats(os, decode_ok);
    }
    else
    {
        ccsoft_engine.print_stats(os, false);
    }
}
    
