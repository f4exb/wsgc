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
#ifndef __CCSOFT_DECISION_BOX_H__
#define __CCSOFT_DECISION_BOX_H__

#include "CCSoft_Engine.h"
#include "WsgcUtils.h"
#include <iostream>

class SourceCodec;

namespace ccsoft
{
    class CC_ReliabilityMatrix;
}

class CCSoft_DecisionBox
{
public:
    CCSoft_DecisionBox(CCSoft_Engine& _ccsoft_engine, const SourceCodec *_source_codec);
    ~CCSoft_DecisionBox();
    
    void run(ccsoft::CC_ReliabilityMatrix& relmat);
    void print_retrieved_message(std::ostream& os);
    void print_stats(const std::vector<unsigned int>& sent_message, const std::vector<unsigned int>& sent_codeword, std::ostream& os);
    
protected:
    CCSoft_Engine& ccsoft_engine;
    const SourceCodec *source_codec;
    std::vector<unsigned int> retrieved_message;
    float retrieved_message_score;
    bool decode_success;
};


#endif // __CCSOFT_DECISION_BOX_H__
