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

     SourceCodec

     Abstract class for all source codecs. Takes a message as a string and encodes it in a
     vector of symbols (unsigned integers)
     
*/
#ifndef __SOURCE_CODEC_H__
#define __SOURCE_CODEC_H__

#include <vector>
#include <string>
#include <sstream>


class SourceCodec
{
public:
    SourceCodec() {}
    virtual ~SourceCodec() {}
    
    /**
     * Coder part 
     * \param in_msg Input textual message
     * \param out_msg Output vector of symbols
     * \return true if successful
     */
    virtual bool encode(const std::string& in_msg, std::vector<unsigned int>& out_msg) const = 0;
    
    /**
     * Decoder part
     * \param in_msg Input vector of symbols
     * \param out_msg Output textual message
     * \return true if successful
     */
    virtual bool decode(const std::vector<unsigned int>& in_msg, std::string& out_msg) const = 0;

    /**
     * Prints the codec data to an output stream
     */
    virtual void print_source_codec_data(std::ostringstream& os) const = 0;
};

#endif //__SOURCE_CODEC_H__
