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
#ifndef __CCSOFT_ENGINE_DEFS_H__
#define __CCSOFT_ENGINE_DEFS_H__

struct CCSoft_Engine_defs
{
public:
    typedef enum AlgoritmType_e
    {
        Algorithm_Stack,
        Algorithm_Fano
    } AlgoritmType;

    typedef enum DecodingMode_e
    {
    	Decoding_normal,
    	Decoding_regex,
    	Decoding_match_str,
    	Decoding_match_msg
    } DecodingMode;
};

#endif // __CCSOFT_ENGINE_DEFS_H__
