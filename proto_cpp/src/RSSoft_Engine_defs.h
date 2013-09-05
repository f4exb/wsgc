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
#ifndef __RSSOFT_ENGINE_DEFS_H__
#define __RSSOFT_ENGINE_DEFS_H__

struct RSSoft_Engine_defs
{
public:
    typedef enum
    {
        MMatrix_retry_arithmetic,            //!< Mul(Mn+1) = Mul(Mn) + inc
        MMatrix_retry_arithmetic_increment,  //!< Mul(Mn+1) = Mul(Mn) + (n+1)*inc
        MMatrix_retry_geometric,             //!< Mul(Mn+1) = Mul(Mn) * inc
        MMatrix_retry_geometric_increment    //!< Mul(Mn+1) = Mul(Mn) * (inc^(n+1))
    } MultiplicityMatrix_RetryStrategy;

    typedef enum
    {
        RSSoft_decoding_all,
        RSSoft_decoding_full,
        RSSoft_decoding_best,
        RSSoft_decoding_first,
        RSSoft_decoding_regex,
        RSSoft_decoding_match,
        RSSoft_decoding_binmatch,
        RSSoft_decoding_relthr
    } RSSoft_decoding_mode;
};

#endif // __RSSOFT_ENGINE_DEFS_H__
