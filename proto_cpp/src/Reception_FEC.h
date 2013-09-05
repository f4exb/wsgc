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
      
     Reception class specialized in FEC processing
*/
#ifndef __RECEPTION_FEC_H__
#define __RECEPTION_FEC_H__

class Options;

#ifdef _RSSOFT
namespace rssoft
{
class RS_ReliabilityMatrix;
}
#endif

#ifdef _CCSOFT
namespace ccsoft
{
class CC_ReliabilityMatrix;
}
#endif


class Reception_FEC
{
public:
#ifdef _RSSOFT
    static void run_rssoft_decoding(Options& options, rssoft::RS_ReliabilityMatrix& relmat);
#endif

#ifdef _CCSOFT
    static void run_ccsoft_decoding(Options& options, ccsoft::CC_ReliabilityMatrix *relmat);
#endif
};

#endif // __RECEPTION_FEC_H__
