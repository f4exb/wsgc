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

     Reception MFSK
*/
#ifndef __RECEPTION_MFSK_H__
#define __RECEPTION_MFSK_H__

#include "Reception.h"

#ifdef _RSSOFT
class RSSoft_Engine;
namespace rssoft
{
    class RS_ReliabilityMatrix;
}
#endif

#ifdef _CCSOFT
class CCSoft_Engine;
namespace ccsoft
{
    class CC_ReliabilityMatrix;
}
#endif

class Reception_MFSK : public Reception
{
public:
    Reception_MFSK(Options& _options, const GoldCodeGenerator& _gc_generator);
    ~Reception_MFSK();

    void message_processing(wsgc_complex *faded_source_samples, unsigned int nb_faded_source_samples);
    void training_processing(wsgc_complex *faded_source_samples, unsigned int nb_faded_source_samples);
};

#endif // __RECEPTION_MFSK_H__

