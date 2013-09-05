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
      
     Reception WSGC: Correlation, BPSK, Pilot
*/
#ifndef __RECEPTION_WSGC_H__
#define __RECEPTION_WSGC_H__

#include "Reception.h"
#include "CodeModulator_BPSK.h"

#ifdef _CUDA
#include "CudaManager.h"
#endif

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

class Reception_WSGC : public Reception
{
public:
    Reception_WSGC(Options& _options, const GoldCodeGenerator& _gc_generator);
    ~Reception_WSGC();
    
    void message_processing(wsgc_complex *faded_source_samples, unsigned int nb_faded_source_samples);
    void training_processing(wsgc_complex *faded_source_samples, unsigned int nb_faded_source_samples);

protected:

#ifdef _RSSOFT
    void run_rssoft_decoding(RSSoft_Engine& rssoft_engine, rssoft::RS_ReliabilityMatrix& rs_relmat);
#endif

#ifdef _CCSOFT
    void run_ccsoft_decoding(CCSoft_Engine& ccsoft_engine, ccsoft::CC_ReliabilityMatrix *relmat);
#endif

    CodeModulator_BPSK localCodeModulator;
#ifdef _CUDA
    CudaManager cuda_manager;
#endif
};


#endif // __RECEPTION_WSGC_H__
