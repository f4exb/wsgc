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
      
     Reception WSGCE: Correlation of all possible PRNs, BPSK
*/
#ifndef __RECEPTION_WSGCE_H__
#define __RECEPTION_WSGCE_H__

#include "Reception.h"
#include "CodeModulator_BPSK.h"

#ifdef _CUDA
#include "CudaManager.h"
#endif

#ifdef _CCSOFT
class CCSoft_Engine;
namespace ccsoft
{
	class CC_ReliabilityMatrix;
}
#endif

class PilotCorrelationRecord;

class Reception_WSGCE : public Reception
{
public:
    Reception_WSGCE(Options& _options, const GoldCodeGenerator& _gc_generator);
    ~Reception_WSGCE();
    
    void message_processing(wsgc_complex *faded_source_samples, unsigned int nb_faded_source_samples);
    void training_processing(wsgc_complex *faded_source_samples, unsigned int nb_faded_source_samples);

protected:
    void print_correlation_records(std::ostream& os, const std::vector<PilotCorrelationRecord>& correlation_records);
    CodeModulator_BPSK localCodeModulator;
#ifdef _CUDA
    CudaManager cuda_manager;
#endif
};


#endif // __RECEPTION_WSGC_H__
