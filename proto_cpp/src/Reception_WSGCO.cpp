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
      
     Reception WSGCO: Correlation, OOK
*/

#include "Reception_FEC.h"
#include "Demodulation.h"
#include "Reception_WSGCO.h"
#include "SampleSequencer.h"
#include "LocalCodesFFT_Host.h"
#include "UnpilotedMultiplePrnCorrelator_Host.h"
#include "WsgcUtils.h"

#include <iostream>
#include <iomanip>

#ifdef _RSSOFT
#include "RS_ReliabilityMatrix.h"
#endif

#ifdef _CCSOFT
#include "CC_ReliabilityMatrix.h"
#endif

//=================================================================================================
Reception_WSGCO::Reception_WSGCO(Options& _options, const GoldCodeGenerator& _gc_generator) :
    Reception(_options, _gc_generator)
{}

//=================================================================================================
Reception_WSGCO::~Reception_WSGCO()
{}

//=================================================================================================
void Reception_WSGCO::message_processing(wsgc_complex *faded_source_samples, unsigned int nb_faded_source_samples)
{
    std::vector<CorrelationRecord> correlation_records;

    std::cout << "Unpiloted correlation with frequency independant modulations" << std::endl;
    unpiloted_message_correlation(localCodeModulator, faded_source_samples, nb_faded_source_samples, correlation_records);

    std::cout << "Do the decoding with the decision box..." << std::endl;
    unpiloted_and_synced_decision(correlation_records);
}

//=================================================================================================
void Reception_WSGCO::training_processing(wsgc_complex *faded_source_samples, unsigned int nb_faded_source_samples)
{
    std::vector<TrainingCorrelationRecord> training_correlation_records;

    unpiloted_training_correlation(localCodeModulator, faded_source_samples, nb_faded_source_samples, training_correlation_records);

    std::ostringstream corr_os;
    corr_os << "--- training correlation records:" << std::endl;
    TrainingCorrelationRecord::dump_banner(corr_os);

    for (std::vector<TrainingCorrelationRecord>::const_iterator it = training_correlation_records.begin(); it != training_correlation_records.end(); ++it)
    {
        it->dump_line(corr_os);
    }

    std::cout << corr_os.str() << std::endl;
}

