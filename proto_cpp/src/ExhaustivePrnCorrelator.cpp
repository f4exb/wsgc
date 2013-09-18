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

     ExhaustivePrnCorrelator

     Search all possible PRNs in frequency and time space

*/

#include "ExhaustivePrnCorrelator.h"
#include "MultiplePrnCorrelator_FreqDep.h"

#ifdef _CCSOFT
#include "CC_ReliabilityMatrix.h"
#endif
#include "PilotCorrelationRecord.h"

//=================================================================================================
ExhaustivePrnCorrelator::ExhaustivePrnCorrelator(const LocalCodes *_local_codes,
		MultiplePrnCorrelator_FreqDep *_ifft_correlator) :
	local_codes(_local_codes),
	ifft_correlator(_ifft_correlator),
    prn_count(0)
{
#ifdef _CCSOFT
    relmat_column = new float[ifft_correlator->get_nb_message_symbols()];
#endif
}

//=================================================================================================
ExhaustivePrnCorrelator::~ExhaustivePrnCorrelator()
{
#ifdef _CCSOFT
    delete[] relmat_column;
#endif
}

//=================================================================================================
#ifdef _CCSOFT
void ExhaustivePrnCorrelator::make_correlation(wsgc_complex *one_prn_samples, ccsoft::CC_ReliabilityMatrix& relmat)
{
    unsigned int batch_number = prn_count / ifft_correlator->get_nb_batch_prns();
    unsigned int prn_position = prn_count % ifft_correlator->get_storage_depth();
    
    ifft_correlator->multiply_and_ifft(one_prn_samples, prn_position);

    // if a batch is ready calculate average. There is one batch delay before average can be calculated so results of batch #0 should be skipped (???)
    //if (((prn_count % ifft_correlator->get_nb_batch_prns()) == (ifft_correlator->get_nb_batch_prns() - 1)) && (batch_number > 0))
    if ((prn_count % ifft_correlator->get_nb_batch_prns()) == (ifft_correlator->get_nb_batch_prns() - 1))
    {
        bool even_batch = (batch_number % 2 == 0);
        ifft_correlator->execute_averaging(!even_batch); // because of the one batch delay, the first half is processed on odd batches
        update_reliability_matrix(relmat); // Update reliability matrix
    }
    
    prn_count++;
}
#endif

//=================================================================================================
void ExhaustivePrnCorrelator::make_correlation(wsgc_complex *one_prn_samples, std::vector<PilotCorrelationRecord>& correlation_records)
{
    unsigned int batch_number = prn_count / ifft_correlator->get_nb_batch_prns();
    unsigned int prn_position = prn_count % ifft_correlator->get_storage_depth();
    
    ifft_correlator->multiply_and_ifft(one_prn_samples, prn_position);

    // if a batch is ready calculate average. There is one batch delay before average can be calculated so results of batch #0 should be skipped (???)
    //if (((prn_count % ifft_correlator->get_nb_batch_prns()) == (ifft_correlator->get_nb_batch_prns() - 1)) && (batch_number > 0))
    if ((prn_count % ifft_correlator->get_nb_batch_prns()) == (ifft_correlator->get_nb_batch_prns() - 1))
    {
        bool even_batch = (batch_number % 2 == 0);
        ifft_correlator->execute_averaging(!even_batch); // Because of the one batch delay, the first half is processed on odd batches
        update_correlation_records(correlation_records); // Update correlation records
    }

    prn_count++;
}

//=================================================================================================
#ifdef _CCSOFT
void ExhaustivePrnCorrelator::update_reliability_matrix(ccsoft::CC_ReliabilityMatrix& relmat)
{
    unsigned int prn_per_symbol = ifft_correlator->get_prn_per_symbol();
    for (unsigned int pi=0; pi < ifft_correlator->get_nb_batch_prns(); pi++)
    {
        if ((pi % prn_per_symbol) == (prn_per_symbol - 1))
        {
            for (unsigned int prni=0; prni<ifft_correlator->get_nb_message_symbols(); prni++)
            {
                relmat_column[prni] = (ifft_correlator->get_batch_magnitudes_max(prni))[pi];
            }
            
            relmat.enter_symbol_data(relmat_column);
        }
    }
}
#endif

//=================================================================================================
void ExhaustivePrnCorrelator::update_correlation_records(std::vector<PilotCorrelationRecord>& correlation_records)
{
    for (unsigned int pi=0; pi < ifft_correlator->get_nb_batch_prns(); pi++)
    {
        wsgc_float max_mag = 0.0;
        unsigned int prni_max = 0;
        
        for (unsigned int prni=0; prni<ifft_correlator->get_nb_message_symbols(); prni++)
        {
            if ((ifft_correlator->get_batch_magnitudes_max(prni))[pi] > max_mag)
            {
                max_mag = (ifft_correlator->get_batch_magnitudes_max(prni))[pi];
                prni_max = prni;
            }
        }
        
        static const PilotCorrelationRecord tmp_correlation_record;
        correlation_records.push_back(tmp_correlation_record);
        PilotCorrelationRecord& correlation_record = correlation_records.back();

        unsigned int ifft_peak_index = (ifft_correlator->get_batch_composite_indexes_max(prni_max))[pi];
        unsigned int batch_index = ifft_correlator->get_batch_index() - 1;
        unsigned int global_prn_index = batch_index*ifft_correlator->get_nb_batch_prns() + pi;
        unsigned int global_averaging_block_index = global_prn_index / ifft_correlator->get_prn_per_symbol();
        unsigned int prn_in_symbol_index = global_prn_index % ifft_correlator->get_prn_per_symbol();
        unsigned int f_index, fs_index, t_index;
        ifft_correlator->calculate_indexes(ifft_peak_index, f_index, fs_index, t_index);
        
        correlation_record.block_count = global_averaging_block_index;
        correlation_record.prn_index = prn_in_symbol_index;
        correlation_record.pilot_index = prni_max;
        correlation_record.magnitude_max = (ifft_correlator->get_batch_magnitudes_max(prni_max))[pi] / ifft_correlator->get_fft_N(); // divide by FFT scaling factor
        correlation_record.phase_at_max = atan2((ifft_correlator->get_batch_complex_values_max(prni_max))[pi].imag(), (ifft_correlator->get_batch_complex_values_max(prni_max))[pi].real());
        correlation_record.f_index_max = f_index * ifft_correlator->get_freq_step_division() + fs_index;
        correlation_record.t_index_max = t_index;
        // frequency shift to get back to zero IF = relative step displacement / step size. Step size is the size of a FFT bin.
        wsgc_float f_bin_size = ifft_correlator->get_f_sampling() / ifft_correlator->get_fft_N();
        correlation_record.delta_f = ((((wsgc_float) correlation_record.f_index_max) / ifft_correlator->get_freq_step_division()) - (ifft_correlator->get_nb_f_bins()/2)) * f_bin_size;
    }
}
