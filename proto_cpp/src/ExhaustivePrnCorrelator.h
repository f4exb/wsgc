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
#ifndef __EXHAUSTIVE_PRN_CORRELATOR_H__
#define __EXHAUSTIVE_PRN_CORRELATOR_H__

#include "WsgcTypes.h"
#include <vector>

class LocalCodes;
class MultiplePrnCorrelator_FreqDep;

#ifdef _CCSOFT
namespace ccsoft
{
    class CC_ReliabilityMatrix;
}
#endif
class PilotCorrelationRecord;

class ExhaustivePrnCorrelator
{
public:
	ExhaustivePrnCorrelator(const LocalCodes *_local_codes,
			MultiplePrnCorrelator_FreqDep *_ifft_correlator);

	~ExhaustivePrnCorrelator();

#ifdef _CCSOFT
	void make_correlation(wsgc_complex *one_prn_samples, ccsoft::CC_ReliabilityMatrix& relmat);
#endif
    void make_correlation(wsgc_complex *one_prn_samples, std::vector<PilotCorrelationRecord>& correlation_records);

protected:
#ifdef _CCSOFT
    void update_reliability_matrix(ccsoft::CC_ReliabilityMatrix& relmat);
#endif
    void update_correlation_records(std::vector<PilotCorrelationRecord>& correlation_records);
    
	const LocalCodes *local_codes; //!< Pointer reference to the local codes
	MultiplePrnCorrelator_FreqDep *ifft_correlator; //!< Frequency and time domain for one PRN correlator
    unsigned int prn_count; //!< Count of PRNs processed. Ever increasing.
#ifdef _CCSOFT
    float *relmat_column;              //!< Reliability matrix column temporary storage
#endif
};


#endif // __EXHAUSTIVE_PRN_CORRELATOR_H__
