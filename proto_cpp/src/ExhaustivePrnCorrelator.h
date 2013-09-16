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
class SinglePrnCorrelator_FreqDep;

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
			const SinglePrnCorrelator_FreqDep *_ifft_correlator);

	~ExhaustivePrnCorrelator();

#ifdef _CCSOFT
	void make_correlation(wsgc_complex *one_prn_samples, ccsoft::CC_ReliabilityMatrix& relmat);
#endif
    void make_correlation(wsgc_complex *one_prn_samples, std::vector<PilotCorrelationRecord>& correlation_records);

protected:
	const LocalCodes *local_codes; //!< Pointer reference to the local codes
	const SinglePrnCorrelator_FreqDep *ifft_correlator; //!< Frequency and time domain for one PRN correlator
};


#endif // __EXHAUSTIVE_PRN_CORRELATOR_H__
