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
#ifndef __RSSOFT_ENGINE_H__
#define __RSSOFT_ENGINE_H__

#include "GFq.h"
#include "GF2_Element.h"
#include "GF2_Polynomial.h"
#include "GF_Utils.h"
#include "EvaluationValues.h"
#include "ReliabilityMatrix.h"
#include "MultiplicityMatrix.h"
#include "GSKV_Interpolation.h"
#include "RR_Factorization.h"
#include "FinalEvaluation.h"
#include "RS_Encoding.h"

class RSSoft_Engine
{
public:
	RSSoft_Engine(unsigned int _m, unsigned int _k);

	~RSSoft_Engine();

	unsigned int get_m() const
	{
		return m;
	}

	unsigned int get_n() const
	{
		return n;
	}

	unsigned int get_q() const
	{
		return q;
	}

	unsigned int get_k() const
	{
		return k;
	}


protected:
	unsigned int m; //!< GF(2^m)
	unsigned int n; //!< 2^m-1
	unsigned int q; //!< 2^m
	unsigned int k; //!< RS(n,k), n = 2^m-1
	unsigned int nb_retries; //!< Number of soft decision retries
	unsigned int M; //!< Global multiplicity for soft decision multiplicity matrix

};



#endif // __RSSOFT_ENGINE_H__
