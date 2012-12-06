// GaussFIR.h: interface for the CGaussFIR class.
//
//    A Gaussian shaped complex FIR LP filter class.
//////////////////////////////////////////////////////////////////////
// Copyright 2000.    Moe Wheatley AE4JY  <ae4jy@mindspring.com>
//
// Edouard Griffiths <f4exb at free dot fr> 2012:
// - from PathSim HF path simulator project see original at 
//   http://www.moetronix.com/ae4jy/pathsim.htm
// - adapted to g++ and WSGC prototype context.
//
//This program is free software; you can redistribute it and/or
//modify it under the terms of the GNU General Public License
//as published by the Free Software Foundation; either version 2
//of the License, or any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program; if not, write to the Free Software
//Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
////////////////////////////////////////////////////////////////////////

#ifndef __GAUSSFIR_H__
#define __GAUSSFIR_H__

#include "WsgcTypes.h"
#include <math.h>


class GaussFIR  
{
public:
	void Init(wsgc_float Fs, wsgc_float F2sig );
	GaussFIR();
	virtual ~GaussFIR();
	wsgc_complex CalcFilter(wsgc_complex in);

private:
	wsgc_float* m_pCoef;
	unsigned int m_FIRlen;
	int m_FirState;
	wsgc_complex* m_pQue;

	wsgc_float dnorm(wsgc_float x, wsgc_float mu, wsgc_float sigma);
};

#endif // __GAUSSFIR_H__
