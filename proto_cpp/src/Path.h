// Path.h: interface for the CPath class.
//
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
//////////////////////////////////////////////////////////////////////

#ifndef __PATH_H__
#define __PATH_H__

#include "WsgcTypes.h"
#include "GaussFIR.h"
#include "FilterTables.h"
#include <math.h>
#include <tr1/random> // c++0x extension (needs -std=c++0x to compile)

class Path  
{
public:
	void CalcPath(const wsgc_complex* pIn, wsgc_complex* pOut);
	void CalcPathSample(const wsgc_complex* sIn, wsgc_complex* sOut);
	void InitPath( wsgc_float Spread, wsgc_float Offset, unsigned int blocksize, unsigned int numpaths, bool active);
	Path();
	~Path();

private:
    typedef std::tr1::ranlux64_base_01 RandomEngine; 
    RandomEngine m_randomEngine;
    std::tr1::uniform_real<wsgc_float> m_unif;
	unsigned m_IIRLength;
	wsgc_complex MakeGaussianDelaySample();
	unsigned m_NoiseSampRate;
	bool m_PathActive;
	int m_inc;
	int m_BlockSize;
	int m_Indx;
	wsgc_float m_Offset;
	wsgc_float m_Spread;
	wsgc_float m_LPGain;
	wsgc_float m_Timeinc;
	wsgc_complex m_pQue0[INTP_QUE_SIZE];
	wsgc_complex m_pQue1[INTP_QUE_SIZE];
	wsgc_complex m_pQue2[INTP_QUE_SIZE];
	wsgc_complex m_pQue3[INTP_QUE_SIZE];
	int m_FirState0;
	int m_FirState1;
	int m_FirState2;
	int m_FirState3;
	GaussFIR* m_pLPFIR;
	bool m_noSpread;
};

#endif // __PATH_H__
