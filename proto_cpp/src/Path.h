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
	/**
	 * Creates a new path object
	 * \param fSampling Main sampling frequency
	 */
	Path(wsgc_float fSampling);
	~Path();

	/**
	 * Calculation for a block of samples (see blocksize parameter of InitPath) through the path
	 */

	void CalcPath(const wsgc_complex* pIn, wsgc_complex* pOut);
	/**
	 * Calculation for one sample through the path
	 */
	void CalcPathSample(const wsgc_complex* sIn, wsgc_complex* sOut);

	/**
	 * Initializes the path
	 * \param Spread Doppler frequency spread
	 * \param Offset Frequency offset
	 * \param blocksize Block size for block processing (see CalcPath)
	 * \param numpaths Total number of paths
	 * \param active True if the path is active
	 */
	void InitPath( wsgc_float Spread, wsgc_float Offset, unsigned int blocksize, unsigned int numpaths, bool active);

private:
    typedef std::tr1::ranlux64_base_01 RandomEngine; 
    RandomEngine m_randomEngine;
    std::tr1::uniform_real<wsgc_float> m_unif;
	wsgc_float m_fSampling; //!< Main sampling frequency
	unsigned m_NoiseSampRate; //!< Doppler Gaussian noise sample rate
	bool m_PathActive; //!< True if this path is active
	int m_inc;
	int m_BlockSize; //!< Block size for samples block processing
	int m_Indx;
	wsgc_float m_Offset; //!< Offset frequency
	wsgc_float m_Spread; //!< Doppler spread frequency
	wsgc_float m_LPGain; //!< LP filter gain
	wsgc_float m_Timeinc;
	wsgc_complex m_pQue0[INTP_QUE_SIZE];
	wsgc_complex m_pQue1[INTP_QUE_SIZE];
	wsgc_complex m_pQue2[INTP_QUE_SIZE];
	wsgc_complex m_pQue3[INTP_QUE_SIZE];
	int m_FirState0;
	int m_FirState1;
	int m_FirState2;
	int m_FirState3;
	GaussFIR* m_pLPFIR; //!< Gaussian LP filter to filter Doppler noise
	bool m_noSpread; //!< True if no Doppler spread is applied
	wsgc_float m_OffsetFreqConst; //!< 2 PI / sampling frequency
	wsgc_float m_SpreadLimit0; //!< Doppler spread frequency ^2 interpolation upper limit (above which an exception is thrown)
	wsgc_float m_SpreadLimit1; //!< Doppler spread frequency ^3 interpolation upper limit
	wsgc_float m_SpreadLimit2; //!< Doppler spread frequency ^4 interpolation upper limit
	wsgc_float m_SpreadLimit3; //!< Doppler spread frequency lower limit (below which no Doppler spread is applied)
	wsgc_float m_NoiseFs_2; //!< Sampling frequency for ^2 interpolation
	wsgc_float m_NoiseFs_3; //!< Sampling frequency for ^3 interpolation
	wsgc_float m_NoiseFs_4; //!< Sampling frequency for ^4 interpolation
	static const wsgc_float rate_2; //!< ^2 interpolation rate
	static const wsgc_float rate_3; //!< ^3 interpolation rate
	static const wsgc_float rate_4; //!< ^4 interpolation rate

	wsgc_complex MakeGaussianDelaySample();
};

#endif // __PATH_H__
