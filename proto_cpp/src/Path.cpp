// Path.cpp: implementation of the CPath class.
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
////////////////////////////////////////////////////////////////////////

#include "Path.h"
#include "WsgcException.h"
#include <math.h>
#include <iostream>
#include <sstream>

#define K_2PI ( 2.0 * M_PI )			// 2 Pi
#define KGNB 0.62665707	//equivalent Noise BW of Gaussian shaped filter

#define RATE_5_5_5_5 0	// Noise interpolation rate = 5^4. Used for 0.01 =< Spread <= 0.4   (based on 8000 Hz sampling frequency)
#define RATE_5_5_5   1	// Noise interpolation rate = 5^3. Used for 0.4   < Spread <= 2.0
#define RATE_5_5     2	// Noise interpolation rate = 5^2. Used for 2.0   < Spread <= 40.0


const wsgc_float Path::rate_2 = 25.0;  // 5^2
const wsgc_float Path::rate_3 = 125.0; // 5^3
const wsgc_float Path::rate_4 = 625.0; // 5^4

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Path::Path(wsgc_float fSampling) :
    m_unif(-1.0,1.0),
    m_fSampling(fSampling),
	m_NoiseSampRate(RATE_5_5),
	m_Indx(0),
	m_pLPFIR(0),
	m_BlockSize(1),
	m_PathActive(false),
	m_inc(0),
	m_Offset(0.0),
	m_Spread(0.0),
	m_LPGain(1.0),
	m_Timeinc(0.0),
	m_FirState0(INTP_QUE_SIZE-1),
	m_FirState1(INTP_QUE_SIZE-1),
	m_FirState2(INTP_QUE_SIZE-1),
	m_FirState3(INTP_QUE_SIZE-1),
	m_noSpread(true),
	m_OffsetFreqConst(K_2PI/fSampling),
	m_SpreadLimit0(fSampling/(10.0*rate_2)), // 200
	m_SpreadLimit1(fSampling/(10.0*rate_3)), // 4000
	m_SpreadLimit2(fSampling/(10.0*rate_4)), // 20000
	m_SpreadLimit3(0.01),
	m_NoiseFs_2(fSampling/rate_2),
	m_NoiseFs_3(fSampling/rate_3),
	m_NoiseFs_4(fSampling/rate_4)
{
    m_randomEngine.seed(time(0));
}

Path::~Path()
{
	if(m_pLPFIR)
		delete m_pLPFIR;
}


//////////////////////////////////////////////////////////////////////
//  Initialize a path
//////////////////////////////////////////////////////////////////////

void Path::InitPath( wsgc_float Spread, wsgc_float Offset, unsigned int blocksize, unsigned int numpaths, bool active)
{
	m_BlockSize = blocksize;
	m_Offset = Offset;
	m_Spread = Spread;
	m_PathActive = active;
	m_FirState0 = INTP_QUE_SIZE-1;
	m_FirState1 = INTP_QUE_SIZE-1;
	m_FirState2 = INTP_QUE_SIZE-1;
	m_FirState3 = INTP_QUE_SIZE-1;
	m_Indx = 0;
	m_inc = 0;
	m_Timeinc = 0.0;
	m_pLPFIR = new GaussFIR;
	m_noSpread = false;
    
	if( (m_Spread > m_SpreadLimit1) && (m_Spread <= m_SpreadLimit0) )
	{
		m_NoiseSampRate = RATE_5_5;
		m_pLPFIR->Init( m_NoiseFs_2, m_Spread );
		m_LPGain = sqrt(m_NoiseFs_2/(4.0*m_Spread*KGNB) );
	}
	else if( (m_Spread > m_SpreadLimit2) && (m_Spread <= m_SpreadLimit1) )
	{
		m_NoiseSampRate = RATE_5_5_5;
		m_pLPFIR->Init( m_NoiseFs_3, m_Spread );
		m_LPGain = sqrt(m_NoiseFs_3/(4.0*m_Spread*KGNB) );
	}
	else if( (m_Spread >= m_SpreadLimit3) && (m_Spread <= m_SpreadLimit2) )
	{
		m_NoiseSampRate = RATE_5_5_5_5;
		m_pLPFIR->Init( m_NoiseFs_4, m_Spread );
		m_LPGain = sqrt(m_NoiseFs_4/(4.0*m_Spread*KGNB) );
	}
	else if( (m_Spread >= 0.0) && (m_Spread < m_SpreadLimit3) )
	{		//here if spread<.01 so will not use any spread just offset
		m_noSpread = true;
		m_NoiseSampRate = RATE_5_5;
		m_LPGain = 1.0;
	}
	else
	{
		std::ostringstream es;
		es << "Path: Error frequency spread over the limit of " << m_SpreadLimit0 << "Hz according to sample rate " << m_fSampling << "Hz";
		throw WsgcException(es.str());
	}
    
	for(unsigned int i=0; i<INTP_QUE_SIZE; i++)
	{
		m_pQue0[i] = wsgc_complex(0.0, 0.0);
		m_pQue1[i] = wsgc_complex(0.0, 0.0);
		m_pQue2[i] = wsgc_complex(0.0, 0.0);
		m_pQue3[i] = wsgc_complex(0.0, 0.0);
	}
    
	m_LPGain = m_LPGain/ sqrt((wsgc_float)numpaths);
    
	for(unsigned int i=0; i<250; i++)
    {
		MakeGaussianDelaySample();		//pre load filter
    }
}

//////////////////////////////////////////////////////////////////////
// Performs a path calculation on pIn and puts it in pOut
//////////////////////////////////////////////////////////////////////
void Path::CalcPath(const wsgc_complex *pIn, wsgc_complex *pOut)
{
	if(m_PathActive)		// if this path is active
	{
		for(unsigned int i=0; i<m_BlockSize; i++)
		{
            CalcPathSample(&pIn[i], &pOut[i]);
        }        
    } 
   	else		// if path is not active just zero the output
	{
		for(unsigned int i=0; i<m_BlockSize; i++)
		{
			pOut[i] = wsgc_complex(0.0, 0.0);
		}
	}

}


//////////////////////////////////////////////////////////////////////
// Performs a path calculation on one sample. Reads it from sIn 
// and puts it in sOut
//
//  Two Low Pass filtered Gaussian random numbers are created at
//	12.8, 64 Hz, or 320 Hz rate.  These form the input to a complex
//	interpolation filter that bumps the sample rate up to 8000Hz.
//
//	Two, three, or four stages of X5 upsampling/interpolation are used.
//	The complex noise is then multiplied by the input I/Q signal
//	to produce the spreading/fading simulation.
//
//  Finally a complex NCO is multiplied by the signal to produce a 
//	Frequency offset.
//////////////////////////////////////////////////////////////////////

void Path::CalcPathSample(const wsgc_complex *sIn, wsgc_complex *sOut)
{
    unsigned int j;
    wsgc_complex acc;
    wsgc_complex tmp;
    const wsgc_float* Kptr;
    wsgc_complex* Firptr;
    wsgc_complex offset;

    if (m_noSpread)
    {
    	acc = wsgc_complex(1.0, 0.0);
    }
    else
    {
		if( m_NoiseSampRate == RATE_5_5_5_5)
		{
			if( m_Indx%(5*5*5*5) == 0 )
			{			//generate noise samples at 12.8Hz rate
				acc = MakeGaussianDelaySample();

				//SweepGenCpx(  &acc, 12.8, 0.0, 6.4, 0.016 );

				j = m_FirState0/INTP_VALUE;
				m_pQue0[j] = acc;
			}
		}
		if( m_NoiseSampRate <= RATE_5_5_5)
		{
			if( m_Indx%(5*5*5) == 0 )
			{
				if( m_NoiseSampRate == RATE_5_5_5)
				{			//generate noise samples at 64Hz rate
					acc = MakeGaussianDelaySample();
				}
				else
				{
					acc = wsgc_complex(0.0, 0.0);
					Firptr = m_pQue0;
					Kptr = X5IntrpFIRCoef+INTP_FIR_SIZE-m_FirState0;

					for(j=0; j<INTP_QUE_SIZE; j++)
					{
						acc.real() += (Firptr->real())*(*Kptr);
						acc.imag() += (Firptr++->imag())*(*Kptr);
						//acc += ((Firptr->real())*(*Kptr), (Firptr++->imag())*(*Kptr));
						Kptr += INTP_VALUE;
					}

					if( --m_FirState0 < 0)
						m_FirState0 = INTP_FIR_SIZE-1;
				}

				//SweepGenCpx(  &acc, 64, 0.0, 32.0, 0.08 );

				j = m_FirState1/INTP_VALUE;
				m_pQue1[j] = acc;
			}
		}
		if( m_Indx%(5*5) == 0 )	//interpolate/upsample x5
		{
			if( m_NoiseSampRate == RATE_5_5)
			{
				acc = MakeGaussianDelaySample();
			}
			else
			{
					acc = wsgc_complex(0.0, 0.0);
					Firptr = m_pQue1;
					Kptr = X5IntrpFIRCoef+INTP_FIR_SIZE-m_FirState1;

					for(j=0; j<INTP_QUE_SIZE; j++)
					{
						acc.real() += (Firptr->real())*(*Kptr);
						acc.imag() += (Firptr++->imag())*(*Kptr);
						//acc += ((Firptr->real())*(*Kptr), (Firptr++->imag())*(*Kptr));
						Kptr += INTP_VALUE;
					}

					if( --m_FirState1 < 0)
						m_FirState1 = INTP_FIR_SIZE-1;
			}

			//SweepGenCpx(  &acc, 320, 0.0, 160.0, 0.4 );

			j = m_FirState2/INTP_VALUE;
			m_pQue2[j] = acc;
		}
		if( m_Indx%(5) == 0 )	//interpolate/upsample x5
		{
			acc = wsgc_complex(0.0, 0.0);
			Firptr = m_pQue2;
			Kptr = X5IntrpFIRCoef+INTP_FIR_SIZE-m_FirState2;

			for(j=0; j<INTP_QUE_SIZE; j++)
			{
				acc.real() += (Firptr->real())*(*Kptr);
				acc.imag() += (Firptr++->imag())*(*Kptr);
				//acc += ((Firptr->real())*(*Kptr), (Firptr++->imag())*(*Kptr));
				Kptr += INTP_VALUE;
			}

			if( --m_FirState2 < 0)
				m_FirState2 = INTP_FIR_SIZE-1;

			//SweepGenCpx(  &acc, 1600, 0.0, 800.0, 2 );

			j = m_FirState3/INTP_VALUE;
			m_pQue3[j] = acc;
		}

		acc = wsgc_complex(0.0, 0.0);
		Firptr = m_pQue3;
		Kptr = X5IntrpFIRCoef+INTP_FIR_SIZE-m_FirState3;

		for(j=0; j<INTP_QUE_SIZE; j++)
		{
			acc.real() += (Firptr->real())*(*Kptr);
			acc.imag() += (Firptr++->imag())*(*Kptr);
			//acc += ((Firptr->real())*(*Kptr), (Firptr++->imag())*(*Kptr));
			Kptr += INTP_VALUE;
		}

		if( --m_FirState3 < 0)
			m_FirState3 = INTP_FIR_SIZE-1;

		//CalcCpxSweepRMS( acc, 8000);
    }

    if (m_Offset == 0)
    {
        *sOut = acc * (*sIn);
    }
    else
    {
        tmp = acc * (*sIn);
        offset = (cos(m_Timeinc), sin(m_Timeinc)); //Cpx multiply by offset frequency
        *sOut = offset * tmp;
        m_Timeinc += (m_OffsetFreqConst*m_Offset);
        double intpart;
        m_Timeinc = K_2PI * modf(m_Timeinc/K_2PI, &intpart); //keep radian counter bounded
    }
    
    if( ++m_Indx > (INTP_VALUE*INTP_VALUE*INTP_VALUE*INTP_VALUE*m_BlockSize) )
        m_Indx = 0;
}


/////////////////////////////////////////////////////////////////
//  Create the complex Rayleigh distributed samples by
//	creating two Gaussian random distributed numbers for the I and Q
//	terms and then passing them through a Gaussian shaped LP IIR.
//	The 2 Sigma bandwidth of the LP filter determines the amount of spread.
/////////////////////////////////////////////////////////////////

wsgc_complex Path::MakeGaussianDelaySample()
{
    wsgc_float u1;
    wsgc_float u2;
    wsgc_float r;
    wsgc_complex val;
    
	if (m_noSpread) //if not using any spread
	{
		val = (m_LPGain, 0);
	}
	else
	{
        // Generate two uniform random numbers between -1 and +1
        // that are inside the unit circle
        
		do {
			//u1 = 1.0 - 2.0 * (wsgc_float)rand()/(wsgc_float)RAND_MAX ;
			//u2 = 1.0 - 2.0 * (wsgc_float)rand()/(wsgc_float)RAND_MAX ;
   			u1 = m_unif(m_randomEngine);
   			u2 = m_unif(m_randomEngine);
			r = u1*u1 + u2*u2;
		} while(r >= 1.0 || r == 0.0);
        
		val.real() = m_LPGain*u1*sqrt(-2.0*log(r)/r);
		val.imag() = m_LPGain*u2*sqrt(-2.0*log(r)/r);

        //SweepGenCpx(  &val, 320, 0.0, 30*5, 30*5/200.0);

        // Now LP filter the Gaussian samples
		val = m_pLPFIR->CalcFilter(val);
	}

    //gDebug1 = CalcCpxRMS( val, 288000);
    //CalcCpxSweepRMS( val, 500);

	return val;
}
