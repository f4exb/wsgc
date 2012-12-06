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
#include <math.h>
#include <iostream>

#define K_2PI ( 2.0 * M_PI )			// 2 Pi
#define OFFSET_FREQ_CONST (K_2PI/8000.0)	//2Pi/8000
#define KGNB 0.62665707	//equivalent Noise BW of Gaussian shaped filter

#define RATE_12_8 0		//Used for 0.1 > Spread >= 0.4
#define RATE_64 1		//Used for 0.4 > Spread >= 2.0
#define RATE_320 2		//Used for 2.0 > Spread >= 10.0


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Path::Path() :
    m_unif(-1.0,1.0),
	m_NoiseSampRate(RATE_320),
	m_Indx(0),
	m_pLPFIR(0)
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
    
	if( (m_Spread > 2.0) && (m_Spread <= 30.0) )
	{
		m_NoiseSampRate = RATE_320;
		m_pLPFIR->Init( 320.0, m_Spread );
		m_LPGain = sqrt(320.0/(4.0*m_Spread*KGNB) );
	}
	else if( (m_Spread > 0.4) && (m_Spread <= 2.0) )
	{
		m_NoiseSampRate = RATE_64;
		m_pLPFIR->Init( 64.0, m_Spread );
		m_LPGain = sqrt(64.0/(4.0*m_Spread*KGNB) );
	}
	else if( (m_Spread >= 0.1) && (m_Spread <= 0.4) )
	{
		m_NoiseSampRate = RATE_12_8;
		m_pLPFIR->Init( 12.8, m_Spread );
		m_LPGain = sqrt(12.8/(4.0*m_Spread*KGNB) );
	}
	else if( (m_Spread >= 0.0) && (m_Spread < 0.1) )
	{		//here if spread<.1 so will not use any spread just offset
		m_NoiseSampRate = RATE_320;
		m_LPGain = 1.0;
	}
    
	for(unsigned int i=0; i<INTP_QUE_SIZE; i++)
	{
		m_pQue0[i] = (0.0, 0.0);
		m_pQue1[i] = (0.0, 0.0);
		m_pQue2[i] = (0.0, 0.0);
		m_pQue3[i] = (0.0, 0.0);
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
			pOut[i] = (0.0, 0.0);
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
    
    if( m_NoiseSampRate == RATE_12_8)
    {
        if( m_Indx%(5*5*5*5) == 0 )	
        {			//generate noise samples at 12.8Hz rate
            acc = MakeGaussianDelaySample();

            //SweepGenCpx(  &acc, 12.8, 0.0, 6.4, 0.016 );

            j = m_FirState0/INTP_VALUE;
            m_pQue0[j] = acc;
        }
    }
    if( m_NoiseSampRate <= RATE_64)
    {
        if( m_Indx%(5*5*5) == 0 )	
        {
            if( m_NoiseSampRate == RATE_64)
            {			//generate noise samples at 64Hz rate
                acc = MakeGaussianDelaySample();
            }
            else
            {
                acc = (0.0, 0.0);
                Firptr = m_pQue0;
                Kptr = X5IntrpFIRCoef+INTP_FIR_SIZE-m_FirState0;
                
                for(j=0; j<INTP_QUE_SIZE; j++)
                {
                    acc += ((Firptr->real())*(*Kptr), (Firptr++->imag())*(*Kptr));
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
        if( m_NoiseSampRate == RATE_320)
        {
            acc = MakeGaussianDelaySample();
        }
        else
        {
                acc = (0.0, 0.0);
                Firptr = m_pQue1;
                Kptr = X5IntrpFIRCoef+INTP_FIR_SIZE-m_FirState1;
                
                for(j=0; j<INTP_QUE_SIZE; j++)
                {
                    acc += ((Firptr->real())*(*Kptr), (Firptr++->imag())*(*Kptr));
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
        acc = (0.0, 0.0);
        Firptr = m_pQue2;
        Kptr = X5IntrpFIRCoef+INTP_FIR_SIZE-m_FirState2;
        
        for(j=0; j<INTP_QUE_SIZE; j++)
        {
            acc += ((Firptr->real())*(*Kptr), (Firptr++->imag())*(*Kptr));
            Kptr += INTP_VALUE;
        }
        
        if( --m_FirState2 < 0)
            m_FirState2 = INTP_FIR_SIZE-1;

        //SweepGenCpx(  &acc, 1600, 0.0, 800.0, 2 );

        j = m_FirState3/INTP_VALUE;
        m_pQue3[j] = acc;
    }
    
    acc = (0.0, 0.0);
    Firptr = m_pQue3;
    Kptr = X5IntrpFIRCoef+INTP_FIR_SIZE-m_FirState3;
    
    for(j=0; j<INTP_QUE_SIZE; j++)
    {
        acc += ((Firptr->real())*(*Kptr), (Firptr++->imag())*(*Kptr));
        Kptr += INTP_VALUE;
    }
    
    if( --m_FirState3 < 0)
        m_FirState3 = INTP_FIR_SIZE-1;

    //CalcCpxSweepRMS( acc, 8000);

    if (m_Offset == 0)
    {
        *sOut = acc * (*sIn);
    }
    else
    {
        tmp = acc * (*sIn);
        offset = (cos(m_Timeinc), sin(m_Timeinc)); //Cpx multiply by offset frequency
        *sOut = offset * tmp;
        m_Timeinc += (OFFSET_FREQ_CONST*m_Offset);
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
    
	if( m_Spread >= 0.1 )
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
        
		val = (m_LPGain*u1*sqrt(-2.0*log(r)/r), m_LPGain*u2*sqrt(-2.0*log(r)/r));

        //SweepGenCpx(  &val, 320, 0.0, 30*5, 30*5/200.0);

        // Now LP filter the Gaussian samples
		val = m_pLPFIR->CalcFilter(val);
	}
	else		//if not using any spread
	{
		val = (m_LPGain, 0);
	}

    //gDebug1 = CalcCpxRMS( val, 288000);
    //CalcCpxSweepRMS( val, 500);

	return val;
}
