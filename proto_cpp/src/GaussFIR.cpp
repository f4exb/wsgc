// GaussFIR.cpp: implementation of the CGaussFIR class.
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

#include "GaussFIR.h"
#include "WsgcException.h"
#include <iostream>

#define SQRT2PI (2.506628275)
#define PI2 ( 2.0 * M_PI )		// 2 Pi
#define SQRT2 (1.414213562)
#define K_GAUSSIAN 1.4	//constant determines accuracy of Gaussian LPFIR.
						// larger number makes FIR longer but more accurate
						// This value is good to about -50 dB

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

GaussFIR::GaussFIR()
{
	m_pCoef = 0;
	m_pQue = 0;
}

GaussFIR::~GaussFIR()
{
	if(m_pCoef)
		delete m_pCoef;
	if(m_pQue)
		delete m_pQue;
}

////////////////////////////////////////////////////////////////////////////
//  This function calculates the length and FIR coefficients for a Gaussian
//		shaped LP filter and also initializes the memory for it.
////////////////////////////////////////////////////////////////////////////
void GaussFIR::Init(wsgc_float Fs, wsgc_float F2sig)
{
    int indx;
    wsgc_float sigma = (Fs*SQRT2)/(PI2*F2sig);		// calculate sigma for coefficients
    
	m_FIRlen = (unsigned int)(K_GAUSSIAN*Fs/F2sig);		//calculate FIR length
    
	if( !(m_FIRlen & 1) ) //make FIR length ODD
    {
		m_FIRlen++;
    }
        
    //allocate buffer and Coefficient memory based on calculated FIR length
	if( (m_pCoef = new wsgc_float[sizeof(wsgc_float)*((m_FIRlen*2)+10)] ) == 0)
    {
		throw WsgcException("GaussFIR Error: cannot allocate coefficient memory");
    }
        
	if( (m_pQue = new wsgc_complex[sizeof(wsgc_complex)*(m_FIRlen+10)] ) == 0)
    {
		throw WsgcException("GaussFIR Error: cannot allocate buffer memory");;
    }
        
    // generate the scaled Gaussian shaped impulse response	to create a 0 dB
    //   passband LP filter with a 2 Sigma frequency bandwidth.
	indx = -((m_FIRlen-1)/2);

	for(unsigned int i=0; i<m_FIRlen;i++, indx++)
	{
		m_pCoef[i] = ( 1.0/(SQRT2PI*sigma) )*dnorm( indx,0.0, sigma)
									/ dnorm( 0.0,0.0, sigma);
		m_pQue[i] = wsgc_complex(0.0, 0.0);
		m_pCoef[i+m_FIRlen] = m_pCoef[i];	//make duplicate for flat FIR
	}
    
	m_FirState = m_FIRlen-1;		//used for flat FIR implementation
}

//////////////////////////////////////////////////////////////////////
//  Calculate complex Gaussian FIR filter iteration for one sample
//////////////////////////////////////////////////////////////////////
wsgc_complex GaussFIR::CalcFilter(wsgc_complex in)
{
    wsgc_complex acc;
    wsgc_complex* Firptr;
    wsgc_float* Kptr;
    
	m_pQue[m_FirState] = in;
	Firptr = m_pQue;
	acc = wsgc_complex(0.0, 0.0);
	Kptr = m_pCoef+m_FIRlen-m_FirState;
    
	for(unsigned int i=0; i<m_FIRlen; i++)
	{
		acc.real() += (Firptr->real())*(*Kptr);
		acc.imag() += (Firptr++->imag())*(*Kptr);
        Kptr++;
	}
    
	if( --m_FirState < 0)
		m_FirState += m_FIRlen;
        
	return acc;
}

//////////////////////////////////////////////////////////////////////
//  implements a Gaussian (Normal) distribution function.
//////////////////////////////////////////////////////////////////////
wsgc_float GaussFIR::dnorm(wsgc_float x, wsgc_float mu, wsgc_float sigma)
{
	return( 1.0/(SQRT2PI*sigma) )*exp( (-1.0/(2.0*sigma*sigma))*(x-mu)*(x-mu) );
}

