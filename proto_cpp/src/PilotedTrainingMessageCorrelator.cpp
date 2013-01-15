/*
     Copyright 2012 Edouard Griffiths <f4exb at free dot fr>

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

     PilotedTrainingMessageCorrelator

     Given the frequency and time displacement of correlation peak given by the
     Pilot Correlator it searches correlation for all PRNs of the training message symbols
     and accumulates results in order to find a peak of magnitude among PRNs.
     The training sequence is the sucession of message PRNs in their PRN index order for the 
     given code set. Accumulation is shifted one PRN place to the right each time (that is one 
     PRN earlier in the index sequence) so that peaks would accumulate always in the same place.
     The PRN index eventually selected gives the relative position in the whole sequence when 
     the process was started. This reveals the point in time when the first PRN in the sequence 
     was sent hence the epoch of the message.
     It uses straightforward time correlation.

*/

#include "PilotedTrainingMessageCorrelator.h"
#include "GoldCodeGenerator.h"
#include "PilotCorrelationAnalyzer.h"
#include "LocalCodes_Host.h"
#include "WsgcUtils.h"
#include <cmath>
#include <cstring>

PilotedTrainingMessageCorrelator::PilotedTrainingMessageCorrelator(
		wsgc_float f_sampling,
		wsgc_float f_chip,
		unsigned int sequence_length) :
	_maxtoavg_max(0.0),
	_mag_max(0.0),
	_prni_max(0),
	_prnai_max(0),
	_f_sampling(f_sampling),
	_f_chip(f_chip),
	_sequence_length(sequence_length),
    _delta_f(0.0)
{
}


PilotedTrainingMessageCorrelator::~PilotedTrainingMessageCorrelator()
{
}


