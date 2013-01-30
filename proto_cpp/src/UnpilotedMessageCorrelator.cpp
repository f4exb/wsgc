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

     UnpilotedMessageCorrelator

	 Does correlation on all possible message PRNs (normally also including noise PRN)
	 to find a correlation peak representing the time delay of the start of sequence.

	 Uses frequency domain correlation to check all possible delays at once which is
	 much more computationnaly efficient

	 This is the super abstract class for Host and CUDA implementations

*/

#include "UnpilotedMessageCorrelator.h"

UnpilotedMessageCorrelator::UnpilotedMessageCorrelator(
			wsgc_float f_sampling,
			wsgc_float f_chip,
			bool simulate_symbol_synchronization,
			unsigned int prn_per_symbol,
			unsigned int nb_batch_prns
			) :
		_f_sampling(f_sampling),
		_f_chip(f_chip),
		_prn_per_symbol(prn_per_symbol),
		_nb_batch_prns((simulate_symbol_synchronization ? nb_batch_prns*prn_per_symbol : nb_batch_prns)),
		_simulate_symbol_synchronization(simulate_symbol_synchronization),
	    _batch_max_magnitudes(_nb_batch_prns,0.0),
	    _batch_max_ifft_indexes(_nb_batch_prns,0),
	    _batch_max_prn_indexes(_nb_batch_prns,0),
	    _batch_complex_values_max(_nb_batch_prns, (1.0, 0.0))
{}


UnpilotedMessageCorrelator::~UnpilotedMessageCorrelator()
{}

