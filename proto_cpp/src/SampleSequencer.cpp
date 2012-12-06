/*
 * SampleSequencer.cpp
 *
 *  Created on: Nov 8, 2012
 *      Author: doudou
 */

#include "SampleSequencer.h"


SampleSequencer::SampleSequencer(wsgc_complex *samples, unsigned int nb_samples, unsigned int nb_code_samples) :
	_samples(samples),
	_nb_samples(nb_samples),
	_nb_code_samples(nb_code_samples),
	_serviced_samples_index(0)
{}


bool SampleSequencer::get_next_code_samples(wsgc_complex **samples)
{
    if (_serviced_samples_index + _nb_code_samples < _nb_samples)
    {
        *samples = &_samples[_serviced_samples_index];
        _serviced_samples_index += _nb_code_samples;
        return true;
    }
    else
    {
        return false;
    }
}
