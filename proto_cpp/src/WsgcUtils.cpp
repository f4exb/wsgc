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
      
     WSGC Utils

     Various utilities
     
*/
#include "WsgcUtils.h"
#include <iostream>
#include <algorithm>
#include <ctime>

// magnitude estimation constants, minimizes average error
const wsgc_float WsgcUtils::magnitude_estimation_alpha = 0.948059448969;
const wsgc_float WsgcUtils::magnitude_estimation_beta  = 0.392699081699;


// time difference in seconds
double WsgcUtils::get_time_difference(const timespec& time1, const timespec& time2)
{
    long long unsigned int time1_ns = time1.tv_sec * 1000000000ull + time1.tv_nsec;
    long long unsigned int time2_ns = time2.tv_sec * 1000000000ull + time2.tv_nsec;
    
    if (time1_ns > time2_ns)
    {
        return ((double) time1_ns - time2_ns) / 1e9;
    }
    else
    {
        return ((double) time2_ns - time1_ns) / 1e9;
    }
}


void WsgcUtils::print_polynomial(unsigned int nb_stages, const std::vector<unsigned int>& polynomial_powers, std::ostringstream& os)
{
    os << "X^" << nb_stages << " + ";
    
    std::vector<unsigned int>::const_iterator ppit = polynomial_powers.begin();
    const std::vector<unsigned int>::const_iterator pp_end = polynomial_powers.end();
    
    for (; ppit != pp_end; ++ppit)
    {
        os << "X^" << *ppit << " + ";
    }
    
    os << "1";
}


bool WsgcUtils::extract_string_vector(std::vector<std::string>& string_elements, std::string scs_string)
{
    std::string element_str;
    
    boost::char_separator<char> sep("/");
    boost::tokenizer<boost::char_separator<char> > tokens(scs_string, sep);
    
    boost::tokenizer<boost::char_separator<char> >::iterator tok_iter = tokens.begin();
    boost::tokenizer<boost::char_separator<char> >::iterator toks_end = tokens.end();
    
    for (; tok_iter != toks_end; ++tok_iter)
    {
        string_elements.push_back(*tok_iter);
    }
    
    return true;
}


unsigned int WsgcUtils::shortest_cyclic_distance(unsigned int i, unsigned int j, unsigned int cycle_length)
{
    unsigned int d1, d2;
    
    if (i<j)
    {
        d1 = j-i;
        d2 = cycle_length+i-j;
    }
    else
    {
        d1 = i-j;
        d2 = cycle_length+j-i;
    }
    
    return std::min(d1,d2);
}

void WsgcUtils::magnitude_estimation(wsgc_complex *in, wsgc_float *magnitude)
{
   /* magnitude ~= alpha * max(|I|, |Q|) + beta * min(|I|, |Q|) */

   wsgc_float abs_inphase = fabs(in->real());
   wsgc_float abs_quadrature = fabs(in->imag());

   if (abs_inphase > abs_quadrature) 
   {
      *magnitude = magnitude_estimation_alpha * abs_inphase + magnitude_estimation_beta * abs_quadrature;
   } 
   else 
   {
      *magnitude = magnitude_estimation_alpha * abs_quadrature + magnitude_estimation_beta * abs_inphase;
   }
}

void WsgcUtils::magnitude_algebraic(wsgc_complex *in, wsgc_float *magnitude)
{
   /* magnitude_algebraic = |I| + |Q| */
	*magnitude = fabs(in->real()) + fabs(in->imag());
}

void WsgcUtils::print_histo(const std::vector<std::pair<unsigned int, unsigned int> >& histo, std::ostringstream& os)
{
    std::vector<std::pair<unsigned int, unsigned int> >::const_iterator histo_it = histo.begin();
    const std::vector<std::pair<unsigned int, unsigned int> >::const_iterator histo_end = histo.end();

    os << "[";

    for (; histo_it != histo_end; ++histo_it)
    {
        if (histo_it != histo.begin())
        {
            os << ", ";
        }

        os << "(" << histo_it->first << "," << histo_it->second << ")";
    }

    os << "]";
}

void WsgcUtils::print_interval(std::ostringstream& os, unsigned int start, unsigned int length, unsigned int wrap_limit)
{
    if (length == 0)
    {
        os << "[]";
    }
    else if (length == 1)
    {
        os << "[" << start << "]";
    }
    else
    {
        if (wrap_limit == 0)
        {
            os << "[" << start << ":" << start + length - 1 << "]";
        }
        else
        {
            os << "[" << start << ":" << (start + length - 1) % wrap_limit << "]";
        }
    }
}


unsigned long WsgcUtils::timenow_usec()
{
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts); 
    return (ts.tv_sec*1000000) + (ts.tv_nsec / 1000);
}


unsigned int WsgcUtils::timenow_usec_hour()
{
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return ((ts.tv_sec % 3600)*1000000) + (ts.tv_nsec / 1000);
}

