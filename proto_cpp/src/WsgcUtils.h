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

#ifndef __WSGC_UTILS__
#define __WSGC_UTILS__

#include "WsgcTypes.h"
#include <time.h>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>


// template to extract a vector of elements from a comma separated string
template<typename TElement> bool extract_vector(std::vector<TElement>& velements, std::string cs_string)
{
    std::string element_str;
    TElement element;
    
    boost::char_separator<char> sep(",");
    boost::tokenizer<boost::char_separator<char> > tokens(cs_string, sep);
    
    boost::tokenizer<boost::char_separator<char> >::iterator tok_iter = tokens.begin();
    boost::tokenizer<boost::char_separator<char> >::iterator toks_end = tokens.end();
    
    try
    {
        for (; tok_iter != toks_end; ++tok_iter)
        {
            element = boost::lexical_cast<TElement>(*tok_iter); 
            velements.push_back(element);
        }
        return true;
    }
    catch (boost::bad_lexical_cast &)
    {
        std::cout << "wrong element in comma separated string argument: " << *tok_iter << std::endl;
        return false;
    }
}


// template to print a vector of printable elements
template<typename TElement, typename TDisplay> void print_vector(const std::vector<TElement>& v, unsigned int width, std::ostringstream& os)
{
    os << "[";
    
    typename std::vector<TElement>::const_iterator it = v.begin();
    const typename std::vector<TElement>::const_iterator v_end = v.end();
        
    try
    {
        for (; it != v_end; ++it)
        {        
            
            os << std::setw(width) << (TDisplay)(*it);
            
            if (it != v.begin()+v.size()-1)
            {
                os << ", ";
            }

        }
    
        os << "]";
    }
    catch (boost::bad_lexical_cast &)
    {
        std::cout << "print_vector error: element cannot be cast for display: " << *it << std::endl;
    }
}


template<typename TElement> void print_array_x(unsigned int lx, TElement *a, std::ostringstream& os, unsigned int precision=6)
{
    for (unsigned int xi=0; xi<lx; xi++)
    {
        os << "[" << xi << "] = " << std::fixed << std::setw(precision+3) << std::setprecision(precision) << a[xi] << std::endl;
    }
}


template<typename TElement> void print_array_xy(unsigned int lx, unsigned int ly, TElement *a, std::ostringstream& os, unsigned int precision=6)
{
    for (unsigned int xi=0; xi<lx; xi++)
    {
        for (unsigned int yi=0; yi<ly; yi++)
        {
            os << "[" << xi << "," << yi << "] = " << std::fixed << std::setw(precision+3) << std::setprecision(precision) << a[xi*ly*+yi] << std::endl;
        }
    }
}

template<typename TElement> void print_array_xzy(unsigned int lx, unsigned int ly, unsigned int lz, TElement *a, std::ostringstream& os, unsigned int precision=6)
{
    for (unsigned int xi=0; xi<lx; xi++)
    {
        for (unsigned int zi=0; zi<lz; zi++)
        {
            for (unsigned int yi=0; yi<ly; yi++)
            {
                os << "[" << xi << "," << yi << "," << zi << "] = " << std::fixed << std::setw(precision+3) << std::setprecision(precision) << a[xi*ly*lz+yi*lz+zi] << std::endl;
            }
        }
    }
}

//void print_array_x(unsigned int lx, wsgc_complex *a, std::ostringstream& os, unsigned int width);
//void print_array_xy(unsigned int lx, unsigned int ly, wsgc_complex *a, std::ostringstream& os, unsigned int width=6);


class WsgcUtils
{
    public:
        // time difference in seconds
        static double get_time_difference(const timespec& time1, const timespec& time2);
        static void print_polynomial(unsigned int nb_stages, const std::vector<unsigned int>& polynomial_powers, std::ostringstream& os);
        static bool extract_string_vector(std::vector<std::string>& string_elements, std::string scs_string);
        static unsigned int shortest_cyclic_distance(unsigned int i, unsigned int j, unsigned int cycle_length);
        static void magnitude_estimation(wsgc_complex *in, wsgc_float *magnitude);
        static void magnitude_algebraic(wsgc_complex *in, wsgc_float *magnitude);
        static void print_histo(const std::vector<std::pair<unsigned int, unsigned int> >& histo, std::ostringstream& os);
        static void print_interval(std::ostringstream& os, unsigned int start, unsigned int length, unsigned int wrap_limit = 0);
        
    private:
        static const wsgc_float magnitude_estimation_alpha;
        static const wsgc_float magnitude_estimation_beta;
};


#endif // __WSGC_UTILS__
