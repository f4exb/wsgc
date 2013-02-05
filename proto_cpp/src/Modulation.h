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

     Modulation

     Represents a modulation its properties, specific methods and attributes
*/
#ifndef __MODULATION_H__
#define __MODULATION_H__

#include <sstream>

class Modulation
{
public:
    typedef enum
    {
        Modulation_BPSK,
        Modulation_DBPSK,
        Modulation_OOK,
        Modulation_CW
    } ModulationScheme_t;

    Modulation(ModulationScheme_t modulation_scheme) : _modulation_scheme(modulation_scheme)
    {}

    void setScheme(ModulationScheme_t modulation_scheme)
    {
        _modulation_scheme = modulation_scheme;
    }

    const ModulationScheme_t getScheme() const
    {
        return _modulation_scheme;
    }

    bool isFrequencyDependant();
    bool isCodeDivisionCapable();
    bool isDifferential();
    bool demodulateBeforeCorrelate();
    void print_modulation_data(std::ostringstream& os);

private:
    ModulationScheme_t _modulation_scheme;
};

#endif // __MODULATION_H__
