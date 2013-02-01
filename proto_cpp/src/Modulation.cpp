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

#include "Modulation.h"

bool Modulation::isFrequencyDependant()
{
    if (_modulation_scheme == Modulation_BPSK)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool Modulation::isCodeDivisionCapable()
{
    if (_modulation_scheme == Modulation_BPSK)
    {
        return true;
    }
    else
    {
        return false;
    }

}

void Modulation::print_modulation_data(std::ostringstream& os)
{
    switch (_modulation_scheme)
    {
        case Modulation_BPSK:
            os << "BPSK";
            break;
        case Modulation_DBPSK:
            os << "DBPSK";
            break;
        case Modulation_OOK:
            os << "OOK";
            break;
        default:
            os << "None";
            break;
    }
}


bool Modulation::isDifferential()
{
    if (_modulation_scheme == Modulation_DBPSK)
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool Modulation::demodulateBeforeCorrelate()
{
    if ((_modulation_scheme == Modulation_DBPSK) || (_modulation_scheme == Modulation_OOK))
    {
        return true;
    }
    else
    {
        return false;
    }
}