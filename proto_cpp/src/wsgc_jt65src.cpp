/*
     Copyright 2012-2013 Edouard Griffiths <f4exb at free dot fr>

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

     wsgc_jt65src

     Test Source codec for JT65 with genuine scheme. 
     Without argument: exercises various hardcoded use cases.
     With string argument: encodes the given string.
     
*/

#include "SourceCodec_JT65.h"

#include <iostream>
#include <sstream>

//=================================================================================================
class InputMessages
{
public:
    InputMessages() {};
    ~InputMessages() {};
    
    const std::string& get_message(unsigned int index) const
    {
        if (index <  messages.size())
        {
            return messages[index];
        }
        else
        {
            return messages[0];
        }
    }
    
    unsigned int nb_messages() const
    {
        return messages.size();
    }
protected:
    static const char *messages_vinit[];
    static const std::vector<std::string> messages;
};

//=================================================================================================
const char *InputMessages::messages_vinit[] = {
    "CQ F4EXB JN33",        // Plain CQ
    "QRZ F4EXB JN33",       // Plain QRZ
    "F4EXB EA3XU JN11",     // Plain 2 calls and loc
    "EA3XU F4EXB R-22",     // Plain 2 calls and report
    "F4EXB EA3XU RRR",      // Plain 2 calls and RRR
    "CQ 112 F4EXB JN33",    // CQ with frequency
    "THIS IS FREEFLOW",     // Freeflow
    "CQ F4EXB/P JN33",      // v.1 suffix
    "F6HTJ 3A/F4EXB RRR",   // v.1 prefix on 2nd call
    "3A/F6HTJ F4EXB RRR",   // v.1 prefix on 1st call
    "CQ 3A2/F4EXB JN33",    // v.2 prefix
    "CQ F4EXB/3A2 JN33",    // v.2 suffix
};
const std::vector<std::string> InputMessages::messages(messages_vinit, messages_vinit + (sizeof messages_vinit / sizeof *messages_vinit));

//=================================================================================================
static void print_msg_symbols(const std::vector<unsigned int>& msg, std::ostream& os);

//=================================================================================================
int main(int argc, char *argv[])
{
	SourceCodec_JT65 src_codec;
	std::vector<unsigned int> out_msg;

	if (argc < 2)
	{
		InputMessages input_messages;

		for (unsigned int ti = 0; ti < input_messages.nb_messages(); ti++)
		{
			std::cout << "=== Test #" << ti+1 << " ===" << std::endl;
			std::cout << "Input: \"" << input_messages.get_message(ti) << "\"" << std::endl;

			if (src_codec.encode(input_messages.get_message(ti), out_msg))
			{
				std::cout << "Output: ";
				print_msg_symbols(out_msg, std::cout);
				std::cout << std::endl;
			}
			else
			{
				std::cout << "Message could not be encoded" << std::endl;
			}

			std::cout << std::endl;
			out_msg.clear();
		}
	}
	else
	{
		std::string input_msg(argv[1]);
		std::cout << "Input: \"" << input_msg << "\"" << std::endl;

		if (src_codec.encode(input_msg, out_msg))
		{
			std::cout << "Output: ";
			print_msg_symbols(out_msg, std::cout);
			std::cout << std::endl;
		}
		else
		{
			std::cout << "Message could not be encoded" << std::endl;
		}
	}
}

//=================================================================================================
void print_msg_symbols(const std::vector<unsigned int>& msg, std::ostream& os)
{
    std::vector<unsigned int>::const_iterator mit = msg.begin();
    const std::vector<unsigned int>::const_iterator mend = msg.end();
    
    os << "[";
    
    for (; mit != mend; ++mit)
    {
        if (mit != msg.begin())
        {
            os << ",";
        }
        
        os << *mit;
    }
    
    os << "]";
}
