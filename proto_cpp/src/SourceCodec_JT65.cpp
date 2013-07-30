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

     SourceCodec_JT65

     Source codec for JT65 genuine scheme. Takes a message as a string and encodes it in a
     vector of symbols (unsigned integers)
     
*/
#include "SourceCodec_JT65.h"
#include "Locator.h"
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

const char *StandardPrefixes::prefix_vinit[] = {
       "1A","1S","3A","3B6","3B8","3B9","3C","3C0",
	   "3D2","3D2C","3D2R","3DA","3V","3W","3X","3Y",
	   "3YB","3YP","4J","4L","4S","4U1I","4U1U","4W",
	   "4X","5A","5B","5H","5N","5R","5T","5U",
	   "5V","5W","5X","5Z","6W","6Y","7O","7P",
	   "7Q","7X","8P","8Q","8R","9A","9G","9H",
	   "9J","9K","9L","9M2","9M6","9N","9Q","9U",
	   "9V","9X","9Y","A2","A3","A4","A5","A6",
	   "A7","A9","AP","BS7","BV","BV9","BY","C2",
	   "C3","C5","C6","C9","CE","CE0X","CE0Y","CE0Z",
	   "CE9","CM","CN","CP","CT","CT3","CU","CX",
	   "CY0","CY9","D2","D4","D6","DL","DU","E3",
	   "E4","EA","EA6","EA8","EA9","EI","EK","EL",
	   "EP","ER","ES","ET","EU","EX","EY","EZ",
	   "F","FG","FH","FJ","FK","FKC","FM","FO",
	   "FOA","FOC","FOM","FP","FR","FRG","FRJ","FRT",
	   "FT5W","FT5X","FT5Z","FW","FY","M","MD","MI",
	   "MJ","MM","MU","MW","H4","H40","HA",
	   "HB","HB0","HC","HC8","HH","HI","HK","HK0A",
	   "HK0M","HL","HM","HP","HR","HS","HV","HZ",
	   "I","IS","IS0","J2","J3","J5","J6",
	   "J7","J8","JA","JDM","JDO","JT","JW",
	   "JX","JY","K","KG4","KH0","KH1","KH2","KH3",
	   "KH4","KH5","KH5K","KH6","KH7","KH8","KH9","KL",
	   "KP1","KP2","KP4","KP5","LA","LU","LX","LY",
	   "LZ","OA","OD","OE","OH","OH0","OJ0","OK",
	   "OM","ON","OX","OY","OZ","P2","P4","PA",
	   "PJ2","PJ7","PY","PY0F","PT0S","PY0T","PZ","R1F",
	   "R1M","S0","S2","S5","S7","S9","SM","SP",
	   "ST","SU","SV","SVA","SV5","SV9","T2","T30",
	   "T31","T32","T33","T5","T7","T8","T9","TA",
	   "TF","TG","TI","TI9","TJ","TK","TL",
	   "TN","TR","TT","TU","TY","TZ","UA","UA2",
	   "UA9","UK","UN","UR","V2","V3","V4","V5",
	   "V6","V7","V8","VE","VK","VK0H","VK0M","VK9C",
	   "VK9L","VK9M","VK9N","VK9W","VK9X","VP2E","VP2M","VP2V",
	   "VP5","VP6","VP6D","VP8","VP8G","VP8H","VP8O","VP8S",
	   "VP9","VQ9","VR","VU","VU4","VU7","XE","XF4",
	   "XT","XU","XW","XX9","XZ","YA","YB","YI",
	   "YJ","YK","YL","YN","YO","YS","YU","YV",
	   "YV0","Z2","Z3","ZA","ZB","ZC4","ZD7","ZD8",
	   "ZD9","ZF","ZK1N","ZK1S","ZK2","ZK3","ZL","ZL7",
	   "ZL8","ZL9","ZP","ZS","ZS8","KC4","E5"
};

const std::vector<std::string> StandardPrefixes::prefixes(prefix_vinit, prefix_vinit + (sizeof prefix_vinit / sizeof *prefix_vinit));

const unsigned int SourceCodec_JT65::call_nbase = 37*36*10*27*27*27;
const unsigned int SourceCodec_JT65::locator_nbase = 180*180;
const unsigned int SourceCodec_JT65::call_DE = 267796945;
const std::string  SourceCodec_JT65::free_text_chars = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ +-./?";
const std::string  SourceCodec_JT65::suffix_chars = "P0123456789A";
const std::string  SourceCodec_JT65::callsign_chars = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ ";
const unsigned int SourceCodec_JT65::v2_prefix_cq_shift  = 262178563; // call_nbase + 1003
const unsigned int SourceCodec_JT65::v2_prefix_qrz_shift = 264002072; 
const unsigned int SourceCodec_JT65::v2_prefix_de_shift  = 265825581;
const unsigned int SourceCodec_JT65::v2_suffix_cq_shift  = 267649090; 
const unsigned int SourceCodec_JT65::v2_suffix_qrz_shift = 267698375;
const unsigned int SourceCodec_JT65::v2_suffix_de_shift  = 267747660; 

//=================================================================================================
StandardPrefixes::StandardPrefixes()
{}

//=================================================================================================
StandardPrefixes::~StandardPrefixes()
{}

//=================================================================================================
const std::string& StandardPrefixes::get_prefix(unsigned int index, bool& found) const
{
	if (index >=  prefixes.size())
	{
		found = false;
		return prefixes[0];
	}
	else
	{
		found = true;
		return prefixes[index];
	}
}

//=================================================================================================
unsigned int StandardPrefixes::get_index(const std::string& prefix, bool& found) const
{
	for (unsigned int index=0; index < prefixes.size(); index++)
	{
		if (prefix == prefixes[index])
		{
			found = true;
			return index;
		}
	}

	found = false;
	return 0;
}

//=================================================================================================
SourceCodec_JT65::SourceCodec_JT65(JT65_Variants _variant) :
    variant(_variant)
{}
   
//=================================================================================================
SourceCodec_JT65::~SourceCodec_JT65()
{}
    
//=================================================================================================
bool SourceCodec_JT65::encode(const std::string& in_msg, std::vector<unsigned int>& out_msg) const
{
	bool success = false;
	unsigned int packed_callsign_1 = 0;
	int pfxsfx_index_1 = -1;
	unsigned int packed_callsign_2 = 0;
	int pfxsfx_index_2 = -1;
	unsigned int packed_locator = 0;
	unsigned int frequency;
	bool seek_you = false;

	std::string msg_str(in_msg); // get a copy
	std::transform(msg_str.begin(), msg_str.end(),msg_str.begin(), ::toupper); // convert to uppercase

	boost::char_separator<char> sep(" ");
	boost::tokenizer<boost::char_separator<char> > tokens(msg_str, sep);
    boost::tokenizer<boost::char_separator<char> >::iterator tok_iter = tokens.begin();
    boost::tokenizer<boost::char_separator<char> >::iterator toks_end = tokens.end();
    unsigned int field_index = 0;
    unsigned int v2_indicator;

    while ((tok_iter != toks_end) && (field_index < 3))
    {
        //std::cout << field_index << "> " << *tok_iter << std::endl;

    	if (field_index == 0) // First call placeholder
    	{
			if (seek_you) // CQ xxx, token on xxx
			{
				try // try to find a frequency
				{
					frequency = boost::lexical_cast<unsigned int>(*tok_iter);

					if (frequency < 1000)
					{
						packed_callsign_1 = call_nbase + 3 + frequency;
						success = true;
					}
					else
					{
						success = false;
					}
				}
				catch (boost::bad_lexical_cast &) // then try to encode the second callsign
				{
					packed_callsign_1 = call_nbase + 1; // encode CQ
					success = pack_callsign(*tok_iter, packed_callsign_2, pfxsfx_index_2, v2_indicator);
					field_index++;
				}

				field_index++;
				++tok_iter;
			}
			else if (*tok_iter == "CQ")
			{
				seek_you = true;
                success = true;
				++tok_iter;
			}
			else if (*tok_iter == "QRZ")
			{
				packed_callsign_1 = call_nbase + 2;
                success = true;
				field_index++;
				++tok_iter;
			}
			else if (*tok_iter == "DE")
			{
				packed_callsign_1 = call_DE;
                success = true;
				field_index++;
				++tok_iter;
			}
			else // try to encode the first callsign
			{
				success = pack_callsign(*tok_iter, packed_callsign_1, pfxsfx_index_1, v2_indicator);
				field_index++;
				++tok_iter;
			}
    	} // field index 0
    	else if (field_index == 1) // Second call placeholder
    	{
    		success = pack_callsign(*tok_iter, packed_callsign_2, pfxsfx_index_2, v2_indicator);
			field_index++;
			++tok_iter;
    	} // field index 1
    	else if (field_index == 2) // Locator placeholder
    	{
    		if (*tok_iter == "RO")
    		{
    			packed_locator = locator_nbase + 62;
    			success = true;
    		}
    		else if (*tok_iter == "RRR")
    		{
    			packed_locator = locator_nbase + 63;
    			success = true;
    		}
    		else if (*tok_iter == "73")
    		{
    			packed_locator = locator_nbase + 64;
    			success = true;
    		}
    		else
    		{
    			std::string loc_rpt(*tok_iter);

    			if (loc_rpt.size() > 2)  
    			{
    				if (loc_rpt[0] == 'R' && loc_rpt[1] == '-') // Report as R-NN
    				{
                        if (loc_rpt.size() == 4)
                        {
                            success = pack_report_jt(true, loc_rpt.substr(2,2), packed_locator);    
                        }
                        else
                        {
                            success = false;
                        }
    				}
    				else if (loc_rpt[0] == '-') // Report as -NN
    				{
                        if (loc_rpt.size() == 3)
                        {
                            success = pack_report_jt(false, loc_rpt.substr(1,2), packed_locator);
                        }
                        else
                        {
                            success = false;
                        }
    				}
    				else // true locator
    				{
    					success = pack_locator_4(loc_rpt, packed_locator);
    				}
    			}
    			else
    			{
    				success = false;
    			}
    		}

			field_index++;
			++tok_iter;
    	} // field index 2

    	if (!success)
    	{
    		break;
    	}
    } // tokens

    if (success) // basic encoding worked now try to overwrite with prefixes or suffixes
    {
		if (pfxsfx_index_1 > 0) // prefix or suffix on 1st callsign
		{
			if (pfxsfx_index_2 > 0)
			{
				success = false;
			}
			else
			{
				success = pack_pfxsfx_v1(pfxsfx_index_1, packed_locator);
			}
		}
		else if (pfxsfx_index_2 > 0) // prefix or suffix on 2nd callsign
		{
			if (v2_indicator > 0)
			{
				if (packed_callsign_1 > call_nbase) // first callsign placeholder is CQ,QRZ or DE
				{
                    success = pack_pfxsfx_v2(v2_indicator, pfxsfx_index_2, packed_callsign_1 - call_nbase, packed_callsign_1); // pack prefix/suffix in v2 version if necessary
				}
				else // first callsign placeholder is a callsign
				{
					success = false; // v.2 prefix/suffix is not supported with two callsigns
				}
			}
			else
			{
				success = pack_pfxsfx_v1(pfxsfx_index_2 + 450, packed_locator);
			}
		}
    }

    if (!success) // if encoding formatted message failed attempt to encode plain text
    {
        //std::cout << "Attempt to code as free text" << std::endl;
        
    	if (!pack_text(msg_str, packed_callsign_1, packed_callsign_2, packed_locator))
        {
            return false;
        }
    }

    //std::cout << packed_callsign_1 << ":" << packed_callsign_2 << ":" << packed_locator << std::endl;
    if (variant == JT65_Classical)
    {
        pack_message(packed_callsign_1, packed_callsign_2, packed_locator, out_msg);
    }
    else if (variant == JT65_256)
    {
        pack_message_256(packed_callsign_1, packed_callsign_2, packed_locator, out_msg);
    }
	return true;
}  

//=================================================================================================
bool SourceCodec_JT65::decode(const std::vector<unsigned int>& in_msg, std::string& out_msg) const
{
    unsigned int packed_callsign_1;
    unsigned int packed_callsign_2;
    unsigned int packed_locator;
    bool is_arbitrary_text;
    std::string field_1;
    std::string field_2;
    std::string field_3;
    std::string pfxsfx_v1_str;
    std::string pfxsfx_v2_str;
    unsigned int pfxsfx_v2 = 0;
    unsigned int pfxsfx_v1 = 0;
    bool success = true;
    
    if (variant == JT65_Classical)
    {
        unpack_message(in_msg, packed_callsign_1, packed_callsign_2, packed_locator, is_arbitrary_text);
    }
    else if (variant == JT65_256)
    {
        unpack_message_256(in_msg, packed_callsign_1, packed_callsign_2, packed_locator, is_arbitrary_text);
    }
    //std::cout << packed_callsign_1 << ":" << packed_callsign_2 << ":" << packed_locator << (is_arbitrary_text ? ":freetext" : "") << std::endl;
    
    if (is_arbitrary_text)
    {
        return unpack_text(out_msg, packed_callsign_1, packed_callsign_2, packed_locator);
    }
    else
    {
        for (unsigned int n_field = 0; (n_field < 3) && success; n_field++)
        {
        	//std::cout << "n_field = " << n_field << std::endl;
            if (n_field == 0)
            {
                success = unpack_callsign_1(packed_callsign_1, field_1, pfxsfx_v2_str, pfxsfx_v2);
            }
            else if (n_field == 1)
            {
                success = unpack_plain_callsign(packed_callsign_2, field_2);
            }
            else if (n_field == 2)
            {
                success = unpack_locator(packed_locator, field_3, pfxsfx_v1_str, pfxsfx_v1);
            }
        }
        
        //std::cout << (success ? "ok" : "ko") << std::endl;

        if (success)
        {
        	//std::cout << "field_1 = " << field_1 << ", field_2 = " << field_2 << ", field_3 = " << field_3 << ", pfxsfx_v1 = " << pfxsfx_v1 << ", pfxsfx_v2 = " << pfxsfx_v2 << ", pfxsfx_v1_str = " << pfxsfx_v1_str << ", pfxsfx_v2_str = " << pfxsfx_v2_str << std::endl;
            compose_formatted_message(out_msg, field_1, field_2, field_3, pfxsfx_v1, pfxsfx_v2, pfxsfx_v1_str, pfxsfx_v2_str);
            return true;
        }
        else
        {
            return false;
        }
    }
}

//=================================================================================================
bool SourceCodec_JT65::pack_pfxsfx_v1(int pfxsfx_index, unsigned int& packed_locator) const
{
    if (pfxsfx_index > 900)
    {
        return false;
    }
    else
    {
        int nlong = 2 * (((pfxsfx_index-1)/5) % 90) - 180;
        
        if (pfxsfx_index > 450)
        {
            nlong += 180;
        }
        
        int nlat = ((pfxsfx_index-1)%5) + 85 + 90;
        
        float dlat  = nlat;
        float dlong = nlong;

        //std::cout << "pfxsfx_index = " << pfxsfx_index << ", pfx_lat = " << dlat << ", pfx_lon = " << dlong << std::endl;
        
        packed_locator = ((180+dlong)/2)*180 + dlat;
        
        return true;
    }
}

//=================================================================================================
bool SourceCodec_JT65::pack_pfxsfx_v2(unsigned int v2_indicator, int pfxsfx_index, unsigned int call_indicator, unsigned int& packed_callsign) const
{
    if (v2_indicator == 1) // prefix
    {
        if (call_indicator == 1) // CQ
        {
            packed_callsign = pfxsfx_index + v2_prefix_cq_shift;
        }
        else if (call_indicator == 2) // QRZ
        {
            packed_callsign = pfxsfx_index + v2_prefix_qrz_shift;
        }
        else if (call_indicator == 3) // DE
        {
            packed_callsign = pfxsfx_index + v2_prefix_qrz_shift;
        }
        else // CQ with frequency => CQ
        {
        	packed_callsign = pfxsfx_index + v2_prefix_cq_shift;
        }
    }
    else if (v2_indicator == 2) // suffix
    {
        if (call_indicator == 1) // CQ
        {
            packed_callsign = pfxsfx_index + v2_suffix_cq_shift;
        }
        else if (call_indicator == 2) // QRZ
        {
            packed_callsign = pfxsfx_index + v2_suffix_qrz_shift;
        }
        else if (call_indicator == 3) // DE
        {
            packed_callsign = pfxsfx_index + v2_suffix_qrz_shift;
        }
        else // CQ with frequency => CQ
        {
        	packed_callsign = pfxsfx_index + v2_suffix_cq_shift;
        }
    }
    else
    {
        return false;
    }
    
    return true;
}

//=================================================================================================
bool SourceCodec_JT65::pack_plain_callsign(const std::string& callsign, unsigned int& packed_callsign) const
{
	std::string padded_callsign("      ");

	if ((callsign.size() > 6) || (callsign.size() < 3))
	{
		return false;
	}
	else
	{
		int callsign_char_index;

		// look for third caracter
		if ((callsign_char_index = callsign_chars.find(callsign[2])) == std::string::npos)
		{
			return false; // invalid callsign
		}

		if (callsign_char_index < 10) // third character numeric => country prefix on 2 characters
		{
			padded_callsign.replace(0, callsign.size(), callsign);
		}
		else // third character alpha => country prefix on 1 character and max 5 characters
		{
			if (callsign.size() > 5)
			{
				return false;
			}
            
            // and second character must be numeric
            if ((callsign_char_index = callsign_chars.find(callsign[1])) == std::string::npos)
            {
                return false; // character invalid
            }
            else
            {
                if (callsign_char_index >= 10) 
                {
                    return false; // not numeric
                }
            }

			padded_callsign.replace(1, callsign.size(), callsign);
		}
        
        //std::cout << "Padded callsign: \"" << padded_callsign << "\"" << std::endl;

		for (unsigned int i=0; i<6; i++)
		{
			if ((callsign_char_index = callsign_chars.find(padded_callsign[i])) == std::string::npos)
			{
				return false;
			}

			if (i==0)
			{
				packed_callsign = callsign_char_index;
			}
            else if (i==1) 
            {
                packed_callsign = packed_callsign*36 + callsign_char_index;
            }
            else if (i==2) // always numeric (already checked)
            {
                packed_callsign = packed_callsign*10 + callsign_char_index;
            }
			else
			{
                if (callsign_char_index < 10)
                {
                    return false; // never numeric in this position
                }
                
				packed_callsign = packed_callsign*27 + callsign_char_index - 10;
			}
		}

		return true;
	}
}

//=================================================================================================
bool SourceCodec_JT65::pack_callsign(const std::string& callsign, unsigned int& packed_callsign, int& pfxsfx_index, unsigned int& pfxsfx_v2) const
{
	int slash_index;
	unsigned int pfxsfx;

	if ((slash_index = callsign.find("/")) == std::string::npos) // no prefix, no suffix, encode normal callsign
	{
		pfxsfx = 0;
		pfxsfx_v2 = 0;
		pfxsfx_index = -1;
		return pack_plain_callsign(callsign, packed_callsign);
	}
	else
	{
        //std::cout << "call size = " << callsign.size() << ", slash index = " << slash_index << std::endl;
        
        if ((slash_index == 0) || (slash_index == callsign.size()-1))
        {
        	return false;
        }

		if (slash_index == callsign.size() - 2) // possible v.1 suffix
		{
			int suffix_char_index;
			pfxsfx = 2;

			if ((suffix_char_index = suffix_chars.find(callsign[callsign.size() - 1])) == std::string::npos)
			{
				pfxsfx_v2 = 2; // v.2 suffix (calculated later)
			}
			else
			{
				pfxsfx_index = 401 + suffix_char_index; // original based on Fortran indexes starting at 1
				pfxsfx_v2 = 0; // v.1 suffix
			}
		}
		else if (slash_index <= 4) // v.1 prefix
		{
			pfxsfx = 1;
			std::string prefix = callsign.substr(0,slash_index);
			bool prefix_index_found;
			pfxsfx_index = standard_prefixes.get_index(prefix, prefix_index_found) + 1; // original based on Fortran indexes

			if (prefix_index_found)
			{
				pfxsfx_v2 = 0;
			}
			else
			{
				pfxsfx_v2 = 1; // v.2 prefix (calculated later)
			}
		}
		else if (slash_index <= callsign.size() - 4)
		{
			pfxsfx = 2;
			pfxsfx_v2 = 2; // v.2 suffix (calculated later)
		}

		if (pfxsfx_v2 == 1) // attempt v.2 prefix coding
		{
			std::string prefix = callsign.substr(0,slash_index);
			int callsign_char_index;

			for (unsigned int i=0; i<4; i++)
			{
				if (i < prefix.size())
				{
					if ((callsign_char_index = callsign_chars.find(prefix[i])) == std::string::npos)
					{
						return false; // invalid character in prefix
					}
				}
				else
				{
					callsign_char_index = callsign_chars.size()-1;
				}

				if (i==0)
				{
					pfxsfx_index = callsign_char_index;
				}
				else
				{
					pfxsfx_index = callsign_chars.size()*pfxsfx_index + callsign_char_index;
				}
			}
		}
		else if (pfxsfx_v2 == 2) // attempt v.2 suffix coding
		{
			std::string suffix = callsign.substr(slash_index+1);
			int callsign_char_index;

			for (unsigned int i=0; i<3; i++)
			{
				if (i < suffix.size())
				{
					if ((callsign_char_index = callsign_chars.find(suffix[i])) == std::string::npos)
					{
						return false; // invalid character in suffix
					}
				}
				else
				{
					callsign_char_index = callsign_chars.size()-1;
				}

				if (i==0)
				{
					pfxsfx_index = callsign_char_index;
				}
				else
				{
					pfxsfx_index = callsign_chars.size()*pfxsfx_index + callsign_char_index;
				}
			}
		}

		// now try to encode callsign
		std::string plain_callsign;

		if (pfxsfx == 1)
		{
			plain_callsign = callsign.substr(slash_index+1);
		}
		else if (pfxsfx == 2)
		{
			plain_callsign = callsign.substr(0, slash_index);
		}

		//std::cout << "v2_ind = " << pfxsfx_v2 << ", pfxsfx_index = " << pfxsfx_index << ", plain call = " << plain_callsign << std::endl;
		return pack_plain_callsign(plain_callsign, packed_callsign);
	} // slash
}

//=================================================================================================
bool SourceCodec_JT65::pack_report_jt(bool r_prefix, const std::string& report_str, unsigned int& packed_locator) const
{
    unsigned int report;
    //std::cout << "report_str = " << report_str << std::endl;
    
    try // try to find a number
    {
        report = boost::lexical_cast<unsigned int>(report_str);

        if (report > 30) // report in [0..30] (31 values)
        {
            report = 30;
        }
        
        if (r_prefix)
        {
            packed_locator = locator_nbase + 1 + 31 + report;
        }
        else
        {
            packed_locator = locator_nbase + 1 + report;
        }
        
        return true;
    }
    catch (boost::bad_lexical_cast &) // report not valid
    {
    }

    return false;
}

//=================================================================================================
bool SourceCodec_JT65::pack_locator_4(const std::string& locator, unsigned int& packed_locator) const
{
    std::string locator_str(locator + "MM");

    try
    {
        Locator locator(locator_str);
        
        if (locator.latitude() < 85.0) // Locators within 5 degrees from the North Pole are used for prefix and suffix information
        {
            int ilong = locator.longitude();
            int ilat  = locator.latitude() + 90.0;
            packed_locator = ((180-ilong)/2)*180 + ilat; // longitude is inverted: positive to the west
            return true;
        }
        else
        {
            return false;
        }
    }
    catch (LocatorInvalidException &)
    {
    }

    return false;
}

//=================================================================================================
bool SourceCodec_JT65::pack_text(const std::string& text,
        unsigned int& packed_callsign_1,
        unsigned int& packed_callsign_2,
        unsigned int& packed_locator) const
{
    if (text.size() < 13)
    {
        return false;
    }
    else
    {
        unsigned int alphabet_size = free_text_chars.size();
        unsigned int i;
        int char_index;
        
        // Characters 1-5 go into packed_callsign_1
        packed_callsign_1 = 0;
        
        for (i=0; i<5; i++) 
        {
            if ((char_index = free_text_chars.find(text[i])) == std::string::npos)
            {
                return false;
            }
            
            packed_callsign_1 = packed_callsign_1*alphabet_size + char_index;
        }
        
        // Characters 6-10 go into packed_callsign_2
        packed_callsign_2 = 0;
        
        for (i=5; i<10; i++) 
        {
            if ((char_index = free_text_chars.find(text[i])) == std::string::npos)
            {
                return false;
            }
            
            packed_callsign_2 = packed_callsign_2*alphabet_size + char_index;
        }
        
        // Characters 11-13 go into packed_locator
        packed_locator = 0;
        
        for (i=10; i<13; i++) 
        {
            if ((char_index = free_text_chars.find(text[i])) == std::string::npos)
            {
                return false;
            }
            
            packed_locator = packed_locator*alphabet_size + char_index;
        }
        
        // We now have used 17 bits in packed_locator.  Must move one each to packed_callsign_1 and packed_callsign_2.
        packed_callsign_1 <<= 1;
        packed_callsign_2 <<= 1;
        
        if (packed_locator & 0x1000)
        {
            packed_callsign_1 += 1;
        }
        
        if (packed_locator & 0x10000)
        {
            packed_callsign_1 += 2;
        }
        
        packed_locator = packed_locator & 0x7FFF;
        packed_locator = packed_locator | 0x8000; // move in text indicator
        
        return true;
    }
}

//=================================================================================================
bool SourceCodec_JT65::unpack_text(std::string& text,
        unsigned int packed_callsign_1,
        unsigned int packed_callsign_2,
        unsigned int packed_locator) const
{
    packed_locator &= 0x7FFF; // Remove arbitrary text indicator bit
    
    if (packed_callsign_1 & 0x01)
    {
        packed_locator |= 0x8000;
    }

    if (packed_callsign_2 & 0x01)
    {
        packed_locator |= 0x10000;
    }
    
    packed_callsign_1 >>= 1;
    packed_callsign_2 >>= 1;
    unsigned int i;
    unsigned int j;
    unsigned int alphabet_size = free_text_chars.size();
    text = "             ";
    
    for (i=5; i>0; i--)
    {
        j = packed_callsign_1 % alphabet_size;
        text[i-1] = free_text_chars[j];
        packed_callsign_1 /= alphabet_size;
    }
    
    for (i=10; i>5; i--)
    {
        j = packed_callsign_2 % alphabet_size;
        text[i-1] = free_text_chars[j];
        packed_callsign_2 /= alphabet_size;
    }
    
    for (i=13; i>10; i--)
    {
        j = packed_locator % alphabet_size;
        text[i-1] = free_text_chars[j];
        packed_locator /= alphabet_size;
    }
    
    return true;
}

//=================================================================================================
bool SourceCodec_JT65::unpack_callsign_1(unsigned int packed_callsign, std::string& field_1, std::string& pfxsfx_str, unsigned int& pfxsfx_v2) const
{
    if (packed_callsign > call_DE) // This is the maximum indicator
    {
        return false; 
    }
    else if (packed_callsign <= call_nbase) // plain callsign
    {
        pfxsfx_v2 = 0;
        return unpack_plain_callsign(packed_callsign, field_1);
    }
    else if (packed_callsign == call_nbase + 1) // plain CQ
    {
        pfxsfx_v2 = 0;
        field_1 = "CQ";
        return true;
    }
    else if (packed_callsign == call_nbase + 2) // plain QRZ
    {
        pfxsfx_v2 = 0;
        field_1 = "QRZ";
        return true;
    }
    else if (packed_callsign < call_nbase + 1003) // CQ + frequency
    {
        pfxsfx_v2 = 0;
        unsigned int frequency = packed_callsign - call_nbase - 3;
        std::ostringstream os;
        os << "CQ " << frequency;
        field_1 = os.str();
        return true;
    }
    else if (packed_callsign <  v2_prefix_qrz_shift) // v.2 prefix CQ range
    {
        pfxsfx_v2 = 1;
        field_1 = "CQ";
        return unpack_pfxsfx_v2(packed_callsign - v2_prefix_cq_shift, pfxsfx_str, 4);
    }
    else if (packed_callsign <  v2_prefix_de_shift)  // v.2 prefix QRZ range
    {
        pfxsfx_v2 = 1;
        field_1 = "QRZ";
        return unpack_pfxsfx_v2(packed_callsign - v2_prefix_qrz_shift, pfxsfx_str, 4);
    }
    else if (packed_callsign <  v2_suffix_cq_shift)  // v.2 prefix DE range
    {
        pfxsfx_v2 = 1;
        field_1 = "DE";
        return unpack_pfxsfx_v2(packed_callsign - v2_prefix_de_shift, pfxsfx_str, 4);
    }
    else if (packed_callsign <  v2_suffix_qrz_shift) // v.2 suffix CQ range
    {
        pfxsfx_v2 = 2;
        field_1 = "CQ";
        return unpack_pfxsfx_v2(packed_callsign - v2_suffix_cq_shift, pfxsfx_str, 3);
    }
    else if (packed_callsign <  v2_suffix_de_shift)  // v.2 suffix QRZ range
    {
        pfxsfx_v2 = 2;
        field_1 = "QRZ";
        return unpack_pfxsfx_v2(packed_callsign - v2_suffix_qrz_shift, pfxsfx_str, 3);
    }
    else if (packed_callsign <  call_DE)             // v.2 suffix DE range
    {
        pfxsfx_v2 = 2;
        field_1 = "DE";
        return unpack_pfxsfx_v2(packed_callsign - v2_suffix_de_shift, pfxsfx_str, 3);
    }
    else if (packed_callsign ==  call_DE) // plain DE
    {
        pfxsfx_v2 = 0;
        field_1 = "DE";
        return true;
    }
    else
    {
    	return false; // unregistered
    }
}

//=================================================================================================
bool SourceCodec_JT65::unpack_pfxsfx_v2(int pfxsfx_index, std::string& pfxsfx_str, unsigned int max_chars) const
{
    unsigned int alphabet_size = callsign_chars.size();
    pfxsfx_str = std::string(max_chars, ' ');
    
    for (unsigned int i = max_chars; i > 0; i--)
    {
        unsigned int j = pfxsfx_index % alphabet_size;
        pfxsfx_str[i-1] = callsign_chars[j];
        pfxsfx_index /= alphabet_size;
    }
    
    boost::algorithm::trim_right(pfxsfx_str); // trim blanks
    return true;
}

//=================================================================================================
bool SourceCodec_JT65::unpack_plain_callsign(unsigned int packed_callsign, std::string& callsign) const
{
    callsign = "      ";
    unsigned int i;
    unsigned int alphabet_size = callsign_chars.size();
    
    i = (packed_callsign % 27) + 10;
    callsign[5] = callsign_chars[i];
    packed_callsign /= 27;
    
    i = (packed_callsign % 27) + 10;
    callsign[4] = callsign_chars[i];
    packed_callsign /= 27;
    
    i = (packed_callsign % 27) + 10;
    callsign[3] = callsign_chars[i];
    packed_callsign /= 27;
    
    i = (packed_callsign % 10);
    callsign[2] = callsign_chars[i];
    packed_callsign /= 10;
    
    i = (packed_callsign % 36);
    callsign[1] = callsign_chars[i];
    packed_callsign /= 36;
    
    if (packed_callsign < alphabet_size)
    {
        callsign[0] = callsign_chars[packed_callsign];
        boost::algorithm::trim(callsign); // trim blanks
        return true;
    }
    else
    {
        return false;
    }
}

//=================================================================================================
bool SourceCodec_JT65::unpack_locator(unsigned int packed_locator, std::string& field_3, std::string& pfxsfx_str, unsigned int& pfxsfx_v1) const
{
    if (packed_locator <= locator_nbase) // true locator
    {
        float dlat = (packed_locator % 180) - 90.0;
        float dlong = (packed_locator/180)*2.0 - 180.0 + 2.0;
        bool pfx_found;
        
        
        if (dlat < 85) // real locator
        {
            //std::cout << "dlat = " << dlat << ", dlong = " << dlong << std::endl;
            Locator locator(dlat, -dlong);
            field_3 = locator.toString().substr(0,4);
            return true;
        }
        else // v.1 prefix or suffix
        {
            dlat += 90.0;
            dlong -= 2.0;
            //std::cout << "dlat = " << dlat << ", dlong = " << dlong << std::endl;
            field_3 = "";
            unsigned int nlong = dlong;
            unsigned int nlat = dlat;
            unsigned int pfxsfx_index = long_lat_to_pfxsfx_index(nlat, nlong);
            unsigned int new_index;
            //std::cout << "pfxsfx_index = " << pfxsfx_index << std::endl;
            
            if (pfxsfx_index < 450) // concerns first callsign
            {
                if (pfxsfx_index < 400) // prefix
                {
                	pfxsfx_v1 = 1;
                    pfxsfx_str = std::string(standard_prefixes.get_prefix(pfxsfx_index, pfx_found));
                    return pfx_found;
                }
                else // suffix
                {
                	pfxsfx_v1 = 3;
                    pfxsfx_index -= 400;
                    
                    if (pfxsfx_index < suffix_chars.size())
                    {
                        pfxsfx_str = suffix_chars[pfxsfx_index];
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
            }
            else // concerns second callsign
            {
                pfxsfx_index -= 450;
                
                if (pfxsfx_index < 400) // prefix
                {
                    pfxsfx_v1 = 2;
                    pfxsfx_str = std::string(standard_prefixes.get_prefix(pfxsfx_index, pfx_found));
                    return pfx_found;
                }
                else // suffix
                {
                    pfxsfx_v1 = 4;
                    pfxsfx_index -= 400;
                    
                    if (pfxsfx_index < suffix_chars.size())
                    {
                        pfxsfx_str = suffix_chars[pfxsfx_index];
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
            }
        }
    }
    else if (packed_locator < locator_nbase + 62) // report
    {
        unsigned int report = packed_locator - locator_nbase - 1;
        std::ostringstream os;

        if (report < 31)
        {
            os << report;
        }
        else
        {
            os << "R-" << report - 31;
        }

        field_3 = os.str();
        pfxsfx_v1 = 0;
        return true;
    }
    else if (packed_locator == locator_nbase + 62)
    {
    	field_3 = "RO";
    	pfxsfx_v1 = 0;
        return true;
    }
    else if (packed_locator == locator_nbase + 63)
    {
    	field_3 = "RRR";
    	pfxsfx_v1 = 0;
        return true;
    }
    else if (packed_locator == locator_nbase + 64)
    {
    	field_3 = "73";
    	pfxsfx_v1 = 0;
        return true;
    }
    else
    {
    	return false; // invalid still
    }
}

//=================================================================================================
unsigned int SourceCodec_JT65::long_lat_to_pfxsfx_index(unsigned int nlat, unsigned int nlong) const
{
    int pfxsfx_q, pfxsfx_r;
    //int nlong = 2 * (((pfxsfx_index-1)/5) % 90) - 180;
    //int nlat = ((pfxsfx_index-1)%5) + 85 + 90;
    pfxsfx_q = (nlong + 180)/2;
    pfxsfx_r = (nlat - 85 - 90);
    //std::cout << "pfxsfx_q = " << pfxsfx_q << ", pfxsfx_r = " << pfxsfx_r << std::endl;
    return pfxsfx_q*5 + pfxsfx_r;
}

//=================================================================================================
void SourceCodec_JT65::compose_formatted_message(std::string& out_msg,
        const std::string& field_1,
        const std::string& field_2,
        const std::string& field_3,
        unsigned int pfxsfx_v1,
        unsigned int pfxsfx_v2,
        const std::string& pfxsfx_str_v1,
        const std::string& pfxsfx_str_v2) const
{
    std::string new_field_1 = field_1;
    std::string new_field_2 = field_2;
    std::string new_field_3 = field_3;

    // v.1 prefix/suffix:
    if (pfxsfx_v1 == 1) // prefix on field 1
    {
        new_field_1 = pfxsfx_str_v1 + "/" + field_1;
        new_field_3 = "";
    }
    else if (pfxsfx_v1 == 2) // prefix on field 2
    {
        new_field_2 = pfxsfx_str_v1 + "/" + field_2;
        new_field_3 = "";
    }
    else if (pfxsfx_v1 == 3) // suffix on field 1
    {
        new_field_1 = field_1 + "/" + pfxsfx_str_v1;
        new_field_3 = "";
    }
    else if (pfxsfx_v1 == 4) // suffix on field 2
    {
        new_field_2 = field_2 + "/" + pfxsfx_str_v1;
        new_field_3 = "";
    }
    
    // v.2 prefix/suffix:
    if (pfxsfx_v2 == 1) // prefix
    {
        new_field_2 = pfxsfx_str_v2 + "/" + field_2;
    }
    else if (pfxsfx_v2 == 2) // suffix
    {
        new_field_2 = field_2 + "/" + pfxsfx_str_v2;
    } 
    
    out_msg = new_field_1 + " " + new_field_2 + " " + new_field_3;
    boost::algorithm::trim(out_msg);
}

//=================================================================================================
bool SourceCodec_JT65::pack_message(unsigned int packed_callsign_1,
        unsigned int packed_callsign_2,
        unsigned int packed_locator,
        std::vector<unsigned int>& message) const
{
    unsigned int code_byte = 0;
    // bits are entered MSB first
    // locator is 16 bit including text indicator bit on MSB
    code_byte = ((packed_callsign_1 >> 22) & (0x3F)); // 6 bits from callsign_1
    message.push_back(code_byte); // symbol 0
    code_byte = ((packed_callsign_1 >> 16) & (0x3F)); // 6 bits from callsign_1
    message.push_back(code_byte); // symbol 1
    code_byte = ((packed_callsign_1 >> 10) & (0x3F)); // 6 bits from callsign_1
    message.push_back(code_byte); // symbol 2
    code_byte = ((packed_callsign_1 >> 4) & (0x3F)); // 6 bits from callsign_1
    message.push_back(code_byte); // symbol 3
    code_byte = ((packed_callsign_1 & 0x0F) << 2) | ((packed_callsign_2 >> 26) & (0x03));  // 4 bits from callsign_1 and 2 bits from callsign_2
    message.push_back(code_byte); // symbol 4
    code_byte = ((packed_callsign_2 >> 20) & (0x3F)); // 6 bits from callsign_2
    message.push_back(code_byte); // symbol 5
    code_byte = ((packed_callsign_2 >> 14) & (0x3F)); // 6 bits from callsign_2
    message.push_back(code_byte); // symbol 6
    code_byte = ((packed_callsign_2 >> 8) & (0x3F)); // 6 bits from callsign_2
    message.push_back(code_byte); // symbol 7
    code_byte = ((packed_callsign_2 >> 2) & (0x3F)); // 6 bits from callsign_2
    message.push_back(code_byte); // symbol 8
    code_byte = ((packed_callsign_2 & 0x03) << 4) | ((packed_locator >> 12) & (0x0F)); // 2 bits from callsign_2 and 4 bits from 16 bit locator
    message.push_back(code_byte); // symbol 9
    code_byte = ((packed_locator >> 6) & (0x3F)); // 6 bits from locator
    message.push_back(code_byte); // symbol 10
    code_byte = (packed_locator & (0x3F)); // 6 bits from locator
    message.push_back(code_byte); // symbol 11
    
    return true;
}

//=================================================================================================
bool SourceCodec_JT65::unpack_message(const std::vector<unsigned int>& message,
        unsigned int& packed_callsign_1,
        unsigned int& packed_callsign_2,
        unsigned int& packed_locator,
        bool& arbitrary_text) const
{
    packed_callsign_1 = 0;
    packed_callsign_2 = 0;
    packed_locator = 0;
    
    if (message.size() < 12)
    {
        return false;
    }
    else
    {
        packed_callsign_1 = message[0] << 22;
        packed_callsign_1 += message[1] << 16;
        packed_callsign_1 += message[2] << 10;
        packed_callsign_1 += message[3] << 4;
        packed_callsign_1 += (message[4] >> 2) & 0x0F;
        
        packed_callsign_2 = (message[4] & 0x03) << 26;
        packed_callsign_2 += message[5] << 20;
        packed_callsign_2 += message[6] << 14;
        packed_callsign_2 += message[7] << 8;
        packed_callsign_2 += message[8] << 2;
        packed_callsign_2 += (message[9] >> 4) & 0x03;
        
        packed_locator = (message[9] & 0x0F) << 12;
        packed_locator += message[10] << 6;
        packed_locator += message[11];
        
        arbitrary_text = (packed_locator > 0x7FFF);
    }
}

//=================================================================================================
bool SourceCodec_JT65::pack_message_256(unsigned int packed_callsign_1,
        unsigned int packed_callsign_2,
        unsigned int packed_locator,
        std::vector<unsigned int>& message) const
{
    // 28 + 28 + 16 bits
    unsigned int code_byte = 0;
    // bits are entered MSB first
    // locator is 16 bit including text indicator bit on MSB
    code_byte = ((packed_callsign_1 >> 20) & (0xFF)); // 8 bits from callsign_1
    message.push_back(code_byte); // symbol 0
    code_byte = ((packed_callsign_1 >> 12) & (0xFF)); // 8 bits from callsign_1
    message.push_back(code_byte); // symbol 1
    code_byte = ((packed_callsign_1 >> 4) & (0xFF)); // 8 bits from callsign_1
    message.push_back(code_byte); // symbol 2
    code_byte = ((packed_callsign_1 & 0x0F) << 4) | ((packed_callsign_2 >> 24) & (0x0F));  // 4 bits from callsign_1 and 4 bits from callsign_2
    message.push_back(code_byte); // symbol 3
    code_byte = ((packed_callsign_2 >> 16) & (0xFF)); // 8 bits from callsign_2
    message.push_back(code_byte); // symbol 4
    code_byte = ((packed_callsign_2 >> 8) & (0xFF)); // 8 bits from callsign_2
    message.push_back(code_byte); // symbol 5
    code_byte = ((packed_callsign_2) & (0xFF)); // 8 bits from callsign_2
    message.push_back(code_byte); // symbol 6
    code_byte = ((packed_locator >> 8) & (0xFF)); // 8 bits from locator
    message.push_back(code_byte); // symbol 7
    code_byte = (packed_locator & (0xFF)); // 8 bits from locator
    message.push_back(code_byte); // symbol 8
    
    return true;
}

//=================================================================================================
bool SourceCodec_JT65::unpack_message_256(const std::vector<unsigned int>& message,
        unsigned int& packed_callsign_1,
        unsigned int& packed_callsign_2,
        unsigned int& packed_locator,
        bool& arbitrary_text) const
{
    packed_callsign_1 = 0;
    packed_callsign_2 = 0;
    packed_locator = 0;
    
    if (message.size() < 9)
    {
        return false;
    }
    else
    {
        packed_callsign_1 = message[0] << 20;
        packed_callsign_1 += message[1] << 12;
        packed_callsign_1 += message[2] << 4;
        packed_callsign_1 += (message[3] >> 4) & 0x0F;
        
        packed_callsign_2 = (message[3] & 0x0F) << 24;
        packed_callsign_2 += message[4] << 16;
        packed_callsign_2 += message[5] << 8;
        packed_callsign_2 += message[6];
        
        packed_locator = message[7] << 8;
        packed_locator += message[8];
        
        arbitrary_text = (packed_locator > 0x7FFF);
    }
}

//=================================================================================================
void SourceCodec_JT65::print_source_codec_data(std::ostringstream& os) const
{
    switch (variant)
    {
        case JT65_Classical:
            os << "JT65 classic with 12 6-bits symbols";
            break;
        case JT65_256:
            os << "JT65 classic with 9 8-bits symbols";
            break;
        default:
            os << "JT65 unrecognized variant";
            break;
    }
}
