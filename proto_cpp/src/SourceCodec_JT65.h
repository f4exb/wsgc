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
     vector of 12 6-bit symbols as unsigned integers.
     
     
     JT65 messages can have one of three basic formats:
        1. Two to four alphanumeric fields with specific contents, as described below
        2. Any other arbitrary text, up to 13 characters
        3. Special shorthand messages RO, RRR, and 73
     
     
     Possible non free-text messages type-1 messages:
        <field_1> <field_2> <field_3> <field_4>
        
        <field_1>:
            {CQ,QRZ,DE} |
            CQ <frequency> |
            <callsign>
            
        <field_2>:
            <callsign>
            
        <field_3>:
            <locator> |
            <report> |
            {RO,RRR,73} |
            none 
            
        <field-4>:
            OOO |
            none
         
        <frequency>: 3 digit usually in kHz
        <callsign>: (<country prefix>/)<country prefix>[0-9]<call body>(/[P,0-9]) 
            - Country prefixes are predefined and statically stored
            - (<country prefix>/) is the optional prefix
            - <call body> has a maximum of 3 alpha only characters
            - (/[P,0-9]) is the optional suffix
        <report>: {R}-NN NN is a numerical value between 01 and 30
        
        Note: if a callsign has a prefix and/or a suffix this information is coded in place of <field_3> which is ignored in this case
     
     
     Type-2 messages:
        Anything not falling into the two other types is treated as arbitrary free text with only 
        the 13 first characters being considered. Valid characters for arbitrary text are: 
        '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ +-./?'
     
     
     Type-3 "shorthand" messages (not treated here):
         {RO,RRR,73}

	 Differences with original JT65:

     1. At least 3 fields should be given even if some may get reused according to the protocol.
     Ex: CQ F4EXB/P JN33 => JN33 is discarded to encode the /P. In original JT65 you may write CQ F4EXB/P or CQ F4EXB.
     
*/
#ifndef __SOURCE_CODEC_JT65_H__
#define __SOURCE_CODEC_JT65_H__

#include "SourceCodec.h"

class StandardPrefixes
{
public:
	StandardPrefixes();
	~StandardPrefixes();

	const std::string& get_prefix(unsigned int index, bool& found) const;
	unsigned int get_index(const std::string& prefix, bool& found) const;

protected:
    static const char *prefix_vinit[];
    static const std::vector<std::string> prefixes;
};

class SourceCodec_JT65 : public SourceCodec
{
public:
    SourceCodec_JT65();
    virtual ~SourceCodec_JT65();
    
    /**
     * Coder part 
     * \param in_msg Input textual message
     * \param out_msg Output vector of symbols
     * \return true if successful
     */
    virtual bool encode(const std::string& in_msg, std::vector<unsigned int>& out_msg) const;
    
    /**
     * Decoder part
     * \param in_msg Input vector of symbols
     * \param out_msg Output textual message
     * \return true if successful
     */
    virtual bool decode(const std::vector<unsigned int>& in_msg, std::string& out_msg) const;

protected:
    /**
     * Pack callsign into 28 bits and process possible prefix or suffix
     * \param callsign Input textual callsign
     * \param packed_callsign Output packed callsign
     * \param pfxsfx_index Prefix/Suffix index
     * \param pfxsfx_v2 Prefix/Suffix in v.2 format: 0: not v.2, 1: v.2 prefix, 2: v.2 suffix
     * \return true if successful
     */
    bool pack_callsign(const std::string& callsign, unsigned int& packed_callsign, int& pfxsfx_index, unsigned int& pfxsfx_v2) const;
    
    /**
     * Pack callsign with no prefix nor suffix into 28 bits
     * \param callsign Input textual callsign
     * \param packed_callsign Output packed callsign
     * \return true if successful
     */
    bool pack_plain_callsign(const std::string& callsign, unsigned int& packed_callsign) const;

    /**
     * Unpack 28 bits first callsign into textual field
     * \param packed_callsign Input packed callsign
     * \param callsign output textual first field 
     * \param pfxsfx_str Prefix or suffix v.2 string
     * \param pfxsfx_v2 Prefix or suffix v.2 indicator
     * \return true if successful
     */
    bool unpack_callsign_1(unsigned int packed_callsign, std::string& field_1, std::string& pfxsfx_str, unsigned int& pfxsfx_v2) const;

    /**
     * Unpack 28 bits plain callsign into textual callsign
     * \param packed_callsign Input packed callsign
     * \param packed_callsign Output packed callsign
     * \return true if successful
     */
    bool unpack_plain_callsign(unsigned int packed_callsign, std::string& callsign) const;

    /**
     * Unpack 15 bits locator into textual field
     * \param packed_locator Input packed locator
     * \param field_3 output textual third field 
     * \param pfxsfx_str Prefix or suffix v.1 string
     * \param pfxsfx_v1 Prefix or suffix v.1 indicator: 0: no prefix nor suffix, 1: prefix for first callsign, 2: prefix for second callsign, 3: suffix for first callsign, 4: suffix for second callsign 
     * \return true if successful
     */
    bool unpack_locator(unsigned int packed_locator, std::string& field_3, std::string& pfxsfx_str, unsigned int& pfxsfx_v1) const;

    /**
     * Pack 4 character locator into 15 bits
     * \param locator Input textual locator
     * \param packed_locator Output packed locator
     * \return true if successful
     */
    bool pack_locator_4(const std::string& locator, unsigned int& packed_locator) const;
    
    /**
     * Pack standard prefix or standard suffix as a fake locator (known as version 1)
     * \param pfxsfx_index Prefix or suffix index
     * \param packed_locator Output packed locator
     * \return true if successful
     */
    bool pack_pfxsfx_v1(int pfxsfx_index, unsigned int& packed_locator) const;

    /**
     * Pack random 4 character prefix or random 3 character suffix as a fake callsign (known as version 2)
     * \param v2_indicator: 0 no v2 prefix/suffix, 1: prefix, 2: suffix
     * \param pfxsfx_index Prefix or suffix index
     * \param call_indicator: 0 for CQ, 1 for QRZ, 2 for DE
     * \param packed_callsign Output packed callsign
     * \return true if successful
     */
    bool pack_pfxsfx_v2(unsigned int v2_indicator, int pfxsfx_index, unsigned int call_indicator, unsigned int& packed_callsign) const;

    /**
     * Unpack v.2 packed prefix or suffix into random 4 character prefix or random 3 character suffix
     * \param pfxsfx_index Prefix or suffix index
     * \param pfxsfx_str Output textual prefix or suffix
     * \param max_chars Maximum number of characters expected
     * \return true if successful
     */
    bool unpack_pfxsfx_v2(int pfxsfx_index, std::string& pfxsfx_str, unsigned int max_chars) const;

    /**
     * Pack report in JT format (R-NN or -NN)
     * \param r_prefix true if "R" prefix was present
     * \param report_str 2 character "NN" string
     * \param packed_locator Output packed fake locator
     * \return true if successful
     */
    bool pack_report_jt(bool r_prefix, const std::string& report_str, unsigned int& packed_locator) const;

    /**
     * Pack 13 character arbitrary text into 2 packed callsigns and a packed locator resulting in 71 bits
     * \param text Input text
     * \param packed_callsign_1 Output first packed callsign
     * \param packed_callsign_2 Output second packed callsign
     * \param packed_locator Output packed locator
     * \return true if successful
     */
    bool pack_text(const std::string& text,
        unsigned int& packed_callsign_1,
        unsigned int& packed_callsign_2,
        unsigned int& packed_locator) const;
    
    /**
     * Unpack 2 packed callsigns and a packed locator into 13 character arbitrary text
     * \param text Output text
     * \param packed_callsign_1 Input first packed callsign
     * \param packed_callsign_2 Input second packed callsign
     * \param packed_locator Input packed locator
     * \return true if successful
     */
    bool unpack_text(std::string& text,
        unsigned int packed_callsign_1,
        unsigned int packed_callsign_2,
        unsigned int packed_locator) const;
        
    /**
     * Pack 2 packed callsigns and a packed locator into 12 6-bit symbols message
     * \param packed_callsign_1 Input first packed callsign
     * \param packed_callsign_2 Input second packed callsign
     * \param packed_locator Input packed locator
     * \param message Output 12 6-bit symbols message
     * \return true if successful
     */
    bool pack_message(unsigned int packed_callsign_1,
        unsigned int packed_callsign_2,
        unsigned int packed_locator,
        std::vector<unsigned int>& message) const;
        
    /**
     * Unpack 12 6-bit symbols message into 2 packed callsigns and a packed locator
     * \param message Input 6-bit symbols message
     * \param packed_callsign_1 Output first packed callsign
     * \param packed_callsign_2 Output second packed callsign
     * \param packed_locator Output packed locator
     * \param arbitrary_text true if this is arbitrary text
     * \return true if successful
     */
    bool unpack_message(const std::vector<unsigned int>& message,
        unsigned int& packed_callsign_1,
        unsigned int& packed_callsign_2,
        unsigned int& packed_locator,
        bool& arbitrary_text) const;
        
    /**
     * Compose the decoded formatted message from its decoded elements
     * \param out_msg formatted message
     * \param field_1 Field 1 string
     * \param field_2 Field 2 string
     * \param field_3 Field 3 string
     * \param pfxsfx_v1 Prefix or suffix v.1 indicator
     * \param pfxsfx_v2 Prefix or suffix v.2 indicator
     * \param pfxsfx_str_v1 Prefix or suffix v.1 string
     * \param pfxsfx_str_v2 Prefix or suffix v.2 string
     */
    void compose_formatted_message(std::string& out_msg,
        const std::string& field_1,
        const std::string& field_2,
        const std::string& field_3,
        unsigned int pfxsfx_v1,
        unsigned int pfxsfx_v2,
        const std::string& pfxsfx_str_v1,
        const std::string& pfxsfx_str_v2) const;
        
        
    unsigned int long_lat_to_pfxsfx_index(unsigned int nlat, unsigned int nlong) const;

    static const unsigned int call_nbase;
    static const unsigned int locator_nbase;
    static const unsigned int call_DE;
    static const std::string free_text_chars;
    static const std::string suffix_chars;
    static const std::string callsign_chars;
    static const unsigned int v2_prefix_cq_shift;
    static const unsigned int v2_prefix_qrz_shift;
    static const unsigned int v2_prefix_de_shift;
    static const unsigned int v2_suffix_cq_shift;
    static const unsigned int v2_suffix_qrz_shift;
    static const unsigned int v2_suffix_de_shift;
    StandardPrefixes standard_prefixes;
};

#endif //__SOURCE_CODEC_JT65_H__
