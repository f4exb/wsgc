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
      
     Exception class for WSGC prototype.
      
*/
#ifndef __WSGC_EXCEPTION_H__
#define __WSGC_EXCEPTION_H__

/**
 * \brief Generic exception class for WSGC prototype.
 */
class WsgcException : public std::exception
{
    public:
        /**
         * Public constructor 
         * \param strError Error message to be returned by the what() method
         */
        WsgcException(std::string strError) : _error_msg(strError) {}
        virtual ~WsgcException() throw() {}
        /**
         * Returns the reason for exception that was specified at exception construction time
         * \return Error message as std::string
         */
        virtual const char* what() const throw()
        {
            return _error_msg.c_str();
        }
        
    protected:
        std::string _error_msg; //!< Error reason shown with what() method
    
    private:
        /**
         *  Default constructor is not meant to be used
         */         
        WsgcException() {}; 
};

#endif // __WSGC_EXCEPTION_H__
