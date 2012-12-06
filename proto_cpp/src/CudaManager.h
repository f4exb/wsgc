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

     CUDA Manager

     Diagnoses CUDA installation and allocates devices for the different functions
*/

#ifndef __CUDA_MANAGER_H__
#define __CUDA_MANAGER_H__

#include <sstream>
#include <vector>
#include <string>
#include <iomanip>

/**
 * \brief Diagnoses CUDA installation and allocates devices for the different functions
 */
class CudaManager
{
public:
	/**
	 * Constructor
	 * \param nb_message_symbols Number of message symbols to explore
	 * \param nb_pilots Number of pilot PRNs to explore (limited to 0,1,2)
	 */
	CudaManager(unsigned int nb_message_symbols, unsigned int nb_pilots);
	virtual ~CudaManager();

	/**
	 * Do CUDA system diagnosis and allocates resources
	 */
	void diagnose();

	/**
	 * Dump CUDA diagnosis and allocation data
	 */
	void dump(std::ostringstream& os) const;

	/**
	 * Get pilot device
	 * \param alternate True for pilot 2
	 * \return CUDA device ID (index start 0)
	 */
	unsigned int get_pilot_device(bool alternate) const
	{
		if ((alternate) && (_nb_pilots > 1))
		{
			return _pilot2_cuda_device;
		}
		else
		{
			return _pilot1_cuda_device;
		}
	}

	/**
	 * Get message device
	 * \param prni Message PRN index in Gold Code sequences
	 * \return CUDA device ID (index start 0)
	 */
	unsigned int get_message_device(unsigned int prni) const
	{
		return _message_prn_allocation[prni];
	}

protected:
    class CudaDeviceProfile
    {
    public:
        unsigned int       _id;          //!< Device sequential ID
        unsigned long long _gmemsize;    //!< Global memory size in bytes
        float              _cpufreq;     //!< Clock rate in MHz
        unsigned int       _nbcores;     //!< Number of cores
        std::string        _name;        //!< CUDA device name
        unsigned int       _pciBusID;    //!< PCI bus ID of the device
        unsigned int       _pciDeviceID; //!< PCI device ID of the device
        unsigned int       _pciDomainID; //!< PCI domain ID of the device
        
        CudaDeviceProfile() : _id(0), _gmemsize(0), _cpufreq(0.0), _nbcores(0) {}
        
        void dump(std::ostringstream& os) const
        {
            os << " " << _id
               << " " << std::left << std::setw(20) << _name
               << " " << std::right  << std::setw(4) << _nbcores
               << " " << std::setw(4) << _cpufreq
               << " " << std::setw(7) << get_bogomips()
               << " " << std::setw(11) << _gmemsize
               << " " << _pciDomainID << "." << _pciBusID << "." << _pciDeviceID
               << std::endl;
        }
        
        static void dump_header(std::ostringstream& os) 
        {
            os << "Id Name................ Cor# fMHz Bogo... Mem(B)..... PCI..." << std::endl;
        }
        
        unsigned int get_bogomips() const
        {
            return _cpufreq * _nbcores;
        }
        
        static bool order(const CudaDeviceProfile& left, const CudaDeviceProfile& right)
        {
            if (left.get_bogomips() == right.get_bogomips())
            {
                return left._id < right._id;
            }
            else
            {
                return left.get_bogomips() > right.get_bogomips();
            }
        }
    };

	unsigned int _nb_message_symbols;
	unsigned int _nb_pilots;
	unsigned int _nb_cuda_devices;
	unsigned int _pilot1_cuda_device;
	unsigned int _pilot2_cuda_device;
	unsigned int _nb_message_devices;
	unsigned int _message_first_device;
    std::vector<CudaDeviceProfile> _device_profiles;
	std::vector<unsigned int> _message_prn_allocation;
    
  	/**
	 * Create device profiles
	 */
	void make_device_profiles();

    /**
     * Prints devices profiles to output stream
     */
    void dump_device_info(std::ostringstream& os) const;
};

#endif /* __CUDA_MANAGER_H__ */
