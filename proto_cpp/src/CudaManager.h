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

#include "cuComplex.h"

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
	 * \param options Processing options
	 */
	CudaManager(
			unsigned int nb_message_symbols,
			unsigned int nb_pilots,
			unsigned int nb_code_samples,
			unsigned int batch_size,
			unsigned int df_steps,
			unsigned int nb_prns_per_symbol,
			unsigned int f_step_division
			);
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
        
        CudaDeviceProfile() : _id(0), _gmemsize(0), _cpufreq(0.0), _nbcores(0), _pciBusID(0), _pciDeviceID(0), _pciDomainID(0) {}
        
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

    class WsgcMemoryProfile
    {
    public:
    	unsigned int _local_codes_matrix;
    	unsigned int _local_codes_code;
    	unsigned int _local_codes_fft_matrix;
    	unsigned int _local_codes_fft_code;
    	unsigned int _pil_msg_corr_corr_in;
    	unsigned int _pil_msg_corr_mul_out;
    	unsigned int _pil_msg_corr_corr_out;
    	unsigned int _pil_msg_corr_corr_out_avg;
    	unsigned int _sprn_corr_fdep_ifft_in;
    	unsigned int _sprn_corr_fdep_ifft_out;
    	unsigned int _source_fft_fft_in;
    	unsigned int _source_fft_fft_out;

    	WsgcMemoryProfile() :
        	_local_codes_matrix(0),
        	_local_codes_code(0),
        	_local_codes_fft_matrix(0),
        	_local_codes_fft_code(0),
        	_pil_msg_corr_corr_in(0),
        	_pil_msg_corr_mul_out(0),
        	_pil_msg_corr_corr_out(0),
        	_pil_msg_corr_corr_out_avg(0),
        	_sprn_corr_fdep_ifft_in(0),
        	_sprn_corr_fdep_ifft_out(0),
        	_source_fft_fft_in(0),
        	_source_fft_fft_out(0)
    	{}

    	virtual ~WsgcMemoryProfile()
    	{}

    	void dump(std::ostringstream& os) const;
    };

	unsigned int _nb_message_symbols;
	unsigned int _nb_pilots;
	unsigned int _nb_code_samples;
	unsigned int _complex_size;
	unsigned int _nb_prns_per_symbol;
	unsigned int _batch_size;
	unsigned int _f_step_division;
	unsigned int _df_steps;
	unsigned int _nb_cuda_devices;
	unsigned int _pilot1_cuda_device;
	unsigned int _pilot2_cuda_device;
	unsigned int _nb_message_devices;
	unsigned int _message_first_device;
    std::vector<CudaDeviceProfile> _device_profiles;
	std::vector<unsigned int> _message_prn_allocation;
	WsgcMemoryProfile _wsgc_memory_profile;
    
  	/**
	 * Create device profiles
	 */
	void make_device_profiles();

    /**
     * Prints devices profiles to output stream
     */
    void dump_device_info(std::ostringstream& os) const;

protected:
	void analyze_memory_profile();
};

#endif /* __CUDA_MANAGER_H__ */
