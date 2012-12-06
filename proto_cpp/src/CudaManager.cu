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

#include "CudaManager.h"
#include "WsgcUtils.h"
#include <iostream>
// CUDA-C includes
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <shrUtils.h>

CudaManager::CudaManager(unsigned int nb_message_symbols, unsigned int nb_pilots) :
	_nb_message_symbols(nb_message_symbols+1), // include noise PRN
	_nb_pilots(nb_pilots),
	_nb_cuda_devices(0),
	_pilot1_cuda_device(0),
	_pilot2_cuda_device(0),
	_nb_message_devices(0),
	_message_first_device(0)
{
	if (_nb_pilots > 2)
	{
		_nb_pilots = 2;
	}
}


CudaManager::~CudaManager()
{
}


void CudaManager::diagnose()
{
    cudaError_t error_id = cudaGetDeviceCount((int *) &_nb_cuda_devices);

    if (error_id != cudaSuccess)
    {
        std::cout << "cudaGetDeviceCount returned %d\n-> %s\n" << (int)error_id << " -> " << cudaGetErrorString(error_id) << std::endl;
        _nb_cuda_devices = 0;
    }

    // TODO: check if device has enough memory
    // TODO: consider more in detail (devices already sorted in bogomips order) each device capabilities to allocate workload more evenly
    make_device_profiles();

    // Pilot allocation
    if (_nb_pilots > 0)
    {
    	if (_nb_cuda_devices > 0)
    	{
    		_pilot1_cuda_device = 0;
    	}
    }

    if (_nb_pilots > 1)
    {
    	if (_nb_cuda_devices > 0)
    	{
    		_pilot2_cuda_device = 0;
    	}
    	else if (_nb_cuda_devices > 1)
    	{
    		_pilot2_cuda_device = 1;
    	}
    }

    // Message allocation

    if (_nb_cuda_devices > _nb_pilots) // route message on devices not processing pilots
    {
    	_nb_message_devices = _nb_cuda_devices - _nb_pilots;
    	_message_first_device = _nb_pilots;
    }
    else // route messages on all devices evenly
    {
    	_nb_message_devices = _nb_cuda_devices;
    	_message_first_device = 0;
    }

    unsigned int orig_block_size = _nb_message_symbols / _nb_message_devices;
    unsigned int nb_extra_prns = _nb_message_symbols % orig_block_size;
    unsigned int message_device = _message_first_device;
    unsigned int nb_messages_for_device = 0;

    for (unsigned int prni=0; prni < _nb_message_symbols; prni++)
    {
    	if (nb_messages_for_device < orig_block_size)
    	{
    		_message_prn_allocation.push_back(message_device);
    		nb_messages_for_device++;
    	}
    	else
    	{
    		if (nb_extra_prns > 0)
    		{
    			_message_prn_allocation.push_back(message_device);
    			message_device++;
    			nb_messages_for_device = 0;
    			nb_extra_prns--;
    		}
    		else
    		{
    			message_device++;
    			nb_messages_for_device = 1;
    			_message_prn_allocation.push_back(message_device);
    		}
    	}
    }
}


void CudaManager::make_device_profiles()
{
    static const CudaDeviceProfile tmp_device_profile;

    for (unsigned int cuda_dev = 0; cuda_dev < _nb_cuda_devices; cuda_dev++)
    {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, cuda_dev);
        
        _device_profiles.push_back(tmp_device_profile);
        _device_profiles.back()._id = cuda_dev;
        _device_profiles.back()._name = deviceProp.name;
        _device_profiles.back()._gmemsize = deviceProp.totalGlobalMem;
        _device_profiles.back()._cpufreq = deviceProp.clockRate * 1e-3f;
        _device_profiles.back()._nbcores = ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) * deviceProp.multiProcessorCount;
        _device_profiles.back()._pciBusID = deviceProp.pciBusID;
        _device_profiles.back()._pciDeviceID = deviceProp.pciDeviceID;
        _device_profiles.back()._pciDomainID = deviceProp.pciDomainID;
    }
    
    std::sort(_device_profiles.begin(), _device_profiles.end(), CudaDeviceProfile::order);
}


void CudaManager::dump(std::ostringstream& os) const
{
	os << "There are " << _nb_cuda_devices << " CUDA devices available" << std::endl;
    dump_device_info(os);
    
	os << "Allocation for Pilots:" << std::endl;
	os << "Pilot1 .....: " << _pilot1_cuda_device << std::endl;

	if (_nb_pilots > 1)
	{
		os << "Pilot2 .....: " << _pilot2_cuda_device << std::endl;
	}

	os << "Allocation for Messages:" << std::endl;

	for (std::vector<unsigned int>::const_iterator it = _message_prn_allocation.begin(); it != _message_prn_allocation.end(); ++it)
	{
		if (it != _message_prn_allocation.begin())
		{
			os << ", ";
		}

		os << "[" << it-_message_prn_allocation.begin() << "]:" << *it;
	}
}


void CudaManager::dump_device_info(std::ostringstream& os) const
{
    std::vector<CudaDeviceProfile>::const_iterator it = _device_profiles.begin();
    const std::vector<CudaDeviceProfile>::const_iterator it_end = _device_profiles.end();
    
    CudaDeviceProfile::dump_header(os);
    
    for (; it != it_end; ++it)
    {
        it->dump(os);
    }
}
