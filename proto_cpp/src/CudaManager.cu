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

CudaManager::CudaManager(
		unsigned int nb_message_symbols,
		unsigned int nb_pilots,
		unsigned int nb_code_samples,
		unsigned int batch_size,
		unsigned int df_steps,
		unsigned int nb_prns_per_symbol,
		unsigned int f_step_division
		) :
	_nb_message_symbols(nb_message_symbols+1), // include noise PRN
	_nb_pilots(nb_pilots),
	_nb_code_samples(nb_code_samples),
	_complex_size(sizeof(cuComplex)),
	_float_size(sizeof(float)),
	_int_size(sizeof(int)),
	_batch_size(batch_size),
	_df_steps(df_steps),
	_nb_prns_per_symbol(nb_prns_per_symbol),
	_f_step_division(f_step_division),
	_nb_cuda_devices(0),
	_pilot_cuda_device(0),
	_nb_message_devices(0),
	_message_first_device(0),
	_gpu_affinity(0),
	_gpu_affinity_specified(false)
{
	if (_nb_pilots > 2)
	{
		_nb_pilots = 2;
	}
}


CudaManager::~CudaManager()
{
}


bool CudaManager::diagnose()
{
    cudaError_t error_id = cudaGetDeviceCount((int *) &_nb_cuda_devices);

    if (error_id != cudaSuccess)
    {
        std::cout << "cudaGetDeviceCount returned %d\n-> %s\n" << (int)error_id << " -> " << cudaGetErrorString(error_id) << std::endl;
        _nb_cuda_devices = 0;
        return false;
    }

   analyze_memory_profile();

    // TODO: check if device has enough memory
    // TODO: consider more in detail (devices already sorted in bogomips order) each device capabilities to allocate workload more evenly
    make_device_profiles();

    if (_gpu_affinity_specified)
    {
        if (_gpu_affinity >= _nb_cuda_devices)
        {
            std::cout << "GPU affinity exceeds the number of GPUs available. GPU affinity is ignored" << std::endl;
            _gpu_affinity_specified = false;
        }
    }
    
	// Pilot allocation
	if (_nb_pilots > 0)
	{
		if (_gpu_affinity_specified)
		{
			_pilot_cuda_device = _gpu_affinity;
		}
		else
		{
			_pilot_cuda_device = 0;
		}
	}

	// Message allocation

	if (_gpu_affinity_specified)
	{
		_nb_message_devices = 1;
		_message_first_device = _gpu_affinity;
	}
	else if ((_nb_pilots > 0) && (_nb_cuda_devices > 1)) // route message on devices not processing pilot
	{
		_nb_message_devices = _nb_cuda_devices - 1;
		_message_first_device = 1;
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

	return true;
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
    
    os << std::endl;
	os << "Allocation for Pilot:" << std::endl;
	os << "Pilot ......: " << _pilot_cuda_device << std::endl;

    os << std::endl;
	os << "Allocation for Messages:" << std::endl;

	for (std::vector<unsigned int>::const_iterator it = _message_prn_allocation.begin(); it != _message_prn_allocation.end(); ++it)
	{
		if (it != _message_prn_allocation.begin())
		{
			os << ", ";
		}

		os << "[" << it-_message_prn_allocation.begin() << "]:" << *it;
	}

	os << std::endl << std::endl;
	os << "Application CUDA memory profile:" << std::endl;
	_wsgc_memory_profile.dump(os);
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

void CudaManager::WsgcMemoryProfile::dump(std::ostringstream& os) const
{
	unsigned int total_pilot, total_message;

	os << "Pilot processing:" << std::endl;
	os << "LocalCodes_FFT matrix ................: " << std::right << std::setw(11) << _local_codes_fft_matrix << std::endl
	   << "LocalCodes_FFT code (t) ..............: " << std::right << std::setw(11) << _local_codes_fft_code << std::endl
	   << "SourceFFT in .........................: " << std::right << std::setw(11) << _source_fft_fft_in << std::endl
	   << "SourceFFT out.........................: " << std::right << std::setw(11) << _source_fft_fft_out << std::endl
	   << "Single PRN correlator fDep IFFT in ...: " << std::right << std::setw(11) << _sprn_corr_fdep_ifft_in << std::endl
	   << "Single PRN correlator fDep IFFT out ..: " << std::right << std::setw(11) << _sprn_corr_fdep_ifft_out << std::endl
	   << "Single PRN correlator fDep avg keys ..: " << std::right << std::setw(11) << _sprn_corr_fdep_avg_keys << std::endl;
	total_pilot = _local_codes_fft_matrix
			+ _local_codes_fft_code
			+ _source_fft_fft_in
			+ _source_fft_fft_out
			+ _sprn_corr_fdep_ifft_in
			+ _sprn_corr_fdep_ifft_out
            + _sprn_corr_fdep_avg_keys;
	os << "Total Pilot ..........................: " << std::right << std::setw(11) << total_pilot << std::endl;

	os << "Message processing:" << std::endl;
	os << "LocalCodes matrix ....................: " << std::right << std::setw(11) << _local_codes_matrix << std::endl
	   << "LocalCodes code (t) ..................: " << std::right << std::setw(11) << _local_codes_code << std::endl
	   << "Piloted Msg correlator in ............: " << std::right << std::setw(11) << _pil_msg_corr_corr_in << std::endl
	   << "Piloted Msg correlator mul ...........: " << std::right << std::setw(11) << _pil_msg_corr_mul_out << std::endl
	   << "Piloted Msg correlator corr ..........: " << std::right << std::setw(11) << _pil_msg_corr_corr_out << std::endl
	   << "Piloted Msg correlator mag ...........: " << std::right << std::setw(11) << _pil_msg_corr_corr_mag << std::endl
	   << "Piloted Msg correlator mag avgsum ....: " << std::right << std::setw(11) << _pil_msg_corr_corr_mag_avgsum << std::endl
	   << "Piloted Msg correlator corr keys .....: " << std::right << std::setw(11) << _pil_msg_corr_keys << std::endl;
	total_message = _local_codes_matrix
			+ _local_codes_code
			+ _pil_msg_corr_corr_in
			+ _pil_msg_corr_mul_out
			+ _pil_msg_corr_corr_out
			+ _pil_msg_corr_corr_mag
			+ _pil_msg_corr_corr_mag_avgsum
			+ _pil_msg_corr_keys;
	os << "Total Message ........................: " << std::right << std::setw(11) << total_message << std::endl;
    os << std::endl;
	os << "Total Pilot + Message ................: " << std::right << std::setw(11) << total_pilot + total_message << std::endl;
    
}

void CudaManager::analyze_memory_profile()
{
	_wsgc_memory_profile._local_codes_matrix = _nb_message_symbols*_nb_code_samples*_complex_size;
	_wsgc_memory_profile._local_codes_code = _nb_code_samples*_complex_size; // transient
	_wsgc_memory_profile._local_codes_fft_matrix = _nb_code_samples*_nb_pilots*_complex_size;
	_wsgc_memory_profile._local_codes_fft_code = _nb_code_samples*_nb_pilots*_complex_size; // transient

	_wsgc_memory_profile._pil_msg_corr_corr_in = _nb_code_samples*_complex_size; // _d_corr_in
	_wsgc_memory_profile._pil_msg_corr_mul_out = _nb_code_samples*_nb_message_symbols*_complex_size; // _d_mul_out
	_wsgc_memory_profile._pil_msg_corr_corr_out = _nb_message_symbols*_nb_prns_per_symbol*_complex_size; // _d_corr_out
	_wsgc_memory_profile._pil_msg_corr_corr_mag = _nb_message_symbols*_nb_prns_per_symbol*_float_size; // _d_corr_mag
	_wsgc_memory_profile._pil_msg_corr_corr_mag_avgsum = _nb_message_symbols*_nb_prns_per_symbol*_float_size; // _d_corr_mag_avgsum
	_wsgc_memory_profile._pil_msg_corr_corr_mag_avgsum_sums = _nb_prns_per_symbol*_float_size; // _d_corr_mag_avgsum_sums
	_wsgc_memory_profile._pil_msg_corr_keys = _nb_message_symbols*_int_size; // _d_keys

	_wsgc_memory_profile._sprn_corr_fdep_ifft_in = 2*_batch_size*_nb_code_samples*_f_step_division*_df_steps*_complex_size;
	_wsgc_memory_profile._sprn_corr_fdep_ifft_out = _wsgc_memory_profile._sprn_corr_fdep_ifft_in;
	_wsgc_memory_profile._sprn_corr_fdep_avg_keys = _batch_size*_nb_code_samples*_f_step_division*_df_steps*_int_size;

	_wsgc_memory_profile._source_fft_fft_in = _nb_code_samples*_f_step_division*_complex_size;
	_wsgc_memory_profile._source_fft_fft_out = _wsgc_memory_profile._source_fft_fft_in;
}


