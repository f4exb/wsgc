#include "WsgcTypes.h"
#include <iostream>
#include <sstream>
#include <vector>

typedef struct options_s
{
	unsigned int nb_message_symbols;
	unsigned int nb_service_symbols;
	unsigned int nb_pilot_prns;
	wsgc_float f_sampling;
	wsgc_float f_chip;
	unsigned int fft_N;
	unsigned int freq_step_division;
	unsigned int gc_nb_stages;
	std::vector<unsigned int> g1_poly_powers;
	std::vector<unsigned int> g2_poly_powers;
	std::vector<unsigned int> prn_list;
	wsgc_float f_tx;
	unsigned int code_shift;
	unsigned int nb_batch_prns;
	unsigned int nb_f_bins;
} options_t;

class CudaManager;
class GoldCodeGenerator;
class CodeModulator_BPSK;

class test_cuda
{
public:
	test_cuda(options_t&);
	virtual ~test_cuda();
	void test1();
	void test2(wsgc_complex *message_samples, GoldCodeGenerator&, CodeModulator_BPSK&);
	void test3(wsgc_complex *message_samples, GoldCodeGenerator&, CodeModulator_BPSK&);
	void test4(wsgc_complex *message_samples, GoldCodeGenerator&, CodeModulator_BPSK&);
protected:
	void decomp_full_index(unsigned int full_index, unsigned int& bi, unsigned int& ffti, unsigned int& fsi, unsigned int& fhi);
	void decomp_strided_index(unsigned int full_index, unsigned int& ffti, unsigned int& fsi, unsigned int& fhi);
	options_t& _options;
	CudaManager* _cuda_manager;
};

