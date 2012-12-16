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

    Collection of operators used in CUDA Thrust transforms

*/

#ifndef __CUDA_OPERATORS_H__
#define __CUDA_OPERATORS_H__

#include <cutil_inline.h>    // includes cuda.h and cuda_runtime_api.h
#include <thrust/tuple.h>
#include <thrust/functional.h>
#include <cuComplex.h>

#include <math.h>
#include <iostream>

/**
 * \brief average sum functor Complex
 */
template <size_t N>
struct avgsum_functor_complex
{
    template <typename Tuple>
    __host__ __device__
    void operator()(Tuple t)
    {
        // do nothing
    }
};

/**
 * \brief average sum (depth 4) functor Complex
 */
template <>
struct avgsum_functor_complex<4>
{
    template <typename Tuple>
    __host__ __device__
    void operator()(Tuple t)
    {
        // Tuple is (A, B, C)
        // C0[i] = C0[i] + C1[i] + C2[i] + C3[i];
        //thrust::get<0>(t).x = thrust::get<0>(t).x + thrust::get<1>(t).x + thrust::get<2>(t).x + thrust::get<3>(t).x;
        //thrust::get<0>(t).y = thrust::get<0>(t).y + thrust::get<1>(t).y + thrust::get<2>(t).y + thrust::get<3>(t).y;
        thrust::get<0>(t).x += thrust::get<1>(t).x + thrust::get<2>(t).x + thrust::get<3>(t).x;
        thrust::get<0>(t).y += thrust::get<1>(t).y + thrust::get<2>(t).y + thrust::get<3>(t).y;
    }
};

/**
 * \brief average sum (depth 6) functor Complex
 */
template <>
struct avgsum_functor_complex<6>
{
    template <typename Tuple>
    __host__ __device__
    void operator()(Tuple t)
    {
        // Tuple is (A, B, C)
        // C0[i] = C0[i] + C1[i] + C2[i] + C3[i] + C4[i] + C5[i];
        thrust::get<0>(t).x = thrust::get<0>(t).x + thrust::get<1>(t).x + thrust::get<2>(t).x + thrust::get<3>(t).x + thrust::get<4>(t).x + thrust::get<5>(t).x;
        thrust::get<0>(t).y = thrust::get<0>(t).y + thrust::get<1>(t).y + thrust::get<2>(t).y + thrust::get<3>(t).y + thrust::get<4>(t).y + thrust::get<5>(t).y;
    }
};

/**
 * \brief average sum (depth 8) functor Complex
 */
template <>
struct avgsum_functor_complex<8>
{
    template <typename Tuple>
    __host__ __device__
    void operator()(Tuple t)
    {
        // Tuple is (A, B, C)
        // C0[i] = C0[i] + C1[i] + C2[i] + C3[i] + C4[i] + C5[i] + C6[i] + C7[i];
        thrust::get<0>(t).x = thrust::get<0>(t).x + thrust::get<1>(t).x + thrust::get<2>(t).x + thrust::get<3>(t).x + thrust::get<4>(t).x + thrust::get<5>(t).x + thrust::get<6>(t).x + thrust::get<7>(t).x;
        thrust::get<0>(t).y = thrust::get<0>(t).y + thrust::get<1>(t).y + thrust::get<2>(t).y + thrust::get<3>(t).y + thrust::get<4>(t).y + thrust::get<5>(t).y + thrust::get<6>(t).y + thrust::get<7>(t).y;
    }
};

/**
 * \brief Complex member to member multiplication
 */
struct cmulc_functor
{
    template <typename Tuple>
    __host__ __device__
    void operator()(Tuple t)
    {
        // Tuple is (A, B, C)
        // C[i] = A[i] * B[i];
        thrust::get<2>(t) = cuCmulf(thrust::get<0>(t), thrust::get<1>(t));
    }
};

/**
 * \brief Complex member to member multiplication
 */
struct cmulc_functor2
{
    __host__ __device__
    cuComplex operator()(const cuComplex& a, const cuComplex& b)
    {
        return cuCmulf(a, b);
    }
};

/**
 * \brief Complex member to member multiplication followed by conjugation
 */
struct cmulconjc_functor
{
    __host__ __device__
    cuComplex operator()(const cuComplex& x, const cuComplex& y)
    {
        return make_cuFloatComplex  ((cuCrealf(x) * cuCrealf(y)) - 
                                     (cuCimagf(x) * cuCimagf(y)),
                                     - (cuCrealf(x) * cuCimagf(y)) - 
                                     (cuCimagf(x) * cuCrealf(y)));
    }
};

/**
 * \brief Complex addition
 */
struct caddc_functor
{
    __host__ __device__
    cuComplex operator()(const cuComplex& a, const cuComplex& b)
    {
        return cuCaddf(a, b);
    }
};

/**
 * \brief Complex conjugate functor
 */
struct conj_functor
{
    __host__ __device__
    cuComplex operator()(cuComplex& z)
    {
        return cuConjf(z);
    }
};

/**
 * \brief Magnitude estimator
 */
struct mag_estim_functor
{
    static const float magnitude_estimation_alpha = 0.948059448969f;
    static const float magnitude_estimation_beta  = 0.392699081699f;
    
    __host__ __device__
    float operator()(cuComplex z)
    {
        /* magnitude ~= alpha * max(|I|, |Q|) + beta * min(|I|, |Q|) */

        float abs_inphase = fabs(z.x);
        float abs_quadrature = fabs(z.y);

        if (abs_inphase > abs_quadrature) 
        {
            return magnitude_estimation_alpha * abs_inphase + magnitude_estimation_beta * abs_quadrature;
        } 
        else 
        {
            return magnitude_estimation_alpha * abs_quadrature + magnitude_estimation_beta * abs_inphase;
        }
    }
};
 
/**
 * \brief Exact squared magnitude
 */
template<typename T2, typename T>
struct mag_squared_functor : thrust::unary_function<T2, T>
{
    __host__ __device__
    T operator()(const T2& z)
    {
        return z.x*z.x + z.y*z.y;
    }
};


/**
 * \brief Algebraic magnitude
 */
struct mag_algebraic_functor
{
    __host__ __device__
    float operator()(cuComplex z)
    {
        return fabs(z.x) + fabs(z.y);
    }
};


/**
 * \brief Compare squared magnitudes (less)
 */
template<typename T2>
struct lesser_mag_squared : thrust::binary_function<T2, T2, bool>
{
    __host__ __device__
    bool operator()(T2 zl, T2 zr)
    {
        return (zl.x*zl.x + zl.y*zl.y) < (zr.x*zr.x + zr.y*zr.y);
    }
};


/**
 * \brief Returns <complex, int> tuple of maximum magnitude using exact squaring
 */
struct bigger_mag_tuple
{
    __host__ __device__
    thrust::tuple<cuComplex, int> operator()(thrust::tuple<cuComplex, int>& t_l, thrust::tuple<cuComplex, int>& t_r)
    {
        if (mag_squared_functor<cuComplex, float>()(thrust::get<0>(t_l)) > mag_squared_functor<cuComplex, float>()(thrust::get<0>(t_r)))
        {
            return t_l;
        }
        else
        {
            return t_r;
        }
    }
};


/**
 * \brief Returns <complex, int> tuple of maximum magnitude using estimated squaring
 */
struct bigger_mag_tuple_estimate
{
    __host__ __device__
    thrust::tuple<cuComplex, int> operator()(thrust::tuple<cuComplex, int>& t_l, thrust::tuple<cuComplex, int>& t_r)
    {
        if (mag_estim_functor()(thrust::get<0>(t_l)) > mag_estim_functor()(thrust::get<0>(t_r)))
        {
            return t_l;
        }
        else
        {
            return t_r;
        }
    }
};


/**
 * Brief null operator
 */
struct null_operator
{
    template <typename Tuple>
    __host__ __device__
    void operator()(Tuple t)
    {
        // do nothing
    }
};


/**
 * Print a cuComplex to output stream
 */
std::ostream& operator<<(std::ostream& os, const cuComplex& z);

#endif // __CUDA_OPERATORS_H__
