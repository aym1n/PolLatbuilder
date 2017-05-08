// This file is part of Lattice Builder.
//
// Copyright (C) 2012-2016  Pierre L'Ecuyer and Universite de Montreal
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

/**
 * \file
 * Tools for streaming and poor man's factorization.
 */

#ifndef POLLATBUILDER__UTIL_H
#define POLLATBUILDER__UTIL_H

#include "PolLatbuilder/Types.h"

#include <map>
#include <vector>

//================================================================================

namespace PolLatBuilder
{

//================================================================================

/**
 * Integer exponentiation.
 *
 * Source:
 * http://en.wikipedia.org/wiki/Exponentiation_by_squaring#Computation_by_powers_of_2
 */
template <typename T>
T intPow(T base, unsigned long exponent)
{
   T result = (T) (1);
   while (exponent) {
      if (exponent % 2 == 1)
         result *= base;
      base *= base;
      exponent /= 2;
   }
   return result;
}


//================================================================================

/**
 * convert Integer to polynomial  
 *
 * the integer $a_{0} + a_{1}2 +... + a_{n}2^{n} is converted to $a_{0} + a_{1}X +... + a_{n}X^{n}$
 * Note that n must be <64
 */
Poly intToPoly(Modulus x);

/**
 * Modular exponentiation.
 *
 * Source:
 * http://en.wikipedia.org/wiki/Modular_exponentiation#Right-to-left_binary_method
 */
Modulus modularPow(Modulus base, Modulus exponent, Modulus modulus);

/**
 * Prime factorization using the naive "trial division" algorithm.
 *
 * Returns a list of prime factors, without their multiplicity if \c raise is
 * \c false, or raised to their multiplicity if it is \c true.
 */
std::vector<Modulus> primeFactors(Modulus n, bool raise = false);

/**
 * Prime factorization using the naive "trial division" algorithm.
 *
 * Returns a map of (factor, multiplicity) pairs.
 */
std::map<Modulus, Modulus> primeFactorsMap(Modulus n);

/**
 * Extended Euclidian algorithm.
 *
 * Source:
 * http://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Iterative_method_2
 */
std::pair<long long, long long> egcd(Modulus a, Modulus b);

}

#endif
