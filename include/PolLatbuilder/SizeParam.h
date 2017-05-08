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

#ifndef POLLATBUILDER__SIZE_PARAM_H
#define POLLATBUILDER__SIZE_PARAM_H

#include "PolLatbuilder/Types.h"

#include <ostream>

namespace PolLatBuilder {

/**
 * Lattice size parameter.
 */
template <LatType> class SizeParam;

/**
 * Abstract CRTP base class for lattice size parameter.
 */
template <class DERIVED>
class BasicSizeParam {
protected:
   BasicSizeParam(const BasicSizeParam&) = default;

public:
   BasicSizeParam(Poly polynomial): m_polynomial(polynomial) { }

   Poly polynomial() const { return m_polynomial; }
   operator Poly() const { return m_polynomial; }

   /**
    * Returns the value of Euler's totient function.
    * It is the number of polynomials over GF2 that have a smaller degree than and coprime
    * to the polynomial.
    */
   size_t totient() const
   { return derived().totient(); }

   /**
    * Return the number of points of the Lattice in the unit cube
    */
   Modulus numPoints() const
   { return derived().numPoints(); }

   template <class D> bool operator== (const BasicSizeParam<D>& other) const { return polynomial() == other.polynomial(); }
   template <class D> bool operator!= (const BasicSizeParam<D>& other) const { return !operator==(other); }
   template <class D> bool operator<  (const BasicSizeParam<D>& other) const { return m_polynomial < other.m_polynomial; }

   /**
    * Divides the merit value \c merit by the number of points.
    */
   void normalize(Real& merit) const
   { derived().normalize(merit); }

  

private:
   Poly m_polynomial;

   template <class D>
   friend std::ostream& operator<<(std::ostream&, const BasicSizeParam<D>&);

   std::ostream& format(std::ostream& os) const
   { return derived().format(os); }

   DERIVED& derived()
   { return static_cast<DERIVED&>(*this); }

   const DERIVED& derived() const
   { return static_cast<const DERIVED&>(*this); }

};

template <class D>
std::ostream& operator<<(std::ostream& os, const BasicSizeParam<D>& sizeParam)
{ return sizeParam.format(os); }

}

#include "PolLatbuilder/SizeParam-ORDINARY.h"

#endif
