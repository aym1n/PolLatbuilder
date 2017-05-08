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

#ifndef PolLATBUILDER__GENSEQ__COPRIME_POLYNOMIALS_H
#define PolLATBUILDER__GENSEQ__COPRIME_POLYNOMIALS_H

#include <NTL/GF2XFactoring.h>
#include "PolLatbuilder/Util.h"
#include "PolLatbuilder/Traversal.h"
#include "PolLatbuilder/CompressTraits.h"

#include <cstdlib>
#include <vector>
#include <limits>

namespace PolLatBuilder { namespace GenSeq {

namespace detail {
   struct CoprimePolynomialsBasisElement {
      Modulus totient;
      Modulus leap;
      Poly irreductible_poly;
      Poly elem; 
      CoprimePolynomialsBasisElement(Modulus t = 0, Modulus l = 0, Poly irr = Poly(0), Poly e = Poly(0)):
         totient(t), leap(l), irreductible_poly(irr), elem(e) {}
   };
}

/**
 * ayman
 * To Do: documentation about chinese theorem and index formula for polynomials
 *
 * 
 *
 * \tparam COMPRESS  Type of compression.
 * \tparam TRAV      Traversal policy.
 *
 * \sa LatBuilder::Compress
 */
template <Compress COMPRESS = Compress::NONE,
         class TRAV = Traversal::Forward>
class CoprimePolynomials :
   public Traversal::Policy<CoprimePolynomials<COMPRESS, TRAV>, TRAV> {

   typedef CoprimePolynomials<COMPRESS, TRAV> self_type;
   typedef Traversal::Policy<self_type, TRAV> TraversalPolicy;
   typedef CompressTraits<COMPRESS> Compress;

public:


   /**
    * value type, polynomials prime with P are taken modulo P
    * ayman : The type name value_type is used by other classes and can not be changed 
    */
   typedef PolyModP value_type;

   /**
    * Size type.
    */
   typedef size_t size_type;

   /**
    * Traversal type.
    */
   typedef TRAV Traversal;

   static std::string name()
   { return std::string("coprime Polynomials / ") + Compress::name() + " / " + Traversal::name(); }

   /**
    * Constructor.
    *
    * \param modulus    Modulus relative to which all numbers in the sequence
    *                   are coprime.
    * \param trav       Traversal instance.
    */
   CoprimePolynomials(Poly polynomial = Poly(1), Traversal trav = Traversal());

   /**
    * Cross-traversal copy-constructor.
    */
   template <class TRAV2>
   CoprimePolynomials(
         const CoprimePolynomials<COMPRESS, TRAV2>& other,
         Traversal trav = Traversal()):
      TraversalPolicy(std::move(trav)),
      m_polynomial(other.m_polynomial),
      m_size(other.m_size),
      m_basis(other.m_basis)
   {}

   /**
    * Rebinds the traversal type.
    */
   template <class TRAV2>
   struct RebindTraversal {
      typedef CoprimePolynomials<COMPRESS, TRAV2> Type;
   };

   /**
    * Returns a copy of this object, but using a different traversal policy.
    */
   template <class TRAV2>
   typename RebindTraversal<TRAV2>::Type rebind(TRAV2 trav) const
   { return typename RebindTraversal<TRAV2>::Type{*this, std::move(trav)}; }

   /**
    * Returns the modulus.
    */
   Poly polynomial() const
   { return m_polynomial; }

   /**
    * Returns the size of the sequence.
    *
    * The size is the value of Euler's totient function.
    */
   size_type size() const
   { return m_size; }

   /**
    * Returns the element at index \c i.
    */
   value_type operator[](size_type i) const;

private:
   template <PolLatBuilder::Compress, class> friend class CoprimePolynomials;

   Poly m_polynomial;
   size_type m_size;
   NTL::vector<detail::CoprimePolynomialsBasisElement> m_basis;
};

}}

//========================================================================
// IMPLEMENTATION
//========================================================================

namespace PolLatBuilder { namespace GenSeq {

template <Compress COMPRESS, class TRAV>
CoprimePolynomials<COMPRESS, TRAV>::CoprimePolynomials(
      Poly polynomial,
      Traversal trav):
   TraversalPolicy(std::move(trav)),
   m_polynomial(polynomial),
   m_size(1)
{
   NTL::vector< NTL::Pair< Poly, long > > factors ;
   CanZass(factors, m_polynomial); // calls "Cantor/Zassenhaus" algorithm from <NTL/GF2XFactoring.h>
   m_basis.resize(factors.size());

   Modulus index = 0;
   for (const auto& b : factors) {
      const auto bk = intPow(b.a, b.b); // b.a is the first element, b.b the second
      const auto m = m_polynomial / bk;
      Modulus totient = intPow(2, b.b * deg(b.a)) / intPow(2,deg(b.a)) * (intPow(2,deg(b.a)) - 1);
      Modulus leap = intPow(2,deg(b.a)) -1;
      detail::CoprimePolynomialsBasisElement e{
         totient,  // totient
         leap,      // leap
         b.a        // irreductible polynomial
      };
      Poly gcd,s,t;
      XGCD(gcd,s,t,bk,m); // bk*s + m*t = gcd = 1
      Poly elem = m * t;
      e.elem = elem % m_polynomial;
      m_size *= e.totient;
      m_basis[index] = std::move(e);
      index ++;
   }

   m_size = Compress::size(m_size + 1) - 1;
}

template <Compress COMPRESS, class TRAV>
auto CoprimePolynomials<COMPRESS, TRAV>::operator[](size_type i) const -> value_type
{
   value_type ret ;
   
   for (const auto& e : m_basis) {
      const ldiv_t qr = ldiv(i, e.totient);
      i = qr.quot;
      Poly Q = intToPoly(qr.rem / e.leap);
      Poly R = intToPoly(qr.rem % e.leap + 1);
      ret += conv <value_type> ((e.irreductible_poly * Q + R) * e.elem);
   }
   
   return Compress::compressIndex(ret , polynomial());
   
}

}}

#endif
// vim: ft=cpp.doxygen
