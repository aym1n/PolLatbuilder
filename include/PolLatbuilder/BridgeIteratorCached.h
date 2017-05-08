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

#ifndef POLLATBUILDER__BRIDGE_ITERATOR_CACHED_H
#define POLLATBUILDER__BRIDGE_ITERATOR_CACHED_H

#include <boost/iterator/iterator_adaptor.hpp>
#include <limits>
#include <stdexcept>

namespace PolLatBuilder {

/**
 * Bridge iterator with cached value.
 */
template <typename SEQ>
class BridgeIteratorCached :
   public boost::iterators::iterator_adaptor<
      BridgeIteratorCached<SEQ>,
      typename SEQ::Base::const_iterator,
      const typename SEQ::value_type
   >
{
   typedef BridgeIteratorCached<SEQ> self_type;

public:
   typedef typename SEQ::value_type value_type;
   typedef typename SEQ::size_type size_type;

   struct end_tag {};

   BridgeIteratorCached():
      self_type::iterator_adaptor_(),
      m_seq(nullptr)
   {}

   explicit BridgeIteratorCached(const SEQ& seq):
      self_type::iterator_adaptor_(seq.base().begin()),
      m_seq(&seq), m_cached(false)
   {}

   BridgeIteratorCached(const SEQ& seq, end_tag):
      self_type::iterator_adaptor_(seq.base().end()),
      m_seq(&seq), m_cached(false)
   {}

   const SEQ& seq() const
   { return *m_seq; }

   size_type index() const
   { return this->base_reference().index(); }

private:
   friend class boost::iterators::iterator_core_access;

   void increment()
   { ++this->base_reference(); m_cached = false; }

   bool equal(const BridgeIteratorCached& other) const
   { return m_seq == other.m_seq and this->base_reference() == other.base_reference(); }

   const value_type& dereference() const
   {
#ifndef NDEBUG
      if (this->base_reference() == m_seq->base().end())
         throw std::runtime_error("BridgeIteratorCached: dereferencing past end of sequence");
#endif
      if (!m_cached) {
         m_value = m_seq->element(this->base_reference());   
         m_cached = true;
      }
      return m_value;
   }

   ptrdiff_t distance_to(const BridgeIteratorCached& other) const
   { return m_seq == other.m_seq ? other.base_reference() - this->base_reference() : std::numeric_limits<ptrdiff_t>::max(); }

private:
   const SEQ* m_seq;
   mutable bool m_cached;
   mutable value_type m_value;
};

}

#endif
