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

#include "PolLatbuilder/SizeParam.h"
#include "PolLatbuilder/Util.h"

namespace PolLatBuilder {

SizeParam<LatType::ORDINARY>::SizeParam(Poly polynomial):
   BasicSizeParam<SizeParam<LatType::ORDINARY>>(polynomial)
{}

/**
    * ayman : not used (at least for the moment) 
    */
size_t
SizeParam<LatType::ORDINARY>::totient() const
{
   return 0;
}

/**
    * TO DO: numPoints below is not efficient, store numPoints sounds better
    */
Modulus
SizeParam<LatType::ORDINARY>::numPoints() const
{
   return intPow(2,deg(polynomial()));
}

void
SizeParam<LatType::ORDINARY>::normalize(Real& merit) const
{ merit /= Real(numPoints()); }



std::ostream&
SizeParam<LatType::ORDINARY>::format(std::ostream& os) const
{ return os << polynomial(); }

}

