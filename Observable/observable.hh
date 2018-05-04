// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef NETKET_OBSERVABLE_HH
#define NETKET_OBSERVABLE_HH

namespace netket{
  class AbstractObservable;
  class CustomObservable;
  class Observables;
}

#include "netket.hh"
#include "abstract_observable.hh"
#include "custom_observable.hh"

namespace netket{

  class Observable:public AbstractObservable{

    using Ptype=std::unique_ptr<AbstractObservable>;
    Ptype o_;

  public:

    using MatType=LocalOperator::MatType;

    Observable(const Hilbert & hilbert,const json & obspars);

    void FindConn(const VectorXd & v,
      vector<std::complex<double>> & mel,
      vector<vector<int>> & connectors,
      vector<vector<double>> & newconfs);

    const Hilbert & GetHilbert()const;

    const std::string Name()const;
  };

}

#include "observables.hh"


#endif
