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

#ifndef NETKET_HILBERT_HH
#define NETKET_HILBERT_HH

namespace netket{
  class AbstractHilbert;
  class Spin;
  class Boson;
  class Qubit;
  class CustomHilbert;
  class LocalOperator;
}

#include "abstract_hilbert.hh"
#include "next_variation.hh"
#include "spins.hh"
#include "bosons.hh"
#include "qubits.hh"
#include "custom_hilbert.hh"

namespace netket{
  class Hilbert:public AbstractHilbert{

    using Ptype=std::unique_ptr<AbstractHilbert>;
    Ptype h_;

  public:

    Hilbert();
    Hilbert(const json & pars);
    void Init(const json & pars);
    bool IsDiscrete()const;
    int LocalSize()const;
    int Size()const;
    vector<double> LocalStates()const;
    void RandomVals(VectorXd & state,netket::default_random_engine & rgen)const;
    void UpdateConf(VectorXd & v,const vector<int>  & tochange,
      const vector<double> & newconf)const;
  };
}

#include "local_operator.hh"

#endif
