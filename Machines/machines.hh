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

#ifndef NETKET_MACHINES_HH
#define NETKET_MACHINES_HH

/** @defgroup machines Machines module
 *
 * Machines module contains implementations of wave-functions.
 */

namespace netket{
  template<class T> class AbstractMachine;
  template<class T> class RbmSpin;
  template<class T> class RbmSpinSymm;
  template<class T> class RbmMultival;
}

#include "netket.hh"
#include "abstract_machine.hh"
#include "rbm_spin.hh"
#include "rbm_spin_symm.hh"
#include "rbm_multival.hh"


namespace netket{

  class Machine:public AbstractMachine<std::complex<double>>{
    using Ptype=std::unique_ptr<AbstractMachine<std::complex<double>>>;

    Ptype m_;

    const Hilbert & hilbert_;
    const Hamiltonian & hamiltonian_;

    int mynode_;

  public:

    using VectorType=typename AbstractMachine<std::complex<double>>::VectorType;
    using MatrixType=typename AbstractMachine<std::complex<double>>::MatrixType;
    using StateType=typename AbstractMachine<std::complex<double>>::StateType;
    using LookupType=typename AbstractMachine<std::complex<double>>::LookupType;


    Machine(const Graph & graph,const Hamiltonian & hamiltonian,const json & pars);

    int Npar()const;
    int Nvisible()const;
    void InitLookup(const VectorXd & v,LookupType & lt);

    void UpdateLookup(const VectorXd & v,const vector<int>  & toflip,
      const vector<double> & newconf,LookupType & lt);

    void UpdateConf(VectorXd & v,const vector<int>  & toflip,const vector<double> & newconf);

    VectorType DerLog(const VectorXd & v);

    MatrixType DerLogDiff(const VectorXd & v,
      const vector<vector<int> >  & toflip,
      const vector<vector<double>> & newconf);

    VectorType GetParameters();

    void SetParameters(const VectorType & pars);

    StateType LogVal(const VectorXd & v);

    StateType LogVal(const VectorXd & v,LookupType & lt);

    VectorType LogValDiff(const VectorXd & v,
      const vector<vector<int> >  & toflip,
      const vector<vector<double>> & newconf);

    StateType LogValDiff(const VectorXd & v,const vector<int>  & toflip,
        const vector<double> & newconf,const LookupType & lt);

    void InitRandomPars(int seed,double sigma);

    const Hilbert& GetHilbert()const;

    const Hamiltonian& GetHamiltonian()const;

    void to_json(json &j)const;

    void from_json(const json&j);
  };

}

#endif
