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

#ifndef NETKET_SAMPLER_CC
#define NETKET_SAMPLER_CC

#include "netket.hh"

namespace netket{

  Sampler::Sampler(Graph & graph,Hamiltonian & hamiltonian,Sampler::MachineType & psi,const json & pars){

    if(!FieldExists(pars,"Sampler")){
      cerr<<"Sampler is not defined in the input"<<endl;
      std::abort();
    }

    if(!FieldExists(pars["Sampler"],"Name")){
      cerr<<"Sampler Name is not defined in the input"<<endl;
      std::abort();
    }

    if(pars["Sampler"]["Name"]=="MetropolisLocal"){
      s_=Ptype(new MetropolisLocal<Sampler::MachineType>(graph,psi,pars));
    }
    else if(pars["Sampler"]["Name"]=="MetropolisLocalPt"){
      s_=Ptype(new MetropolisLocalPt<Sampler::MachineType>(graph,psi,pars));
    }
    else if(pars["Sampler"]["Name"]=="MetropolisExchange"){
      s_=Ptype(new MetropolisExchange<Sampler::MachineType>(graph,psi,pars));
    }
    else if(pars["Sampler"]["Name"]=="MetropolisExchangePt"){
      s_=Ptype(new MetropolisExchangePt<Sampler::MachineType>(graph,psi,pars));
    }
    else if(pars["Sampler"]["Name"]=="MetropolisHop"){
      s_=Ptype(new MetropolisHop<Sampler::MachineType>(graph,psi,pars));
    }
    else if(pars["Sampler"]["Name"]=="MetropolisHamiltonian"){
      s_=Ptype(new MetropolisHamiltonian<Sampler::MachineType,Hamiltonian>(graph,psi,hamiltonian,pars));
    }
    else if(pars["Sampler"]["Name"]=="MetropolisHamiltonianPt"){
      s_=Ptype(new MetropolisHamiltonianPt<Sampler::MachineType,Hamiltonian>(graph,psi,hamiltonian,pars));
    }
    else{
      cout<<"Sampler not found"<<endl;
      std::abort();
    }
  }

  void Sampler::Reset(bool initrandom){
    return s_->Reset(initrandom);
  }
  void Sampler::Sweep(){
    return s_->Sweep();
  }
  VectorXd Sampler::Visible(){
    return s_->Visible();
  }
  void Sampler::SetVisible(const VectorXd & v){
    return s_->SetVisible(v);
  }
  Sampler::MachineType & Sampler::Psi(){
    return s_->Psi();
  }
  VectorXd Sampler::Acceptance()const{
    return s_->Acceptance();
  }
}

#endif
