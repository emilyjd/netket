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

#ifndef NETKET_ADADELTA_HH
#define NETKET_ADADELTA_HH

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

namespace netket{

using namespace std;
using namespace Eigen;

class AdaDelta: public AbstractStepper{

  int npar_;

  //decay constant
  double rho_;

  //small parameter
  double eps_;

  VectorXd Eg2_;
  VectorXd Edx2_;

  int mynode_;

  const std::complex<double> I_;

public:

  AdaDelta(double rho=0.95,double eps=1.0e-6):
  rho_(rho),eps_(eps),I_(0,1){
    npar_=-1;

    PrintParameters();
  }

  //Json constructor
  AdaDelta(const json & pars):
    rho_(FieldOrDefaultVal(pars["Learning"],"rho",0.95)),
    eps_(FieldOrDefaultVal(pars["Learning"],"eps",1.0e-6)),
    I_(0,1)
  {
    npar_=-1;

    PrintParameters();
  }

  void PrintParameters(){
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode_);
    if(mynode_==0){
      cout<<"# AdaDelta stepper initialized with these parameters : "<<endl;
      cout<<"# Rho = "<<rho_<<endl;
      cout<<"# Eps = "<<eps_<<endl;
    }
  }

  void Init(const VectorXd & pars){

    npar_=pars.size();
    Eg2_.setZero(npar_);
    Edx2_.setZero(npar_);

  }

  void Init(const VectorXcd & pars){

    npar_=2*pars.size();
    Eg2_.setZero(npar_);
    Edx2_.setZero(npar_);

  }

  void Update(const VectorXd & grad,VectorXd & pars){
    assert(npar_>0);

    Eg2_=rho_*Eg2_+(1.-rho_)*grad.cwiseAbs2();

    VectorXd Dx(npar_);

    for(int i=0;i<npar_;i++){
      Dx(i)=-std::sqrt(Edx2_(i)+eps_)*grad(i);
      Dx(i)/=std::sqrt(Eg2_(i)+eps_);
      pars(i)+=Dx(i);
    }

    Edx2_=rho_*Edx2_+(1.-rho_)*Dx.cwiseAbs2();

  }

  void Update(const VectorXcd & grad,VectorXd & pars){
    Update(VectorXd(grad.real()),pars);
  }

  void Update(const VectorXcd & grad,VectorXcd & pars){

    assert(npar_==2*pars.size());

    VectorXd Dx(npar_);

    for(int i=0;i<pars.size();i++){
      Eg2_(2*i)=rho_*Eg2_(2*i)+(1.-rho_)*std::pow(grad(i).real(),2);
      Eg2_(2*i+1)=rho_*Eg2_(2*i+1)+(1.-rho_)*std::pow(grad(i).imag(),2);

      Dx(2*i)=-std::sqrt(Edx2_(2*i)+eps_)*grad(i).real();
      Dx(2*i+1)=-std::sqrt(Edx2_(2*i+1)+eps_)*grad(i).imag();
      Dx(2*i)/=std::sqrt(Eg2_(2*i)+eps_);
      Dx(2*i+1)/=std::sqrt(Eg2_(2*i+1)+eps_);

      pars(i)+=Dx(2*i);
      pars(i)+=I_*Dx(2*i+1);

      Edx2_(2*i)=rho_*Edx2_(2*i)+(1.-rho_)*std::pow(Dx(2*i),2);
      Edx2_(2*i+1)=rho_*Edx2_(2*i+1)+(1.-rho_)*std::pow(Dx(2*i+1),2);
    }
  }

  void Reset(){
    Eg2_=VectorXd::Zero(npar_);
    Edx2_=VectorXd::Zero(npar_);
  }
};

}

#endif
