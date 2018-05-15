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

#ifndef NETKET_AMSGRAD_HH
#define NETKET_AMSGRAD_HH

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <complex>

namespace netket{

using namespace std;
using namespace Eigen;

class AMSGrad: public AbstractStepper{

  int npar_;

  double eta_;
  double beta1_;
  double beta2_;

  VectorXd mt_;
  VectorXd vt_;

  double epscut_;

  int mynode_;

  const std::complex<double> I_;

public:

  AMSGrad(double eta=0.001,double beta1=0.9,double beta2=0.999,double epscut=1.0e-7):
    eta_(eta),beta1_(beta1),beta2_(beta2),epscut_(epscut),I_(0,1)
  {
    npar_=-1;

    PrintParameters();
  }

  //Json constructor
  AMSGrad(const json & pars):
    eta_(FieldOrDefaultVal(pars["Learning"],"Learning Rate",0.001)),
    beta1_(FieldOrDefaultVal(pars["Learning"],"Beta1",0.9)),
    beta2_(FieldOrDefaultVal(pars["Learning"],"Beta2",0.999)),
    epscut_(FieldOrDefaultVal(pars["Learning"],"Epscut",1.0e-7)),
    I_(0,1)
  {
    npar_=-1;

    PrintParameters();
  }

  void PrintParameters(){
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode_);
    if(mynode_==0){
      cout<<"# AMSGrad stepper initialized with these parameters : "<<endl;
      cout<<"# Eta = "<<eta_<<endl;
      cout<<"# Beta1 = "<<beta1_<<endl;
      cout<<"# Beta2 = "<<beta2_<<endl;
      cout<<"# Epscut = "<<epscut_<<endl;
    }
  }

  void Init(const VectorXd & pars){

    npar_=pars.size();
    vt_.setZero(npar_);
    mt_.setZero(npar_);
  }

  void Init(const VectorXcd & pars){

    npar_=2*pars.size();
    vt_.setZero(npar_);
    mt_.setZero(npar_);
  }

  void Update(const VectorXd & grad,VectorXd & pars){

    assert(npar_>0);

    mt_=beta1_*mt_+(1.-beta1_)*grad;

    for(int i=0;i<npar_;i++){
      vt_(i)=std::max(vt_(i),beta2_*vt_(i)+(1-beta2_)*std::pow(grad(i),2));
    }

    for(int i=0;i<npar_;i++){
      pars(i)-=eta_*mt_(i)/(std::sqrt(vt_(i))+epscut_);
    }
  }

  void Update(const VectorXcd & grad,VectorXd & pars){
    Update(VectorXd(grad.real()),pars);
  }

  void Update(const VectorXcd & grad,VectorXcd & pars){

    assert(npar_==2*pars.size());

    for(int i=0;i<pars.size();i++){
      mt_(2*i)=beta1_*mt_(2*i)+(1.-beta1_)*grad(i).real();
      mt_(2*i+1)=beta1_*mt_(2*i+1)+(1.-beta1_)*grad(i).imag();
    }

    for(int i=0;i<pars.size();i++){
      vt_(2*i)=std::max(vt_(2*i),beta2_*vt_(2*i)+(1-beta2_)*std::pow(grad(i).real(),2));
      vt_(2*i+1)=std::max(vt_(2*i+1),beta2_*vt_(2*i+1)+(1-beta2_)*std::pow(grad(i).imag(),2));
    }

    for(int i=0;i<pars.size();i++){
      pars(i)-=eta_*mt_(2*i)/(std::sqrt(vt_(2*i))+epscut_);
      pars(i)-=eta_*I_*mt_(2*i+1)/(std::sqrt(vt_(2*i+1))+epscut_);
    }
  }

  void Reset(){
    vt_=VectorXd::Zero(npar_);
    mt_=VectorXd::Zero(npar_);
  }

};


}

#endif
