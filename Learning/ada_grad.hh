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

#ifndef NETKET_ADAGRAD_HH
#define NETKET_ADAGRAD_HH

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <complex>

namespace netket{

using namespace std;
using namespace Eigen;

class AdaGrad: public AbstractStepper{

  int npar_;

  double eta_;

  VectorXd Gt_;

  double epscut_;

  int mynode_;

  const std::complex<double> I_;

public:

  AdaGrad(double eta=0.001,double epscut=1.0e-7):
    eta_(eta),epscut_(epscut),I_(0,1)
  {
    npar_=-1;
    PrintParameters();
  }

  //Json constructor
  AdaGrad(const json & pars):
    eta_(FieldOrDefaultVal(pars["Learning"],"Learning Rate",0.001)),
    epscut_(FieldOrDefaultVal(pars["Learning"],"Epscut",1.0e-7)),
    I_(0,1)
  {
    npar_=-1;

    PrintParameters();
  }

  void PrintParameters(){
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode_);
    if(mynode_==0){
      cout<<"# Adagrad stepper initialized with these parameters : "<<endl;
      cout<<"# Learning Rate = "<<eta_<<endl;
      cout<<"# Epscut = "<<epscut_<<endl;
    }
  }

  void Init(const VectorXd & pars){

    npar_=pars.size();
    Gt_.setZero(npar_);
  }

  void Init(const VectorXcd & pars){

    npar_=2*pars.size();
    Gt_.setZero(npar_);
  }

  void Update(const VectorXd & grad,VectorXd & pars){

    assert(npar_>0);

    Gt_+=grad.cwiseAbs2();

    for(int i=0;i<npar_;i++){
      pars(i)-=eta_*grad(i)/std::sqrt(Gt(i)+epscut);
    }
  }

  void Update(const VectorXcd & grad,VectorXd & pars){
    Update(VectorXd(grad.real()),pars);
  }

  void Update(const VectorXcd & grad,VectorXcd & pars){

    assert(npar_==2*pars.size());

    for(int i=0;i<pars.size();i++){
      Gt_(2*i)+=std::pow(grad(i).real(),2);
      Gt_(2*i+1)+=std::pow(grad(i).imag(),2);
    }

    for(int i=0;i<pars.size();i++){
      pars(i)-=eta_*grad(i)/std::sqrt(Gt(2*i)+epscut);
      pars(i)-=eta_*I_*grad(i)/std::sqrt(Gt(2*i+1)+epscut);
    }
  }

  void Reset(){
    Gt_=VectorXd::Zero(npar_);
  }
};


}

#endif
