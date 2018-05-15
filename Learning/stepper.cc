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

#ifndef NETKET_STEPPER_CC
#define NETKET_STEPPER_CC

namespace netket{

class Stepper:public AbstractStepper{

  AbstractStepper * s_;

public:

  Stepper(const json & pars){

    if(!FieldExists(pars,"Learning")){
      cerr<<"Learning is not defined in the input"<<endl;
      std::abort();
    }

    if(!FieldExists(pars["Learning"],"StepperType")){
      cerr<<"Stepper Type is not defined in the input"<<endl;
      std::abort();
    }

    if(pars["Learning"]["StepperType"]=="Sgd"){
      s_=new Sgd(pars);
    }
    else if(pars["Learning"]["StepperType"]=="AdaMax"){
      s_=new AdaMax(pars);
    }
    else if(pars["Learning"]["StepperType"]=="Momentum"){
      s_=new Momentum(pars);
    }
    else if(pars["Learning"]["StepperType"]=="AdaGrad"){
          s_=new AMSGrad(pars);
    }
    else if(pars["Learning"]["StepperType"]=="AdaDelta"){
          s_=new AdaDelta(pars);
    }
    else if(pars["Learning"]["StepperType"]=="RMSProp"){
          s_=new RMSProp(pars);
    }
    else if(pars["Learning"]["StepperType"]=="AMSGrad"){
          s_=new AMSGrad(pars);
    }

    else{
      cout<<"StepperType not found"<<endl;
      std::abort();
    }
  }

  void Init(const VectorXd & pars){
    return s_->Init(pars);
  }

  void Init(const VectorXcd & pars){
    return s_->Init(pars);
  }

  void Update(const VectorXd & grad,VectorXd & pars){
    return s_->Update(grad,pars);
  }

  void Update(const VectorXcd & grad,VectorXd & pars){
    return s_->Update(grad,pars);
  }

  void Update(const VectorXcd & grad,VectorXcd & pars){
    return s_->Update(grad,pars);
  }

  void Reset(){
    return s_->Reset();
  }

};
}
#endif
