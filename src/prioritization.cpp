/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2013, Willow Garage, Inc.
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of the Willow Garage nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

/* Author: Roland Philippsen */

#include <kinematic_objectives/prioritization.h>
#include <kinematic_objectives/compound_objective.h>
#include <kinematic_objectives/kinematic_model.h>
#include <kinematic_objectives/pseudo_inverse.h>
#include <kinematic_objectives/print.h>
#include <kinematic_objectives/objective.h>


namespace kinematic_objectives {
  
  
  Prioritization::
  Prioritization()
    : dbgos_(0),
      dbgpre_("  ")
  {
  }
  
  
  void Prioritization::
  process(Matrix const & N_in,
	  Vector const & bias_in,
	  Objective const * objective,
	  Matrix * opt_N_up,
	  Vector & bias_out)
  {
    objective->bias_comp_ = objective->getBias() - objective->getJacobian() * bias_in;
    objective->jbar_ = objective->getJacobian() * N_in;
    
    pseudo_inverse_moore_penrose(objective->jbar_,
				 objective->jbar_inv_,
				 opt_N_up,
				 &objective->jbar_svd_);
    
    bias_out = objective->jbar_inv_ * objective->bias_comp_;
    
    if (dbgos_) {
      *dbgos_ << dbgpre_ << "objective " << objective->name_ << ":\n";
      string pre(dbgpre_);
      pre += "  ";
      print(objective->getBias(), *dbgos_, "bias", pre);
      print(objective->getJacobian(), *dbgos_, "Jacobian", pre);
      Eigen::JacobiSVD<Matrix> svd;
      svd.compute(objective->getJacobian(), Eigen::ComputeThinU | Eigen::ComputeThinV);
      print(svd.singularValues(), *dbgos_, "singular value of Jacobian", pre);
      print(objective->bias_comp_, *dbgos_, "bias_comp", pre);
      print(objective->jbar_, *dbgos_, "J_bar", pre);
      print(objective->jbar_svd_.original_sigma, *dbgos_, "original J_bar sigma", pre);
      print(objective->jbar_svd_.regularized_sigma, *dbgos_, "regularized J_bar sigma", pre);
      print(objective->jbar_inv_, *dbgos_, "J_bar pseudo inverse", pre);
      if (opt_N_up) {
	print(*opt_N_up, *dbgos_, "nullspace projector update", pre);
      }
      print(bias_out, *dbgos_, "bias update", pre);
    }
  }
  
  
  void Prioritization::
  projectObjectives(Matrix const & N_in,
		    Vector const & bias_in,
		    vector<Objective*> const & objectives,
		    Vector & bias_out,
		    Matrix & N_out)
  {
    if (&bias_in != &bias_out) {
      bias_out = bias_in;
    }
    if (&N_out != &N_in) {
      N_out = N_in;
    }
    
    Vector bias_up;
    Matrix N_up;
    
    for (size_t ii(0); ii < objectives.size(); ++ii) {
      Objective const * obj(objectives[ii]);
      
      if ( ! obj->isActive()) {
	if (dbgos_) {
	  *dbgos_ << dbgpre_ << "objective " << ii << " \"" << obj->name_ << "\": INACTIVE\n";
	}
	obj->clearFeedback();
	continue;
      }
      
      process(N_out, bias_out, obj, &N_up, bias_up);
      
      bias_out += bias_up;
      N_out -= N_up;
      
      if (dbgos_) {
        string pre (dbgpre_);
        pre += "  ";
	print(bias_out, *dbgos_, "accumulated bias", pre);
	print(N_out, *dbgos_, "accumulated nullspace projector", pre);
	Eigen::JacobiSVD<Matrix> svd;
	svd.compute(N_out, Eigen::ComputeThinU | Eigen::ComputeThinV);
	print(svd.singularValues(), *dbgos_, "singular values of nullspace projector", pre);
      }
    }
  }
  
  
  void Prioritization::
  addUpObjectives(Matrix const & N_in,
		  ////		  Vector const & bias_in,
		  vector<Objective*> const & objectives,
		  Vector & bias_out)
  {
    Vector accu = Vector::Zero(N_in.rows());
    for (size_t ii(0); ii < objectives.size(); ++ii) {
      Objective const * obj(objectives[ii]);
      obj->clearFeedback();
      if ( ! obj->isActive()) {
        continue;
      }
      Matrix Jinv;
      pseudo_inverse_moore_penrose(obj->getJacobian(), Jinv);
      accu += Jinv * obj->getBias();
    }
    bias_out = N_in * accu;
    //
    //
    // Vector accu = Vector::Zero(N_in.rows());
    // for (size_t ii(0); ii < objectives.size(); ++ii) {
    //   Objective const * obj(objectives[ii]);
    //   obj->clearFeedback();
    //   if ( ! obj->isActive()) {
    //     continue;
    //   }
    //   obj->jbar_ = obj->getJacobian() * N_in;
    //   pseudo_inverse_moore_penrose(obj->jbar_,
    // 				   obj->jbar_inv_,
    // 				   0,
    // 				   &obj->jbar_svd_);
    //   accu += obj->jbar_inv_ * obj->getBias();
    // }
    // bias_out = N_in * accu;
    //
    //
    // Vector accu(bias_in);	// needed in case &bias_in == &bias_out
    // Vector bias_up;
    // for (size_t ii(0); ii < objectives.size(); ++ii) {
    //   Objective const * obj(objectives[ii]);
    //   if ( ! obj->isActive()) {
    // 	obj->clearFeedback();
    //     continue;
    //   }
    //   process(N_in, bias_in, obj, 0, bias_up);
    //   accu += bias_up;
    // }
    // bias_out = accu;
    //
    //
    // Vector zero = Vector::Zero(bias_in.size());
    // Vector accu = zero;
    // Vector bias_up;
    // for (size_t ii(0); ii < objectives.size(); ++ii) {
    //   Objective const * obj(objectives[ii]);
    //   if ( ! obj->isActive()) {
    // 	obj->clearFeedback();
    //     continue;
    //   }
    //   process(N_in, zero, obj, 0, bias_up);
    //   accu += bias_up;
    // }
    // bias_out = accu - bias_in;
  }
  
}
