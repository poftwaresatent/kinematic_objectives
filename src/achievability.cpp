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

#include <kinematic_objectives/achievability.h>
#include <kinematic_objectives/objective.h>
#include <kinematic_objectives/compound_objective.h>
#include <kinematic_objectives/kinematic_model.h>
#include <kinematic_objectives/print.h>
#include <stdio.h>


namespace kinematic_objectives {
  
  
  static void compute(Achievability::objective_tag_t tag,
		      vector<Objective*> const & objectives,
		      Vector const & joint_velocity,
		      vector<Achievability> & information)
  {
    Achievability info;
    Matrix Jinv;
    PseudoInverseFeedback extra;
    Vector null_residual;
    
    for (size_t ii(0); ii < objectives.size(); ++ii) {
      
      if ( ! objectives[ii]->isActive()) {
	continue;
      }
      
      info.tag_                 = tag;
      info.objective_           = objectives[ii];
      info.available_dimension_ = info.objective_->jbar_svd_.truncated_range;
      info.required_dimension_  = info.objective_->jbar_svd_.original_range;
      info.residual_error_
	= info.objective_->getBias()
	- info.objective_->getJacobian() * joint_velocity;
      info.residual_error_magnitude_ = info.objective_->computeResidualErrorMagnitude(info.residual_error_);
      
      pseudo_inverse_moore_penrose(objectives[ii]->getJacobian(), Jinv, 0, &extra);
      info.nullspace_residuals_.resize(extra.input_space.rows() - extra.truncated_range);
      for (ssize_t jj(extra.truncated_range); jj < joint_velocity.size(); ++jj) {
	null_residual
	  = info.objective_->getJacobian()
	  * extra.input_space.block(0, jj, joint_velocity.size(), 1);
        info.nullspace_residuals_[jj - extra.truncated_range]
	  = info.objective_->computeResidualErrorMagnitude(null_residual);
      }
      
      information.push_back(info);
    }
  }
  
  
  void Achievability::
  compute(CompoundObjective const & co,
	  vector<Achievability> & information)
  {
    // Compute some soft-objective data that does not get computed by
    // blender (because it is not required).
    
    Matrix Nct;
    if (co.hard_objectives_.empty()) {
      Nct = Matrix::Identity(co.model_.getJointPosition().size(), co.model_.getJointPosition().size());
    }
    else {
      Nct = co.fb_.hard_objective_nullspace_projector_;
    }
    Matrix Jbinv;
    for (size_t ii(0); ii < co.soft_objectives_.size(); ++ii) {
      if ( ! co.soft_objectives_[ii]->isActive()) {
	continue;
      }
      pseudo_inverse_moore_penrose(co.soft_objectives_[ii]->getJacobian() * Nct, Jbinv, 0,
				   &co.soft_objectives_[ii]->jbar_svd_);
    }
    
    // OK now we can compute all the achievability information...
    
    information.clear();
    
    compute(UNILATERAL_CONSTRAINT, co.unilateral_constraints_, co.model_.getJointVelocity(), information);
    compute(HARD_OBJECTIVE,        co.hard_objectives_,        co.model_.getJointVelocity(), information);
    compute(SOFT_OBJECTIVE,        co.soft_objectives_,        co.model_.getJointVelocity(), information);
  }
  
  
  void Achievability::
  print(vector<Achievability> const & information,
	ostream & os, string const & pfx)
  {
    if (information.empty()) {
      return;
    }
    
    int namlen(information[0].objective_->name_.size());
    for (size_t ii(1); ii < information.size(); ++ii) {
      int nl(information[ii].objective_->name_.size());
      if (nl > namlen) {
	namlen = nl;
      }
    }
    
    size_t const buflen(128+namlen);
    char buf[buflen];
    for (size_t ii(0); ii < information.size(); ++ii) {
      snprintf(buf, buflen, "%s %-*s   %3zu / %3zu   ",
	       information[ii].tag_ == UNILATERAL_CONSTRAINT ? "CONS"
	       : (information[ii].tag_ == HARD_OBJECTIVE ? "HARD" : "SOFT"),
	       namlen, information[ii].objective_->name_.c_str(),
	       information[ii].available_dimension_,
	       information[ii].required_dimension_);
      os << pfx << buf << pstring(information[ii].residual_error_magnitude_)
	 << "   (" << pstring(information[ii].residual_error_)
	 << ")   (" << pstring(information[ii].nullspace_residuals_) << ")\n";
    }
  }
  
}
