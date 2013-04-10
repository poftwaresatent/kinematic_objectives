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
  
  
  void Achievability::
  compute(CompoundObjective const & co,
	  vector<Achievability> & information)
  {
    information.clear();
    Achievability info;
    
    for (size_t ii(0); ii < co.unilateral_constraints_.size(); ++ii) {
      if ( ! co.unilateral_constraints_[ii]->isActive()) {
	continue;
      }
      info.tag_                 = UNILATERAL_CONSTRAINT;
      info.objective_           = co.unilateral_constraints_[ii];
      info.available_dimension_ = info.objective_->jbar_svd_.truncated_range;
      info.required_dimension_  = info.objective_->jbar_svd_.original_range;
      info.residual_error_      = info.objective_->getBias()
	- info.objective_->getJacobian() * co.model_.getJointVelocity();
      ////      info.score_;
      information.push_back(info);
    }
    
    for (size_t ii(0); ii < co.hard_objectives_.size(); ++ii) {
      if ( ! co.hard_objectives_[ii]->isActive()) {
	continue;
      }
      info.tag_                 = HARD_OBJECTIVE;
      info.objective_           = co.hard_objectives_[ii];
      info.available_dimension_ = info.objective_->jbar_svd_.truncated_range;
      info.required_dimension_  = info.objective_->jbar_svd_.original_range;
      info.residual_error_      = info.objective_->getBias()
	- info.objective_->getJacobian() * co.model_.getJointVelocity();
      ////      info.score_;
      information.push_back(info);
    }
    
    /**< \note Soft objectives are a bit special because they do not
       individually get projected. So here we compute a few spurious
       SVDs just in order to find out how many dimensions each soft
       objective wants and gets. Assuming that this info is not
       computed too often, that should not result in an appreciable
       slowdown.
    */
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
      info.tag_                 = SOFT_OBJECTIVE;
      info.objective_           = co.soft_objectives_[ii];
      pseudo_inverse_moore_penrose(co.soft_objectives_[ii]->getJacobian() * Nct, Jbinv, 0,
				   &co.soft_objectives_[ii]->jbar_svd_);
      info.available_dimension_ = info.objective_->jbar_svd_.truncated_range;
      info.required_dimension_  = info.objective_->jbar_svd_.original_range;
      info.residual_error_      = info.objective_->getBias()
	- info.objective_->getJacobian() * co.model_.getJointVelocity();
      ////      info.score_;
      information.push_back(info);
    }
    
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
      os << pfx << buf << pstring(information[ii].residual_error_) << "\n";
    }
  }
  
}
