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
  
  
  static size_t compute_dimension(Vector const & sigma)
  {
    ssize_t dim;
    for (dim = 0; dim < sigma.size(); ++dim) {
      if (sigma[dim] == 0.0) {
	break;
      }
    }
    return dim;
  }
  
  
  static void helper(Achievability::objective_tag_t tag,
		     vector<Objective*> const & objectives,
		     Vector const & joint_velocity,
		     vector<Achievability> & information)
  {
    Achievability info;
    
    for (size_t ii(0); ii < objectives.size(); ++ii) {
      
      if ( ! objectives[ii]->isActive()) {
	continue;
      }
      
      info.tag_                    = tag;
      info.objective_              = objectives[ii];
      
      Matrix Jinv;
      pseudo_inverse_moore_penrose(objectives[ii]->getJacobian(), Jinv, 0, &info.j_svd_);
      
      info.jbar_svd_                   = info.objective_->jbar_svd_;
      info.original_j_dimension_       = compute_dimension(info.j_svd_.original_sigma);
      info.original_jbar_dimension_    = compute_dimension(info.jbar_svd_.original_sigma);
      info.regularized_j_dimension_    = compute_dimension(info.j_svd_.regularized_sigma);
      info.regularized_jbar_dimension_ = compute_dimension(info.jbar_svd_.regularized_sigma);
      
      info.residual_error_
	= info.objective_->getBias()
	- info.objective_->getJacobian() * joint_velocity;
      info.residual_error_magnitude_
	= info.objective_->computeResidualErrorMagnitude(info.residual_error_);
      
      if (0 == info.j_svd_.input_space.rows()) {
	info.jbar_nullspace_residuals_.resize(0);
      }
      else {
	info.j_nullspace_residuals_.resize(info.j_svd_.input_space.rows() - info.regularized_j_dimension_);
	for (ssize_t jj(info.regularized_j_dimension_); jj < joint_velocity.size(); ++jj) {
	  Vector null_residual
	    = info.objective_->getJacobian()
	    * info.j_svd_.input_space.block(0, jj, joint_velocity.size(), 1);
	  info.j_nullspace_residuals_[jj - info.regularized_j_dimension_]
	    = info.objective_->computeResidualErrorMagnitude(null_residual);
	}
      }
      
      if (0 == info.jbar_svd_.input_space.rows()) {
	info.jbar_nullspace_residuals_.resize(0);
	info.cross_nullspace_residuals_.resize(0);
      }
      else {
	info.jbar_nullspace_residuals_.resize(info.jbar_svd_.input_space.rows() - info.regularized_jbar_dimension_);
	for (ssize_t jj(info.regularized_jbar_dimension_); jj < joint_velocity.size(); ++jj) {
	  Vector null_residual
	    = info.objective_->jbar_
	    * info.jbar_svd_.input_space.block(0, jj, joint_velocity.size(), 1);
	  info.jbar_nullspace_residuals_[jj - info.regularized_jbar_dimension_]
	    = info.objective_->computeResidualErrorMagnitude(null_residual);
	}
	info.cross_nullspace_residuals_.resize(info.jbar_svd_.input_space.rows() - info.regularized_jbar_dimension_);
	for (ssize_t jj(info.regularized_jbar_dimension_); jj < joint_velocity.size(); ++jj) {
	  Vector null_residual
	    = info.objective_->getJacobian()
	    * info.jbar_svd_.input_space.block(0, jj, joint_velocity.size(), 1);
	  info.cross_nullspace_residuals_[jj - info.regularized_jbar_dimension_]
	    = info.objective_->computeResidualErrorMagnitude(null_residual);
	}
      }
      
      information.push_back(info);
    }
  }
  
  
  void Achievability::
  compute(KinematicModel & model,
	  CompoundObjective const & co,
	  vector<Achievability> & information)
  {
    information.clear();
    
    helper(UNILATERAL_CONSTRAINT, co.unilateral_constraints_, model.getJointVelocity(), information);
    helper(HARD_OBJECTIVE,        co.hard_objectives_,        model.getJointVelocity(), information);
    helper(SOFT_OBJECTIVE,        co.soft_objectives_,        model.getJointVelocity(), information);
  }
  
  
  void Achievability::
  print(vector<Achievability> const & information,
	ostream & os, string const & pfx)
  {
    if (information.empty()) {
      return;
    }
    
    int namlen(strlen("objective"));
    for (size_t ii(0); ii < information.size(); ++ii) {
      int nl(information[ii].objective_->name_.size());
      if (nl > namlen) {
	namlen = nl;
      }
    }
    
    size_t const buflen(128+namlen);
    char buf[buflen];
    
    snprintf(buf, buflen, "  [pos/idx] %-*s   J have/need   Jbar have/need magnitude\n", namlen, "objective");
    os << pfx << "==========================================================\n"
       << pfx << "SUMMARY: dimensions and magnitude of residual error\n"
       << pfx << buf
       << pfx << "----------------------------------------------------------\n";
    
    for (size_t ii(0), jj(0); ii < information.size(); ++ii, ++jj) {
      if ((0 == ii) || information[ii].tag_ != information[ii-1].tag_) {
	jj = 0;
	if (UNILATERAL_CONSTRAINT == information[ii].tag_) {
	  os << pfx << "UNILATERAL_CONSTRAINT\n";
	}
	else if (HARD_OBJECTIVE == information[ii].tag_) {
	  os << pfx << "HARD_OBJECTIVE\n";
	}
	else {
	  os << pfx << "SOFT_OBJECTIVE\n";
	}
      }
      snprintf(buf, buflen, "  [%3zu/%3zu] %-*s   %3zu / %3zu  %3zu / %3zu   ", ii, jj, namlen,
	       information[ii].objective_->name_.c_str(),
	       information[ii].regularized_j_dimension_,
	       information[ii].original_j_dimension_,
	       information[ii].regularized_jbar_dimension_,
	       information[ii].original_jbar_dimension_);
      os << pfx << buf << pstring(information[ii].residual_error_magnitude_) << "\n";
    }
    
    os << pfx << "----------------------------------------------------------\n"
       << pfx << "joint-space residual error\n"
       << pfx << "----------------------------------------------------------\n";
    
    for (size_t ii(0), jj(0); ii < information.size(); ++ii, ++jj) {
      if ((0 == ii) || information[ii].tag_ != information[ii-1].tag_) {
	jj = 0;
      }
      snprintf(buf, buflen, "  [%3zu/%3zu] ", ii, jj);
      os << pfx << buf << pstring(information[ii].residual_error_) << "\n";
    }
    
    os << pfx << "----------------------------------------------------------\n"
       << pfx << "error magnitudes of Jacobian nullspace basis vectors\n"
       << pfx << "----------------------------------------------------------\n";
    
    for (size_t ii(0), jj(0); ii < information.size(); ++ii, ++jj) {
      if ((0 == ii) || information[ii].tag_ != information[ii-1].tag_) {
    	jj = 0;
      }
      snprintf(buf, buflen, "  [%3zu/%3zu] ", ii, jj);
      os << pfx << buf << pstring(information[ii].j_nullspace_residuals_) << "\n";
    }
    
    os << pfx << "----------------------------------------------------------\n"
       << pfx << "error magnitudes of J_bar nullspace basis vectors\n"
       << pfx << "----------------------------------------------------------\n";
    
    for (size_t ii(0), jj(0); ii < information.size(); ++ii, ++jj) {
      if ((0 == ii) || information[ii].tag_ != information[ii-1].tag_) {
    	jj = 0;
      }
      snprintf(buf, buflen, "  [%3zu/%3zu] ", ii, jj);
      os << pfx << buf << pstring(information[ii].jbar_nullspace_residuals_) << "\n";
    }
    
    os << pfx << "----------------------------------------------------------\n"
       << pfx << "magnitudes of J_bar nullspace vectors mapped through J\n"
       << pfx << "----------------------------------------------------------\n";
    
    for (size_t ii(0), jj(0); ii < information.size(); ++ii, ++jj) {
      if ((0 == ii) || information[ii].tag_ != information[ii-1].tag_) {
    	jj = 0;
      }
      snprintf(buf, buflen, "  [%3zu/%3zu] ", ii, jj);
      os << pfx << buf << pstring(information[ii].cross_nullspace_residuals_) << "\n";
    }
    
    os << pfx << "==========================================================\n";
  }
  
}
