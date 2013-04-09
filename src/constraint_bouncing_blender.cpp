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

#include <kinematic_objectives/constraint_bouncing_blender.h>
#include <kinematic_objectives/compound_objective.h>
#include <kinematic_objectives/prioritization.h>
#include <kinematic_objectives/kinematic_model.h>
#include <kinematic_objectives/pseudo_inverse.h>
#include <kinematic_objectives/print.h>
#include <kinematic_objectives/objective.h>


namespace kinematic_objectives {
  
  
  ConstraintBouncingBlender::
  ConstraintBouncingBlender(double timestep)
    : Blender("ConstraintBouncingBlender"),
      timestep_(timestep)
  {
  }
  
  
  void ConstraintBouncingBlender::
  update(CompoundObjective * wpt)
  {
    wpt->preUpdateHook();
    
    for (size_t ii(0); ii < wpt->unilateral_constraints_.size(); ++ii) {
      wpt->unilateral_constraints_[ii]->update(wpt->model_);
    }
    for (size_t ii(0); ii < wpt->hard_objectives_.size(); ++ii) {
      wpt->hard_objectives_[ii]->update(wpt->model_);
    }
    for (size_t ii(0); ii < wpt->soft_objectives_.size(); ++ii) {
      wpt->soft_objectives_[ii]->update(wpt->model_);
    }
    
    ssize_t const ndof(wpt->model_.getJointPosition().size());
    Vector & qdd_c(wpt->fb_.constraint_bias_);
    Matrix & N_c(wpt->fb_.constraint_nullspace_projector_);
    prioritization_siciliano1991(Matrix::Identity(ndof, ndof),
				 wpt->unilateral_constraints_,
				 qdd_c,
				 N_c,
				 0,
				 "");
    static double const kp_c(100.0);
    static double const kd_c(20.0);
    qdd_c = kp_c * qdd_c + kd_c * wpt->model_.getJointVelocity();
    
    Vector & qdd_t(wpt->fb_.hard_objective_bias_);
    Matrix & N_t(wpt->fb_.hard_objective_nullspace_projector_);
    prioritization_siciliano1991(N_c,
				 wpt->hard_objectives_,
				 qdd_t,
				 N_t,
				 0,
				 "");
    
    Vector & qdd_o(wpt->fb_.soft_objective_bias_);
    qdd_o = Vector::Zero(wpt->model_.getJointPosition().size());
    for (size_t ii(0); ii < wpt->soft_objectives_.size(); ++ii) {
      if (wpt->soft_objectives_[ii]->isActive()) {
	Matrix Jinv;
	pseudo_inverse_moore_penrose(wpt->soft_objectives_[ii]->getJacobian(), Jinv);
	qdd_o += Jinv * wpt->soft_objectives_[ii]->getBias();
      }
    }
    qdd_o = N_t * qdd_o;
    
    Vector qdd_res(qdd_c + qdd_t + qdd_o);
    Vector qd_res(wpt->model_.getJointVelocity() + timestep_ * qdd_res);
    Vector q_res(wpt->model_.getJointPosition() + timestep_ * qd_res);
    
    wpt->model_.update(q_res, qd_res);
  }
  
}
