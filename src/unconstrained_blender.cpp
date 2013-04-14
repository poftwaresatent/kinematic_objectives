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

#include <kinematic_objectives/unconstrained_blender.h>
#include <kinematic_objectives/compound_objective.h>
#include <kinematic_objectives/prioritization.h>
#include <kinematic_objectives/kinematic_model.h>
#include <kinematic_objectives/pseudo_inverse.h>
#include <kinematic_objectives/print.h>
#include <kinematic_objectives/objective.h>


namespace kinematic_objectives {
  
  
  UnconstrainedBlender::
  UnconstrainedBlender(double timestep)
    : Blender("UnconstrainedBlender"),
      timestep_(timestep)
  {
  }
  
  
  void UnconstrainedBlender::
  update(KinematicModel & model, CompoundObjective * wpt)
  {
    for (size_t ii(0); ii < wpt->hard_objectives_.size(); ++ii) {
      wpt->hard_objectives_[ii]->update(model);
    }
    for (size_t ii(0); ii < wpt->soft_objectives_.size(); ++ii) {
      wpt->soft_objectives_[ii]->update(model);
    }
    
    ssize_t const ndof(model.getJointPosition().size());
    Vector & qdd_t(wpt->fb_.hard_objective_bias_);
    Matrix & N_t(wpt->fb_.hard_objective_nullspace_projector_);
    prioritization_.projectObjectives(Matrix::Identity(ndof, ndof),
				      Vector::Zero(ndof),
				      wpt->hard_objectives_,
				      qdd_t,
				      N_t);
    
    Vector & qdd_o(wpt->fb_.soft_objective_bias_);
    qdd_o = Vector::Zero(model.getJointPosition().size());
    for (size_t ii(0); ii < wpt->soft_objectives_.size(); ++ii) {
      if (wpt->soft_objectives_[ii]->isActive()) {
	Matrix Jinv;
	pseudo_inverse_moore_penrose(wpt->soft_objectives_[ii]->getJacobian(), Jinv);
	qdd_o += Jinv * wpt->soft_objectives_[ii]->getBias();
      }
    }
    qdd_o = N_t * qdd_o;
    
    Vector qdd_res(qdd_t + qdd_o);
    Vector qd_res(model.getJointVelocity() + timestep_ * qdd_res);
    Vector q_res(model.getJointPosition() + timestep_ * qd_res);
    
    model.update(q_res, qd_res);
    
    wpt->fb_.constraint_bias_.resize(0);
    wpt->fb_.constraint_nullspace_projector_.resize(0, 0);
  }
  
}
