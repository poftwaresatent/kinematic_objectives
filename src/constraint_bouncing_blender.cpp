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
  ConstraintBouncingBlender(double stepsize, ostream * dbgos, string const & dbgpre)
    : Blender("ConstraintBouncingBlender", stepsize, dbgos, dbgpre)
  {
  }
  
  
  void ConstraintBouncingBlender::
  update(KinematicModel & model, CompoundObjective * wpt)
  {
    for (size_t ii(0); ii < wpt->unilateral_constraints_.size(); ++ii) {
      wpt->unilateral_constraints_[ii]->update(model);
    }
    for (size_t ii(0); ii < wpt->hard_objectives_.size(); ++ii) {
      wpt->hard_objectives_[ii]->update(model);
    }
    for (size_t ii(0); ii < wpt->soft_objectives_.size(); ++ii) {
      wpt->soft_objectives_[ii]->update(model);
    }
    
    ssize_t const ndof(model.getJointPosition().size());
    Vector & qd_c(wpt->fb_.constraint_bias_);
    Matrix & N_c(wpt->fb_.constraint_nullspace_projector_);
    prioritization_.projectObjectives(Matrix::Identity(ndof, ndof),
				      Vector::Zero(ndof),
				      wpt->unilateral_constraints_,
				      qd_c,
				      N_c);
    
    Vector & qd_t(wpt->fb_.hard_objective_bias_);
    Matrix & N_t(wpt->fb_.hard_objective_nullspace_projector_);
    prioritization_.projectObjectives(N_c,
				      Vector::Zero(ndof), // XXXX likewise wrong, but try to not change everything at the same time...
				      wpt->hard_objectives_,
				      qd_t,
				      N_t);
    
    Vector & qd_o(wpt->fb_.soft_objective_bias_);
    qd_o = Vector::Zero(model.getJointPosition().size());
    for (size_t ii(0); ii < wpt->soft_objectives_.size(); ++ii) {
      if (wpt->soft_objectives_[ii]->isActive()) {
	Matrix Jinv;
	pseudo_inverse_moore_penrose(wpt->soft_objectives_[ii]->getJacobian(), Jinv);
	qd_o += Jinv * wpt->soft_objectives_[ii]->getBias();
      }
    }
    qd_o = N_t * qd_o;
    
    Vector qd_res(qd_c + qd_t + qd_o);
    Vector q_res(model.getJointPosition() + stepsize_ * qd_res);
    
    model.update(q_res, qd_res);
  }
  
}
