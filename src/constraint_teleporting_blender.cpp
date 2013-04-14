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

#include <kinematic_objectives/constraint_teleporting_blender.h>
#include <kinematic_objectives/compound_objective.h>
#include <kinematic_objectives/prioritization.h>
#include <kinematic_objectives/kinematic_model.h>
#include <kinematic_objectives/pseudo_inverse.h>
#include <kinematic_objectives/print.h>
#include <kinematic_objectives/objective.h>


namespace kinematic_objectives {
  
  
  ConstraintTeleportingBlender::
  ConstraintTeleportingBlender(double stepsize, ostream * dbgos, string const & dbgpre)
    : Blender("ConstraintTeleportingBlender", stepsize, dbgos, dbgpre)
  {
  }
  
  
  void ConstraintTeleportingBlender::
  update(KinematicModel & model, CompoundObjective * wpt)
  {
    if (dbgos_) {
      *dbgos_ << dbgpre_ << "==================================================\n"
	      << dbgpre_ << "ConstraintTeleportingBlender::updateCompoundObjective()\n";
      print(model.getJointPosition(), *dbgos_, "current position", dbgpre_ + "  ");
      print(model.getJointVelocity(), *dbgos_, "current velocity", dbgpre_ + "  ");
    }
    
    for (size_t ii(0); ii < wpt->hard_objectives_.size(); ++ii) {
      wpt->hard_objectives_[ii]->update(model);
    }
    for (size_t ii(0); ii < wpt->soft_objectives_.size(); ++ii) {
      wpt->soft_objectives_[ii]->update(model);
    }
    
    if (dbgos_) {
      *dbgos_ << dbgpre_ << "--------------------------------------------------\n"
	      << dbgpre_ << "trying without constraints first\n";
    }
    
    ssize_t const ndof(model.getJointPosition().size());
    Vector & qd_t(wpt->fb_.hard_objective_bias_);
    Matrix & N_t(wpt->fb_.hard_objective_nullspace_projector_);
    prioritization_.dbgos_ = dbgos_;
    prioritization_.dbgpre_ = dbgpre_ + "  ";
    prioritization_.projectObjectives(Matrix::Identity(ndof, ndof),
				      Vector::Zero(ndof),
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
    
    Vector qd_res(qd_t + qd_o);
    Vector q_res(model.getJointPosition() + stepsize_ * qd_res);
    
    if (dbgos_) {
      print(qd_res, *dbgos_, "unconstrained velocity", dbgpre_ + "  ");
      print(q_res, *dbgos_, "resulting unconstrained position", dbgpre_ + "  ");
    }
    
    Vector const oldpos(model.getJointPosition()); // will need this in case of constraints
    Vector const oldvel(model.getJointVelocity()); // will need this in case of constraints
    model.update(q_res, qd_res);
    
    bool need_constraints(false);
    for (size_t ii(0); ii < wpt->unilateral_constraints_.size(); ++ii) {
      wpt->unilateral_constraints_[ii]->update(model);
      if (wpt->unilateral_constraints_[ii]->isActive()) {
	if (dbgos_) {
	  *dbgos_ << dbgpre_ << "constraint [" << ii << "] is active\n";
	}
	need_constraints = true;
      }
    }
    
    if ( ! need_constraints) {
      if (dbgos_) {
	*dbgos_ << dbgpre_ << "all constraints are inactive\n";
      }
      wpt->fb_.constraint_bias_.resize(0);
      wpt->fb_.constraint_nullspace_projector_.resize(0, 0);
      return;
    }
    
    if (dbgos_) {
      *dbgos_ << dbgpre_ << "--------------------------------------------------\n"
	      << dbgpre_ << "recomputing with constraints enabled\n";
    }
    
    Vector & qd_c(wpt->fb_.constraint_bias_);
    Matrix & N_c(wpt->fb_.constraint_nullspace_projector_);
    prioritization_.projectObjectives(Matrix::Identity(ndof, ndof),
				      Vector::Zero(ndof),
				      wpt->unilateral_constraints_,
				      qd_c,
				      N_c);
    
    // Semi-open question: after repairing the position and velocity
    // to something consistent with the constraints, do we then then
    // re-run the hard and soft objectives? I think yes, to give soft
    // objectives a chance to get fulfilled within the
    // constraints. But that does not seem to work quite yet...
    //
    // Also, wouldn't the position and velocity change for the
    // constraints then influence the constraints themselves? More
    // importantly, shouldn't the progress monitors be aware of what
    // we've done to achieve constraints? But what if the constraints
    // that triggered this then get switched off by the new
    // update... that would make the whole picture rather
    // inconsistent. Hm...
    //
    // Note that the non-constrained velocity would be (q_res + qd_c -
    // oldpos) / but we're pre-multiplying with N_c and qd_c is
    // perpendicular to that so we don't need to add it.
    //
    model.update(q_res + stepsize_ * qd_c, N_c * (q_res - oldpos));
    
    if (dbgos_) {
      print(qd_c, *dbgos_, "velocity to satisfy constraints", dbgpre_ + "  ");
      print(N_c, *dbgos_, "nullspace of constrains", dbgpre_ + "  ");
      print(model.getJointPosition(), *dbgos_, "resulting position", dbgpre_ + "  ");
      print(model.getJointVelocity(), *dbgos_, "resulting velocity", dbgpre_ + "  ");
    }
    
    for (size_t ii(0); ii < wpt->hard_objectives_.size(); ++ii) {
      wpt->hard_objectives_[ii]->update(model);
    }
    for (size_t ii(0); ii < wpt->soft_objectives_.size(); ++ii) {
      wpt->soft_objectives_[ii]->update(model);
    }
    
    // Re-run objective priority scheme, but seed it with the constraint nullspace this time.
    
    prioritization_.projectObjectives(N_c,
				      // XXXX this here is wrong, just
				      // shows an old bug (hopefully)
				      // that will now be removed
				      // (hopefully)
				      Vector::Zero(ndof),
				      wpt->hard_objectives_,
				      qd_t,
				      N_t);
    
    qd_o = Vector::Zero(model.getJointPosition().size());
    for (size_t ii(0); ii < wpt->soft_objectives_.size(); ++ii) {
      if (wpt->soft_objectives_[ii]->isActive()) {
	Matrix Jinv;
	pseudo_inverse_moore_penrose(wpt->soft_objectives_[ii]->getJacobian(), Jinv);
	qd_o += Jinv * wpt->soft_objectives_[ii]->getBias();
      }
    }
    qd_o = N_t * qd_o;
    
    qd_res = qd_t + qd_o;
    q_res = model.getJointPosition() + stepsize_ * qd_res;
    
    if (dbgos_) {
      print(qd_res, *dbgos_, "constrained velocity", dbgpre_ + "  ");
      print(q_res, *dbgos_, "resulting constrained position", dbgpre_ + "  ");
    }
    
    model.update(q_res, qd_res);
  }
  
}
