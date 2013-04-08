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

#include <kinematic_objectives/blender.h>
#include <kinematic_objectives/compound_objective.h>
#include <kinematic_objectives/kinematic_model.h>
#include <kinematic_objectives/pseudo_inverse.h>
#include <kinematic_objectives/print.h>
#include <kinematic_objectives/objective.h>


namespace kinematic_objectives {
  
  
  static void perform_prioritization(/** initial nullspace projector */
				     Matrix const & N_init,
				     /** hierarchy of objectives to fuse/blend */
				     vector<Objective*> const & objectives,
				     /** resulting (fused/blended) bias */
				     Vector & bias_res,
				     /** accumulated nullspace projector at the end of the fusion */
				     Matrix & N_res,
				     /** stream for debug output (use 0 for silent operation) */
				     ostream * dbgos,
				     /** prefix for debug output (prepended to each line) */
				     string const & dbgpre)
  {
    bias_res = Vector::Zero(N_init.rows());
    N_res = N_init;
    
    Matrix Jbinv;		// pseudo-inverse of J_bar (which is J * N)
    Matrix Nup;			// nullspace updater: N -= N_up at each hierarchy level
    
    if (dbgos) {
      string pre (dbgpre);
      pre += "  ";
      print(N_res, *dbgos, "initial nullspace", pre);
      Eigen::JacobiSVD<Matrix> svd;
      svd.compute(N_res, Eigen::ComputeThinU | Eigen::ComputeThinV);
      print(svd.singularValues(), *dbgos, "eigenvalues of nullspace", pre);
    }
    
    for (size_t ii(0); ii < objectives.size(); ++ii) {
      
      if ( ! objectives[ii]->isActive()) {
	if (dbgos) {
	  *dbgos << dbgpre << "objective " << ii << " is INACTIVE\n";
	}	
	continue;
      }
      
      if (dbgos) {
	*dbgos << dbgpre << "objective " << ii << " \"" << objectives[ii]->name_ << "\":\n";
	string pre (dbgpre);
	pre += "  ";
	print(objectives[ii]->getBias(), *dbgos, "objective bias", pre);
	print(objectives[ii]->getJacobian(), *dbgos, "objective Jacobian", pre);
	Eigen::JacobiSVD<Matrix> svd;
	svd.compute(objectives[ii]->getJacobian(), Eigen::ComputeThinU | Eigen::ComputeThinV);
	print(svd.singularValues(), *dbgos, "eigenvalues of objective Jacobian", pre);
	Vector vtmp;
	vtmp = objectives[ii]->getBias() - objectives[ii]->getJacobian() * bias_res;
	print(vtmp, *dbgos, "bias_comp", pre);
	Matrix mtmp;
	mtmp = objectives[ii]->getJacobian() * N_res;
	print(mtmp, *dbgos, "J_bar", pre);
      }
      
      pseudo_inverse_moore_penrose(objectives[ii]->getJacobian() * N_res, Jbinv, &Nup,
				   &objectives[ii]->jbar_svd_);
      
      if (dbgos) {
        string pre (dbgpre);
        pre += "  ";
	print(objectives[ii]->jbar_svd_.singular_values, *dbgos, "eigenvalues of J_bar", pre);
	print(Jbinv, *dbgos, "J_bar pseudo inverse", pre);
	print(Nup, *dbgos, "nullspace projector update", pre);
	Vector vtmp;
        vtmp = Jbinv * (objectives[ii]->getBias() - objectives[ii]->getJacobian() * bias_res);
	print(vtmp, *dbgos, "bias update", pre);
        vtmp = objectives[ii]->getJacobian() * vtmp;
	print(vtmp, *dbgos, "biased objective update", pre);
        vtmp = objectives[ii]->getJacobian() * Jbinv * objectives[ii]->getBias();
	print(vtmp, *dbgos, "unbiased objective update", pre);
      }
      
      bias_res += Jbinv * (objectives[ii]->getBias() - objectives[ii]->getJacobian() * bias_res);
      
      if (Nup.cols() == 0) {
	if (dbgos) {
	  *dbgos << dbgpre << "NO DEGREES OF FREEDOM LEFT\n";
	}
	break;
      }
      
      N_res -= Nup;
      
      if (dbgos) {
        string pre (dbgpre);
        pre += "  ";
	print(bias_res, *dbgos, "accumulated bias", pre);
	print(N_res, *dbgos, "accumulated nullspace projector", pre);
	Eigen::JacobiSVD<Matrix> svd;
	svd.compute(N_res, Eigen::ComputeThinU | Eigen::ComputeThinV);
	print(svd.singularValues(), *dbgos, "eigenvalues of nullspace projector", pre);
      }
      
    }
    
    if (dbgos) {
      string pre (dbgpre);
      pre += "  ";
      print(N_res, *dbgos, "final nullspace projectir", pre);
      Eigen::JacobiSVD<Matrix> svd;
      svd.compute(N_res, Eigen::ComputeThinU | Eigen::ComputeThinV);
      print(svd.singularValues(), *dbgos, "eigenvalues of final nullspace projector", pre);
    }

  }
  
  
  Blender::
  Blender(double timestep, ostream * dbgos, string const & dbgpre)
    : dbgos_(dbgos),
      dbgpre_(dbgpre),
      dbgpre2_(dbgpre + "  "),
      timestep_(timestep)
  {
  }
  
  
  Blender::
  ~Blender()
  {
    clear();
  }
  
  
  void Blender::
  clear()
  {
    for (path_t::iterator ii(path_.begin()); ii != path_.end(); ++ii) {
      delete *ii;
    }
    path_.clear();
  }
  
  
  void Blender::
  update()
  {
    if (dbgos_) {
      *dbgos_ << dbgpre_ << "\n"
	      << dbgpre_ << "**************************************************\n";
    }
    for (path_t::iterator ii(path_.begin()); ii != path_.end(); ++ii) {
      updateCompoundObjective(*ii);
    }
  }
  
  
  void Blender::
  updateCompoundObjective(CompoundObjective * wpt)
  {
    wpt->preUpdateHook();
    
    if (dbgos_) {
      *dbgos_ << dbgpre_ << "==================================================\n"
	      << dbgpre_ << "Blender::updateCompoundObjective()\n";
      print(wpt->model_.getJointPosition(), *dbgos_, "current position", dbgpre2_);
      print(wpt->model_.getJointVelocity(), *dbgos_, "current velocity", dbgpre2_);
    }
    
    for (size_t ii(0); ii < wpt->hard_objectives_.size(); ++ii) {
      wpt->hard_objectives_[ii]->update(wpt->model_);
    }
    for (size_t ii(0); ii < wpt->soft_objectives_.size(); ++ii) {
      wpt->soft_objectives_[ii]->update(wpt->model_);
    }
    
    if (dbgos_) {
      *dbgos_ << dbgpre_ << "--------------------------------------------------\n"
	      << dbgpre_ << "trying without constraints first\n";
    }
    
    ssize_t const ndof(wpt->model_.getJointPosition().size());
    Vector & qdd_t(wpt->fb_.hard_objective_bias_);
    Matrix & N_t(wpt->fb_.hard_objective_nullspace_projector_);
    perform_prioritization(Matrix::Identity(ndof, ndof),
			   wpt->hard_objectives_,
			   qdd_t,
			   N_t,
			   dbgos_,
			   dbgpre_ + "  ");
    
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
    
    Vector qdd_res(qdd_t + qdd_o);
    Vector qd_res(wpt->model_.getJointVelocity() + timestep_ * qdd_res);
    Vector q_res(wpt->model_.getJointPosition() + timestep_ * qd_res);
    
    if (dbgos_) {
      print(qdd_res, *dbgos_, "unconstrained acceleration", dbgpre2_);
      print(qd_res, *dbgos_, "resulting unconstrained velocity", dbgpre2_);
      print(q_res, *dbgos_, "resulting unconstrained position", dbgpre2_);
    }
    
    Vector const oldpos(wpt->model_.getJointPosition()); // will need this in case of constraints
    Vector const oldvel(wpt->model_.getJointVelocity()); // will need this in case of constraints
    wpt->model_.update(q_res, qd_res);
    
    bool need_constraints(false);
    for (size_t ii(0); ii < wpt->constraints_.size(); ++ii) {
      wpt->constraints_[ii]->update(wpt->model_);
      if (wpt->constraints_[ii]->isActive()) {
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
    
    Vector & dq_c(wpt->fb_.constraint_bias_);
    Matrix & N_c(wpt->fb_.constraint_nullspace_projector_);
    perform_prioritization(Matrix::Identity(ndof, ndof),
			   wpt->constraints_,
			   dq_c,
			   N_c,
			   dbgos_,
			   dbgpre_ + "  ");
    
    // The constraints cheat with the robot state: they directly work
    // on the positions that would have been achieved without
    // constraints.
    //
    // Semi-open question: after repairing the position and velocity
    // to something consistent with the constraints, do we then then
    // re-run the hard and soft objectives? I think yes, to give soft objectives a
    // chance to get fulfilled within the constraints. But that does
    // not seem to work quite yet...
    //
    // Also, wouldn't the position and velocity change for the
    // constraints then influence the constraints themselves? We
    // should thus re-run the constraints as well, possibly leading to
    // another correction and so forth ad infinitum. But the nullspace
    // of the constraints at least should not change too much, so we
    // can probably skip the chicken-and-egg constraint update
    // problem.
    //
    // Note that the non-constrained velocity would be (q_res + dq_c -
    // oldpos) / timestep_ but we're pre-multiplying with N_c and dq_c
    // is perpendicular to that so we don't need to add it.
    //
    wpt->model_.update(q_res + dq_c, N_c * (q_res - oldpos) / timestep_);
    
    if (dbgos_) {
      print(dq_c, *dbgos_, "position correction to satisfy constraints", dbgpre2_);
      print(N_c, *dbgos_, "nullspace of constrains", dbgpre2_);
      print(wpt->model_.getJointVelocity(), *dbgos_, "resulting constrained velocity", dbgpre2_);
      print(wpt->model_.getJointPosition(), *dbgos_, "resulting constrained position", dbgpre2_);
    }
    
    for (size_t ii(0); ii < wpt->hard_objectives_.size(); ++ii) {
      wpt->hard_objectives_[ii]->update(wpt->model_);
    }
    for (size_t ii(0); ii < wpt->soft_objectives_.size(); ++ii) {
      wpt->soft_objectives_[ii]->update(wpt->model_);
    }
    
    // Re-run objective priority scheme, but seed it with the constraint nullspace this time.
    
    perform_prioritization(N_c,
			   wpt->hard_objectives_,
			   qdd_t,
			   N_t,
			   dbgos_,
			   dbgpre_ + "objective   ");
    
    qdd_o = Vector::Zero(wpt->model_.getJointPosition().size());
    for (size_t ii(0); ii < wpt->soft_objectives_.size(); ++ii) {
      if (wpt->soft_objectives_[ii]->isActive()) {
	Matrix Jinv;
	pseudo_inverse_moore_penrose(wpt->soft_objectives_[ii]->getJacobian(), Jinv);
	qdd_o += Jinv * wpt->soft_objectives_[ii]->getBias();
      }
    }
    
    qdd_res = qdd_t + N_t * qdd_o;
    qd_res = wpt->model_.getJointVelocity() + timestep_ * qdd_res;
    q_res = wpt->model_.getJointPosition() + timestep_ * qd_res;
    
    if (dbgos_) {
      print(qdd_res, *dbgos_, "constrained acceleration", dbgpre2_);
      print(qd_res, *dbgos_, "resulting constrained velocity", dbgpre2_);
      print(q_res, *dbgos_, "resulting constrained position", dbgpre2_);
    }
    
    wpt->model_.update(q_res, qd_res);
  }
  
}
