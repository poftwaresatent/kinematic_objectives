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

#include "elastic.hpp"
#include "waypoint.hpp"
#include "model.hpp"
#include "pseudo_inverse.hpp"
#include "print.hpp"
#include "task.hpp"


namespace kinematic_elastic {
  
  
  static void perform_prioritization(Matrix const & N_init,
				     vector<Task*> const & tasks,
				     Vector & delta_res,
				     Matrix & N_res,
				     ostream * dbgos,
				     string const & dbgpre)
  {
    delta_res = Vector::Zero(N_init.rows());
    N_res = N_init;
    
    Matrix Jbinv;
    Matrix Nup;
    
    Vector sigma;
    Vector * dbgsigma(0);
    if (dbgos) {
      dbgsigma = &sigma;
      
      string pre (dbgpre);
      pre += "  ";
      print(N_res, *dbgos, "initial nullspace", pre);
      Eigen::JacobiSVD<Matrix> svd;
      svd.compute(N_res, Eigen::ComputeThinU | Eigen::ComputeThinV);
      print(svd.singularValues(), *dbgos, "eigenvalues of nullspace", pre);
    }
    
    for (size_t ii(0); ii < tasks.size(); ++ii) {
      
      if ( ! tasks[ii]->isActive()) {
	if (dbgos) {
	  *dbgos << dbgpre << "task " << ii << " is INACTIVE\n";
	}	
	continue;
      }
      
      if (dbgos) {
	*dbgos << dbgpre << "task " << ii << ":\n";
	string pre (dbgpre);
	pre += "  ";
	print(tasks[ii]->delta_, *dbgos, "task delta", pre);
	print(tasks[ii]->Jacobian_, *dbgos, "task Jacobian", pre);
	Eigen::JacobiSVD<Matrix> svd;
	svd.compute(tasks[ii]->Jacobian_, Eigen::ComputeThinU | Eigen::ComputeThinV);
	print(svd.singularValues(), *dbgos, "eigenvalues of task Jacobian", pre);
	Vector vtmp;
	vtmp = tasks[ii]->delta_ - tasks[ii]->Jacobian_ * delta_res;
	print(vtmp, *dbgos, "delta_comp", pre);
	Matrix mtmp;
	mtmp = tasks[ii]->Jacobian_ * N_res;
	print(mtmp, *dbgos, "J_bar", pre);
      }
      
      pseudo_inverse_moore_penrose(tasks[ii]->Jacobian_ * N_res, Jbinv, &Nup, dbgsigma);
      
      if (dbgos) {
        string pre (dbgpre);
        pre += "  ";
	print(sigma, *dbgos, "eigenvalues of J_bar", pre);
	print(Jbinv, *dbgos, "J_bar pseudo inverse", pre);
	print(Nup, *dbgos, "nullspace update", pre);
	Vector vtmp;
        vtmp = Jbinv * (tasks[ii]->delta_ - tasks[ii]->Jacobian_ * delta_res);
	print(vtmp, *dbgos, "delta update", pre);
        vtmp = tasks[ii]->Jacobian_ * vtmp;
	print(vtmp, *dbgos, "biased task update", pre);
        vtmp = tasks[ii]->Jacobian_ * Jbinv * tasks[ii]->delta_;
	print(vtmp, *dbgos, "unbiased task update", pre);
      }
      
      delta_res += Jbinv * (tasks[ii]->delta_ - tasks[ii]->Jacobian_ * delta_res);
      
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
	print(delta_res, *dbgos, "accumulated delta", pre);
	print(N_res, *dbgos, "accumulated nullspace", pre);
	Eigen::JacobiSVD<Matrix> svd;
	svd.compute(N_res, Eigen::ComputeThinU | Eigen::ComputeThinV);
	print(svd.singularValues(), *dbgos, "eigenvalues of nullspace", pre);
      }
      
    }
    
    if (dbgos) {
      string pre (dbgpre);
      pre += "  ";
      print(N_res, *dbgos, "final nullspace", pre);
      Eigen::JacobiSVD<Matrix> svd;
      svd.compute(N_res, Eigen::ComputeThinU | Eigen::ComputeThinV);
      print(svd.singularValues(), *dbgos, "eigenvalues of nullspace", pre);
    }

  }
  
  
  Elastic::
  Elastic(double timestep, ostream * dbgos, string const & dbgpre)
    : dbgos_(dbgos),
      dbgpre_(dbgpre),
      dbgpre2_(dbgpre + "  "),
      timestep_(timestep)
  {
  }
  
  
  Elastic::
  ~Elastic()
  {
    clear();
  }
  
  
  void Elastic::
  clear()
  {
    for (path_t::iterator ii(path_.begin()); ii != path_.end(); ++ii) {
      delete *ii;
    }
    path_.clear();
  }
  
  
  void Elastic::
  update()
  {
    if (dbgos_) {
      *dbgos_ << dbgpre_ << "\n"
	      << dbgpre_ << "**************************************************\n";
    }
    for (path_t::iterator ii(path_.begin()); ii != path_.end(); ++ii) {
      updateWaypoint(*ii);
    }
  }
  
  
  void Elastic::
  updateWaypoint(Waypoint * wpt)
  {
    wpt->preUpdateHook();
    
    if (dbgos_) {
      *dbgos_ << dbgpre_ << "==================================================\n"
	      << dbgpre_ << "Elastic::updateWaypoint()\n";
      print(wpt->model_.getPosition(), *dbgos_, "current position", dbgpre2_);
      print(wpt->model_.getVelocity(), *dbgos_, "current velocity", dbgpre2_);
    }
    
    for (size_t ii(0); ii < wpt->tasks_.size(); ++ii) {
      wpt->tasks_[ii]->update(wpt->model_);
    }
    for (size_t ii(0); ii < wpt->objectives_.size(); ++ii) {
      wpt->objectives_[ii]->update(wpt->model_);
    }
    
    if (dbgos_) {
      *dbgos_ << dbgpre_ << "--------------------------------------------------\n"
	      << dbgpre_ << "trying without constraints first\n";
    }
    
    ssize_t const ndof(wpt->model_.getPosition().size());
    Vector qdd_t;
    Matrix N_t;
    perform_prioritization(Matrix::Identity(ndof, ndof),
			   wpt->tasks_,
			   qdd_t,
			   N_t,
			   dbgos_,
			   dbgpre_ + "task   ");
    
    Vector qdd_o(Vector::Zero(wpt->model_.getPosition().size()));
    for (size_t ii(0); ii < wpt->objectives_.size(); ++ii) {
      if (wpt->objectives_[ii]->isActive()) {
	Matrix Jinv;
	pseudo_inverse_moore_penrose(wpt->objectives_[ii]->Jacobian_, Jinv);
	qdd_o += Jinv * wpt->objectives_[ii]->delta_;
      }
    }
    
    Vector qdd_res(qdd_t + N_t * qdd_o);
    Vector qd_res(wpt->model_.getVelocity() + timestep_ * qdd_res);
    Vector q_res(wpt->model_.getPosition() + timestep_ * qd_res);
    
    if (dbgos_) {
      print(qdd_res, *dbgos_, "unconstrained acceleration", dbgpre2_);
      print(qd_res, *dbgos_, "resulting unconstrained velocity", dbgpre2_);
      print(q_res, *dbgos_, "resulting unconstrained position", dbgpre2_);
    }
    
    Vector const oldpos(wpt->model_.getPosition()); // will need this in case of constraints
    Vector const oldvel(wpt->model_.getVelocity()); // will need this in case of constraints
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
      return;
    }
    
    if (dbgos_) {
      *dbgos_ << dbgpre_ << "--------------------------------------------------\n"
	      << dbgpre_ << "recomputing with constraints enabled\n";
    }
    
    Vector dq_c;
    Matrix N_c;
    perform_prioritization(Matrix::Identity(ndof, ndof),
			   wpt->constraints_,
			   dq_c,
			   N_c,
			   dbgos_,
			   dbgpre_ + "constr ");
    
    // The constraints cheat with the robot state: they directly work
    // on the positions that would have been achieved without
    // constraints.
    //
    // Semi-open question: after repairing the position and velocity
    // to something consistent with the constraints, do we then then
    // re-run the tasks and objectives? I think yes, to give tasks a
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
      print(wpt->model_.getVelocity(), *dbgos_, "resulting constrained velocity", dbgpre2_);
      print(wpt->model_.getPosition(), *dbgos_, "resulting constrained position", dbgpre2_);
    }
    
    for (size_t ii(0); ii < wpt->tasks_.size(); ++ii) {
      wpt->tasks_[ii]->update(wpt->model_);
    }
    for (size_t ii(0); ii < wpt->objectives_.size(); ++ii) {
      wpt->objectives_[ii]->update(wpt->model_);
    }
    
    // Re-run task priority scheme, but seed it with the constraint nullspace this time.
    
    perform_prioritization(N_c,
			   wpt->tasks_,
			   qdd_t,
			   N_t,
			   dbgos_,
			   dbgpre_ + "task   ");
    
    qdd_o = Vector::Zero(wpt->model_.getPosition().size());
    for (size_t ii(0); ii < wpt->objectives_.size(); ++ii) {
      if (wpt->objectives_[ii]->isActive()) {
	Matrix Jinv;
	pseudo_inverse_moore_penrose(wpt->objectives_[ii]->Jacobian_, Jinv);
	qdd_o += Jinv * wpt->objectives_[ii]->delta_;
      }
    }
    
    qdd_res = qdd_t + N_t * qdd_o;
    qd_res = wpt->model_.getVelocity() + timestep_ * qdd_res;
    q_res = wpt->model_.getPosition() + timestep_ * qd_res;
    
    if (dbgos_) {
      print(qdd_res, *dbgos_, "constrained acceleration", dbgpre2_);
      print(qd_res, *dbgos_, "resulting constrained velocity", dbgpre2_);
      print(q_res, *dbgos_, "resulting constrained position", dbgpre2_);
    }
    
    wpt->model_.update(q_res, qd_res);
  }
  
}
