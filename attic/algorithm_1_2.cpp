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

#include "algorithm.hpp"
#include "pseudo_inverse.hpp"
#include "print.hpp"
#include "task.hpp"
#include "joint_limits.hpp"
#include <Eigen/LU>
#include <iostream>


namespace kinematic_elastic {
  
  
  /**
     \pre tl must not be empty, and tl[0] must be non-singular.
  */
  static Vector compute_dq (TaskData * alpha,
			    TaskData * beta,
			    TaskData * gamma,
			    ostream * dbgos,
			    char const * dbgpre)
  {
    Matrix Ja_inv;
    pseudo_inverse_nonsingular (alpha->Jacobian_, Ja_inv);
    Vector dq(Ja_inv * alpha->delta_);
    
    if (dbgos) {
      *dbgos << dbgpre << "primary task:\n";
      string pre (dbgpre);
      pre += "  ";
      print (alpha->delta_, *dbgos, "dxa",    pre);
      print (alpha->Jacobian_, *dbgos, "Ja",     pre);
      print (Ja_inv,    *dbgos, "Ja_inv", pre);
      print (dq,        *dbgos, "dq",     pre);
    }
    
    if (0 == beta) {
      if (dbgos) {
	*dbgos << dbgpre << "DONE (single task)\n";
      }
      return dq;
    }
    
    Matrix Jb_inv;
    pseudo_inverse_damped (beta->Jacobian_, 1.0 / beta->step_hint_, Jb_inv);
    size_t const ndof(alpha->Jacobian_.cols());
    Matrix Na(Matrix::Identity(ndof, ndof) - Ja_inv * alpha->Jacobian_);
    dq +=  Na * Jb_inv * beta->delta_;
    
    if (dbgos) {
      Matrix Na_times_Jb_inv(Na * Jb_inv);
      Vector dqb(Na_times_Jb_inv * beta->delta_);
      
      *dbgos << dbgpre << "secondary task:\n";
      string pre (dbgpre);
      pre += "  ";
      print (Na,              *dbgos, "Na", pre);
      print (beta->delta_,       *dbgos, "dxb", pre);
      print (beta->Jacobian_,       *dbgos, "Jb", pre);
      print (Jb_inv,          *dbgos, "Jb_inv", pre);
      print (Na_times_Jb_inv, *dbgos, "Na * Jb_inv", pre);
      print (dqb,             *dbgos, "Na * Jb_inv * dxb", pre);
      print (dq,              *dbgos, "dq", pre);
    }
    
    if (0 == gamma) {
      if (dbgos) {
	*dbgos << dbgpre << "DONE (two tasks)\n";
      }
      return dq;
    }
    
    Matrix Nb(Matrix::Identity(ndof, ndof) - Jb_inv * beta->Jacobian_);
    Matrix Jc_inv;
    pseudo_inverse_damped (gamma->Jacobian_, 1.0 / gamma->step_hint_, Jc_inv);
    dq += Na * Nb * Jc_inv * gamma->delta_;
    
    if (dbgos) {
      *dbgos << dbgpre << "remainder:\n";
      string pre (dbgpre);
      pre += "  ";
      print (Nb,        *dbgos, "Nb",     pre);
      print (gamma->delta_, *dbgos, "dxc",    pre);
      print (gamma->Jacobian_, *dbgos, "Jc",     pre);
      print (Jc_inv,    *dbgos, "Jc_inv", pre);
      print (dq,        *dbgos, "dq",     pre);
      *dbgos << dbgpre << "DONE\n";
    }
    
    return dq;
  }


  static inline Vector compute_dq(vector<TaskData*> const & tl,
				  ostream * dbgos,
				  char const * dbgpre)
  {
    return compute_dq (tl.size() > 0 ? tl[0] : 0,
		       tl.size() > 1 ? tl[1] : 0,
		       tl.size() > 2 ? tl[2] : 0,
		       dbgos,
		       dbgpre);
  }
  
  
  Vector algorithm1 (JointLimits & joint_limits,
		     Vector const & state,
		     vector<TaskData*> const & tasklist,
		     ostream * dbgos,
		     char const * dbgpre)
  {
    size_t const ndof(state.size());
    Vector dq(Vector::Zero(ndof));
    
    if ( ! tasklist.empty()) {
      TaskData * alpha(0);
      TaskData * beta(0);
      TaskData * gamma(0);
      TaskData tmp1, tmp2;
      
      alpha = tasklist[0];
      if (tasklist.size() > 1) {
	tmp1.stack(*tasklist[0], *tasklist[1]);
	beta = &tmp1;
      }
      if (tasklist.size() > 3) {
	tmp2.stack(tasklist.begin() + 2, tasklist.end());
	gamma = &tmp2;
      }
      else if (tasklist.size() == 3) {
	gamma = tasklist[2];
      }
      dq = compute_dq(alpha, beta, gamma, dbgos, dbgpre);
    }
    
    Vector const next_state(state + dq);
    joint_limits.update(next_state);
    if ( ! joint_limits.isActive()) {
      if (dbgos) {
	*dbgos << dbgpre << "joint limit check passed\n";
      }
      return dq;
    }
    
    if (tasklist.empty()) {
      dq = compute_dq(&joint_limits, 0, 0, dbgos, dbgpre);
      if (dbgos) {
	*dbgos << dbgpre << "joint limits adjusted, but empty task list\n";
      }
      return dq;
    }
    
    // Try stacking the joint limits on top of the primary task. That
    // means knocking out the locked joints from the primary task's
    // Jacobian, too.
    
    if (dbgos) {
      string pre (dbgpre);
      pre += "  ";
      *dbgos << dbgpre << "try removing joints from the primary task...\n";
      print (joint_limits.Jacobian_, *dbgos, "J_lim", pre);
    }
    
    Matrix J_try(tasklist[0]->Jacobian_);
    // XXXX it should also work without knocking the columns to zero,
    // just because of the nullspace of the joint_limits.Jacobian_.
    for (size_t ii(0); ii < joint_limits.locked_joints_.size(); ++ii) {
      J_try.block(0, joint_limits.locked_joints_[ii], J_try.rows(), 1) = Vector::Zero(J_try.rows());
      if (dbgos) {
	*dbgos << dbgpre << "  joint " << joint_limits.locked_joints_[ii] << " is locked!\n";
      }
    }
    Eigen::FullPivLU<Matrix> plu(J_try * J_try.transpose());
    if (plu.isInvertible()) {
      if (dbgos) {
	*dbgos << dbgpre << "primary task remains achievable with locked joints\n";
      }
      
      TaskData * alpha(0);
      TaskData * beta(0);
      TaskData * gamma(0);
      TaskData tmp1, tmp2;
      
      tmp1.stack(joint_limits, *tasklist[0]);
      // grr, spurious extra work... and maybe not even necessary
      stackMatrix(joint_limits.Jacobian_, J_try, tmp1.Jacobian_);
      alpha = &tmp1;
      
      if (tasklist.size() > 1) {
	// Now this is a bit of an open question. Ideally we should
	// select a set of tasks that are somewhat likely to be not
	// overly conflicting. Reasonable heuristics are either:
	//  A: joint_limits, tasklist[0], tasklist[1]
	//     but that is likely to kill tasklist[1] entirely.
	//  B: tasklist[0], tasklist[1]
	//     i.e. try to achieve both (ignoring joint limits, which
	//     get respected due to the nullspace of tl[0]).
	//  C: just tasklist[1]
	//     This would appear to heighten our chances of achieving
	//     the secondary task, which is likely to be something that
	//     users actually care about (whereas ternary etc are
	//     conceivably just meant as "nice to have" anyway).
	//
	// Tried (B) and (C) but saw no difference, so going with (C)
	// because that avoids creating yet another tmp.
	//
	beta = tasklist[1];
      }
      if (tasklist.size() > 3) {
	// At the time of writing, this is the same as tmp2 in the
	// first trial
	tmp2.stack(tasklist.begin() + 2, tasklist.end());
	gamma = &tmp2;
      }
      else if (tasklist.size() == 3) {
	gamma = tasklist[2];
      }
      // We could also knock the locked joints from the other tasks,
      // but the null-space projection "should" take care of it anyway
      // from tl[1] onward.
      dq = compute_dq(alpha, beta, gamma, dbgos, dbgpre);
      if (dbgos) {
	*dbgos << dbgpre << "joint limits merged into primary task\n";
      }
      return dq;
    }
    
    // Well, that didn't work, which means that the primary task
    // cannot be achieved due to joint limits. The best we can do is
    // prepend the joint limits to the task list.
    //
    // Again, some open questions regarding what to stack together for
    // the second task that gets passed to compute_dq(). But it seems
    // natural to now push the (original) secondary into the "nice to
    // have" category and keep the (original) primary "please try as
    // hard as you can."
    
    TaskData * alpha(&joint_limits);
    TaskData * beta(tasklist[0]);
    TaskData * gamma(0);
    TaskData tmp;
    if (tasklist.size() > 1) {
      if (tasklist.size() > 2) {
	tmp.stack(tasklist.begin() + 1, tasklist.end());
	gamma = &tmp;
      }
      else {
	gamma = tasklist[1];
      }
    }
    dq = compute_dq(alpha, beta, gamma, dbgos, dbgpre);
    if (dbgos) {
      *dbgos << dbgpre << "joint limits take precedence over primary task\n";
    }
    return dq;
  }
  
  
  Vector algorithm2 (JointLimits & joint_limits,
		     Vector const & state,
		     vector<TaskData*> const & tasklist,
		     ostream * dbgos,
		     char const * dbgpre)
  {
    size_t const ndof(state.size());
    Vector dq(Vector::Zero(ndof));
    
    if ( ! tasklist.empty()) {
      vector<TaskData*> tl;
      TaskData tmp;
      tl.push_back(tasklist[0]);
      if (tasklist.size() > 1) {
	tl.push_back(tasklist[1]);
      }
      if (tasklist.size() > 2) {
	if (tasklist.size() == 3) {
	  tl.push_back(tasklist[2]);
	}
	else {
	  tmp.stack(tasklist.begin() + 2, tasklist.end());
	  tl.push_back(&tmp);
	}
      }
      dq = compute_dq(tl, dbgos, dbgpre);
    }
    
    Vector const next_state(state + dq);
    joint_limits.update(next_state);
    if ( ! joint_limits.isActive()) {
      if (dbgos) {
	*dbgos << dbgpre << "joint limit check passed\n";
      }
      return dq;
    }
    
    if (tasklist.empty()) {
      vector<TaskData*> tl;
      tl.push_back(&joint_limits);
      dq = compute_dq(tl, dbgos, dbgpre);
      if (dbgos) {
	*dbgos << dbgpre << "joint limits adjusted, but empty task list\n";
      }
      return dq;
    }
    
    // Try stacking the joint limits on top of the primary task. That
    // means knocking out the locked joints from the primary task's
    // Jacobian, too.
    
    if (dbgos) {
      string pre (dbgpre);
      pre += "  ";
      *dbgos << dbgpre << "try removing joints from the primary task...\n";
      print (joint_limits.Jacobian_, *dbgos, "J_lim", pre);
    }
    
    Matrix J_try(tasklist[0]->Jacobian_);
    for (size_t ii(0); ii < joint_limits.locked_joints_.size(); ++ii) {
      J_try.block(0, joint_limits.locked_joints_[ii], J_try.rows(), 1) = Vector::Zero(J_try.rows());
      if (dbgos) {
	*dbgos << dbgpre << "  joint " << joint_limits.locked_joints_[ii] << " is locked!\n";
      }
    }
    Eigen::FullPivLU<Matrix> plu(J_try * J_try.transpose());
    if (plu.isInvertible()) {
      if (dbgos) {
	*dbgos << dbgpre << "primary task remains achievable with locked joints\n";
      }
      vector<TaskData*> tl;
      TaskData tmp1, tmp2;
      tmp1.stack(joint_limits, *tasklist[0]);
      // grr, spurious extra work... and maybe not even necessary
      stackMatrix(joint_limits.Jacobian_, J_try, tmp1.Jacobian_);
      tl.push_back(&tmp1);
      if (tasklist.size() > 1) {
	tl.push_back(tasklist[1]);
      }
      if (tasklist.size() > 2) {
	if (tasklist.size() == 3) {
	  tl.push_back(tasklist[2]);
	}
	else {
	  tmp2.stack(tasklist.begin() + 2, tasklist.end());
	  tl.push_back(&tmp2);
	}
      }
      // We could also knock the locked joints from the other tasks,
      // but the null-space projection "should" take care of it anyway
      // from tl[1] onward.
      dq = compute_dq(tl, dbgos, dbgpre);
      if (dbgos) {
	*dbgos << dbgpre << "joint limits merged into primary task\n";
      }
      return dq;
    }
    
    // Well, that didn't work, which means that the primary task
    // cannot be achieved due to joint limits. The best we can do is
    // prepend the joint limits to the task list.
    
    vector<TaskData*> tl;
    TaskData tmp;
    tl.push_back(&joint_limits);
    tl.push_back(tasklist[0]);
    if (tasklist.size() > 1) {
      if (tasklist.size() == 2) {
	tl.push_back(tasklist[1]);
      }
      else {
	tmp.stack(tasklist.begin() + 1, tasklist.end());
	tl.push_back(&tmp);
      }
    }
    dq = compute_dq(tl, dbgos, dbgpre);
    if (dbgos) {
      *dbgos << dbgpre << "joint limits take precedence over primary task\n";
    }
    return dq;
  }
  
}
