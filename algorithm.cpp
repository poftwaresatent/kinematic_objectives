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
#include <Eigen/LU>
#include <iostream>


namespace kinematic_elastic {
  
  
  /**
     \pre tl must not be empty, and tl[0] must be non-singular.
  */
  static Vector compute_dq (tasklist_t const & tl,
			    ostream * dbgos,
			    char const * dbgpre)
  {
    Matrix Ja_inv;
    pseudo_inverse_nonsingular (tl[0]->Jx, Ja_inv);
    Vector dq(Ja_inv * tl[0]->dx);
    
    if (dbgos) {
      *dbgos << dbgpre << "primary task:\n";
      string pre (dbgpre);
      pre += "  ";
      print (tl[0]->dx, *dbgos, "dxa",    pre);
      print (tl[0]->Jx, *dbgos, "Ja",     pre);
      print (Ja_inv,    *dbgos, "Ja_inv", pre);
      print (dq,        *dbgos, "dq",     pre);
    }
    
    if (tl.size() == 1) {
      if (dbgos) {
	*dbgos << dbgpre << "DONE (single task)\n";
      }
      return dq;
    }
    
    Matrix Jb_inv;
    pseudo_inverse_damped (tl[1]->Jx, 1.0 / tl[1]->b_max, Jb_inv);
    size_t const ndof(tl[0]->Jx.cols());
    Matrix Na(Matrix::Identity(ndof, ndof) - Ja_inv * tl[0]->Jx);
    dq +=  Na * Jb_inv * tl[1]->dx;
    
    if (dbgos) {
      Matrix Na_times_Jb_inv(Na * Jb_inv);
      Vector dqb(Na_times_Jb_inv * tl[1]->dx);
      
      *dbgos << dbgpre << "secondary task:\n";
      string pre (dbgpre);
      pre += "  ";
      print (Na,              *dbgos, "Na", pre);
      print (tl[1]->dx,       *dbgos, "dxb", pre);
      print (tl[1]->Jx,       *dbgos, "Jb", pre);
      print (Jb_inv,          *dbgos, "Jb_inv", pre);
      print (Na_times_Jb_inv, *dbgos, "Na * Jb_inv", pre);
      print (dqb,             *dbgos, "Na * Jb_inv * dxb", pre);
      print (dq,              *dbgos, "dq", pre);
    }
    
    if (tl.size() == 2) {
      if (dbgos) {
	*dbgos << dbgpre << "DONE (two tasks)\n";
      }
      return dq;
    }
    
    Matrix Nb(Matrix::Identity(ndof, ndof) - Jb_inv * tl[1]->Jx);
    Matrix Jc_inv;
    pseudo_inverse_damped (tl[2]->Jx, 1.0 / tl[2]->b_max, Jc_inv);
    dq += Na * Nb * Jc_inv * tl[2]->dx;
    
    if (dbgos) {
      *dbgos << dbgpre << "remainder:\n";
      string pre (dbgpre);
      pre += "  ";
      print (Nb,        *dbgos, "Nb",     pre);
      print (tl[2]->dx, *dbgos, "dxc",    pre);
      print (tl[2]->Jx, *dbgos, "Jc",     pre);
      print (Jc_inv,    *dbgos, "Jc_inv", pre);
      print (dq,        *dbgos, "dq",     pre);
      *dbgos << dbgpre << "DONE\n";
    }
    
    return dq;
  }
  
  
  Vector algorithm (Model const & model,
		    Vector const & state,
		    tasklist_t const & tasklist,
		    ostream * dbgos,
		    char const * dbgpre)
  {
    size_t const ndof(state.size());
    Vector dq(Vector::Zero(ndof));
    
    if ( ! tasklist.empty()) {
      tasklist_t tl;
      task_s tmp1, tmp2;
      tl.push_back(tasklist[0]);
      if (tasklist.size() > 1) {
	tmp1 = stack(tasklist, 0, 2);
	tl.push_back(&tmp1);
      }
      if (tasklist.size() > 3) {
	tmp2 = stack(tasklist, 2, tasklist.size() - 2);
	tl.push_back(&tmp2);
      }
      else if (tasklist.size() == 3) {
	tl.push_back(tasklist[2]);
      }
      dq = compute_dq(tl, dbgos, dbgpre);
    }
    
    Vector const next_state(state + dq);
    if (model.checkJointLimits(next_state)) {
      if (dbgos) {
	*dbgos << dbgpre << "joint limit check passed\n";
      }
      return dq;
    }
    
    task_s t_lim;
    vector<size_t> locked;
    // Now here's an intricate point: need to lock based on what we
    // would get AFTER this iteration.
    model.createJointLimitTask(next_state, t_lim, locked);
    
    if (tasklist.empty()) {
      tasklist_t tl;
      tl.push_back(&t_lim);
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
      print (t_lim.Jx, *dbgos, "J_lim", pre);
    }
    
    Matrix J_try(tasklist[0]->Jx);
    for (size_t ii(0); ii < locked.size(); ++ii) {
      J_try.block(0, locked[ii], J_try.rows(), 1) = Vector::Zero(J_try.rows());
      if (dbgos) {
	*dbgos << dbgpre << "  joint " << locked[ii] << " is locked!\n";
      }
    }
    Eigen::FullPivLU<Matrix> plu(J_try * J_try.transpose());
    if (plu.isInvertible()) {
      if (dbgos) {
	*dbgos << dbgpre << "primary task remains achievable with locked joints\n";
      }
      tasklist_t tl;
      task_s tmp1, tmp2, tmp3;
      tmp1 = stack(t_lim, *tasklist[0]);
      // grr, spurious extra work...
      stack(t_lim.Jx, J_try, tmp1.Jx);
      tl.push_back(&tmp1);
      if (tasklist.size() > 1) {
	// Now this is a bit of an open question. Ideally we should
	// select a set of tasks that are somewhat likely to be not
	// overly conflicting. Reasonable heuristics are either:
	//  A: t_lim, tasklist[0], tasklist[1]
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
	// Tried (C) but that lead to no improvement (although there
	// were tiny little dq that appeared to be in the right kind
	// of direction).
	//
	// So, now let's try (B): I figure at least the secondary will
	// be aware of the primary, and that combo works well in the
	// absence of joint limits.
	//
	// BTW yes, the tmpX scheme is not very pretty... and yes,
	// tmp3 is the same as tmp1 in the first trial. This is
	// prototyping code after all.
	//
	tmp3 = stack(tasklist, 0, 2);
	tl.push_back(&tmp3);
      }
      if (tasklist.size() > 3) {
	// At the time of writing, this is the same as tmp2 in the
	// first trial, so why not just reuse that? But beware of
	// future changes, this code is highly experimental at the
	// moment.
	tmp2 = stack(tasklist, 2, tasklist.size() - 2);
	tl.push_back(&tmp2);
      }
      else if (tasklist.size() == 3) {
	tl.push_back(tasklist[2]);
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
    //
    // Again, some open questions regarding what to stack together for
    // the second task that gets passed to compute_dq(). But it seems
    // natural to now push the (original) secondary into the "nice to
    // have" category and keep the (original) primary "please try as
    // hard as you can."
    
    tasklist_t tl;
    task_s tmp;
    tl.push_back(&t_lim);
    tl.push_back(tasklist[0]);
    if (tasklist.size() > 1) {
      if (tasklist.size() > 2) {
	tmp = stack(tasklist, 1, tasklist.size() - 1);
	tl.push_back(&tmp);
      }
      else {
	tl.push_back(tasklist[1]);
      }
    }
    dq = compute_dq(tl, dbgos, dbgpre);
    if (dbgos) {
      *dbgos << dbgpre << "joint limits take precedence over primary task\n";
    }
    return dq;
  }
  
}
