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
#include <deque>


namespace kinematic_elastic {
  
  
  typedef deque<task_s> taskdeque_t; // we need a local modifiable copy of the task list given by the user
  
  
  static inline void stack (Vector const & v1,
			    Vector const & v2,
			    Vector & vv)
  {
    vv.resize(v1.size() + v2.size());
    vv.block(0,         0, v1.size(), 1) = v1;
    vv.block(v1.size(), 0, v2.size(), 1) = v2;
  }
  
  
  static inline void stack (Matrix const & m1,
			    Matrix const & m2,
			    Matrix & mm)
  {
    mm.resize(m1.rows() + m2.rows(), m1.cols());
    mm.block(0, 0, m1.rows(), m1.cols()) = m1;
    mm.block(m1.rows(), 0, m2.rows(), m2.cols()) = m2; // this should abort automatically if m2.cols() != m1.cols()
  }
  
  
  static Vector compute_dq (size_t ndof,
			    taskdeque_t const & tasks,
			    //// sensible speedup: Matrix const * precomputed_Ja_inv,
			    ostream * dbgos,
			    char const * dbgpre)
  {
    Vector dq(Vector::Zero(ndof));
    
    if (tasks.empty()) {
      if (dbgos) {
	*dbgos << dbgpre << "DONE (no tasks):\n";
	print (dq, *dbgos, "dq", string(dbgpre) + "  ");
      }
      return dq;
    }
    
    Vector dxa(tasks[0].desired - tasks[0].current);
    Matrix Ja_inv;
    pseudo_inverse_nonsingular (tasks[0].Jacobian, Ja_inv);
    dq += Ja_inv * dxa;
    
    if (dbgos) {
      *dbgos << dbgpre << "primary task:\n";
      string pre (dbgpre);
      pre += "  ";
      print (dxa,               *dbgos, "dxa",    pre);
      print (tasks[0].Jacobian, *dbgos, "Ja",     pre);
      print (Ja_inv,            *dbgos, "Ja_inv", pre);
      print (dq,                *dbgos, "dq",     pre);
    }
    
    if (tasks.size() == 1) {
      if (dbgos) {
	*dbgos << dbgpre << "DONE (single task):\n";
      }
      return dq;
    }
    
    Vector dxb(tasks[0].ndim + tasks[1].ndim);
    dxb.block(0,              0, tasks[0].ndim, 1) = dxa;
    dxb.block(tasks[0].ndim, 0, tasks[1].ndim, 1) = tasks[1].desired - tasks[1].current;
    Matrix Jb;
    stack(tasks[0].Jacobian, tasks[1].Jacobian, Jb);
    Matrix Jb_inv;
    double tmp(1.0 / (tasks[0].b_max < tasks[1].b_max ? tasks[0].b_max : tasks[1].b_max));
    pseudo_inverse_damped (Jb, tmp, Jb_inv);
    Matrix Na(Matrix::Identity(ndof, ndof) - Ja_inv * tasks[0].Jacobian);
    dq +=  Na * Jb_inv * dxb;
    
    if (dbgos) {
      *dbgos << dbgpre << "secondary task:\n";
      string pre (dbgpre);
      pre += "  ";
      print (Na,     *dbgos, "Na",     pre);
      print (dxb,    *dbgos, "dxb",    pre);
      print (Jb,     *dbgos, "Jb",     pre);
      print (Jb_inv, *dbgos, "Jb_inv", pre);
      print (dq,     *dbgos, "dq",     pre);
    }
    
    if (tasks.size() == 2) {
      if (dbgos) {
	*dbgos << dbgpre << "DONE (two tasks):\n";
      }
      return dq;
    }
    
    Matrix Nb(Matrix::Identity(ndof, ndof) - Jb_inv * Jb);
    size_t nrest(0);
    for (size_t it(2); it < tasks.size(); ++it) {
      nrest += tasks[it].ndim;
    }
    Vector dxc(nrest);
    Matrix Jc(nrest, ndof);
    tmp = numeric_limits<double>::max();
    for (size_t it(2), ir(0); it < tasks.size(); ir += tasks[it++].ndim) { // obscure update expression, grr
      dxc.block(ir, 0, tasks[it].ndim, 1) = tasks[it].desired - tasks[it].current;
      Jc.block( ir, 0, tasks[it].ndim, ndof) = tasks[it].Jacobian;
      if (tasks[it].b_max < tmp) {
	tmp = tasks[it].b_max; // beware: tmp is 1.0/lambda here (whereas it was lambda above)
      }
    }
    Matrix Jc_inv;
    pseudo_inverse_damped (Jc, 1.0 / tmp, Jc_inv);
    dq += Na * Nb * Jc_inv * dxc;
    
    if (dbgos) {
      *dbgos << dbgpre << "remainder:\n";
      string pre (dbgpre);
      pre += "  ";
      print (Nb,     *dbgos, "Nb",     pre);
      print (dxc,    *dbgos, "dxc",    pre);
      print (Jc,     *dbgos, "Jc",     pre);
      print (Jc_inv, *dbgos, "Jc_inv", pre);
      print (dq,     *dbgos, "dq",     pre);
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
    // XXXX Make this copy lazy. But that requires templatizing
    // compute_dq on the container type. Do that later, once
    // everything works.
    taskdeque_t tasks(tasklist.size());
    for (size_t ii(0); ii < tasklist.size(); ++ii) {
      tasks[ii] = tasklist[ii];
    }
    
    size_t const ndof(state.size());
    Vector dq(compute_dq(ndof, tasks, dbgos, dbgpre));
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
    if (tasks.empty()) {
      tasks.push_back(t_lim);	// XXXX this ends up copy-constructing I guess
      dq = compute_dq(ndof, tasks, dbgos, dbgpre);
      if (dbgos) {
	*dbgos << dbgpre << "joint limits adjusted, but empty task list\n";
      }
      return dq;
    }
    
    // Try stacking the joint limits on top of the primary task. That
    // means knocking out the locked joints from the primary task's
    // Jacobian, too.
    //
    // XXXX There may be some caching opportunities here, at the very
    // least all deltas except for the first set would remain
    // untouched and should just be reused instead of
    // recomputed.
    
    if (dbgos) {
      string pre (dbgpre);
      pre += "  ";
      *dbgos << dbgpre << "knock locked joints out from primary task...\n";
      print (t_lim.Jacobian, *dbgos, "J_lim", pre);
    }
    
    Matrix J_try(tasks[0].Jacobian);
    for (size_t ii(0); ii < locked.size(); ++ii) {
      J_try.block(0, locked[ii], J_try.rows(), 1) = Vector::Zero(J_try.rows());
      if (dbgos) {
	*dbgos << dbgpre << "  joint " << locked[ii] << " is locked!\n";
      }
    }
    Eigen::FullPivLU<Matrix> plu(J_try * J_try.transpose());
    if (plu.isInvertible()) {
      if (dbgos) {
	*dbgos << dbgpre << "great: primary task should still be achievable even with locked joints\n";
      }
      stack(t_lim.Jacobian, J_try, tasks[0].Jacobian);
      tasks[0].ndim += t_lim.ndim;
      Vector tmp;
      stack(t_lim.current, tasks[0].current, tmp);
      tasks[0].current = tmp;
      stack(t_lim.desired, tasks[0].desired, tmp);
      tasks[0].desired = tmp;
      // We could also knock the locked joints from the other
      // tasks. Maybe we should, but the null-space projection
      // "should" take care of it anyway from tasks[1] onward.
      dq = compute_dq(ndof, tasks, dbgos, dbgpre);
      if (dbgos) {
	*dbgos << dbgpre << "joint limits merged into primary task\n";
      }
      return dq;
    }
    
    // Well, that didn't work, which means that the primary task
    // cannot be achieved due to joint limits. The best we can do is
    // prepend the joint limits to the task list.
    
    tasks.push_front(t_lim);	// XXXX another copy-construction?
    dq = compute_dq(ndof, tasks, dbgos, dbgpre);
    if (dbgos) {
      *dbgos << dbgpre << "joint limits take precedence over primary task\n";
    }
    return dq;
  }
  
}
