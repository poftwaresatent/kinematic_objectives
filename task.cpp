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

#include "task.hpp"
#include <iostream>


namespace kinematic_elastic {
  
  
  task_s::
  task_s()
    : b_max(numeric_limits<double>::max()),
      ndim(0)
  {
  }
  
  
  task_s::
  task_s(size_t ndof, size_t ndim_, double b_max_)
    : b_max(b_max_),
      ndim(ndim_)
  {
    xcur = Vector::Zero(ndim_);
    xdes = Vector::Zero(ndim_);
    dx = Vector::Zero(ndim_);
    Jx = Matrix::Zero(ndim_, ndof);
  }
  
  
  task_s stack(task_s const & t1, task_s const & t2)
  {
    task_s tt;
    stack(t1.xcur, t2.xcur, tt.xcur);
    stack(t1.xdes, t2.xdes, tt.xdes);
    stack(t1.dx, t2.dx, tt.dx);
    stack(t1.Jx, t2.Jx, tt.Jx);
    if (t1.b_max < t2.b_max) {
      tt.b_max = t1.b_max;
    }
    else {
      tt.b_max = t2.b_max;
    }
    tt.ndim = t1.ndim + t2.ndim;
    return tt;
  }
  
  
  task_s stack(tasklist_t const & tl, size_t first, size_t count)
  {
    count += first;
    
    size_t ndim(0);
    double b_max(numeric_limits<double>::max());
    for (size_t ii(first); ii < count; ++ii) {
      ndim += tl[ii]->ndim;
      if (b_max > tl[ii]->b_max) {
	b_max = tl[ii]->b_max;
      }
    }
    size_t const ndof(tl[first]->Jx.cols());
    task_s tt(ndof, ndim, b_max);
    
    // Beware of obscure update expression... it has two effects.
    for (size_t it(first), ir(0); it < count; ir += tl[it++]->ndim) {
      tt.xcur.block(ir, 0, tl[it]->ndim, 1) = tl[it]->xcur;
      tt.xdes.block(ir, 0, tl[it]->ndim, 1) = tl[it]->xdes;
      tt.dx.block(ir, 0, tl[it]->ndim, 1) = tl[it]->dx;
      tt.Jx.block(ir, 0, tl[it]->ndim, ndof) = tl[it]->Jx;
    }
    
    return tt;
  }
  
  
  void dump (Vector const & state,
	     tasklist_t const & tasklist,
	     Vector const & dq)
  {
    for (ssize_t ii(0); ii < state.size(); ++ii) {
      cout << state[ii] << "  ";
    }
    for (size_t ii(0); ii < tasklist.size(); ++ii) {
      cout << "  ";
      for (size_t jj(0); jj < tasklist[ii]->ndim; ++jj) {
	cout << tasklist[ii]->xcur[jj] << "  ";
      }
    }
    cout << "  ";
    for (ssize_t ii(0); ii < dq.size(); ++ii) {
      cout << dq[ii] << "  ";
    }
    cout << "\n";
  }
  
  
  void dbg (Vector const & state,
	    tasklist_t const & tasklist,
	    Vector const & dq)
  {
    cout << "==================================================\n"
	 << "state:";
    for (ssize_t ii(0); ii < state.size(); ++ii) {
      cout << "\t" << state[ii];
    }
    for (size_t ii(0); ii < tasklist.size(); ++ii) {
      cout << "\ntask " << ii << "\n";
      cout << "  current:";
      for (size_t jj(0); jj < tasklist[ii]->ndim; ++jj) {
	cout << "\t" << tasklist[ii]->xcur[jj];
      }
      cout << "\n  desired:";
      for (size_t jj(0); jj < tasklist[ii]->ndim; ++jj) {
	cout << "\t" << tasklist[ii]->xdes[jj];
      }
      cout << "\n  delta:";
      for (size_t jj(0); jj < tasklist[ii]->ndim; ++jj) {
	cout << "\t" << tasklist[ii]->dx[jj];
      }
      cout << "\n  Jacobian:";	// hardcoded for 1xN matrices
      for (ssize_t jj(0); jj < state.size(); ++jj) {
	cout << "\t" << tasklist[ii]->Jx(0, jj);
      }
    }
    cout << "\ndelta_q:";
    for (ssize_t ii(0); ii < dq.size(); ++ii) {
      cout << "\t" << dq[ii];
    }
    cout << "\n";
  }

}
