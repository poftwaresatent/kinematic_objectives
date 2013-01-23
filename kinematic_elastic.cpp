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

#include "baerlocher_algorithm.hpp"
#include "print.hpp"
#include <iostream>


namespace kinematic_elastic {
  
  
  task_s::
  task_s(size_t ndof, size_t ndim_, double b_max_)
    : b_max(b_max_),
      ndim(ndim_)
  {
    current = Vector::Zero(ndim_);
    desired = Vector::Zero(ndim_);
    Jacobian = Matrix::Zero(ndim_, ndof);
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
      for (size_t jj(0); jj < tasklist[ii].ndim; ++jj) {
	cout << tasklist[ii].current[jj] << "  ";
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
      for (size_t jj(0); jj < tasklist[ii].ndim; ++jj) {
	cout << "\t" << tasklist[ii].current[jj];
      }
      cout << "\n  desired:";
      for (size_t jj(0); jj < tasklist[ii].ndim; ++jj) {
	cout << "\t" << tasklist[ii].desired[jj];
      }
      cout << "\n  Jacobian:";	// hardcoded for 1xN matrices
      for (ssize_t jj(0); jj < state.size(); ++jj) {
	cout << "\t" << tasklist[ii].Jacobian(0, jj);
      }
    }
    cout << "\ndelta_q:";
    for (ssize_t ii(0); ii < dq.size(); ++ii) {
      cout << "\t" << dq[ii];
    }
    cout << "\n";
  }

}
