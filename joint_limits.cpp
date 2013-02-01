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

#include "joint_limits.hpp"
#include "task.hpp"
#include <limits>


namespace kinematic_elastic {
  
  
  void JointLimits::
  init(size_t ndof)
  {
    limits_.resize(ndof, 4);
    for (size_t ii(0); ii < ndof; ++ii) {
      size_t jj(0);
      for (; jj < 2; ++jj) {
	limits_(ii, jj) = -numeric_limits<double>::max();
      }
      for (; jj < 4; ++jj) {
	limits_(ii, jj) =  numeric_limits<double>::max();
      }
    }
  }
  
  
  bool JointLimits::
  check(Vector const & state)
    const
  {
    for (ssize_t ii(0); ii < state.size(); ++ii) {
      if (state[ii] < limits_(ii, 1)) {
	return false;
      }
      if (state[ii] > limits_(ii, 2)) {
	return false;
      }
    }
    return true;
  }
  
  
  void JointLimits::
  createTask(Vector const & state,
	     TaskData & jl,
	     vector<size_t> & locked)
    const
  {
    jl.step_hint_ = numeric_limits<double>::max();
    
    vector<double> delta;
    locked.clear();
    for (ssize_t ii(0); ii < state.size(); ++ii) {
      if (state[ii] < limits_(ii, 1)) {
	locked.push_back(ii);
	delta.push_back(limits_(ii, 0) - state[ii]);
      }
      else if (state[ii] > limits_(ii, 2)) {
	locked.push_back(ii);
	delta.push_back(limits_(ii, 3) - state[ii]);
      }
    }
    
    jl.delta_ = Vector::Map(&delta[0], delta.size());
    jl.Jacobian_ = Matrix::Zero(delta.size(), state.size());
    for (size_t ii(0); ii < delta.size(); ++ii) {
      jl.Jacobian_(ii, locked[ii]) = 1.0;
    }
  }
  
}
