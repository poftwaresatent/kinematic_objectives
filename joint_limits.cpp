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
#include "model.hpp"
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
    locked_joints_.clear();
    delta_.resize(0);
    Jacobian_.resize(0, 0);
    step_hint_ = numeric_limits<double>::max();
  }
  
  
  void JointLimits::
  update(Vector const & position, Vector const & velocity)
  {
    vector<double> delta;
    locked_joints_.clear();
    
    for (ssize_t ii(0); ii < position.size(); ++ii) {
      if (position[ii] < limits_(ii, 1)) {
	locked_joints_.push_back(ii);
	delta.push_back(limits_(ii, 0) - position[ii]);
      }
      else if (position[ii] > limits_(ii, 2)) {
	locked_joints_.push_back(ii);
	delta.push_back(limits_(ii, 3) - position[ii]);
      }
    }
    
    delta_ = Vector::Map(&delta[0], delta.size());
    Jacobian_ = Matrix::Zero(delta.size(), position.size());
    for (size_t ii(0); ii < delta.size(); ++ii) {
      Jacobian_(ii, locked_joints_[ii]) = 1.0;
    }
  }
  
  
  bool JointLimits::
  isActive()
    const
  {
    return ! locked_joints_.empty();
  }
  
}
