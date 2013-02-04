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

#include "position_control.hpp"
#include "model.hpp"


namespace kinematic_elastic {
  
  PositionControl::
  PositionControl(size_t node,
		  Vector const & point)
    : kp_(100.0),
      kd_(20.0),
      node_(node),
      point_(point)
  {
  }
  
  
  void PositionControl::
  init(Model const & model)
  {
    Eigen::Vector3d tmp1, tmp2;
    tmp1 << point_[0], point_[1], 0.0;
    tmp2 = model.frame(node_) * tmp1;
    gpoint_.resize(2);
    gpoint_ << tmp2[0], tmp2[1];
    goal_ = gpoint_;
    delta_ = Vector::Zero(point_.size());
    Matrix tmp3(model.computeJx(node_, gpoint_));
    Jacobian_ = tmp3.block(0, 0, 2, tmp3.cols());
  }
  
  
  void PositionControl::
  update(Model const & model)
  {
    Eigen::Vector3d tmp1, tmp2;
    tmp1 << point_[0], point_[1], 0.0;
    tmp2 = model.frame(node_) * tmp1;
    gpoint_ << tmp2[0], tmp2[1];
    Matrix tmp3(model.computeJx(node_, gpoint_));
    Jacobian_ = tmp3.block(0, 0, 2, tmp3.cols());
    delta_ = kp_ * (goal_ - gpoint_) - kd_ * Jacobian_ * model.getVelocity();
  }

}
