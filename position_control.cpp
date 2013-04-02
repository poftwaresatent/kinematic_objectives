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


namespace kinematic_objectives {
  
  LinkPositionObjective::
  LinkPositionObjective(size_t node,
		  double px,
		  double py,
		  double pz,
		  double kp,
		  double kd)
    : kp_(kp),
      kd_(kd),
      node_(node),
      point_(3)
  {
    point_ << px, py, pz;
  }
  
  
  void LinkPositionObjective::
  init(KinematicModel const & model)
  {
    gpoint_ = model.getLinkFrame(node_) * point_.homogeneous();
    goal_ = gpoint_;
    jacobian_ = model.getLinkJacobian(node_, gpoint_).block(0, 0, 3, model.getJointPosition().size());
    bias_ = Vector::Zero(point_.size());
  }
  
  
  void LinkPositionObjective::
  update(KinematicModel const & model)
  {
    gpoint_ = model.getLinkFrame(node_) * point_.homogeneous();
    jacobian_ = model.getLinkJacobian(node_, gpoint_).block(0, 0, 3, model.getJointPosition().size());
    bias_ = kp_ * (goal_ - gpoint_) - kd_ * jacobian_ * model.getJointVelocity();
  }

}
