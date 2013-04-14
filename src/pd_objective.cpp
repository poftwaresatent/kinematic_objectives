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

#include <kinematic_objectives/pd_objective.h>
#include <kinematic_objectives/kinematic_model.h>

namespace kinematic_objectives {
  
  
  PDObjective::
  PDObjective(string const & name, Objective * src, bool own_src, double kp, double kd)
    : Objective(name),
      src_(src),
      own_src_(own_src),
      kp_(kp),
      kd_(kd)
  {
  }
  
  
  PDObjective::
  ~PDObjective()
  {
    if (own_src_) {
      delete src_;
    }
  }
  
  
  void PDObjective::
  init(KinematicModel const & model)
  {
    src_->init(model);
  }
  
  
  bool PDObjective::
  isActive() const
  {
    return src_->isActive();
  }
  
  
  void PDObjective::
  update(KinematicModel const & model)
  {
    src_->update(model);
    jacobian_ = src_->getJacobian();
    if (src_->isActive()) {
      bias_ = kp_ * src_->getBias() - kd_ * jacobian_ * model.getJointVelocity();
    }
  }
  
  
  double PDObjective::
  computeResidualErrorMagnitude(Vector const & ee) const
  {
    return src_->computeResidualErrorMagnitude(ee);
  }
  
}
