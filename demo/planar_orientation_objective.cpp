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

#include <kinematic_objectives/kinematic_model.h>
#include "planar_orientation_objective.h"


namespace kinematic_objectives {
  
  namespace demo {
    
    
    PlanarOrientationObjective::
    PlanarOrientationObjective(string const & name,
			       size_t node,
			       double kp,
			       double kd)
      : Objective(name),
	kp_(kp),
	kd_(kd),
	node_(node)
    {
    }
  
  
    void PlanarOrientationObjective::
    init(KinematicModel const & model)
    {
      Vector const ex(model.getLinkFrame(node_).linear().block(0, 0, 3, 1));
      angle_ = atan2(ex[1], ex[0]);
      goal_ = angle_;
      bias_ = Vector::Zero(1);
      jacobian_ = model.getLinkJacobian(node_, Vector::Zero(3)).block(5, 0, 1, model.getJointPosition().size());
    }
  
  
    void PlanarOrientationObjective::
    update(KinematicModel const & model)
    {
      Vector const ex(model.getLinkFrame(node_).linear().block(0, 0, 3, 1));
      angle_ = atan2(ex[1], ex[0]);
      jacobian_ = model.getLinkJacobian(node_, Vector::Zero(3)).block(5, 0, 1, model.getJointPosition().size());
      bias_[0] = kp_ * (goal_ - angle_);
      bias_ -= kd_ * jacobian_ * model.getJointVelocity();
    }

  }

}
