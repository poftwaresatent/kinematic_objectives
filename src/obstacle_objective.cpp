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

#include <kinematic_objectives/obstacle_objective.h>
#include <kinematic_objectives/kinematic_model.h>
#include <kinematic_objectives/distance_model.h>


namespace kinematic_objectives {
  
  
  ObstacleObjective::
  ObstacleObjective(DistanceModel const & distance_model,
		    size_t node,
		    double mindist)
    : distance_model_(distance_model),
      node_(node),
      mindist_(mindist)
  {
  }
  
  
  void ObstacleObjective::
  init(KinematicModel const & model)
  {
    bias_ = Vector::Zero(1);
    jacobian_.resize(0, 0);
  }
  
  
  void ObstacleObjective::
  update(KinematicModel const & model)
  {
    double const dist(distance_model_.computeMinimumSeparation(node_, gpoint_, obstacle_));
    if (dist >= mindist_) {
      jacobian_.resize(0, 0);
      return;
    }
    Vector tmp(gpoint_ - obstacle_);
    tmp /= dist;
    jacobian_
      = tmp.transpose()
      * model.getLinkJacobian(node_, obstacle_).block(0, 0, 3, model.getJointPosition().size());
    bias_[0] = mindist_ - dist;
  }
  
  
  bool ObstacleObjective::
  isActive() const
  {
    return jacobian_.rows() > 0;
  }

}
