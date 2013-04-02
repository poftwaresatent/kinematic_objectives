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

#include <kinematic_objectives/point_attraction_objective.h>
#include <kinematic_objectives/kinematic_model.h>


namespace kinematic_objectives {
  
  
  PointAttractionObjective::
  PointAttractionObjective(size_t node,
		  double gain,
		  double distance)
  {
    construct(node, Vector::Zero(3), gain, distance);
  }
  
  
  PointAttractionObjective::
  PointAttractionObjective(size_t node,
		  double px,
		  double py,
		  double pz,
		  double gain,
		  double distance)
  {
    Vector silly(3);
    silly << px, py, pz;
    construct(node, silly, gain, distance);
  }
  
  
  void PointAttractionObjective::
  construct(size_t node,
		  Vector const & point,
		  double gain,
		  double distance)
  {
    gain_ = gain;
    distance_ = distance;
    node_ = node;
    point_ = point;
    attractor_.resize(0);
  }
  
  
  void PointAttractionObjective::
  init(KinematicModel const & model)
  {
    gpoint_.resize(point_.size());
    update(model);
  }
  
  
  void PointAttractionObjective::
  update(KinematicModel const & model)
  {
    if (0 == attractor_.size()) {
      jacobian_.resize(0, 0);
      return;
    }
    gpoint_ = model.getLinkFrame(node_) * point_.homogeneous();
    bias_ = attractor_ - gpoint_;
    double const dist(bias_.norm());
    if (dist < 1e-9) {
      jacobian_.resize(0, 0);
      return;
    }
    if (distance_ < 0.0) {
      // no saturation
      bias_ *= - gain_ / distance_;
    }
    else {
      // saturate at the given distance
      if (dist < distance_) {
	bias_ *= gain_ / distance_;
      }
      else {
	bias_ *= gain_ / dist;
      }
    }
    jacobian_ = model.getLinkJacobian(node_, gpoint_).block(0, 0, 3, model.getJointPosition().size());
  }
  
  
  bool PointAttractionObjective::
  isActive() const
  {
    return jacobian_.rows() > 0;
  }

}
