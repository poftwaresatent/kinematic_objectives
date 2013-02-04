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

#include "point_attraction.hpp"
#include "model.hpp"


namespace kinematic_elastic {
  
  
  PointAttraction::
  PointAttraction(size_t node,
		  Vector const & point)
    : gain_(100.0),
      distance_(2.0),
      node_(node),
      point_(point)
  {
    attractor_.resize(0);
  }
  
  
  void PointAttraction::
  init(Model const & model)
  {
    gpoint_.resize(point_.size());
    update(model);
  }
  
  
  void PointAttraction::
  update(Model const & model)
  {
    if (0 == attractor_.size()) {
      Jacobian_.resize(0, 0);
      return;
    }
    gpoint_ = model.frame(node_) * point_.homogeneous();
    delta_ = attractor_ - gpoint_;
    double const dist(delta_.norm());
    if (dist < 1e-9) {
      Jacobian_.resize(0, 0);
      return;
    }
    if (dist < distance_) {
      delta_ *= gain_ / distance_;
    }
    else {
      delta_ *= gain_ / dist;
    }
    Jacobian_ = model.computeJxo(node_, gpoint_).block(0, 0, 3, model.getPosition().size());
  }
  
  
  bool PointAttraction::
  isActive() const
  {
    return Jacobian_.rows() > 0;
  }

}