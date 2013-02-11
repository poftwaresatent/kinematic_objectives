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

#ifndef KINEMATIC_ELASTIC_BASE_WAYPOINT_HPP
#define KINEMATIC_ELASTIC_BASE_WAYPOINT_HPP

#include "waypoint.hpp"

#include "joint_limit_constraint.hpp"
#include "point_mindist_constraint.hpp"
#include "position_control.hpp"
#include "point_attraction.hpp"
#include "point_repulsion.hpp"
#include "posture_damping.hpp"
#include "qh_ori_z_control.hpp"

#include "example_robot.hpp"	// rfct
#include <cairo/cairo.h>	// rfct


namespace kinematic_elastic {
  
  class BaseWaypoint
    : public Waypoint
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    BaseWaypoint(double qh_obstacle_radius,
		 double qh_repulsor_radius);
    
    virtual void draw(cairo_t * cr, double weight, double pixelsize);
    virtual void update(Vector const & qh_obstacle_point,
			Vector const & qh_repulsor_point,
			double qh_zangle);
    
    //// XXXX protected or so...
    
    double timestep_;
    ExampleRobot robot_; // XXXX keep this before any constraints so we can use its values for initializing them
    
    JointLimitConstraint joint_limits_;
    
    PointMindistConstraint avoid_base_;
    PointMindistConstraint avoid_ellbow_;
    PointMindistConstraint avoid_wrist_;
    PointMindistConstraint avoid_ee_;
    
    OriZControl orient_ee_;
    
    PointRepulsion repulse_base_;
    PointRepulsion repulse_ellbow_;
    PointRepulsion repulse_wrist_;
    PointRepulsion repulse_ee_;
    
    PostureDamping joint_damping_;
  };
  
}

#endif // KINEMATIC_ELASTIC_BASE_WAYPOINT_HPP
