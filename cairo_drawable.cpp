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

#include "cairo_drawable.hpp"
#include "example_robot.hpp"
#include "base_waypoint.hpp"


namespace kinematic_elastic {
  
  
  void drawExampleRobot(ExampleRobot const & robot, cairo_t * cr, double weight, double pixelsize)
  {
    cairo_save(cr);
    
    // translucent disk for base
    cairo_set_source_rgba(cr, 0.7, 0.7, 0.7, 0.5);
    cairo_arc(cr, robot.position_[0], robot.position_[1], robot.radius_, 0., 2. * M_PI);
    cairo_fill(cr);
    
    // thick circle outline for base
    cairo_set_source_rgb(cr, 0.2, 0.2, 0.2);
    cairo_set_line_width(cr, weight * 3.0 / pixelsize);
    cairo_arc(cr, robot.position_[0], robot.position_[1], robot.radius_, 0., 2. * M_PI);
    cairo_stroke(cr);
    
    // thick line for arms
    cairo_set_source_rgb(cr, 0.2, 0.2, 0.2);
    cairo_set_line_width(cr, weight * 3.0 / pixelsize);
    cairo_move_to(cr, robot.position_[0], robot.position_[1]);
    cairo_line_to(cr, robot.pos_a_[0], robot.pos_a_[1]);
    cairo_line_to(cr, robot.pos_b_[0], robot.pos_b_[1]);
    cairo_line_to(cr, robot.pos_c_[0], robot.pos_c_[1]);
    cairo_stroke(cr);
    
    cairo_restore(cr);
  }
  
  
  void drawBaseWaypoint(BaseWaypoint const & wpt, cairo_t * cr, double weight, double pixelsize)
  {
    drawExampleRobot(wpt.robot_, cr, weight, pixelsize);
    
    cairo_save(cr);
    
    // orientation task
    
    cairo_set_source_rgba(cr, 0.0, 1.0, 0.5, 0.3);
    cairo_set_line_width(cr, weight * 6.0 / pixelsize);
    static double const len(0.5);
    double const dx(len * cos(wpt.orient_ee_.goal_));
    double const dy(len * sin(wpt.orient_ee_.goal_));
    cairo_move_to(cr, wpt.robot_.pos_b_[0], wpt.robot_.pos_b_[1]);
    cairo_line_to(cr, wpt.robot_.pos_b_[0] + dx, wpt.robot_.pos_b_[1] + dy);
    cairo_stroke(cr);
    
    // joint limits
    
    if (wpt.joint_limits_.isActive()) {
      cairo_set_source_rgba(cr, 1.0, 0.2, 0.8, 0.8);
      cairo_set_line_width(cr, weight * 1.0 / pixelsize);
      for (ssize_t ii(0); ii < wpt.joint_limits_.Jacobian_.rows(); ++ii) {
	if (0.0 < wpt.joint_limits_.Jacobian_(ii, 3)) {
	  cairo_move_to(cr, wpt.robot_.pos_a_[0], wpt.robot_.pos_a_[1]);
	  cairo_arc(cr, wpt.robot_.pos_a_[0], wpt.robot_.pos_a_[1], 0.1,
		    normangle(normangle(wpt.robot_.position_[2]) + wpt.joint_limits_.limits_(3, 0)),
		    normangle(normangle(wpt.robot_.position_[2]) + wpt.joint_limits_.limits_(3, 3)));
	  cairo_line_to(cr, wpt.robot_.pos_a_[0], wpt.robot_.pos_a_[1]);
	  cairo_fill(cr);
	}
	if (0.0 < wpt.joint_limits_.Jacobian_(ii, 4)) {
	  cairo_move_to(cr, wpt.robot_.pos_b_[0], wpt.robot_.pos_b_[1]);
	  cairo_arc(cr, wpt.robot_.pos_b_[0], wpt.robot_.pos_b_[1], 0.1,
		    normangle(normangle(wpt.robot_.q23_) + wpt.joint_limits_.limits_(4, 0)),
		    normangle(normangle(wpt.robot_.q23_) + wpt.joint_limits_.limits_(4, 3)));
	  cairo_line_to(cr, wpt.robot_.pos_b_[0], wpt.robot_.pos_b_[1]);
	  cairo_fill(cr);
	}
      }
    }
    
    // avoidance points
    
    cairo_set_source_rgb(cr, 1.0, 0.4, 1.0);
    cairo_set_line_width(cr, weight * 5.0 / pixelsize);
    
    if (wpt.avoid_base_.isActive()) {
      cairo_move_to(cr, wpt.avoid_base_.gpoint_[0], wpt.avoid_base_.gpoint_[1]);
      cairo_line_to(cr, wpt.avoid_base_.gpoint_[0], wpt.avoid_base_.gpoint_[1]);
      cairo_stroke(cr);
    }
    if (wpt.avoid_ellbow_.isActive()) {
      cairo_move_to(cr, wpt.avoid_ellbow_.gpoint_[0], wpt.avoid_ellbow_.gpoint_[1]);
      cairo_line_to(cr, wpt.avoid_ellbow_.gpoint_[0], wpt.avoid_ellbow_.gpoint_[1]);
      cairo_stroke(cr);
    }
    if (wpt.avoid_wrist_.isActive()) {
      cairo_move_to(cr, wpt.avoid_wrist_.gpoint_[0], wpt.avoid_wrist_.gpoint_[1]);
      cairo_line_to(cr, wpt.avoid_wrist_.gpoint_[0], wpt.avoid_wrist_.gpoint_[1]);
      cairo_stroke(cr);
    }
    if (wpt.avoid_ee_.isActive()) {
      cairo_move_to(cr, wpt.avoid_ee_.gpoint_[0], wpt.avoid_ee_.gpoint_[1]);
      cairo_line_to(cr, wpt.avoid_ee_.gpoint_[0], wpt.avoid_ee_.gpoint_[1]);
      cairo_stroke(cr);
    }
    
    // repulsion vectors
    
    cairo_set_source_rgb(cr, 0.4, 0.4, 1.0);
    cairo_set_line_width(cr, weight * 1.0 / pixelsize);
    if (wpt.repulse_base_.isActive()) {
	cairo_move_to(cr, wpt.repulse_base_.gpoint_[0], wpt.repulse_base_.gpoint_[1]);
	cairo_line_to(cr,
		      wpt.repulse_base_.gpoint_[0] + wpt.repulse_base_.delta_[0] / wpt.repulse_base_.gain_,
		      wpt.repulse_base_.gpoint_[1] + wpt.repulse_base_.delta_[1] / wpt.repulse_base_.gain_);
	cairo_stroke(cr);
    }
    if (wpt.repulse_ellbow_.isActive()) {
      cairo_move_to(cr, wpt.repulse_ellbow_.gpoint_[0], wpt.repulse_ellbow_.gpoint_[1]);
      cairo_line_to(cr,
		    wpt.repulse_ellbow_.gpoint_[0] + wpt.repulse_ellbow_.delta_[0] / wpt.repulse_ellbow_.gain_,
		    wpt.repulse_ellbow_.gpoint_[1] + wpt.repulse_ellbow_.delta_[1] / wpt.repulse_ellbow_.gain_);
      cairo_stroke(cr);
    }
    if (wpt.repulse_wrist_.isActive()) {
      cairo_move_to(cr, wpt.repulse_wrist_.gpoint_[0], wpt.repulse_wrist_.gpoint_[1]);
      cairo_line_to(cr,
		    wpt.repulse_wrist_.gpoint_[0] + wpt.repulse_wrist_.delta_[0] / wpt.repulse_wrist_.gain_,
		    wpt.repulse_wrist_.gpoint_[1] + wpt.repulse_wrist_.delta_[1] / wpt.repulse_wrist_.gain_);
      cairo_stroke(cr);
    }
    if (wpt.repulse_ee_.isActive()) {
      cairo_move_to(cr, wpt.repulse_ee_.gpoint_[0], wpt.repulse_ee_.gpoint_[1]);
      cairo_line_to(cr,
		    wpt.repulse_ee_.gpoint_[0] + wpt.repulse_ee_.delta_[0] / wpt.repulse_ee_.gain_,
		    wpt.repulse_ee_.gpoint_[1] + wpt.repulse_ee_.delta_[1] / wpt.repulse_ee_.gain_);
      cairo_stroke(cr);
    }
    
    cairo_restore(cr);
  }
  
}
