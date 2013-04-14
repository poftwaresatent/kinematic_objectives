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

#include <kinematic_objectives/print.h>
#include <kinematic_objectives/util.h>
#include "interactive_compounds.h"

#include <err.h>


namespace kinematic_objectives {

  namespace demo {
    
    
    InteractiveCompound::
    InteractiveCompound(PlanarRobot & robot)
      : h_ee_      (0.2, 0.0, 1.0, 0.0, 0.5),
	h_ee_ori_  (0.1, 0.0, 1.0, 0.0, 0.3),
	h_base_    (0.2, 0.0, 1.0, 0.5, 0.5),
	h_repulsor_(1.5, 1.0, 0.5, 0.0, 0.2),
	h_obstacle_(1.5, 0.7, 0.0, 0.2, 0.5),
	robot_(robot),
	distance_api_(robot_, &h_obstacle_),
	joint_limits_  ("joint_limits"),
	avoid_base_    ("avoid_base",   distance_api_, 0, 0.0),
	avoid_ellbow_  ("avoid_ellbow", distance_api_, 1, 0.0),
	avoid_wrist_   ("avoid_wrist",  distance_api_, 2, 0.0),
	avoid_ee_      ("avoid_ee",     distance_api_, 3, 0.0),
	orient_ee_     ("orient_ee", 3, 1.0),
	repulse_base_  ("repulse_base",   0,             0.0, 0.0, 0.0, h_repulsor_.radius_),
	repulse_ellbow_("repulse_ellbow", 1,   robot_.len_a_, 0.0, 0.0, h_repulsor_.radius_),
	repulse_wrist_ ("repulse_wrist",  2,   robot_.len_b_, 0.0, 0.0, h_repulsor_.radius_),
	repulse_ee_    ("repulse_ee",     3,   robot_.len_c_, 0.0, 0.0, h_repulsor_.radius_)
    {
      handles_.push_back(&h_ee_);
      handles_.push_back(&h_ee_ori_);
      handles_.push_back(&h_base_);
      handles_.push_back(&h_repulsor_);
      handles_.push_back(&h_obstacle_);
      
      joint_limits_.init(5);
      joint_limits_.limits_(3, 0) = -120.0 * deg;
      joint_limits_.limits_(3, 1) = -119.999 * deg;
      joint_limits_.limits_(3, 2) =  119.999 * deg;
      joint_limits_.limits_(3, 3) =  120.0 * deg;
      joint_limits_.limits_(4, 0) = -120.0 * deg;
      joint_limits_.limits_(4, 1) = -119.999 * deg;
      joint_limits_.limits_(4, 2) =  119.999 * deg;
      joint_limits_.limits_(4, 3) =  120.0 * deg;
    
      compound_objective_.unilateral_constraints_.push_back(&joint_limits_);
      compound_objective_.unilateral_constraints_.push_back(&avoid_ee_);
      compound_objective_.unilateral_constraints_.push_back(&avoid_wrist_);
      compound_objective_.unilateral_constraints_.push_back(&avoid_ellbow_);
      compound_objective_.unilateral_constraints_.push_back(&avoid_base_);
    
      compound_objective_.hard_objectives_.push_back(&orient_ee_);
      
      compound_objective_.soft_objectives_.push_back(&repulse_base_);
      compound_objective_.soft_objectives_.push_back(&repulse_ellbow_);
      compound_objective_.soft_objectives_.push_back(&repulse_wrist_);
      compound_objective_.soft_objectives_.push_back(&repulse_ee_);
    }
  
  
    void InteractiveCompound::
    init(double gui_dimx, double gui_dimy)
    {
      h_ee_.point_         <<             1.0, gui_dimy / 2.0      ,     0.0;
      h_ee_ori_.point_     <<             2.0, gui_dimy / 2.0 + 1.0,     0.0;
      h_base_.point_       <<             1.0,                  1.0,     0.0;
      h_repulsor_.point_   <<  gui_dimx / 2.0,                  1.0,     0.0;
      h_obstacle_.point_   <<  gui_dimx / 2.0,       gui_dimy - 1.0,     0.0;
      
      Vector posture(5);
      posture <<
	gui_dimx / 2.0,
	gui_dimy / 2.0,
	80.0 * deg,
	- 40.0 * deg,
	25.0 * deg;
      robot_.update(posture, Vector::Zero(posture.size()));
      compound_objective_.init(robot_);
    }
    
    
    void InteractiveCompound::
    draw(cairo_t * cr, double weight, double pixelsize)
      const
    {
      cairo_save(cr);
      
      // handles
      
      for (size_t ii(0); ii < handles_.size(); ++ii) {
	handles_[ii]->draw(cr, weight, pixelsize);
      }
      
      cairo_set_source_rgba(cr, 0.0, 1.0, 0.0, 0.5);
      cairo_set_line_width(cr, weight * 1.0 / pixelsize);
      cairo_move_to(cr, h_ee_.point_[0], h_ee_.point_[1]);
      cairo_line_to(cr, h_ee_ori_.point_[0], h_ee_ori_.point_[1]);
      cairo_stroke(cr);
      
      // orientation objective
    
      cairo_set_source_rgba(cr, 0.0, 1.0, 0.5, 0.3);
      cairo_set_line_width(cr, weight * 6.0 / pixelsize);
      static double const len(0.5);
      double const dx(len * cos(orient_ee_.goal_));
      double const dy(len * sin(orient_ee_.goal_));
      cairo_move_to(cr, robot_.pos_b_[0], robot_.pos_b_[1]);
      cairo_line_to(cr, robot_.pos_b_[0] + dx, robot_.pos_b_[1] + dy);
      cairo_stroke(cr);
    
      // joint limits
    
      if (joint_limits_.isActive()) {
	cairo_set_source_rgba(cr, 1.0, 0.2, 0.8, 0.8);
	cairo_set_line_width(cr, weight * 1.0 / pixelsize);
	for (ssize_t ii(0); ii < joint_limits_.getJacobian().rows(); ++ii) {
	  if (0.0 < joint_limits_.getJacobian()(ii, 3)) {
	    cairo_move_to(cr, robot_.pos_a_[0], robot_.pos_a_[1]);
	    cairo_arc(cr, robot_.pos_a_[0], robot_.pos_a_[1], 0.1,
		      normangle(normangle(robot_.position_[2]) + joint_limits_.limits_(3, 0)),
		      normangle(normangle(robot_.position_[2]) + joint_limits_.limits_(3, 3)));
	    cairo_line_to(cr, robot_.pos_a_[0], robot_.pos_a_[1]);
	    cairo_fill(cr);
	  }
	  if (0.0 < joint_limits_.getJacobian()(ii, 4)) {
	    cairo_move_to(cr, robot_.pos_b_[0], robot_.pos_b_[1]);
	    cairo_arc(cr, robot_.pos_b_[0], robot_.pos_b_[1], 0.1,
		      normangle(normangle(robot_.q23_) + joint_limits_.limits_(4, 0)),
		      normangle(normangle(robot_.q23_) + joint_limits_.limits_(4, 3)));
	    cairo_line_to(cr, robot_.pos_b_[0], robot_.pos_b_[1]);
	    cairo_fill(cr);
	  }
	}
      }
    
      // avoidance points
    
      cairo_set_source_rgb(cr, 1.0, 0.4, 1.0);
      cairo_set_line_width(cr, weight * 5.0 / pixelsize);
    
      if (avoid_base_.isActive()) {
	cairo_move_to(cr, avoid_base_.gpoint_[0], avoid_base_.gpoint_[1]);
	cairo_line_to(cr, avoid_base_.gpoint_[0], avoid_base_.gpoint_[1]);
	cairo_stroke(cr);
      }
      if (avoid_ellbow_.isActive()) {
	cairo_move_to(cr, avoid_ellbow_.gpoint_[0], avoid_ellbow_.gpoint_[1]);
	cairo_line_to(cr, avoid_ellbow_.gpoint_[0], avoid_ellbow_.gpoint_[1]);
	cairo_stroke(cr);
      }
      if (avoid_wrist_.isActive()) {
	cairo_move_to(cr, avoid_wrist_.gpoint_[0], avoid_wrist_.gpoint_[1]);
	cairo_line_to(cr, avoid_wrist_.gpoint_[0], avoid_wrist_.gpoint_[1]);
	cairo_stroke(cr);
      }
      if (avoid_ee_.isActive()) {
	cairo_move_to(cr, avoid_ee_.gpoint_[0], avoid_ee_.gpoint_[1]);
	cairo_line_to(cr, avoid_ee_.gpoint_[0], avoid_ee_.gpoint_[1]);
	cairo_stroke(cr);
      }
    
      // repulsion vectors
    
      cairo_set_source_rgb(cr, 0.4, 0.4, 1.0);
      cairo_set_line_width(cr, weight * 1.0 / pixelsize);
      if (repulse_base_.isActive()) {
	cairo_move_to(cr, repulse_base_.gpoint_[0], repulse_base_.gpoint_[1]);
	cairo_line_to(cr,
		      repulse_base_.gpoint_[0] + repulse_base_.getBias()[0],
		      repulse_base_.gpoint_[1] + repulse_base_.getBias()[1]);
	cairo_stroke(cr);
      }
      if (repulse_ellbow_.isActive()) {
	cairo_move_to(cr, repulse_ellbow_.gpoint_[0], repulse_ellbow_.gpoint_[1]);
	cairo_line_to(cr,
		      repulse_ellbow_.gpoint_[0] + repulse_ellbow_.getBias()[0],
		      repulse_ellbow_.gpoint_[1] + repulse_ellbow_.getBias()[1]);
	cairo_stroke(cr);
      }
      if (repulse_wrist_.isActive()) {
	cairo_move_to(cr, repulse_wrist_.gpoint_[0], repulse_wrist_.gpoint_[1]);
	cairo_line_to(cr,
		      repulse_wrist_.gpoint_[0] + repulse_wrist_.getBias()[0],
		      repulse_wrist_.gpoint_[1] + repulse_wrist_.getBias()[1]);
	cairo_stroke(cr);
      }
      if (repulse_ee_.isActive()) {
	cairo_move_to(cr, repulse_ee_.gpoint_[0], repulse_ee_.gpoint_[1]);
	cairo_line_to(cr,
		      repulse_ee_.gpoint_[0] + repulse_ee_.getBias()[0],
		      repulse_ee_.gpoint_[1] + repulse_ee_.getBias()[1]);
	cairo_stroke(cr);
      }
    
      cairo_restore(cr);
    }
  
  
    void InteractiveCompound::
    update()
    {
      orient_ee_.goal_ = atan2(h_ee_ori_.point_[1] - h_ee_.point_[1], h_ee_ori_.point_[0] - h_ee_.point_[0]);
      
      repulse_base_.repulsor_ = h_repulsor_.point_;
      repulse_ellbow_.repulsor_ = h_repulsor_.point_;
      repulse_wrist_.repulsor_ = h_repulsor_.point_;
      repulse_ee_.repulsor_ = h_repulsor_.point_;
    }
    
    
    ElasticLinksCompound::
    ElasticLinksCompound(PlanarRobot & robot)
      : InteractiveCompound(robot),
	h2_ee_    (0.2, 0.0, 0.6, 0.0, 0.5),
	h1_wrist_ (0.2, 0.0, 0.5, 0.5, 0.5),
	h2_wrist_ (0.2, 0.0, 0.3, 0.3, 0.5),
	h1_ellbow_(0.2, 0.0, 0.8, 0.8, 0.5),
	h2_ellbow_(0.2, 0.0, 0.4, 0.4, 0.5),
	h2_base_  (0.2, 0.0, 0.6, 0.3, 0.5),
	ee_left_     ("ee_left",      3, robot_.len_c_, 0.0, 0.0, -10.0),
	ee_right_    ("ee_right",     3, robot_.len_c_, 0.0, 0.0, -10.0),
	wrist_left_  ("wrist_left",   2, robot_.len_b_, 0.0, 0.0, -10.0),
	wrist_right_ ("wrist_right",  2, robot_.len_b_, 0.0, 0.0, -10.0),
	ellbow_left_ ("ellbow_left",  1, robot_.len_a_, 0.0, 0.0, -10.0),
	ellbow_right_("ellbow_right", 1, robot_.len_a_, 0.0, 0.0, -10.0),
	base_left_   ("base_left",    0,           0.0, 0.0, 0.0, -10.0),
	base_right_  ("base_right",   0,           0.0, 0.0, 0.0, -10.0)
    {
      handles_.push_back(&h2_ee_);
      handles_.push_back(&h1_wrist_);
      handles_.push_back(&h2_wrist_);
      handles_.push_back(&h1_ellbow_);
      handles_.push_back(&h2_ellbow_);
      handles_.push_back(&h2_base_);

      compound_objective_.soft_objectives_.push_back(&ee_left_);
      compound_objective_.soft_objectives_.push_back(&ee_right_);
      compound_objective_.soft_objectives_.push_back(&wrist_left_);
      compound_objective_.soft_objectives_.push_back(&wrist_right_);
      compound_objective_.soft_objectives_.push_back(&ellbow_left_);
      compound_objective_.soft_objectives_.push_back(&ellbow_right_);
      compound_objective_.soft_objectives_.push_back(&base_left_);
      compound_objective_.soft_objectives_.push_back(&base_right_);
    }
    
    
    void ElasticLinksCompound::
    init(double gui_dimx, double gui_dimy)
    {
      h2_ee_.point_     <<             3.0, gui_dimy / 2.0      ,     0.0;
      h1_wrist_.point_  <<             1.0, gui_dimy / 2.0 - 0.5,     0.0;
      h2_wrist_.point_  <<             3.0, gui_dimy / 2.0 - 0.5,     0.0;
      h1_ellbow_.point_ <<             1.0, gui_dimy / 2.0 - 1.0,     0.0;
      h2_ellbow_.point_ <<             3.0, gui_dimy / 2.0 - 1.0,     0.0;
      h2_base_.point_   <<             3.0,                  1.0,     0.0;
      
      InteractiveCompound::init(gui_dimx, gui_dimy);
    }
    
    
    static void draw_elastic(cairo_t * cr, double weight, double pixelsize,
			     InteractionHandle const & handle,
			     PointAttractionObjective const & objective)
    {
      cairo_set_source_rgba(cr, handle.red_, handle.green_, handle.blue_, 0.5);
      cairo_set_line_width(cr, weight * 1.0 / pixelsize);
      cairo_move_to(cr,
		    objective.gpoint_[0],
		    objective.gpoint_[1]);
      cairo_line_to(cr,
		    handle.point_[0],
		    handle.point_[1]);
      cairo_stroke(cr);
      
      if (objective.isActive()) {
	cairo_set_source_rgba(cr, handle.red_, handle.green_, handle.blue_, 0.8);
	cairo_set_line_width(cr, weight * 3.0 / pixelsize);
	cairo_move_to(cr,
		      objective.gpoint_[0],
		      objective.gpoint_[1]);
	cairo_line_to(cr,
		      objective.gpoint_[0] + objective.getBias()[0],
		      objective.gpoint_[1] + objective.getBias()[1]);
	cairo_stroke(cr);
      }
    }
    
    
    void ElasticLinksCompound::
    draw(cairo_t * cr, double weight, double pixelsize)
      const
    {
      InteractiveCompound::draw(cr, weight, pixelsize);
      
      draw_elastic(cr, weight, pixelsize,  h_ee_,     ee_left_);
      draw_elastic(cr, weight, pixelsize, h2_ee_,     ee_right_);
      draw_elastic(cr, weight, pixelsize, h1_wrist_,  wrist_left_);
      draw_elastic(cr, weight, pixelsize, h2_wrist_,  wrist_right_);
      draw_elastic(cr, weight, pixelsize, h1_ellbow_, ellbow_left_);
      draw_elastic(cr, weight, pixelsize, h2_ellbow_, ellbow_right_);
      draw_elastic(cr, weight, pixelsize,  h_base_,   base_left_);
      draw_elastic(cr, weight, pixelsize, h2_base_,   base_right_);
    }
    
    
    void ElasticLinksCompound::
    update()
    {
      ee_left_.attractor_ =       h_ee_.point_;
      ee_right_.attractor_ =     h2_ee_.point_;
      wrist_left_.attractor_ =   h1_wrist_.point_;
      wrist_right_.attractor_ =  h2_wrist_.point_;
      ellbow_left_.attractor_ =  h1_ellbow_.point_;
      ellbow_right_.attractor_ = h2_ellbow_.point_;
      base_left_.attractor_ =     h_base_.point_;
      base_right_.attractor_ =   h2_base_.point_;
      
      InteractiveCompound::update();
    }
    
    
    EEGoalCompound::
    EEGoalCompound(PlanarRobot & robot)
      : InteractiveCompound(robot),
	attract_ee_       ("attract_ee",   3, robot_.len_c_, 0.0, 0.0, -1.0),
	attract_base_     ("attract_base", 0,           0.0, 0.0, 0.0, -1.0)
    {
      compound_objective_.hard_objectives_.push_back(&attract_ee_);
      compound_objective_.soft_objectives_.push_back(&attract_base_);
    }
    
    
    void EEGoalCompound::
    draw(cairo_t * cr, double weight, double pixelsize)
      const
    {
      InteractiveCompound::draw(cr, weight, pixelsize);
      
      // thin line for end effector objective
      cairo_set_source_rgb(cr, 1.0, 0.4, 0.4);
      cairo_set_line_width(cr, weight * 1.0 / pixelsize);
      cairo_move_to(cr, attract_ee_.gpoint_[0], attract_ee_.gpoint_[1]);
      cairo_line_to(cr, attract_ee_.attractor_[0], attract_ee_.attractor_[1]);
      cairo_stroke(cr);
      
      // base attraction
      if (attract_base_.isActive()) {
	cairo_set_source_rgb(cr, 0.4, 1.0, 0.4);
	cairo_set_line_width(cr, weight * 1.0 / pixelsize);
	cairo_move_to(cr, attract_base_.gpoint_[0], attract_base_.gpoint_[1]);
	cairo_line_to(cr,
		      attract_base_.gpoint_[0] + attract_base_.getBias()[0],
		      attract_base_.gpoint_[1] + attract_base_.getBias()[1]);
	cairo_stroke(cr);
      }
    }
    
    
    void EEGoalCompound::
    update()
    {
      attract_ee_.attractor_ = h_ee_.point_;
      attract_base_.attractor_ = h_base_.point_;
      InteractiveCompound::update();
    }
    
  }

}
