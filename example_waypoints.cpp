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

#include "example_waypoints.hpp"
#include "print.hpp"

#include <err.h>


namespace kinematic_elastic {

  namespace example {
    
    
    InteractionHandle::
    InteractionHandle(double radius, double red, double green, double blue, double alpha)
      : point_(3),
	radius_(radius),
	red_(red),
	green_(green),
	blue_(blue),
	alpha_(alpha)
    {
    }
    
    
    void InteractionHandle::
    draw(cairo_t * cr, double weight, double pixelsize)
      const
    {
      cairo_set_line_width(cr, weight * 1.0 / pixelsize);
      cairo_set_source_rgba(cr, red_, green_, blue_, alpha_);
      cairo_arc(cr, point_[0], point_[1], radius_, 0.0, 2.0 * M_PI);
      cairo_fill(cr);
    }
    
    
    BaseWaypoint::
    BaseWaypoint(InteractiveElastic const & elastic,
		 InteractionHandle const & repulsor,
		 double const & z_angle)
      : Waypoint(robot_),
	distance_api_(robot_, elastic),
	repulsor_(repulsor),
	z_angle_(z_angle),
	avoid_base_    (distance_api_, 0, 0.0),
	avoid_ellbow_  (distance_api_, 1, 0.0),
	avoid_wrist_   (distance_api_, 2, 0.0),
	avoid_ee_      (distance_api_, 3, 0.0),
	orient_ee_     (3, 100.0, 20.0),
	repulse_base_  (0,                 0.0, 0.0, 0.0, 100.0, repulsor.radius_),
	repulse_ellbow_(1,       robot_.len_a_, 0.0, 0.0, 100.0, repulsor.radius_),
	repulse_wrist_ (2,       robot_.len_b_, 0.0, 0.0, 100.0, repulsor.radius_),
	repulse_ee_    (3,       robot_.len_c_, 0.0, 0.0, 100.0, repulsor.radius_),
	joint_damping_ (10.0)
    {
      joint_limits_.init(5);
      joint_limits_.limits_(3, 0) = -120.0 * deg;
      joint_limits_.limits_(3, 1) = -119.999 * deg;
      joint_limits_.limits_(3, 2) =  119.999 * deg;
      joint_limits_.limits_(3, 3) =  120.0 * deg;
      joint_limits_.limits_(4, 0) = -120.0 * deg;
      joint_limits_.limits_(4, 1) = -119.999 * deg;
      joint_limits_.limits_(4, 2) =  119.999 * deg;
      joint_limits_.limits_(4, 3) =  120.0 * deg;
    
      constraints_.push_back(&joint_limits_);
      constraints_.push_back(&avoid_ee_);
      constraints_.push_back(&avoid_wrist_);
      constraints_.push_back(&avoid_ellbow_);
      constraints_.push_back(&avoid_base_);
    
      tasks_.push_back(&orient_ee_);
    
      objectives_.push_back(&repulse_base_);
      objectives_.push_back(&repulse_ellbow_);
      objectives_.push_back(&repulse_wrist_);
      objectives_.push_back(&repulse_ee_);
      objectives_.push_back(&joint_damping_);
    }
  
  
    void BaseWaypoint::
    draw(cairo_t * cr, double weight, double pixelsize)
      const
    {
      robot_.draw(cr, weight, pixelsize);
      
      cairo_save(cr);
    
      // orientation task
    
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
	for (ssize_t ii(0); ii < joint_limits_.Jacobian_.rows(); ++ii) {
	  if (0.0 < joint_limits_.Jacobian_(ii, 3)) {
	    cairo_move_to(cr, robot_.pos_a_[0], robot_.pos_a_[1]);
	    cairo_arc(cr, robot_.pos_a_[0], robot_.pos_a_[1], 0.1,
		      normangle(normangle(robot_.position_[2]) + joint_limits_.limits_(3, 0)),
		      normangle(normangle(robot_.position_[2]) + joint_limits_.limits_(3, 3)));
	    cairo_line_to(cr, robot_.pos_a_[0], robot_.pos_a_[1]);
	    cairo_fill(cr);
	  }
	  if (0.0 < joint_limits_.Jacobian_(ii, 4)) {
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
		      repulse_base_.gpoint_[0] + repulse_base_.delta_[0] / repulse_base_.gain_,
		      repulse_base_.gpoint_[1] + repulse_base_.delta_[1] / repulse_base_.gain_);
	cairo_stroke(cr);
      }
      if (repulse_ellbow_.isActive()) {
	cairo_move_to(cr, repulse_ellbow_.gpoint_[0], repulse_ellbow_.gpoint_[1]);
	cairo_line_to(cr,
		      repulse_ellbow_.gpoint_[0] + repulse_ellbow_.delta_[0] / repulse_ellbow_.gain_,
		      repulse_ellbow_.gpoint_[1] + repulse_ellbow_.delta_[1] / repulse_ellbow_.gain_);
	cairo_stroke(cr);
      }
      if (repulse_wrist_.isActive()) {
	cairo_move_to(cr, repulse_wrist_.gpoint_[0], repulse_wrist_.gpoint_[1]);
	cairo_line_to(cr,
		      repulse_wrist_.gpoint_[0] + repulse_wrist_.delta_[0] / repulse_wrist_.gain_,
		      repulse_wrist_.gpoint_[1] + repulse_wrist_.delta_[1] / repulse_wrist_.gain_);
	cairo_stroke(cr);
      }
      if (repulse_ee_.isActive()) {
	cairo_move_to(cr, repulse_ee_.gpoint_[0], repulse_ee_.gpoint_[1]);
	cairo_line_to(cr,
		      repulse_ee_.gpoint_[0] + repulse_ee_.delta_[0] / repulse_ee_.gain_,
		      repulse_ee_.gpoint_[1] + repulse_ee_.delta_[1] / repulse_ee_.gain_);
	cairo_stroke(cr);
      }
    
      cairo_restore(cr);
    }
  
  
    void BaseWaypoint::
    preUpdateHook()
    {
      orient_ee_.goal_ = z_angle_;
    
      repulse_base_.repulsor_ = repulsor_.point_;
      repulse_ellbow_.repulsor_ = repulsor_.point_;
      repulse_wrist_.repulsor_ = repulsor_.point_;
      repulse_ee_.repulsor_ = repulsor_.point_;
    }
    
    
    NormalWaypoint::
    NormalWaypoint(InteractiveElastic const & elastic,
		   InteractionHandle const & repulsor,
		   double const & z_angle)
      : BaseWaypoint(elastic, repulsor, z_angle)
    {
    }
    
    
    NormalWaypoint::
    ~NormalWaypoint()
    {
      for (size_t ii(0); ii < attract_prev_.size(); ++ii) {
	delete attract_prev_[ii];
      }
      for (size_t ii(0); ii < attract_next_.size(); ++ii) {
	delete attract_next_[ii];
      }
    }
    
    
    void NormalWaypoint::
    init(Vector const & position, Vector const & velocity)
    {
      if (attract_prev_.empty()) {
	errx(EXIT_FAILURE, "please call NormalWaypoint::setNeighbors exactly once on every waypoint");
      }
      BaseWaypoint::init(position, velocity);
    }
    
    
    // virtual void draw(cairo_t * cr, double weight, double pixelsize)
    // {
    //   BaseWaypoint::draw(cr, pixelsize);
    
    //   cairo_set_source_rgb(cr, 0.4, 1.0, 0.4);
    //   cairo_set_line_width(cr, weight * 1.0 / pixelsize);
    
    //   for (size_t ii(0); ii < attract_prev_.size(); ++ii) {
    //     if (attract_prev_[ii]->isActive()) {
    // 	cairo_move_to(cr, attract_prev_[ii]->gpoint_[0], attract_prev_[ii]->gpoint_[1]);
    // 	cairo_line_to(cr, attract_prev_[ii]->gpoint_[0] + attract_prev_[ii]->delta_[0] / attract_prev_[ii]->gain_, attract_prev_[ii]->gpoint_[1] + attract_prev_[ii]->delta_[1] / attract_prev_[ii]->gain_);
    // 	cairo_stroke(cr);
    //     }
    //   }
    //
    //   for (size_t ii(0); ii < attract_next_.size(); ++ii) {
    //     if (attract_next_[ii]->isActive()) {
    // 	cairo_move_to(cr, attract_next_[ii]->gpoint_[0], attract_next_[ii]->gpoint_[1]);
    // 	cairo_line_to(cr, attract_next_[ii]->gpoint_[0] + attract_next_[ii]->delta_[0] / attract_next_[ii]->gain_, attract_next_[ii]->gpoint_[1] + attract_next_[ii]->delta_[1] / attract_next_[ii]->gain_);
    // 	cairo_stroke(cr);
    //     }
    //   }
    // }
    
    
    void NormalWaypoint::
    preUpdateHook()
    {
      for (size_t ii(0); ii < attract_prev_.size(); ++ii) {
	attract_prev_[ii]->attractor_
	  = prev_->robot_.frame(attract_prev_[ii]->node_)
	  * attract_prev_[ii]->point_.homogeneous();
      }
      
      for (size_t ii(0); ii < attract_next_.size(); ++ii) {
	attract_next_[ii]->attractor_
	  = next_->robot_.frame(attract_next_[ii]->node_)
	  * attract_next_[ii]->point_.homogeneous();
      }
      
      BaseWaypoint::preUpdateHook();
    }
    
    
    void NormalWaypoint::
    setNeighbors(BaseWaypoint const * prev,
		 BaseWaypoint const * next)
    {
      if ( ! attract_prev_.empty()) {
	errx(EXIT_FAILURE, "please do not call NormalWaypoint::setNeighbors multiple times");
      }
      
      prev_ = prev;
      next_ = next;
      
      PointAttraction * pa;
      
      pa = new PointAttraction(0,           0.0, 0.0, 0.0, 500.0, -10.0);
      attract_prev_.push_back(pa);
      objectives_.push_back(pa);
      pa = new PointAttraction(1, robot_.len_a_, 0.0, 0.0, 500.0, -10.0);
      attract_prev_.push_back(pa);
      objectives_.push_back(pa);
      pa = new PointAttraction(2, robot_.len_b_, 0.0, 0.0, 500.0, -10.0);
      attract_prev_.push_back(pa);
      objectives_.push_back(pa);
      pa = new PointAttraction(3, robot_.len_c_, 0.0, 0.0, 500.0, -10.0);
      attract_prev_.push_back(pa);
      objectives_.push_back(pa);
      
      pa = new PointAttraction(0,           0.0, 0.0, 0.0, 500.0, -10.0);
      attract_next_.push_back(pa);
      objectives_.push_back(pa);
      pa = new PointAttraction(1, robot_.len_a_, 0.0, 0.0, 500.0, -10.0);
      attract_next_.push_back(pa);
      objectives_.push_back(pa);
      pa = new PointAttraction(2, robot_.len_b_, 0.0, 0.0, 500.0, -10.0);
      attract_next_.push_back(pa);
      objectives_.push_back(pa);
      pa = new PointAttraction(3, robot_.len_c_, 0.0, 0.0, 500.0, -10.0);
      attract_next_.push_back(pa);
      objectives_.push_back(pa);
    }
    
    
    BoundaryWaypoint::
    BoundaryWaypoint(InteractiveElastic const & elastic,
		     InteractionHandle const & repulsor,
		     double const & z_angle,
		     Vector const * eegoal,
		     Vector const * baseattractor)
      : BaseWaypoint(elastic, repulsor, z_angle),
	eetask_      (3, robot_.len_c_, 0.0, 0.0, 100.0, 20.0),
	attract_base_(0,           0.0, 0.0, 0.0, 100.0, 2.0),
	eegoal_(eegoal),
	baseattractor_(baseattractor)
    {
      tasks_.push_back(&eetask_);
      objectives_.push_back(&attract_base_);
    }
    
    
    void BoundaryWaypoint::
    draw(cairo_t * cr, double weight, double pixelsize)
      const
    {
      BaseWaypoint::draw(cr, weight, pixelsize);
      
      // thin line for end effector task
      cairo_set_source_rgb(cr, 1.0, 0.4, 0.4);
      cairo_set_line_width(cr, weight * 1.0 / pixelsize);
      cairo_move_to(cr, eetask_.gpoint_[0], eetask_.gpoint_[1]);
      cairo_line_to(cr, eetask_.goal_[0], eetask_.goal_[1]);
      cairo_stroke(cr);
      
      // base attraction
      if (attract_base_.isActive()) {
	cairo_set_source_rgb(cr, 0.4, 1.0, 0.4);
	cairo_set_line_width(cr, weight * 1.0 / pixelsize);
	cairo_move_to(cr, attract_base_.gpoint_[0], attract_base_.gpoint_[1]);
	cairo_line_to(cr, attract_base_.gpoint_[0] + attract_base_.delta_[0] / attract_base_.gain_, attract_base_.gpoint_[1] + attract_base_.delta_[1] / attract_base_.gain_);
	cairo_stroke(cr);
      }
    }
    
    
    void BoundaryWaypoint::
    preUpdateHook()
    {
      eetask_.goal_ = *eegoal_;
      attract_base_.attractor_ = *baseattractor_;
      BaseWaypoint::preUpdateHook();
    }
    
  }

}
