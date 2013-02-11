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

#include "algorithm.hpp"	// rfct
#include "pseudo_inverse.hpp"	// rfct

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
    BaseWaypoint(InteractionHandle const & obstacle,
		 InteractionHandle const & repulsor,
		 double const & z_angle)
      : Waypoint(robot_),
	timestep_(1e-2),
	obstacle_(obstacle),
	repulsor_(repulsor),
	z_angle_(z_angle),
	avoid_base_    (0,                 0.0, 0.0, 0.0, obstacle.radius_ + robot_.radius_),
	avoid_ellbow_  (1,       robot_.len_a_, 0.0, 0.0, obstacle.radius_),
	avoid_wrist_   (2,       robot_.len_b_, 0.0, 0.0, obstacle.radius_),
	avoid_ee_      (3, robot_.len_c_ / 2.0, 0.0, 0.0, obstacle.radius_),
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
    update()
    {
      avoid_ee_.obstacle_ = obstacle_.point_;
      avoid_wrist_.obstacle_ = obstacle_.point_;
      avoid_ellbow_.obstacle_ = obstacle_.point_;
      avoid_base_.obstacle_ = obstacle_.point_;
    
      orient_ee_.goal_ = z_angle_;
    
      repulse_base_.repulsor_ = repulsor_.point_;
      repulse_ellbow_.repulsor_ = repulsor_.point_;
      repulse_wrist_.repulsor_ = repulsor_.point_;
      repulse_ee_.repulsor_ = repulsor_.point_;
    
      if (dbgos_) {
	*dbgos_ << dbgpre_ << "==================================================\n"
		<< dbgpre_ << "BaseWaypoint::update()\n";
	print(robot_.getPosition(), *dbgos_, "current position", dbgpre2_);
	print(robot_.getVelocity(), *dbgos_, "current velocity", dbgpre2_);
      }
    
      for (size_t ii(0); ii < tasks_.size(); ++ii) {
	tasks_[ii]->update(robot_);
      }
      for (size_t ii(0); ii < objectives_.size(); ++ii) {
	objectives_[ii]->update(robot_);
      }
    
      if (dbgos_) {
	*dbgos_ << dbgpre_ << "--------------------------------------------------\n"
		<< dbgpre_ << "trying without constraints first\n";
      }
    
      ssize_t const ndof(robot_.getPosition().size());
      Vector qdd_t;
      Matrix N_t;
      perform_prioritization(Matrix::Identity(ndof, ndof),
			     tasks_,
			     qdd_t,
			     N_t,
			     dbgos_,
			     dbgpre_ + "task   ");
    
      Vector qdd_o(Vector::Zero(robot_.getPosition().size()));
      for (size_t ii(0); ii < objectives_.size(); ++ii) {
	if (objectives_[ii]->isActive()) {
	  Matrix Jinv;
	  pseudo_inverse_moore_penrose(objectives_[ii]->Jacobian_, Jinv);
	  qdd_o += Jinv * objectives_[ii]->delta_;
	}
      }
    
      Vector qdd_res(qdd_t + N_t * qdd_o);
      Vector qd_res(robot_.getVelocity() + timestep_ * qdd_res);
      Vector q_res(robot_.getPosition() + timestep_ * qd_res);
    
      if (dbgos_) {
	print(qdd_res, *dbgos_, "unconstrained acceleration", dbgpre2_);
	print(qd_res, *dbgos_, "resulting unconstrained velocity", dbgpre2_);
	print(q_res, *dbgos_, "resulting unconstrained position", dbgpre2_);
      }
    
      Vector const oldpos(robot_.getPosition()); // will need this in case of constraints
      Vector const oldvel(robot_.getVelocity()); // will need this in case of constraints
      robot_.update(q_res, qd_res);
    
      bool need_constraints(false);
      for (size_t ii(0); ii < constraints_.size(); ++ii) {
	constraints_[ii]->update(robot_);
	if (constraints_[ii]->isActive()) {
	  if (dbgos_) {
	    *dbgos_ << dbgpre_ << "constraint [" << ii << "] is active\n";
	  }
	  need_constraints = true;
	}
      }
    
      if ( ! need_constraints) {
	if (dbgos_) {
	  *dbgos_ << dbgpre_ << "all constraints are inactive\n";
	}
	return;
      }
    
      if (dbgos_) {
	*dbgos_ << dbgpre_ << "--------------------------------------------------\n"
		<< dbgpre_ << "recomputing with constraints enabled\n";
      }
    
      Vector dq_c;
      Matrix N_c;
      perform_prioritization(Matrix::Identity(ndof, ndof),
			     constraints_,
			     dq_c,
			     N_c,
			     dbgos_,
			     dbgpre_ + "constr ");
    
      // The constraints cheat with the robot state: they directly work
      // on the positions that would have been achieved without
      // constraints.
      //
      // Semi-open question: after repairing the position and velocity
      // to something consistent with the constraints, do we then then
      // re-run the tasks and objectives? I think yes, to give tasks a
      // chance to get fulfilled within the constraints. But that does
      // not seem to work quite yet...
      //
      // Also, wouldn't the position and velocity change for the
      // constraints then influence the constraints themselves? We
      // should thus re-run the constraints as well, possibly leading to
      // another correction and so forth ad infinitum. But the nullspace
      // of the constraints at least should not change too much, so we
      // can probably skip the chicken-and-egg constraint update
      // problem.
      //
      // Note that the non-constrained velocity would be (q_res + dq_c -
      // oldpos) / timestep_ but we're pre-multiplying with N_c and dq_c
      // is perpendicular to that so we don't need to add it.
      //
      robot_.update(q_res + dq_c, N_c * (q_res - oldpos) / timestep_);
    
      if (dbgos_) {
	print(dq_c, *dbgos_, "position correction to satisfy constraints", dbgpre2_);
	print(N_c, *dbgos_, "nullspace of constrains", dbgpre2_);
	print(robot_.getVelocity(), *dbgos_, "resulting constrained velocity", dbgpre2_);
	print(robot_.getPosition(), *dbgos_, "resulting constrained position", dbgpre2_);
      }
    
      for (size_t ii(0); ii < tasks_.size(); ++ii) {
	tasks_[ii]->update(robot_);
      }
      for (size_t ii(0); ii < objectives_.size(); ++ii) {
	objectives_[ii]->update(robot_);
      }
    
      // Re-run task priority scheme, but seed it with the constraint nullspace this time.
    
      perform_prioritization(N_c,
			     tasks_,
			     qdd_t,
			     N_t,
			     dbgos_,
			     dbgpre_ + "task   ");
    
      qdd_o = Vector::Zero(robot_.getPosition().size());
      for (size_t ii(0); ii < objectives_.size(); ++ii) {
	if (objectives_[ii]->isActive()) {
	  Matrix Jinv;
	  pseudo_inverse_moore_penrose(objectives_[ii]->Jacobian_, Jinv);
	  qdd_o += Jinv * objectives_[ii]->delta_;
	}
      }
    
      qdd_res = qdd_t + N_t * qdd_o;
      qd_res = robot_.getVelocity() + timestep_ * qdd_res;
      q_res = robot_.getPosition() + timestep_ * qd_res;
    
      if (dbgos_) {
	print(qdd_res, *dbgos_, "constrained acceleration", dbgpre2_);
	print(qd_res, *dbgos_, "resulting constrained velocity", dbgpre2_);
	print(q_res, *dbgos_, "resulting constrained position", dbgpre2_);
      }
    
      robot_.update(q_res, qd_res);
    }
    
    
    NormalWaypoint::
    NormalWaypoint(InteractionHandle const & obstacle,
		   InteractionHandle const & repulsor,
		   double const & z_angle)
      : BaseWaypoint(obstacle, repulsor, z_angle)
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
    update()
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
      
      BaseWaypoint::update();
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
    BoundaryWaypoint(InteractionHandle const & obstacle,
		     InteractionHandle const & repulsor,
		     double const & z_angle,
		     Vector const * eegoal,
		     Vector const * baseattractor)
      : BaseWaypoint(obstacle, repulsor, z_angle),
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
    update()
    {
      eetask_.goal_ = *eegoal_;
      attract_base_.attractor_ = *baseattractor_;
      BaseWaypoint::update();
    }
    
  }

}
