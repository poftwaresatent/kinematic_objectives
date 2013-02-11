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

#include "base_waypoint.hpp"
#include "print.hpp"

#include "cairo_drawable.hpp"	// rfct
#include "algorithm.hpp"	// rfct
#include "pseudo_inverse.hpp"	// rfct


namespace kinematic_elastic {
  
  
  BaseWaypoint::
  BaseWaypoint(double qh_obstacle_radius,
	       double qh_repulsor_radius)
    : Waypoint(robot_),		// XXXX rfct?
      timestep_(1e-2),
      avoid_base_    (0,                 0.0, 0.0, 0.0, qh_obstacle_radius + robot_.radius_),
      avoid_ellbow_  (1,       robot_.len_a_, 0.0, 0.0, qh_obstacle_radius),
      avoid_wrist_   (2,       robot_.len_b_, 0.0, 0.0, qh_obstacle_radius),
      avoid_ee_      (3, robot_.len_c_ / 2.0, 0.0, 0.0, qh_obstacle_radius),
      orient_ee_     (3, 100.0, 20.0),
      repulse_base_  (0,                 0.0, 0.0, 0.0, 100.0, qh_repulsor_radius),
      repulse_ellbow_(1,       robot_.len_a_, 0.0, 0.0, 100.0, qh_repulsor_radius),
      repulse_wrist_ (2,       robot_.len_b_, 0.0, 0.0, 100.0, qh_repulsor_radius),
      repulse_ee_    (3,       robot_.len_c_, 0.0, 0.0, 100.0, qh_repulsor_radius),
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
  {
    drawBaseWaypoint(*this, cr, weight, pixelsize);
  }
  
  
  void BaseWaypoint::
  update(Vector const & qh_obstacle_point,
	 Vector const & qh_repulsor_point,
	 double qh_zangle)
  {
    avoid_ee_.obstacle_ = qh_obstacle_point;
    avoid_wrist_.obstacle_ = qh_obstacle_point;
    avoid_ellbow_.obstacle_ = qh_obstacle_point;
    avoid_base_.obstacle_ = qh_obstacle_point;
    
    orient_ee_.goal_ = qh_zangle;
    
    repulse_base_.repulsor_ = qh_repulsor_point;
    repulse_ellbow_.repulsor_ = qh_repulsor_point;
    repulse_wrist_.repulsor_ = qh_repulsor_point;
    repulse_ee_.repulsor_ = qh_repulsor_point;
    
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
  
}
