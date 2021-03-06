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

#include "algorithm.hpp"
#include "task.hpp"
#include "model.hpp"
#include "joint_limit_constraint.hpp"
#include "point_mindist_constraint.hpp"
#include "position_control.hpp"
#include "point_attraction.hpp"
#include "point_repulsion.hpp"
#include "posture_damping.hpp"
#include "print.hpp"
#include "pseudo_inverse.hpp"
#include <gtk/gtk.h>
#include <cmath>
#include <iostream>
#include <list>
#include <err.h>

using namespace kinematic_elastic;


static double const dimx(10.);
static double const dimy(8.);
static double const lwscale(5.0);

static GtkWidget * gw(0);
static gint gw_width(800), gw_height(640);
static gint gw_sx, gw_sy, gw_x0, gw_y0;

static bool verbose(false);
static int play(0);

struct handle_s {
  handle_s(double radius, double red, double green, double blue, double alpha)
    : point_(3),
      radius_(radius),
      red_(red),
      green_(green),
      blue_(blue),
      alpha_(alpha)
  {
  }
  
  Vector point_;
  double radius_, red_, green_, blue_, alpha_;
};

static handle_s eeleft    (0.2, 0.0, 1.0, 0.0, 0.5);
static handle_s baseleft  (0.2, 0.0, 1.0, 0.5, 0.5);
static handle_s eeright   (0.2, 1.0, 0.0, 0.0, 0.5);
static handle_s baseright (0.2, 1.0, 0.0, 0.5, 0.5);
static handle_s repulsor1 (1.5, 0.0, 0.0, 1.0, 0.2);
static handle_s repulsor2 (1.5, 0.0, 0.0, 1.0, 0.2);

static handle_s * handle[] = { &eeleft, &baseleft, &eeright, &baseright, &repulsor1, &repulsor2, 0 };
static handle_s * grabbed(0);
static Vector grab_offset(3);


class Robot
  : public Model
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  Robot()
    : radius_(0.5),
      len_a_(0.8),
      len_b_(0.6),
      len_c_(0.3),
      pos_a_(3),
      pos_b_(3),
      pos_c_(3)
  {
  }
  
  
  virtual Vector const & getPosition() const
  {
    return position_;
  }
  
  
  virtual Vector const & getVelocity() const
  {
    return velocity_;
  }
  
  
  virtual Transform frame(size_t node) const
  {
    Transform tf(Transform::Identity());
    switch (node) {
    case 0:
      tf.translation() << position_[0], position_[1], 0.0;
      break;
    case 1:
      tf.translation() << position_[0], position_[1], 0.0;
      tf.linear() << c2_, -s2_, 0.0, s2_, c2_, 0.0, 0.0, 0.0, 1.0;
      break;
    case 2:
      tf.translation() << pos_a_[0], pos_a_[1], 0.0;
      tf.linear() << c23_, -s23_, 0.0, s23_, c23_, 0.0, 0.0, 0.0, 1.0;
      break;
    case 3:
      tf.translation() << pos_b_[0], pos_b_[1], 0.0;
      tf.linear() << c234_, -s234_, 0.0, s234_, c234_, 0.0, 0.0, 0.0, 1.0;
      break;
    default:
      errx (EXIT_FAILURE, "Robot::frame() called on invalid node %zu", node);
    }
    return tf;
  }
  
  
  virtual Matrix computeJxo(size_t node, Vector const & gpoint) const
  {
    Matrix Jxo(Matrix::Zero(6, 5));
    switch (node) {
    case 3:
      Jxo(0, 4) = pos_b_[1] - gpoint[1];
      Jxo(1, 4) = gpoint[0] - pos_b_[0];
      Jxo(5, 4) = 1.0;
    case 2:
      Jxo(0, 3) = pos_a_[1] - gpoint[1];
      Jxo(1, 3) = gpoint[0] - pos_a_[0];
      Jxo(5, 3) = 1.0;
    case 1:
      Jxo(0, 2) = position_[1] - gpoint[1];
      Jxo(1, 2) = gpoint[0]    - position_[0];
      Jxo(5, 2) = 1.0;
    case 0:
      Jxo(0, 0) = 1.0;
      Jxo(1, 1) = 1.0;
      break;
    default:
      errx (EXIT_FAILURE, "Robot::computeJxo() called on invalid node %zu", node);
    }
    // if (verbose) {
    //   cout << "Robot::computeJxo()\n"
    // 	   << "  node " << node << "\n";
    //   print(gpoint, cout, "  gpoint", "    ");
    //   print(Jxo, cout, "  Jxo", "    ");
    // }
    return Jxo;
  }
  
  
  virtual void update(Vector const & position, Vector const & velocity)
  {
    if (position.size() != 5) {
      errx (EXIT_FAILURE, "Robot::update(): position has %zu DOF (but needs 5)", (size_t) position.size());
    }
    if (velocity.size() != 5) {
      errx (EXIT_FAILURE, "Robot::update(): velocity has %zu DOF (but needs 5)", (size_t) velocity.size());
    }
    position_ = position;
    velocity_ = velocity;
    
    c2_ = cos(position_[2]);
    s2_ = sin(position_[2]);
    ac2_ = len_a_ * c2_;
    as2_ = len_a_ * s2_;
    
    q23_ = position_[2] + position_[3];
    c23_ = cos(q23_);
    s23_ = sin(q23_);
    bc23_ = len_b_ * c23_;
    bs23_ = len_b_ * s23_;
    
    q234_ = q23_ + position_[4];
    c234_ = cos(q234_);
    s234_ = sin(q234_);
    cc234_ = len_c_ * c234_;
    cs234_ = len_c_ * s234_;
    
    pos_a_ <<
      position_[0] + ac2_,
      position_[1] + as2_,
      0.0;
    pos_b_ <<
      pos_a_[0] + bc23_,
      pos_a_[1] + bs23_,
      0.0;
    pos_c_ <<
      pos_b_[0] + cc234_,
      pos_b_[1] + cs234_,
      0.0;
  }
  
  
  void draw(cairo_t * cr, double pixelsize) const
  {
    cairo_save(cr);
    
    // translucent disk for base
    cairo_set_source_rgba(cr, 0.7, 0.7, 0.7, 0.5);
    cairo_arc(cr, position_[0], position_[1], radius_, 0., 2. * M_PI);
    cairo_fill(cr);
    
    // thick circle outline for base
    cairo_set_source_rgb(cr, 0.2, 0.2, 0.2);
    cairo_set_line_width(cr, lwscale * 3.0 / pixelsize);
    cairo_arc(cr, position_[0], position_[1], radius_, 0., 2. * M_PI);
    cairo_stroke(cr);
    
    // thick line for arms
    cairo_set_source_rgb(cr, 0.2, 0.2, 0.2);
    cairo_set_line_width(cr, lwscale * 3.0 / pixelsize);
    cairo_move_to(cr, position_[0], position_[1]);
    cairo_line_to(cr, pos_a_[0], pos_a_[1]);
    cairo_line_to(cr, pos_b_[0], pos_b_[1]);
    cairo_line_to(cr, pos_c_[0], pos_c_[1]);
    cairo_stroke(cr);
    
    cairo_restore(cr);
  }

  //protected: or whatnot...
  
  double const radius_;
  double const len_a_;
  double const len_b_;
  double const len_c_;
  
  Vector position_;
  Vector velocity_;
  Vector pos_a_;
  Vector pos_b_;
  Vector pos_c_;
  
  double c2_;
  double s2_;
  double c23_;
  double s23_;
  double c234_;
  double s234_;
  double q23_;
  double q234_;
  double ac2_;
  double as2_;
  double bc23_;
  double bs23_;
  double cc234_;
  double cs234_;
};


class OriZControl
  : public Task
{
public:
  OriZControl(size_t node,
	      double kp,
	      double kd)
    : kp_(kp),
      kd_(kd),
      node_(node)
  {
  }
  
  virtual void init(Model const & model)
  {
    Vector const ex(model.frame(node_).linear().block(0, 0, 3, 1));
    angle_ = atan2(ex[1], ex[0]);
    goal_ = angle_;
    delta_ = Vector::Zero(1);
    Jacobian_ = model.computeJxo(node_, Vector::Zero(3)).block(5, 0, 1, model.getPosition().size());
  }
  
  virtual void update(Model const & model)
  {
    Vector const ex(model.frame(node_).linear().block(0, 0, 3, 1));
    angle_ = atan2(ex[1], ex[0]);
    Jacobian_ = model.computeJxo(node_, Vector::Zero(3)).block(5, 0, 1, model.getPosition().size());
    delta_[0] = kp_ * (goal_ - angle_);
    delta_ -= kd_ * Jacobian_ * model.getVelocity();
  }
  
  double kp_;
  double kd_;
  size_t node_;
  double angle_;
  double goal_;
};


class BaseWaypoint
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  
  BaseWaypoint()
    : timestep_(1e-2),
      repulse1_base_  (0,                 0.0, 0.0, 0.0, 100.0, repulsor1.radius_),
      repulse1_ellbow_(1,       robot_.len_a_, 0.0, 0.0, 100.0, repulsor1.radius_),
      repulse1_wrist_ (2,       robot_.len_b_, 0.0, 0.0, 100.0, repulsor1.radius_),
      repulse1_ee_    (3,       robot_.len_c_, 0.0, 0.0, 100.0, repulsor1.radius_),
      repulse2_base_  (0,                 0.0, 0.0, 0.0, 100.0, repulsor2.radius_),
      repulse2_ellbow_(1,       robot_.len_a_, 0.0, 0.0, 100.0, repulsor2.radius_),
      repulse2_wrist_ (2,       robot_.len_b_, 0.0, 0.0, 100.0, repulsor2.radius_),
      repulse2_ee_    (3,       robot_.len_c_, 0.0, 0.0, 100.0, repulsor2.radius_),
      joint_damping_ (10.0)
  {
    objectives_.push_back(&repulse1_base_);
    objectives_.push_back(&repulse1_ellbow_);
    objectives_.push_back(&repulse1_wrist_);
    objectives_.push_back(&repulse1_ee_);
    objectives_.push_back(&repulse2_base_);
    objectives_.push_back(&repulse2_ellbow_);
    objectives_.push_back(&repulse2_wrist_);
    objectives_.push_back(&repulse2_ee_);
    objectives_.push_back(&joint_damping_);
  }
  
  
  virtual ~BaseWaypoint()
  {
  }
  
  
  virtual void draw(cairo_t * cr, double pixelsize)
  {
    robot_.draw(cr, pixelsize);
    
    cairo_save(cr);
    
    // repulsion vectors
    
    cairo_set_source_rgb(cr, 0.4, 0.4, 1.0);
    cairo_set_line_width(cr, lwscale * 1.0 / pixelsize);
    if (repulse1_base_.isActive()) {
      cairo_move_to(cr, repulse1_base_.gpoint_[0], repulse1_base_.gpoint_[1]);
      cairo_line_to(cr, repulse1_base_.gpoint_[0] + repulse1_base_.delta_[0] / repulse1_base_.gain_, repulse1_base_.gpoint_[1] + repulse1_base_.delta_[1] / repulse1_base_.gain_);
      cairo_stroke(cr);
    }
    if (repulse1_ellbow_.isActive()) {
      cairo_move_to(cr, repulse1_ellbow_.gpoint_[0], repulse1_ellbow_.gpoint_[1]);
      cairo_line_to(cr, repulse1_ellbow_.gpoint_[0] + repulse1_ellbow_.delta_[0] / repulse1_ellbow_.gain_, repulse1_ellbow_.gpoint_[1] + repulse1_ellbow_.delta_[1] / repulse1_ellbow_.gain_);
      cairo_stroke(cr);
    }
    if (repulse1_wrist_.isActive()) {
      cairo_move_to(cr, repulse1_wrist_.gpoint_[0], repulse1_wrist_.gpoint_[1]);
      cairo_line_to(cr, repulse1_wrist_.gpoint_[0] + repulse1_wrist_.delta_[0] / repulse1_wrist_.gain_, repulse1_wrist_.gpoint_[1] + repulse1_wrist_.delta_[1] / repulse1_wrist_.gain_);
      cairo_stroke(cr);
    }
    if (repulse1_ee_.isActive()) {
      cairo_move_to(cr, repulse1_ee_.gpoint_[0], repulse1_ee_.gpoint_[1]);
      cairo_line_to(cr, repulse1_ee_.gpoint_[0] + repulse1_ee_.delta_[0] / repulse1_ee_.gain_, repulse1_ee_.gpoint_[1] + repulse1_ee_.delta_[1] / repulse1_ee_.gain_);
      cairo_stroke(cr);
    }
    if (repulse2_base_.isActive()) {
      cairo_move_to(cr, repulse2_base_.gpoint_[0], repulse2_base_.gpoint_[1]);
      cairo_line_to(cr, repulse2_base_.gpoint_[0] + repulse2_base_.delta_[0] / repulse2_base_.gain_, repulse2_base_.gpoint_[1] + repulse2_base_.delta_[1] / repulse2_base_.gain_);
      cairo_stroke(cr);
    }
    if (repulse2_ellbow_.isActive()) {
      cairo_move_to(cr, repulse2_ellbow_.gpoint_[0], repulse2_ellbow_.gpoint_[1]);
      cairo_line_to(cr, repulse2_ellbow_.gpoint_[0] + repulse2_ellbow_.delta_[0] / repulse2_ellbow_.gain_, repulse2_ellbow_.gpoint_[1] + repulse2_ellbow_.delta_[1] / repulse2_ellbow_.gain_);
      cairo_stroke(cr);
    }
    if (repulse2_wrist_.isActive()) {
      cairo_move_to(cr, repulse2_wrist_.gpoint_[0], repulse2_wrist_.gpoint_[1]);
      cairo_line_to(cr, repulse2_wrist_.gpoint_[0] + repulse2_wrist_.delta_[0] / repulse2_wrist_.gain_, repulse2_wrist_.gpoint_[1] + repulse2_wrist_.delta_[1] / repulse2_wrist_.gain_);
      cairo_stroke(cr);
    }
    if (repulse2_ee_.isActive()) {
      cairo_move_to(cr, repulse2_ee_.gpoint_[0], repulse2_ee_.gpoint_[1]);
      cairo_line_to(cr, repulse2_ee_.gpoint_[0] + repulse2_ee_.delta_[0] / repulse2_ee_.gain_, repulse2_ee_.gpoint_[1] + repulse2_ee_.delta_[1] / repulse2_ee_.gain_);
      cairo_stroke(cr);
    }
    
    cairo_restore(cr);
  }
  
  
  virtual void init(Vector const & position, Vector const & velocity)
  {
    robot_.update(position, velocity);
    
    for (size_t ii(0); ii < constraints_.size(); ++ii) {
      constraints_[ii]->init(robot_);
    }
    for (size_t ii(0); ii < tasks_.size(); ++ii) {
      tasks_[ii]->init(robot_);
    }
    for (size_t ii(0); ii < objectives_.size(); ++ii) {
      objectives_[ii]->init(robot_);
    }
  }
  
  
  virtual void update()
  {
    repulse1_base_.repulsor_ = repulsor1.point_;
    repulse1_ellbow_.repulsor_ = repulsor1.point_;
    repulse1_wrist_.repulsor_ = repulsor1.point_;
    repulse1_ee_.repulsor_ = repulsor1.point_;
    repulse2_base_.repulsor_ = repulsor2.point_;
    repulse2_ellbow_.repulsor_ = repulsor2.point_;
    repulse2_wrist_.repulsor_ = repulsor2.point_;
    repulse2_ee_.repulsor_ = repulsor2.point_;
    
    ostream * dbgos(0);
    if (verbose) {
      dbgos = &cout;
      cout << "==================================================\n"
	   << "BaseWaypoint::update()\n";
      print(robot_.getPosition(), cout, "current position", "  ");
      print(robot_.getVelocity(), cout, "current velocity", "  ");
    }
    
    for (size_t ii(0); ii < tasks_.size(); ++ii) {
      tasks_[ii]->update(robot_);
    }
    for (size_t ii(0); ii < objectives_.size(); ++ii) {
      objectives_[ii]->update(robot_);
    }
    
    if (verbose) {
      cout << "--------------------------------------------------\n"
	   << "trying without constraints first\n";
    }
    
    ssize_t const ndof(robot_.getPosition().size());
    Vector qdd_t;
    Matrix N_t;
    perform_prioritization(Matrix::Identity(ndof, ndof),
			   tasks_,
			   qdd_t,
			   N_t,
			   dbgos, "task   ");
    
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
    
    if (verbose) {
      print(qdd_res, cout, "unconstrained acceleration", "  ");
      print(qd_res, cout, "resulting unconstrained velocity", "  ");
      print(q_res, cout, "resulting unconstrained position", "  ");
    }
    
    Vector const oldpos(robot_.getPosition()); // will need this in case of constraints
    Vector const oldvel(robot_.getVelocity()); // will need this in case of constraints
    robot_.update(q_res, qd_res);
    
    bool need_constraints(false);
    for (size_t ii(0); ii < constraints_.size(); ++ii) {
      constraints_[ii]->update(robot_);
      if (constraints_[ii]->isActive()) {
    	if (verbose) {
    	  cout << "constraint [" << ii << "] is active\n";
    	}
    	need_constraints = true;
      }
    }
    
    if ( ! need_constraints) {
      if (verbose) {
    	cout << "all constraints are inactive\n";
      }
      return;
    }
    
    if (verbose) {
      cout << "--------------------------------------------------\n"
    	   << "recomputing with constraints enabled\n";
    }
    
    Vector dq_c;
    Matrix N_c;
    perform_prioritization(Matrix::Identity(ndof, ndof),
			   constraints_,
			   dq_c,
			   N_c,
			   dbgos, "constr ");
    
    // The constraints cheat with the robot state: they directly work
    // on the positions that would have been achieved without
    // constraints.
    //
    // Open question: after repairing the position and velocity to
    // something consistent with the constraints, do we then then
    // re-run the tasks and objectives? Or is that taken care of
    // automatically in the next iteration? Also, if we re-run
    // everything here, wouldn't the changed position and velocity
    // also influence the constraints themselves? In the latter case,
    // we'd need to re-run the constraints as well, possibly leading
    // to another correction and so forth ad infinitum. But the
    // nullspace of the constraints at least should not change (too
    // much anyhow) so we can probably skip the chicken-and-egg
    // constraint update problem.
    //
    // Note that the non-constrained velocity would be (q_res + dq_c -
    // oldpos) / timestep_ but we're pre-multiplying with N_c and dq_c
    // is perpendicular to that so we don't need to add it.
    //
    robot_.update(q_res + dq_c, N_c * (q_res - oldpos) / timestep_);
    
    if (verbose) {
      print(dq_c, cout, "position correction to satisfy constraints", "  ");
      print(N_c, cout, "nullspace of constrains", "  ");
      print(robot_.getVelocity(), cout, "resulting constrained velocity", "  ");
      print(robot_.getPosition(), cout, "resulting constrained position", "  ");
    }
    
    // // for (size_t ii(0); ii < tasks_.size(); ++ii) {
    // //   tasks_[ii]->update(robot_);
    // // }
    // // for (size_t ii(0); ii < objectives_.size(); ++ii) {
    // //   objectives_[ii]->update(robot_);
    // // }
    
    // // // Re-run task priority scheme, but seed it with the constraint nullspace this time.
    
    // // perform_prioritization(N_c,
    // // 			   tasks_,
    // // 			   qdd_t,
    // // 			   N_t,
    // // 			   dbgos, "task   ");
    
    // // qdd_o = Vector::Zero(robot_.getPosition().size());
    // // for (size_t ii(0); ii < objectives_.size(); ++ii) {
    // //   if (objectives_[ii]->isActive()) {
    // // 	Matrix Jinv;
    // // 	pseudo_inverse_moore_penrose(objectives_[ii]->Jacobian_, Jinv);
    // // 	qdd_o += Jinv * objectives_[ii]->delta_;
    // //   }
    // // }
    
    // // qdd_res = qdd_t + N_t * qdd_o;
    // // qd_res = robot_.getVelocity() + timestep_ * qdd_res;
    // // q_res = robot_.getPosition() + timestep_ * qd_res;
    
    // // if (verbose) {
    // //   print(qdd_res, cout, "constrained acceleration", "  ");
    // //   print(qd_res, cout, "resulting constrained velocity", "  ");
    // //   print(q_res, cout, "resulting constrained position", "  ");
    // // }
    
    // // robot_.update(q_res, qd_res);
  }

  ////protected:
  double timestep_;
  Robot robot_; // XXXX keep this before any constraints so we can use its values for initializing them
  
  PointRepulsion repulse1_base_;
  PointRepulsion repulse1_ellbow_;
  PointRepulsion repulse1_wrist_;
  PointRepulsion repulse1_ee_;
  PointRepulsion repulse2_base_;
  PointRepulsion repulse2_ellbow_;
  PointRepulsion repulse2_wrist_;
  PointRepulsion repulse2_ee_;
  
  PostureDamping joint_damping_;
  
  vector<Task *> constraints_;
  vector<Task *> tasks_;
  vector<Task *> objectives_;
};  


class LeftWaypoint
  : public BaseWaypoint
{
public:
  LeftWaypoint()
    : attract_ee_  (3, robot_.len_c_, 0.0, 0.0, 100.0, 2.0),
      attract_base_(0,           0.0, 0.0, 0.0, 100.0, 2.0)
  {
    objectives_.push_back(&attract_ee_);
    objectives_.push_back(&attract_base_);
  }
  
  
  virtual void draw(cairo_t * cr, double pixelsize)
  {
    BaseWaypoint::draw(cr, pixelsize);
    
    cairo_set_source_rgb(cr, 0.4, 1.0, 0.4);
    cairo_set_line_width(cr, lwscale * 1.0 / pixelsize);
    
    if (attract_ee_.isActive()) {
      cairo_move_to(cr, attract_ee_.gpoint_[0], attract_ee_.gpoint_[1]);
      cairo_line_to(cr, attract_ee_.gpoint_[0] + attract_ee_.delta_[0] / attract_ee_.gain_, attract_ee_.gpoint_[1] + attract_ee_.delta_[1] / attract_ee_.gain_);
      cairo_stroke(cr);
    }
    
    if (attract_base_.isActive()) {
      cairo_move_to(cr, attract_base_.gpoint_[0], attract_base_.gpoint_[1]);
      cairo_line_to(cr, attract_base_.gpoint_[0] + attract_base_.delta_[0] / attract_base_.gain_, attract_base_.gpoint_[1] + attract_base_.delta_[1] / attract_base_.gain_);
      cairo_stroke(cr);
    }
  }
  
  
  virtual void update()
  {
    attract_ee_.attractor_ = eeleft.point_;
    attract_base_.attractor_ = baseleft.point_;
    BaseWaypoint::update();
  }
  
  PointAttraction attract_ee_;
  PointAttraction attract_base_;
};


class RightWaypoint
  : public BaseWaypoint
{
public:
  RightWaypoint()
    : control_ee_  (3, robot_.len_c_, 0.0, 0.0, 100.0, 20.0),
      attract_base_(0,           0.0, 0.0, 0.0, 100.0, 2.0)
  {
    tasks_.push_back(&control_ee_);
    objectives_.push_back(&attract_base_);
  }
  
  
  virtual void draw(cairo_t * cr, double pixelsize)
  {
    BaseWaypoint::draw(cr, pixelsize);

    // thin line for end effector task
    cairo_set_source_rgb(cr, 1.0, 0.4, 0.4);
    cairo_set_line_width(cr, lwscale * 1.0 / pixelsize);
    cairo_move_to(cr, control_ee_.gpoint_[0], control_ee_.gpoint_[1]);
    cairo_line_to(cr, control_ee_.goal_[0], control_ee_.goal_[1]);
    cairo_stroke(cr);
    
    // base attraction
    if (attract_base_.isActive()) {
      cairo_set_source_rgb(cr, 0.4, 1.0, 0.4);
      cairo_set_line_width(cr, lwscale * 1.0 / pixelsize);
      cairo_move_to(cr, attract_base_.gpoint_[0], attract_base_.gpoint_[1]);
      cairo_line_to(cr, attract_base_.gpoint_[0] + attract_base_.delta_[0] / attract_base_.gain_, attract_base_.gpoint_[1] + attract_base_.delta_[1] / attract_base_.gain_);
      cairo_stroke(cr);
    }
  }
  
  
  virtual void update()
  {
    control_ee_.goal_ = eeright.point_;
    attract_base_.attractor_ = baseright.point_;
    BaseWaypoint::update();
  }
  
  PositionControl control_ee_;
  PointAttraction attract_base_;
};


// could templatize on Waypoint or pass in a factory or something along
// those lines to make it robot agnostic...
class Elastic
{
public:
  typedef list<BaseWaypoint *> path_t;

  ~Elastic()
  {
    clear();
  }
  
  void clear()
  {
    for (path_t::iterator ii(path_.begin()); ii != path_.end(); ++ii) {
      delete *ii;
    }
    path_.clear();
  }
  
  
  void init(Vector const & state)
  {
    clear();
    
    path_.push_back(new LeftWaypoint());
    path_.push_back(new RightWaypoint());
    
    for (path_t::iterator ii(path_.begin()); ii != path_.end(); ++ii) {
      (*ii)->init(state, Vector::Zero(state.size()));
    }
  }
  
  
  void update()
  {
    if (verbose) {
      cout << "\n**************************************************\n";
    }
    for (path_t::iterator ii(path_.begin()); ii != path_.end(); ++ii) {
      (*ii)->update();
    }
  }
  
  void draw(cairo_t * cr, double pixelsize)
  {
    for (path_t::reverse_iterator ii(path_.rbegin()); ii != path_.rend(); ++ii) {
      (*ii)->draw(cr, pixelsize);
    }
  }
  
  ////private:
  path_t path_;
};


static Elastic elastic;


static void update()
{
  elastic.update();
  gtk_widget_queue_draw(gw);
}


static void cb_play(GtkWidget * ww, gpointer data)
{
  if (play) {
    play = 0;
  }
  else {
    play = 1;
  }
}


static void cb_normalize(GtkWidget * ww, gpointer data)
{
  for (Elastic::path_t::iterator ii(elastic.path_.begin()); ii != elastic.path_.end(); ++ii) {
    if (verbose) {
      if (fabs((*ii)->robot_.position_[2]) > M_PI) {
	cout << "normalize " << (*ii)->robot_.position_[2]
	     << " to " << normangle((*ii)->robot_.position_[2]) << "\n";
      }
    }
    (*ii)->robot_.position_[2] = normangle((*ii)->robot_.position_[2]);
    (*ii)->robot_.update((*ii)->robot_.position_, (*ii)->robot_.velocity_);
  }
}


static void cb_next(GtkWidget * ww, gpointer data)
{
  if (play) {
    play = 0;
  }
  else {
    update();
  }    
}


static void cb_quit(GtkWidget * ww, gpointer data)
{
  gtk_main_quit();
}


static gint cb_expose(GtkWidget * ww,
		      GdkEventExpose * ee,
		      gpointer data)
{
  cairo_t * cr = gdk_cairo_create(ee->window);
  
  cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
  cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
  cairo_rectangle(cr, 0, 0, gw_width, gw_height);
  cairo_fill(cr);
  
  cairo_translate(cr, gw_x0, gw_y0);
  cairo_scale(cr, gw_sx, gw_sy);
  
  elastic.draw(cr, gw_sx);
  
  cairo_set_line_width(cr, lwscale * 1.0 / gw_sx);
  for (handle_s ** hh(handle); *hh != 0; ++hh) {
    cairo_set_source_rgba(cr, (*hh)->red_, (*hh)->green_, (*hh)->blue_, (*hh)->alpha_);
    cairo_arc(cr, (*hh)->point_[0], (*hh)->point_[1], (*hh)->radius_, 0.0, 2.0 * M_PI);
    cairo_fill(cr);
  }
  
  cairo_destroy(cr);
  
  return TRUE;
}


static gint cb_size_allocate(GtkWidget * ww,
			     GtkAllocation * aa,
			     gpointer data)
{
  gw_width = aa->width;
  gw_height = aa->height;
  
  gw_sx = gw_width / dimx;
  if (gw_sx < 1) {
    gw_sx = 1;
  }
  gw_sy = - gw_height / dimy;
  if ( - gw_sy < 1) {
    gw_sy = -1;
  }
  if (gw_sx > - gw_sy) {
    gw_sx = - gw_sy;
  }
  else {
    gw_sy = - gw_sx;
  }
  gw_x0 = (gw_width - dimx * gw_sx) / 2;
  gw_y0 = gw_height - (gw_height + dimy * gw_sy) / 2;
  
  return TRUE;
}


static gint cb_click(GtkWidget * ww,
		     GdkEventButton * bb,
		     gpointer data)
{
  if (bb->type == GDK_BUTTON_PRESS) {
    Vector point(3);
    point << (bb->x - gw_x0) / (double) gw_sx, (bb->y - gw_y0) / (double) gw_sy, 0.0;
    for (handle_s ** hh(handle); *hh != 0; ++hh) {
      Vector offset = (*hh)->point_ - point;
      if (offset.norm() <= (*hh)->radius_) {
    	grab_offset = offset;
    	grabbed = *hh;
    	break;
      }
    }
  }
  else if (bb->type == GDK_BUTTON_RELEASE) {
    grabbed = 0;
    cout << "Why don't I get GDK_BUTTON_RELEASE (under OSX)?\n";
  }
  
  return TRUE;
}


static gint cb_motion(GtkWidget * ww,
		      GdkEventMotion * ee)
{
  int mx, my;
  GdkModifierType modifier;
  gdk_window_get_pointer(ww->window, &mx, &my, &modifier);
  
  if (0 != grabbed) {
    Vector point(3);
    point << (mx - gw_x0) / (double) gw_sx, (my - gw_y0) / (double) gw_sy, 0.0;
    grabbed->point_ = point + grab_offset;
    gtk_widget_queue_draw(gw);
  }
  
  return TRUE;
}


static gint idle(gpointer data)
{
  if (play) {
    update();
  }
  return TRUE;
}


static void init_gui(int * argc, char *** argv)
{
  GtkWidget *window, *vbox, *hbox, *btn;
  
  gtk_init(argc, argv);
  
  window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  
  vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER (window), vbox);
  gtk_widget_show(vbox);
  
  gw = gtk_drawing_area_new();
  g_signal_connect(gw, "expose_event", G_CALLBACK (cb_expose), NULL);
  g_signal_connect(gw, "size_allocate", G_CALLBACK (cb_size_allocate), NULL);
  g_signal_connect(gw, "button_press_event", G_CALLBACK (cb_click), NULL);
  g_signal_connect(gw, "motion_notify_event", G_CALLBACK (cb_motion), NULL);
  gtk_widget_set_events(gw,
			GDK_BUTTON_PRESS_MASK |
			GDK_BUTTON_RELEASE_MASK |
			GDK_BUTTON_MOTION_MASK);
  
  gtk_widget_show(gw);
  
  gtk_widget_set_size_request(gw, gw_width, gw_height);
  gtk_box_pack_start(GTK_BOX (vbox), gw, TRUE, TRUE, 0);
  
  hbox = gtk_hbox_new(TRUE, 3);
  gtk_box_pack_start(GTK_BOX (vbox), hbox, FALSE, TRUE, 0);
  gtk_widget_show(hbox);
  
  btn = gtk_button_new_with_label("play");
  g_signal_connect(btn, "clicked", G_CALLBACK (cb_play), NULL);
  gtk_box_pack_start(GTK_BOX (hbox), btn, TRUE, TRUE, 0);
  gtk_widget_show(btn);
  
  btn = gtk_button_new_with_label("next");
  g_signal_connect(btn, "clicked", G_CALLBACK (cb_next), NULL);
  gtk_box_pack_start(GTK_BOX (hbox), btn, TRUE, TRUE, 0);
  gtk_widget_show(btn);
  
  btn = gtk_button_new_with_label("normalize");
  g_signal_connect(btn, "clicked", G_CALLBACK (cb_normalize), NULL);
  gtk_box_pack_start(GTK_BOX (hbox), btn, TRUE, TRUE, 0);
  gtk_widget_show(btn);
  
  btn = gtk_button_new_with_label("quit");
  g_signal_connect(btn, "clicked", G_CALLBACK (cb_quit), NULL);
  gtk_box_pack_start(GTK_BOX (hbox), btn, TRUE, TRUE, 0);
  gtk_widget_show(btn);
  
  gtk_idle_add(idle, 0);
  
  gtk_widget_show(window);
}


int main(int argc, char ** argv)
{
  init_gui(&argc, &argv);
  
  if ((argc > 1) && (0 == strcmp("-v", argv[1]))) {
    verbose = true;
  }
  
  eeleft.point_    <<                 2.0, dimy / 2.0 + 1.0,     0.0;
  baseleft.point_  <<                 2.0,              1.0,     0.0;
  eeright.point_   <<          dimx - 2.0, dimy / 2.0 + 1.0,     0.0;
  baseright.point_ <<          dimx - 2.0,              1.0,     0.0;
  repulsor1.point_ <<                 3.0, dimy / 2.0 - 2.0,     0.0;
  repulsor2.point_ <<          dimx - 3.0, dimy / 2.0 - 2.0,     0.0;
  
  Vector posture(5);
  posture <<
    dimx / 2.0,
    dimy / 2.0,
    80.0 * deg,
    - 40.0 * deg,
    25.0 * deg;
  elastic.init(posture);
  
  gtk_main();
  
  return 0;
}
