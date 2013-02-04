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


static double const deg(M_PI / 180.);
static double const dimx(10.);
static double const dimy(8.);

static GtkWidget * gw(0);
static gint gw_width(800), gw_height(640);
static gint gw_sx, gw_sy, gw_x0, gw_y0;

static bool verbose(false);
static int play(0);
static Vector eegoal(3);
static Vector repulsor(3);
static Vector attractor(3);
static Vector * handle[] = { &eegoal, &repulsor, &attractor, 0 };
static Vector * grabbed(0);
static double grab_radius(0.2);
static Vector grab_offset(3);


static inline double bound(double lower, double value, double upper)
{
  if (value < lower) {
    value = lower;
  }
  else if (value > upper) {
    value = upper;
  }
  return value;
}


static inline double normangle(double phi)
{
  phi = fmod(phi, 2.0 * M_PI);
  if (phi > M_PI) {
    phi -= M_PI;
  }
  else if (phi < -M_PI) {
    phi += M_PI;
  }
  return phi;
}


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
    if (verbose) {
      cout << "Robot::computeJxo()\n"
	   << "  node " << node << "\n";
      print(gpoint, cout, "  gpoint", "    ");
      print(Jxo, cout, "  Jxo", "    ");
    }
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
    cairo_set_line_width(cr, 3.0 / pixelsize);
    cairo_arc(cr, position_[0], position_[1], radius_, 0., 2. * M_PI);
    cairo_stroke(cr);
    
    // thick line for arms
    cairo_set_source_rgb(cr, 0.2, 0.2, 0.2);
    cairo_set_line_width(cr, 3.0 / pixelsize);
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


class Waypoint
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  Waypoint()
    : timestep_(1e-2),
      eetask_(3, Vector::Zero(3)),
      attract_base_(0, Vector::Zero(3)),
      repulse_ellbow_(1, Vector::Zero(3)),
      repulse_wrist_(2, Vector::Zero(3))
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
    
    eetask_.point_ << robot_.len_c_, 0.0, 0.0;
    
    tasks_.push_back(&eetask_);
    
    repulse_ellbow_.point_ << robot_.len_a_, 0.0, 0.0;
    repulse_wrist_.point_ << robot_.len_b_, 0.0, 0.0;
    
    objectives_.push_back(&attract_base_);
    objectives_.push_back(&repulse_ellbow_);
    objectives_.push_back(&repulse_wrist_);
    objectives_.push_back(&joint_damping_);
  }
  
  
  ~Waypoint()
  {
  }
  
  
  void draw(cairo_t * cr, double pixelsize)
  {
    robot_.draw(cr, pixelsize);
    
    cairo_save(cr);
    
    // thin arcs for arm joint limits
    cairo_set_source_rgb(cr, 0.5, 0.5, 1.0);
    cairo_set_line_width(cr, 1.0 / pixelsize);
    cairo_move_to(cr, robot_.pos_a_[0], robot_.pos_a_[1]);
    cairo_arc(cr, robot_.pos_a_[0], robot_.pos_a_[1], 0.1,
	      bound(-2.0*M_PI, robot_.position_[2] + joint_limits_.limits_(3, 0), 2.0*M_PI),
	      bound(-2.0*M_PI, robot_.position_[2] + joint_limits_.limits_(3, 3), 2.0*M_PI));
    cairo_line_to(cr, robot_.pos_a_[0], robot_.pos_a_[1]);
    cairo_stroke(cr);
    cairo_move_to(cr, robot_.pos_b_[0], robot_.pos_b_[1]);
    cairo_arc(cr, robot_.pos_b_[0], robot_.pos_b_[1], 0.1,
	      bound(-2.0*M_PI, robot_.q23_ + joint_limits_.limits_(4, 0), 2.0*M_PI),
	      bound(-2.0*M_PI, robot_.q23_ + joint_limits_.limits_(4, 3), 2.0*M_PI));
    cairo_line_to(cr, robot_.pos_b_[0], robot_.pos_b_[1]);
    cairo_stroke(cr);
    
    // thin line for end effector task
    cairo_set_source_rgb(cr, 1.0, 0.4, 0.4);
    cairo_set_line_width(cr, 1.0 / pixelsize);
    cairo_move_to(cr, eetask_.gpoint_[0], eetask_.gpoint_[1]);
    cairo_line_to(cr, eetask_.goal_[0], eetask_.goal_[1]);
    cairo_stroke(cr);
    
    // ellbow repulsion
    if (repulse_ellbow_.isActive()) {
      cairo_set_source_rgb(cr, 0.4, 0.4, 1.0);
      cairo_set_line_width(cr, 1.0 / pixelsize);
      cairo_move_to(cr, repulse_ellbow_.gpoint_[0], repulse_ellbow_.gpoint_[1]);
      cairo_line_to(cr, repulse_ellbow_.gpoint_[0] + repulse_ellbow_.delta_[0] / repulse_ellbow_.gain_, repulse_ellbow_.gpoint_[1] + repulse_ellbow_.delta_[1] / repulse_ellbow_.gain_);
      cairo_stroke(cr);
    }
    
    // wrist repulsion
    if (repulse_wrist_.isActive()) {
      cairo_set_source_rgb(cr, 0.4, 0.4, 1.0);
      cairo_set_line_width(cr, 1.0 / pixelsize);
      cairo_move_to(cr, repulse_wrist_.gpoint_[0], repulse_wrist_.gpoint_[1]);
      cairo_line_to(cr, repulse_wrist_.gpoint_[0] + repulse_wrist_.delta_[0] / repulse_wrist_.gain_, repulse_wrist_.gpoint_[1] + repulse_wrist_.delta_[1] / repulse_wrist_.gain_);
      cairo_stroke(cr);
    }
    
    // base attraction
    if (attract_base_.isActive()) {
      cairo_set_source_rgb(cr, 0.4, 1.0, 0.4);
      cairo_set_line_width(cr, 1.0 / pixelsize);
      cairo_move_to(cr, attract_base_.gpoint_[0], attract_base_.gpoint_[1]);
      cairo_line_to(cr, attract_base_.gpoint_[0] + attract_base_.delta_[0] / attract_base_.gain_, attract_base_.gpoint_[1] + attract_base_.delta_[1] / attract_base_.gain_);
      cairo_stroke(cr);
    }
    
    cairo_restore(cr);
  }
  
  
  void setEEGoal(Vector const & goal)
  {
    eetask_.goal_ = goal;
  }
  
  
  void setBaseAttractor(Vector const & attractor)
  {
    attract_base_.attractor_ = attractor;
  }
  
  
  void setEllbowRepulsor(Vector const & repulsor)
  {
    repulse_ellbow_.repulsor_ = repulsor;
  }
  
  
  void setWristRepulsor(Vector const & repulsor)
  {
    repulse_wrist_.repulsor_ = repulsor;
  }
  
  
  void init(Vector const & position, Vector const & velocity)
  {
    robot_.update(position, velocity);
    
    for (size_t ii(0); ii < constraints_.size(); ++ii) {
      constraints_[ii];
    }
    for (size_t ii(0); ii < tasks_.size(); ++ii) {
      tasks_[ii]->init(robot_);
    }
    for (size_t ii(0); ii < objectives_.size(); ++ii) {
      objectives_[ii]->init(robot_);
    }
  }
  
  
  void update()
  {
    ostream * dbgos(0);
    if (verbose) {
      dbgos = &cout;
      cout << "==================================================\n"
	   << "Waypoint::update()\n";
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
			   dbgos, "");
    
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
			   dbgos, "");
    
    // The constraints cheat with the robot state: they directly work
    // on the positions (and velocities) that would have been achieved
    // without constraints. So, we then need to repair the position
    // and velocity to something consistent with the constraints, then
    // re-run the tasks and objectives.
    
    robot_.update(oldpos + dq_c, N_c * dq_c / timestep_);
    
    if (verbose) {
      print(dq_c, cout, "position correction to satisfy constraints", "  ");
      print(N_c, cout, "nullspace of constrains", "  ");
      print(robot_.getVelocity(), cout, "resulting constrained velocity", "  ");
      print(robot_.getPosition(), cout, "resulting constrained position", "  ");
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
			   dbgos, "");
    
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
    
    if (verbose) {
      print(qdd_res, cout, "constrained acceleration", "  ");
      print(qd_res, cout, "resulting constrained velocity", "  ");
      print(q_res, cout, "resulting constrained position", "  ");
    }
    
    robot_.update(q_res, qd_res);
  }
  
protected:
  double timestep_;
  Robot robot_;
  
  JointLimitConstraint joint_limits_;
  
  PositionControl eetask_;
  
  PointAttraction attract_base_;
  PointRepulsion repulse_ellbow_;
  PointRepulsion repulse_wrist_;
  PostureDamping joint_damping_;
  
  vector<Task *> constraints_;
  vector<Task *> tasks_;
  vector<Task *> objectives_;
};


// could templatize on Waypoint or pass in a factory or something along
// those lines to make it robot agnostic...
class Elastic
{
public:
  typedef list<Waypoint *> path_t;

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
    
    wpt_ = new Waypoint();
    wpt_->init(state, Vector::Zero(state.size()));
    
    eegoal <<
      1.0,
      dimy - 1.0,
      0.0;
    repulsor <<
      dimx - 1.0,
      dimy - 1.0,
      0.0;
    attractor <<
      1.0,
      1.0,
      0.0;
    
    wpt_->setEEGoal(eegoal);
    wpt_->setEllbowRepulsor(repulsor);
    wpt_->setWristRepulsor(repulsor);
    wpt_->setBaseAttractor(attractor);
    path_.push_back(wpt_);
  }
  
  
  void update()
  {
    if (verbose) {
      cout << "\n**************************************************\n";
    }
    for (path_t::iterator ii(path_.begin()); ii != path_.end(); ++ii) {
      (*ii)->setEEGoal(eegoal);
      (*ii)->setEllbowRepulsor(repulsor);
      (*ii)->setWristRepulsor(repulsor);
      (*ii)->setBaseAttractor(attractor);
      (*ii)->update();
    }
  }
  
  void draw(cairo_t * cr, double pixelsize)
  {
    for (path_t::reverse_iterator ii(path_.rbegin()); ii != path_.rend(); ++ii) {
      (*ii)->draw(cr, pixelsize);
    }
  }
  
private:
  path_t path_;
  Waypoint * wpt_;
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
  
  cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
  cairo_rectangle(cr, 0, 0, gw_width, gw_height);
  cairo_fill(cr);
  
  cairo_translate(cr, gw_x0, gw_y0);
  cairo_scale(cr, gw_sx, gw_sy);
  
  elastic.draw(cr, gw_sx);
  
  cairo_set_source_rgba(cr, 0.6, 0.0, 0.0, 0.5);
  cairo_arc(cr, eegoal[0], eegoal[1], grab_radius, 0., 2. * M_PI);
  cairo_fill(cr);
  
  cairo_set_source_rgba(cr, 0.0, 0.0, 0.6, 0.5);
  cairo_arc(cr, repulsor[0], repulsor[1], grab_radius, 0., 2. * M_PI);
  cairo_fill(cr);
  
  cairo_set_source_rgba(cr, 0.0, 0.6, 0.0, 0.5);
  cairo_arc(cr, attractor[0], attractor[1], grab_radius, 0., 2. * M_PI);
  cairo_fill(cr);
  
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
    for (Vector ** hh(handle); *hh != 0; ++hh) {
      Vector offset = **hh - point;
      if (offset.norm() <= grab_radius) {
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
    *grabbed = point + grab_offset;
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
  
  // btn = gtk_button_new_with_label("reset");
  // g_signal_connect(btn, "clicked", G_CALLBACK (cb_reset), NULL);
  // gtk_box_pack_start(GTK_BOX (hbox), btn, TRUE, TRUE, 0);
  // gtk_widget_show(btn);
  
  btn = gtk_button_new_with_label("play");
  g_signal_connect(btn, "clicked", G_CALLBACK (cb_play), NULL);
  gtk_box_pack_start(GTK_BOX (hbox), btn, TRUE, TRUE, 0);
  gtk_widget_show(btn);
  
  btn = gtk_button_new_with_label("next");
  g_signal_connect(btn, "clicked", G_CALLBACK (cb_next), NULL);
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
