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
#include "joint_limits.hpp"
#include "print.hpp"
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
static Vector eegoal(2);
static Vector ellbowgoal(2);
static Vector basegoal(2);
static Vector * handle[] = { &eegoal, &ellbowgoal, &basegoal, 0 };
static Vector * grabbed(0);
static double grab_radius(0.2);
static Vector grab_offset(2);


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


class PositionTask
  : public Task
{
public:
  PositionTask(size_t node,
	       Vector const & point)
    : node_(node),
      point_(point)
  {
  }
  
  
  virtual bool init(Model const & model)
  {
    step_hint_ = 0.1;
    Eigen::Vector3d tmp1, tmp2;
    tmp1 << point_[0], point_[1], 0.0;
    tmp2 = model.frame(node_) * tmp1;
    gpoint_.resize(2);
    gpoint_ << tmp2[0], tmp2[1];
    goal_ = gpoint_;
    delta_ = Vector::Zero(point_.size());
    Matrix tmp3(model.computeJx(node_, gpoint_));
    Jacobian_ = tmp3.block(0, 0, 2, tmp3.cols());
    return true;
  }
  
  
  virtual bool update(Model const & model)
  {
    Eigen::Vector3d tmp1, tmp2;
    tmp1 << point_[0], point_[1], 0.0;
    tmp2 = model.frame(node_) * tmp1;
    gpoint_ << tmp2[0], tmp2[1];
    delta_ = goal_ - gpoint_;
    Matrix tmp3(model.computeJx(node_, gpoint_));
    Jacobian_ = tmp3.block(0, 0, 2, tmp3.cols());
    return true;
  }
  
  size_t node_;
  Vector point_;
  Vector gpoint_;
  Vector goal_;
};

  
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
      pos_a_(2),
      pos_b_(2),
      pos_c_(2)
  {
  }
  
  
  virtual Vector const & getState() const
  {
    return state_;
  }
  
  
  virtual Transform frame(size_t node) const
  {
    Transform tf(Transform::Identity());
    switch (node) {
    case 0:
      tf.translation() << state_[0], state_[1], 0.0;
      break;
    case 1:
      tf.translation() << state_[0], state_[1], 0.0;
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
  
  
  virtual Matrix computeJx(size_t node, Vector const & gpoint) const
  {
    Matrix Jx(Matrix::Zero(3, 5));
    Vector delta;
    switch (node) {
    case 3:
      delta = gpoint - pos_b_;
      Jx.block(0, 4, 3, 1) << -delta[1], delta[0], 1.0;
    case 2:
      delta = gpoint - pos_a_;
      Jx.block(0, 3, 3, 1) << -delta[1], delta[0], 1.0;
    case 1:
      delta = gpoint - state_.block(0, 0, 2, 1);
      Jx.block(0, 2, 3, 1) << -delta[1], delta[0], 1.0;
    case 0:
      Jx.block(0, 0, 2, 2) = Matrix::Identity(2, 2);
      break;
    default:
      errx (EXIT_FAILURE, "Robot::computeJx() called on invalid node %zu", node);
    }
    if (verbose) {
      cout << "Robot::computeJx(" << node << ", [" << gpoint[0] << "  " << gpoint[1] << "])\n";
      print(Jx, cout, "", "  ");
    }
    return Jx;
  }
  
  
  virtual void update(Vector const & state)
  {
    if (state.size() != 5) {
      errx (EXIT_FAILURE, "Robot::update(): state has %zu DOF (but needs 5)", (size_t) state.size());
    }
    state_ = state;
    
    c2_ = cos(state_[2]);
    s2_ = sin(state_[2]);
    ac2_ = len_a_ * c2_;
    as2_ = len_a_ * s2_;
    
    q23_ = state_[2] + state_[3];
    c23_ = cos(q23_);
    s23_ = sin(q23_);
    bc23_ = len_b_ * c23_;
    bs23_ = len_b_ * s23_;
    
    q234_ = q23_ + state_[4];
    c234_ = cos(q234_);
    s234_ = sin(q234_);
    cc234_ = len_c_ * c234_;
    cs234_ = len_c_ * s234_;
    
    pos_a_ <<
      state_[0] + ac2_,
      state_[1] + as2_;
    pos_b_ <<
      pos_a_[0] + bc23_,
      pos_a_[1] + bs23_;
    pos_c_ <<
      pos_b_[0] + cc234_,
      pos_b_[1] + cs234_;
  }
  
  
  void draw(cairo_t * cr, double pixelsize) const
  {
    cairo_save(cr);
    
    // translucent disk for base
    cairo_set_source_rgba(cr, 0.7, 0.7, 0.7, 0.5);
    cairo_arc(cr, state_[0], state_[1], radius_, 0., 2. * M_PI);
    cairo_fill(cr);
    
    // thick circle outline for base
    cairo_set_source_rgb(cr, 0.2, 0.2, 0.2);
    cairo_set_line_width(cr, 3.0 / pixelsize);
    cairo_arc(cr, state_[0], state_[1], radius_, 0., 2. * M_PI);
    cairo_stroke(cr);
    
    // thick line for arms
    cairo_set_source_rgb(cr, 0.2, 0.2, 0.2);
    cairo_set_line_width(cr, 3.0 / pixelsize);
    cairo_move_to(cr, state_[0], state_[1]);
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
  
  Vector state_;
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
    : eetask_(3, Vector::Zero(2)),
      ellbowtask_(1, Vector::Zero(2)),
      basetask_(0, Vector::Zero(2))
  {
    eetask_.point_ << robot_.len_c_, 0.0;
    ellbowtask_.point_ << robot_.len_a_, 0.0;
    tasks_.push_back(&eetask_);
    tasks_.push_back(&ellbowtask_);
    tasks_.push_back(&basetask_);
    
    joint_limits_.init(5);
    
    // yet another subtlety: soft limits must not be too close to hard
    // limits, otherwise we get jitter from the joint-limit avoidance
    // algorithm.
    
    joint_limits_.limits_(3, 0) = -120.0 * deg;
    joint_limits_.limits_(3, 1) = -119.0 * deg;
    joint_limits_.limits_(3, 2) =  119.0 * deg;
    joint_limits_.limits_(3, 3) =  120.0 * deg;
    
    joint_limits_.limits_(4, 0) = -120.0 * deg;
    joint_limits_.limits_(4, 1) = -119.0 * deg;
    joint_limits_.limits_(4, 2) =  119.0 * deg;
    joint_limits_.limits_(4, 3) =  120.0 * deg;
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
	      bound(-2.0*M_PI, robot_.state_[2] + joint_limits_.limits_(3, 0), 2.0*M_PI),
	      bound(-2.0*M_PI, robot_.state_[2] + joint_limits_.limits_(3, 3), 2.0*M_PI));
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
    
    // thin line for ellbow task
    cairo_set_source_rgb(cr, 0.4, 0.4, 1.0);
    cairo_set_line_width(cr, 1.0 / pixelsize);
    cairo_move_to(cr, ellbowtask_.gpoint_[0], ellbowtask_.gpoint_[1]);
    cairo_line_to(cr, ellbowtask_.goal_[0], ellbowtask_.goal_[1]);
    cairo_stroke(cr);
    
    // thin line for base task
    cairo_set_source_rgb(cr, 0.4, 1.0, 0.4);
    cairo_move_to(cr, basetask_.gpoint_[0], basetask_.gpoint_[1]);
    cairo_line_to(cr, basetask_.goal_[0], basetask_.goal_[1]);
    cairo_stroke(cr);
    
    // // thin line for laser task
    // cairo_set_source_rgb(cr, 0.4, 0.4, 1.0);
    // cairo_move_to(cr, model_.pos_c_[0], model_.pos_c_[1]);
    // double const dx(3.0 * cos(lasertask_.goal_[0]));
    // double const dy(3.0 * sin(lasertask_.goal_[0]));
    // cairo_line_to(cr, model_.pos_c_[0] + dx, model_.pos_c_[1] + dy);
    // cairo_stroke(cr);
    
    cairo_restore(cr);
  }
  
  
  void setEEGoal(Vector const & goal)
  {
    eetask_.goal_ = goal;
  }
  
  
  void setEllbowGoal(Vector const & goal)
  {
    ellbowtask_.goal_ = goal;
  }
  
  
  void setBaseGoal(Vector const & goal)
  {
    basetask_.goal_ = goal;
  }
  
  
  // void setLaserGoal(Vector const & goalpoint)
  // {
  //   Vector dg(goalpoint - model_.pos_c_);
  //   lasertask_.goal_ << atan2(dg[1], dg[0]);
  //   lasertask_.dx << normangle(lasertask_.goal_[0] - lasertask_.xcur[0]);
  // }
  
  
  bool init(Vector const & state)
  {
    robot_.update(state);
    
    for (size_t ii(0); ii < tasks_.size(); ++ii) {
      if ( ! ((Task*)tasks_[ii])->init(robot_)) {
	cerr << "Waypoint::init(): tasks_[" << ii << "]->init() failed\n";
	return false;
      }
    }
    next_state_ = state;
    
    return true;
  }
  
  
  bool update()
  {
    ostream * dbgos(0);
    if (verbose) {
      dbgos = &cout;
      cout << "--------------------------------------------------\n"
	   << "Waypoint::update()\n";
    }
    
    robot_.update(next_state_);
    
    for (size_t ii(0); ii < tasks_.size(); ++ii) {
      if ( ! ((Task*)tasks_[ii])->update(robot_)) {
	cerr << "Waypoint::update(): tasks_[" << ii << "]->update() failed\n";
	return false;
      }
    }
    
    Vector const dq(algorithm2(joint_limits_,
			       next_state_, // that's now the current state, btw
			       tasks_,
			       dbgos,
			       "  "));
    next_state_ += dq;
    
    return true;
    
    // eetask_.Jx <<
    //   1,
    //   0,
    //   -model_.as2_ - model_.bs23_ - model_.cs234_,
    //   -model_.bs23_ - model_.cs234_,
    //   -model_.cs234_,
    //   0,
    //   1,
    //   model_.ac2_ + model_.bc23_ + model_.cc234_,
    //   model_.bc23_ + model_.cc234_,
    //   model_.cc234_;
    
    // lasertask_.xcur << model_.q234_;
    // lasertask_.dx << normangle(lasertask_.goal_[0] - lasertask_.xcur[0]);
    // // strictly speaking, the Jacobian of the laser task also has
    // // entries for the base... but for now just treat it as a moving
    // // goal when the base moves, even if the laser target point has
    // // not moved.
    // lasertask_.Jx <<
    //   0, 0, 1, 1, 1;
  }
  
protected:
  Robot robot_;
  Vector next_state_;
  
  JointLimits joint_limits_;
  
private:  
  PositionTask eetask_;
  PositionTask ellbowtask_;
  PositionTask basetask_;
  
  // Don't you wish C++ templates allowed polymorphism on the
  // collection level based on polymorphism at the item level?
  vector<TaskData*> tasks_;
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
  
  
  bool init(Vector const & state)
  {
    clear();
    
    wpt_ = new Waypoint();
    if ( ! wpt_->init(state)) {
      delete wpt_;
      wpt_ = 0;
      return false;
    }
    
    eegoal <<
      1.0,
      dimy - 1.0;
    ellbowgoal <<
      dimx * 0.5,
      dimy * 0.5;
    basegoal <<
      dimx - 1.0,
      1.0;
    
    wpt_->setEEGoal(eegoal);
    wpt_->setEllbowGoal(ellbowgoal);
    wpt_->setBaseGoal(basegoal);
    path_.push_back(wpt_);
    
    return true;
  }
  
  
  bool update()
  {
    if (verbose) {
      cout << "\n**************************************************\n";
    }
    for (path_t::iterator ii(path_.begin()); ii != path_.end(); ++ii) {
      (*ii)->setEEGoal(eegoal);
      (*ii)->setEllbowGoal(ellbowgoal);
      (*ii)->setBaseGoal(basegoal);
      if ( ! (*ii)->update()) {
	return false;
      }
    }
    return true;
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
  cairo_arc(cr, ellbowgoal[0], ellbowgoal[1], grab_radius, 0., 2. * M_PI);
  cairo_fill(cr);
  
  cairo_set_source_rgba(cr, 0.0, 0.6, 0.0, 0.5);
  cairo_arc(cr, basegoal[0], basegoal[1], grab_radius, 0., 2. * M_PI);
  cairo_fill(cr);
  
  // cairo_set_source_rgba(cr, 0.0, 0.0, 0.6, 0.5);
  // cairo_arc(cr, lasergoal[0], lasergoal[1], grab_radius, 0., 2. * M_PI);
  // cairo_fill(cr);
  
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
    Vector point(2);
    point << (bb->x - gw_x0) / (double) gw_sx, (bb->y - gw_y0) / (double) gw_sy;
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
    Vector point(2);
    point << (mx - gw_x0) / (double) gw_sx, (my - gw_y0) / (double) gw_sy;
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
