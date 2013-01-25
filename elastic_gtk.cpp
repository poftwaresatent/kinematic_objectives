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

#include "mistry_algorithm.hpp"
#include "algorithm.hpp"
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
static Vector basegoal(2);
static Vector * handle[] = { &eegoal, &basegoal, 0 };
static Vector * grabbed(0);
static double grab_radius(0.2);
static Vector grab_offset(2);


static inline double bound (double lower, double value, double upper)
{
  if (value < lower) {
    value = lower;
  }
  else if (value > upper) {
    value = upper;
  }
  return value;
}


class Robot
  : public Model
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  explicit Robot (Vector const & state)
    : Model(5),
      radius_(0.5),
      len_a_(0.8),
      len_b_(0.6),
      len_c_(0.3),
      pos_a_(2),
      pos_b_(2),
      pos_c_(2)
  {
    // yet another subtlety: soft limits must not be too close to hard
    // limits, otherwise we get jitter from the joint-limit avoidance
    // algorithm.
    
    joint_limits_(3, 0) = -120.0 * deg;
    joint_limits_(3, 1) = -119.0 * deg;
    joint_limits_(3, 2) =  119.0 * deg;
    joint_limits_(3, 3) =  120.0 * deg;

    joint_limits_(4, 0) = -120.0 * deg;
    joint_limits_(4, 1) = -119.0 * deg;
    joint_limits_(4, 2) =  119.0 * deg;
    joint_limits_(4, 3) =  120.0 * deg;
    
    update(state);
  }
  
  void update (Vector const & state)
  {
    if (state.size() != 5) {
      errx (EXIT_FAILURE, "only NDOF=5 allowed for now...");
    }
    state_ = state;
    
    ac2_ = len_a_ * cos(state_[2]);
    as2_ = len_a_ * sin(state_[2]);
    q23_ = state_[2] + state_[3];
    bc23_ = len_b_ * cos(q23_);
    bs23_ = len_b_ * sin(q23_);
    q234_ = q23_ + state_[4];
    cc234_ = len_c_ * cos(q234_);
    cs234_ = len_c_ * sin(q234_);
    
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
  
  void draw (cairo_t * cr, double pixelsize)
  {
    cairo_save (cr);
    
    // translucent disk for base
    cairo_set_source_rgba (cr, 0.7, 0.7, 0.7, 0.5);
    cairo_arc (cr, state_[0], state_[1], radius_, 0., 2. * M_PI);
    cairo_fill (cr);
    
    // thick circle outline for base
    cairo_set_source_rgb (cr, 0.2, 0.2, 0.2);
    cairo_set_line_width (cr, 3.0 / pixelsize);
    cairo_arc (cr, state_[0], state_[1], radius_, 0., 2. * M_PI);
    cairo_stroke (cr);
    
    // thin arcs for arm joint limits
    cairo_set_source_rgb (cr, 0.5, 0.5, 1.0);
    cairo_set_line_width (cr, 1.0 / pixelsize);
    cairo_move_to (cr, pos_a_[0], pos_a_[1]);
    cairo_arc (cr, pos_a_[0], pos_a_[1], 0.1,
	       bound(-2.0*M_PI, state_[2] + joint_limits_(3, 0), 2.0*M_PI),
	       bound(-2.0*M_PI, state_[2] + joint_limits_(3, 3), 2.0*M_PI));
    cairo_line_to (cr, pos_a_[0], pos_a_[1]);
    cairo_stroke (cr);
    cairo_move_to (cr, pos_b_[0], pos_b_[1]);
    cairo_arc (cr, pos_b_[0], pos_b_[1], 0.1,
	       bound(-2.0*M_PI, q23_ + joint_limits_(4, 0), 2.0*M_PI),
	       bound(-2.0*M_PI, q23_ + joint_limits_(4, 3), 2.0*M_PI));
    cairo_line_to (cr, pos_b_[0], pos_b_[1]);
    cairo_stroke (cr);
    
    // thick line for arms
    cairo_set_source_rgb (cr, 0.2, 0.2, 0.2);
    cairo_set_line_width (cr, 3.0 / pixelsize);
    cairo_move_to (cr, state_[0], state_[1]);
    cairo_line_to (cr, pos_a_[0], pos_a_[1]);
    cairo_line_to (cr, pos_b_[0], pos_b_[1]);
    cairo_line_to (cr, pos_c_[0], pos_c_[1]);
    cairo_stroke (cr);
    
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
  
  explicit Waypoint (Vector const & state)
    : model_(state),
      eetask_(5, 2, 0.1),
      basetask_(5, 2, 0.1)
  {
    task_.push_back(&eetask_);
    task_.push_back(&basetask_);
    
    basetask_.Jx <<
      1, 0, 0, 0, 0,
      0, 1, 0, 0, 0;
    
    update(state);
    for (size_t ii(0); ii < task_.size(); ++ii) {
      task_[ii]->xdes = task_[ii]->xcur;
      task_[ii]->dx = Vector::Zero(4);
    }
  }
  
  void draw (cairo_t * cr, double pixelsize)
  {
    model_.draw(cr, pixelsize);
    
    cairo_save (cr);
    
    // thin line for end effector task
    cairo_set_source_rgb (cr, 1.0, 0.4, 0.4);
    cairo_set_line_width (cr, 1.0 / pixelsize);
    cairo_move_to (cr, eetask_.xcur[0], eetask_.xcur[1]);
    cairo_line_to (cr, eetask_.xdes[0], eetask_.xdes[1]);
    cairo_stroke (cr);
    
    // thin line for base task
    cairo_set_source_rgb (cr, 0.4, 1.0, 0.4);
    cairo_move_to (cr, basetask_.xcur[0], basetask_.xcur[1]);
    cairo_line_to (cr, basetask_.xdes[0], basetask_.xdes[1]);
    cairo_stroke (cr);
    
    cairo_restore (cr);
  }
  
  void setEEGoal (Vector const & goal)
  {
    eetask_.xdes = goal;
    eetask_.dx = eetask_.xdes - eetask_.xcur;
  }
  
  void setEEGoal (double gx, double gy)
  {
    eetask_.xdes << gx, gy;
    eetask_.dx = eetask_.xdes - eetask_.xcur;
  }
  
  void setBaseGoal (Vector const & goal)
  {
    basetask_.xdes = goal;
    basetask_.dx = basetask_.xdes - basetask_.xcur;
  }
  
  void setBaseGoal (double gx, double gy)
  {
    basetask_.xdes << gx, gy;
    basetask_.dx = basetask_.xdes - basetask_.xcur;
  }
  
  void update (Vector const & state)
  {
    model_.update(state);
    
    eetask_.xcur = model_.pos_c_;
    eetask_.dx = eetask_.xdes - eetask_.xcur;
    eetask_.Jx <<
      1,
      0,
      -model_.as2_ - model_.bs23_ - model_.cs234_,
      -model_.bs23_ - model_.cs234_,
      -model_.cs234_,
      0,
      1,
      model_.ac2_ + model_.bc23_ + model_.cc234_,
      model_.bc23_ + model_.cc234_,
      model_.cc234_;
    
    basetask_.xcur << state[0], state[1];
    basetask_.dx = basetask_.xdes - basetask_.xcur;
  }
  
  Model const & getModel () const { return model_; }
  
  Vector const & getState () const { return model_.state_; }
  
  tasklist_t const getTasks () const { return task_; }
  
  
private:  
  Robot model_;
  task_s eetask_;
  task_s basetask_;
  tasklist_t task_;
};


// could templatize on Waypoint or pass in a factory or something along
// those lines to make it robot agnostic...
class Elastic
{
public:
  typedef list<Waypoint *> path_t;

  ~Elastic ()
  {
    clear();
  }
  
  void clear ()
  {
    for (path_t::iterator ii(path_.begin()); ii != path_.end(); ++ii) {
      delete *ii;
    }
    path_.clear();
  }
  
  void init ()
  {
    clear();
    
    eegoal <<
      1.0,
      dimy - 1.0;
    basegoal <<
      dimx - 1.0, 1.0;
    
    Vector posture(5);
    posture <<
      dimx / 2.0,
      dimy / 2.0,
      80.0 * deg,
      - 40.0 * deg,
      25.0 * deg;
    
    wpt_ = new Waypoint(posture);
    wpt_->setEEGoal(eegoal);
    wpt_->setBaseGoal(basegoal);
    path_.push_back (wpt_);
  }
  
  void update ()
  {
    ostream * dbgos(0);
    if (verbose) {
      dbgos = &cout;
      cout << "\n**************************************************\n";
    }
    for (path_t::iterator ii(path_.begin()); ii != path_.end(); ++ii) {
      if (verbose) {
	cout << "--------------------------------------------------\n";
      }
      (*ii)->setEEGoal(eegoal);
      (*ii)->setBaseGoal(basegoal);
      Vector dq(algorithm((*ii)->getModel(), (*ii)->getState(), (*ii)->getTasks(), dbgos, "  "));
      (*ii)->update ((*ii)->getState() + dq);
    }
  }
  
  void draw (cairo_t * cr, double pixelsize)
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


static void update ()
{
  elastic.update ();
  gtk_widget_queue_draw (gw);
}


static void cb_play (GtkWidget * ww, gpointer data)
{
  if (play) {
    play = 0;
  }
  else {
    play = 1;
  }
}


static void cb_next (GtkWidget * ww, gpointer data)
{
  if (play) {
    play = 0;
  }
  else {
    update ();
  }    
}


static void cb_quit (GtkWidget * ww, gpointer data)
{
  gtk_main_quit();
}


static gint cb_expose (GtkWidget * ww,
		       GdkEventExpose * ee,
		       gpointer data)
{
  cairo_t * cr = gdk_cairo_create (ee->window);
  
  cairo_set_source_rgb (cr, 1.0, 1.0, 1.0);
  cairo_rectangle (cr, 0, 0, gw_width, gw_height);
  cairo_fill (cr);
  
  cairo_translate (cr, gw_x0, gw_y0);
  cairo_scale (cr, gw_sx, gw_sy);
  
  elastic.draw(cr, gw_sx);
  
  cairo_set_source_rgba (cr, 0.6, 0.0, 0.0, 0.5);
  cairo_arc (cr, eegoal[0], eegoal[1], grab_radius, 0., 2. * M_PI);
  cairo_fill (cr);
  
  cairo_set_source_rgba (cr, 0.0, 0.6, 0.0, 0.5);
  cairo_arc (cr, basegoal[0], basegoal[1], grab_radius, 0., 2. * M_PI);
  cairo_fill (cr);
  
  cairo_destroy (cr);
  
  return TRUE;
}


static gint cb_size_allocate (GtkWidget * ww,
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


static gint cb_click (GtkWidget * ww,
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


static gint cb_motion (GtkWidget * ww,
		       GdkEventMotion * ee)
{
  int mx, my;
  GdkModifierType modifier;
  gdk_window_get_pointer (ww->window, &mx, &my, &modifier);
  
  if (0 != grabbed) {
    Vector point(2);
    point << (mx - gw_x0) / (double) gw_sx, (my - gw_y0) / (double) gw_sy;
    *grabbed = point + grab_offset;
    gtk_widget_queue_draw (gw);
  }
  
  return TRUE;
}


static gint idle (gpointer data)
{
  if (play) {
    update ();
  }
  return TRUE;
}


static void init_gui (int * argc, char *** argv)
{
  GtkWidget *window, *vbox, *hbox, *btn;
  
  gtk_init (argc, argv);
  
  window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
  
  vbox = gtk_vbox_new (FALSE, 0);
  gtk_container_add (GTK_CONTAINER (window), vbox);
  gtk_widget_show (vbox);
  
  gw = gtk_drawing_area_new ();
  g_signal_connect (gw, "expose_event", G_CALLBACK (cb_expose), NULL);
  g_signal_connect (gw, "size_allocate", G_CALLBACK (cb_size_allocate), NULL);
  g_signal_connect (gw, "button_press_event", G_CALLBACK (cb_click), NULL);
  g_signal_connect (gw, "motion_notify_event", G_CALLBACK (cb_motion), NULL);
  gtk_widget_set_events (gw,
			 GDK_BUTTON_PRESS_MASK |
			 GDK_BUTTON_RELEASE_MASK |
			 GDK_BUTTON_MOTION_MASK);
  
  gtk_widget_show (gw);
  
  gtk_widget_set_size_request (gw, gw_width, gw_height);
  gtk_box_pack_start (GTK_BOX (vbox), gw, TRUE, TRUE, 0);
  
  hbox = gtk_hbox_new (TRUE, 3);
  gtk_box_pack_start (GTK_BOX (vbox), hbox, FALSE, TRUE, 0);
  gtk_widget_show (hbox);
  
  // btn = gtk_button_new_with_label ("reset");
  // g_signal_connect (btn, "clicked", G_CALLBACK (cb_reset), NULL);
  // gtk_box_pack_start (GTK_BOX (hbox), btn, TRUE, TRUE, 0);
  // gtk_widget_show (btn);
  
  btn = gtk_button_new_with_label ("play");
  g_signal_connect (btn, "clicked", G_CALLBACK (cb_play), NULL);
  gtk_box_pack_start (GTK_BOX (hbox), btn, TRUE, TRUE, 0);
  gtk_widget_show (btn);
  
  btn = gtk_button_new_with_label ("next");
  g_signal_connect (btn, "clicked", G_CALLBACK (cb_next), NULL);
  gtk_box_pack_start (GTK_BOX (hbox), btn, TRUE, TRUE, 0);
  gtk_widget_show (btn);
  
  btn = gtk_button_new_with_label ("quit");
  g_signal_connect (btn, "clicked", G_CALLBACK (cb_quit), NULL);
  gtk_box_pack_start (GTK_BOX (hbox), btn, TRUE, TRUE, 0);
  gtk_widget_show (btn);
  
  gtk_idle_add (idle, 0);
  
  gtk_widget_show (window);
}


int main (int argc, char ** argv)
{
  init_gui (&argc, &argv);
  
  if ((argc > 1) && (0 == strcmp("-v", argv[1]))) {
    verbose = true;
  }
  
  elastic.init ();
  
  gtk_main ();
  
  return 0;
}
