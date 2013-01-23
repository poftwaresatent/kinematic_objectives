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

#include "pbmockup.hpp"

#include <gtk/gtk.h>
#include <cmath>
#include <iostream>
#include <list>

#include <err.h>

using namespace pbmockup;


class Waypoint
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  explicit Waypoint (Vector const & state)
    : radius_(0.5),
      len_a_(0.8),
      len_b_(0.6),
      state_(state),
      pos_a_(2),
      pos_b_(2)
  {
    task_.push_back (task_s(4, 2, 0.1));
    task_.push_back (task_s(4, 2, 0.1));
    task_[1].Jacobian << 1, 0, 0, 0, 0, 1, 0, 0;
    setState (state);
    task_[0].desired = task_[0].current;
    task_[1].desired = task_[1].current;
  }
  
  void draw (cairo_t * cr, double pixelsize)
  {
    cairo_save (cr);
    
    cairo_set_source_rgba (cr, 0.7, 0.7, 0.7, 0.5);
    cairo_arc (cr, state_[0], state_[1], radius_, 0., 2. * M_PI);
    cairo_fill (cr);
    
    cairo_set_source_rgb (cr, 0.2, 0.2, 0.2);
    cairo_set_line_width (cr, 3.0 / pixelsize);
    cairo_arc (cr, state_[0], state_[1], radius_, 0., 2. * M_PI);
    cairo_stroke (cr);
    
    cairo_move_to (cr, state_[0], state_[1]);
    cairo_line_to (cr, pos_a_[0], pos_a_[1]);
    cairo_line_to (cr, pos_b_[0], pos_b_[1]);
    cairo_stroke (cr);
    
    cairo_set_source_rgb (cr, 1.0, 0.2, 0.2);
    cairo_set_line_width (cr, 1.0 / pixelsize);
    cairo_move_to (cr, task_[0].current[0], task_[0].current[1]);
    cairo_line_to (cr, task_[0].desired[0], task_[0].desired[1]);
    cairo_stroke (cr);
    
    cairo_set_source_rgb (cr, 0.2, 1.0, 0.2);
    cairo_move_to (cr, task_[1].current[0], task_[1].current[1]);
    cairo_line_to (cr, task_[1].desired[0], task_[1].desired[1]);
    cairo_stroke (cr);
    
    cairo_restore (cr);
  }
  
  void setEEGoal (Vector const & goal)
  {
    task_[0].desired = goal;
    cout << "setEEGoal (v): " << task_[0].desired[0] << "  " << task_[0].desired[1] << "\n";
  }
  
  void setEEGoal (double gx, double gy)
  {
    task_[0].desired << gx, gy;
    cout << "setEEGoal (dd): " << task_[0].desired[0] << "  " << task_[0].desired[1] << "\n";
  }
  
  void setBaseGoal (Vector const & goal)
  {
    task_[1].desired = goal;
    cout << "setBaseGoal (v): " << task_[1].desired[0] << "  " << task_[1].desired[1] << "\n";
  }
  
  void setBaseGoal (double gx, double gy)
  {
    task_[1].desired << gx, gy;
    cout << "setBaseGoal (dd): " << task_[1].desired[0] << "  " << task_[1].desired[1] << "\n";
  }
  
  void setState (Vector const & state)
  {
    if (state.size() != 4) {
      errx (EXIT_FAILURE, "only NDOF=4 allowed for now...");
    }
    state_ = state;
    
    double const ac2(len_a_ * cos(state_[2]));
    double const as2(len_a_ * sin(state_[2]));
    double const bc23(len_b_ * cos(state_[2] + state_[3]));
    double const bs23(len_b_ * sin(state_[2] + state_[3]));
    
    pos_a_ << state_[0] + ac2, state_[1] + as2;
    pos_b_ << pos_a_[0] + bc23, pos_a_[1] + bs23;
    
    task_[0].current = pos_b_;
    task_[0].Jacobian <<
      1, 0, -as2 - bs23, -bs23,
      0, 1,  ac2 + bc23,  bc23;
    
    task_[1].current << state_[0], state_[1];
  }
  
  Vector const & getState () const { return state_; }
  
  tasklist_t const getTasks () const { return task_; }
  
  
private:  
  double const radius_;
  double const len_a_;
  double const len_b_;
  
  Vector state_;
  Vector pos_a_;
  Vector pos_b_;
  
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
  
  void init (Vector const & start, Vector const & dest, size_t nsteps)
  {
    if (0 == nsteps) {
      nsteps = 1;
    }
    clear();
    
    Waypoint * start_wpt(new Waypoint(start));
    Vector start_eegoal(start_wpt->getTasks()[0].current);
    Vector start_basegoal(start_wpt->getTasks()[1].current);
    
    Waypoint * dest_wpt(new Waypoint(dest));
    Vector dest_eegoal(dest_wpt->getTasks()[0].current);
    Vector dest_basegoal(dest_wpt->getTasks()[1].current);
    
    Vector delta_eegoal = (dest_eegoal - start_eegoal) / nsteps;
    Vector delta_basegoal = (dest_basegoal - start_basegoal) / nsteps;
    
    Vector eegoal = start_eegoal + delta_eegoal;
    Vector basegoal = start_basegoal + delta_basegoal;
    
    Vector middle = (start + dest) / 2.0;
    
    path_.push_back (start_wpt);
    for (size_t ii(1); ii < nsteps; ++ii) {
      Waypoint * wpt(new Waypoint(middle));
      wpt->setEEGoal(eegoal);
      wpt->setBaseGoal(basegoal);
      path_.push_back (wpt);
      eegoal += delta_eegoal;
      basegoal += delta_basegoal;
    }
    path_.push_back (dest_wpt);
  }
  
  void update ()
  {
    for (path_t::iterator ii(path_.begin()); ii != path_.end(); ++ii) {
      Vector dq = recursive_task_priority_algorithm (4, (*ii)->getTasks());
      (*ii)->setState ((*ii)->getState() + dq);
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
};


static double const dimx(10.);
static double const dimy(8.);
static double const deg(M_PI / 180.);

static GtkWidget * gw(0);
static gint gw_width(800), gw_height(640);
static gint gw_sx, gw_sy, gw_x0, gw_y0;
static int play(0);

static Elastic elastic;


static void update ()
{
  static size_t tick(0);
  cout << "tick " << tick++ << "\n";
  elastic.update();
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
  gdouble const cx = (bb->x - gw_x0) / gw_sx;
  gdouble const cy = (bb->y - gw_y0) / gw_sy;
  
  cout << "click " << cx << "  " << cy << "\n";
  
  ////  gtk_widget_queue_draw (w_phi);
  
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
  gtk_widget_set_events (gw, GDK_BUTTON_PRESS_MASK);
  
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
  
  Vector start(4), dest(4);
  start << 1.0, dimy - 1.0, 80*deg, -40*deg;
  dest  << dimx -  1.0, 1.0, -90*deg, 60*deg;
  elastic.init (start, dest, 7);
  
  gtk_main ();
  
  return 0;
}
