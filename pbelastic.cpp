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


struct robot_s : public system_s {
  
  double const radius;
  double const len_a;
  double const len_b;
  
  Vector pos_a;
  Vector pos_b;
  

  // robot_s ()
  // : system_s(4),
  //   radius(0.5),
  //   len_a(0.8),
  //   len_b(0.6),
  //   pos_a(2),
  //   pos_b(2)
  // {
  //   state << 2.5, 5.7, 15 * deg, 32 * deg;
  //   update();
  // }
  
  explicit robot_s (Vector const & state_)
    : system_s(4),
      radius(0.5),
      len_a(0.8),
      len_b(0.6),
      pos_a(2),
      pos_b(2)
  {
    if (state_.size() != 4) {
      errx (EXIT_FAILURE, "only NDOF=4 allowed for now...");
    }
    state = state_;
    update();
  }
  
  void draw (cairo_t * cr, double pixelsize)
  {
    cairo_save (cr);
    
    cairo_set_source_rgba (cr, 0.7, 0.7, 0.7, 0.5);
    cairo_arc (cr, state[0], state[1], radius, 0., 2. * M_PI);
    cairo_fill (cr);
    
    cairo_set_source_rgb (cr, 0.2, 0.2, 0.2);
    cairo_set_line_width (cr, 3.0 / pixelsize);
    cairo_arc (cr, state[0], state[1], radius, 0., 2. * M_PI);
    cairo_stroke (cr);
    
    cairo_move_to (cr, state[0], state[1]);
    cairo_line_to (cr, pos_a[0], pos_a[1]);
    cairo_line_to (cr, pos_b[0], pos_b[1]);
    cairo_stroke (cr);
    
    cairo_restore (cr);
  }
  
  void update ()
  {
    pos_a <<
      state[0] + len_a * cos(state[2]),
      state[1] + len_a * sin(state[2]);
    pos_b <<
      len_b * cos(state[2] + state[3]),
      len_b * sin(state[2] + state[3]);
    pos_b += pos_a;
  }
};


// could templatize on robot_s or pass in a factory or something along
// those lines to make it robot agnostic...
class Elastic
{
public:
  typedef list<robot_s *> path_t;
  
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
    clear();
    
    if (0 == nsteps) {
      nsteps = 1;
    }
    Vector delta = (dest - start) / nsteps;
    Vector state = start;
    for (size_t ii(0); ii < nsteps; ++ii) {
      path_.push_back (new robot_s(state));
      state += delta;
    }
    path_.push_back (new robot_s(dest));    
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

////static robot_s robot;
static Elastic elastic;


static void update ()
{
  static size_t tick(0);
  cout << "tick " << tick++ << "\n";
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
