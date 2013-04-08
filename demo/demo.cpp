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

#include <kinematic_objectives/util.h>
#include <kinematic_objectives/print.h>
#include <kinematic_objectives/constraint_teleporting_blender.h>
#include "interactive_blender.h"
#include <gtk/gtk.h>
#include <cmath>
#include <iostream>
#include <list>
#include <err.h>

using namespace kinematic_objectives::demo;
using namespace kinematic_objectives;


static double const dimx(10.);
static double const dimy(8.);
static double const lwscale(1.0);

static GtkWidget * gw(0);
static gint gw_width(800), gw_height(640);
static gint gw_sx, gw_sy, gw_x0, gw_y0;

static bool verbose(false);
static int play(0);

static InteractiveBlender * blender;
static InteractionHandle * grabbed(0);
static Vector grab_offset(3);



static void dump_achievability(string const & type, size_t index, Objective const * obj)
{
  cout << type << " #" << index << " \"" << obj->name_ << "\"";
  if (obj->isActive()) {
    PseudoInverseFeedback const & fb(obj->jbar_svd_);
    if (0 == fb.truncated_range) {
      cout << " zero range (original " << fb.original_range << ")\n";
    }
    else {
      cout << "\n  range " << fb.truncated_range << " (original " << fb.original_range << ")\n";
      print(fb.singular_values, cout, "singular values", "  ");
      print(fb.input_space, cout, "input_space", "  ");
      print(fb.output_space, cout, "output_space", "  ");
    }
  }
  else {
    cout << " inactive\n";
  }
  fflush(stdout);
}


static void update()
{
  blender->update();
  gtk_widget_queue_draw(gw);
  
  CompoundObjective const & co(blender->robot_);
  for (size_t ii(0); ii < co.constraints_.size(); ++ii) {
    dump_achievability("constraint", ii, co.constraints_[ii]);
  }
  for (size_t ii(0); ii < co.hard_objectives_.size(); ++ii) {
    dump_achievability("hard_objective", ii, co.hard_objectives_[ii]);
  }
  for (size_t ii(0); ii < co.soft_objectives_.size(); ++ii) {
    dump_achievability("soft_objective", ii, co.soft_objectives_[ii]);
  }
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
  //  for (Blender::path_t::iterator ii(blender->path_.begin()); ii != blender->path_.end(); ++ii) {
  BaseCompoundObjective * wpt(&blender->robot_);//dynamic_cast<BaseCompoundObjective*>(*ii));
  // if ( ! wpt) {
  //   continue;
  // }
  if (verbose) {
    if (fabs(wpt->robot_.position_[2]) > M_PI) {
      cout << "normalize " << wpt->robot_.position_[2]
	   << " to " << normangle(wpt->robot_.position_[2]) << "\n";
    }
  }
  wpt->robot_.position_[2] = normangle(wpt->robot_.position_[2]);
  wpt->robot_.update(wpt->robot_.position_, wpt->robot_.velocity_);
  //  }
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
  
  blender->draw(cr, lwscale, gw_sx);
  
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
    for (size_t ii(0); ii < blender->handles_.size(); ++ii) {
      Vector offset = blender->handles_[ii]->point_ - point;
      if (offset.norm() <= blender->handles_[ii]->radius_) {
    	grab_offset = offset;
    	grabbed = blender->handles_[ii];
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
  
  Blender * blender_imp; // XXXX make this runtime configurable (timestep is another good candidate for that)
  if ((argc > 1) && (0 == strcmp("-v", argv[1]))) {
    verbose = true;
    blender_imp = new ConstraintTeleportingBlender(1.0e-2, &cout, "ctb  ");
    blender = new InteractiveBlender(blender_imp, &cout, "");
  }
  else {
    blender_imp = new ConstraintTeleportingBlender(1.0e-2, 0, "");
    blender = new InteractiveBlender(blender_imp, 0, "");
  }
  
  blender->ee_.point_         <<         1.0, dimy / 2.0      ,     0.0;
  blender->ee_ori_.point_     <<         2.0, dimy / 2.0 + 1.0,     0.0;
  blender->base_.point_       <<         1.0,              1.0,     0.0;
  blender->repulsor_.point_   <<  dimx / 2.0,              1.0,     0.0;
  blender->obstacle_.point_   <<  dimx / 2.0,       dimy - 1.0,     0.0;
  //  blender->eegoal_.point_     <<  dimx - 1.0, dimy / 2.0      ,     0.0;
  //  blender->basegoal_.point_   <<  dimx - 1.0,              1.0,     0.0;
  
  Vector posture(5);
  posture <<
    dimx / 2.0,
    dimy / 2.0,
    80.0 * deg,
    - 40.0 * deg,
    25.0 * deg;
  blender->init(posture);
  
  gtk_main();
  
  delete blender;
  delete blender_imp;
  
  return 0;
}
