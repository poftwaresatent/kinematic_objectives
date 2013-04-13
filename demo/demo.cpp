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

#include "interactive_compound_objectives.h"
#include <kinematic_objectives/util.h>
#include <kinematic_objectives/print.h>
#include <kinematic_objectives/unconstrained_blender.h>
#include <kinematic_objectives/constraint_teleporting_blender.h>
#include <kinematic_objectives/constraint_bouncing_blender.h>
#include <kinematic_objectives/achievability.h>
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

static PlanarRobot planar_robot;
static InteractiveCompoundObjective * interactive_compound(0);
static Blender * blender(0);
static InteractionHandle * grabbed(0);
static Vector grab_offset(3);



static void update()
{
  interactive_compound->update();
  blender->update(planar_robot, &interactive_compound->compound_objective_);
  gtk_widget_queue_draw(gw);

  vector<Achievability> info;
  Achievability::compute(planar_robot, interactive_compound->compound_objective_, info);
  Achievability::print(info, cout, "");
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
  for (size_t ii(2); ii < 5; ++ii) {
    if (verbose) {
      if (fabs(planar_robot.position_[ii]) > M_PI) {
	cout << "normalize " << planar_robot.position_[ii]
	     << " to " << normangle(planar_robot.position_[ii]) << "\n";
      }
    }
    planar_robot.position_[ii] = normangle(planar_robot.position_[ii]);
    planar_robot.update(planar_robot.position_, planar_robot.velocity_);
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
  
  planar_robot.draw(cr, lwscale, gw_sx);
  interactive_compound->draw(cr, lwscale, gw_sx);
  
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
    for (size_t ii(0); ii < interactive_compound->handles_.size(); ++ii) {
      Vector offset = interactive_compound->handles_[ii]->point_ - point;
      if (offset.norm() <= interactive_compound->handles_[ii]->radius_) {
    	grab_offset = offset;
    	grabbed = interactive_compound->handles_[ii];
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


void parse_options(int argc, char ** argv)
{
  printf("\n"
	 "  Demo for kinematic constrained optimization.\n"
	 "\n"
	 "  Copyright (c) 2013, Willow Garage, Inc. All rights reserved.\n"
	 "                      Released under the 3-clause BSD licence.\n"
	 "                      Written by Roland Philippsen.\n"
	 "\n");
  
  string opt_blender("teleporting");
  string opt_compound("eegoal");
  double opt_timestep(10.0);	// milliseconds
  
  for (int iopt(1); iopt < argc; ++iopt) {
    
    if (0 == strcmp("-h", argv[iopt])) {
      printf("usage [-vh] [-b blender] [-c compound] [-t milliseconds]\n"
	     "\n"
	     "  -h             print this message\n"
	     "  -v             enable verbose messages\n"
	     "  -b  blender    name of the blender class to use (default '%s')\n"
	     "                 specify an invalid name to obtain a list\n"
	     "  -c  compound   name of the objectives compound to run (default '%s')\n"
	     "                 specify an invalid name to obtain a list\n"
	     "  -t  ms         blender integration timestep in milliseconds (default %3g)\n"
	     // "  -x  dimx              x-dimension of the world, in meters (default 10)\n"
	     // "  -y  dimy              y-dimension of the world, in meters (default 8)\n"
	     , opt_blender.c_str(), opt_compound.c_str(), opt_timestep);
      exit(EXIT_SUCCESS);
    }
    
    else if (0 == strcmp("-v", argv[iopt])) {
      verbose = true;
    }
    
    else if (0 == strcmp("-b", argv[iopt])) {
      if (++iopt >= argc) {
	errx(EXIT_FAILURE, "%s requires an argument (use -h for help)", argv[iopt-1]);
      }
      opt_blender = argv[iopt];
    }
    
    else if (0 == strcmp("-c", argv[iopt])) {
      if (++iopt >= argc) {
	errx(EXIT_FAILURE, "%s requires an argument (use -h for help)", argv[iopt-1]);
      }
      opt_compound = argv[iopt];
    }
    
    else if (0 == strcmp("-t", argv[iopt])) {
      if (++iopt >= argc) {
	errx(EXIT_FAILURE, "%s requires an argument (use -h for help)", argv[iopt-1]);
      }
      if (1 != sscanf(argv[iopt], "%lf", &opt_timestep)) {
	err(EXIT_FAILURE, "sscanf('%s'...)", argv[iopt]);
      }
    }
    
    else {
      errx(EXIT_FAILURE, "invalid option %s (use -h for help)", argv[iopt]);
    }
    
  }
  
  if ("teleporting" == opt_blender) {
    blender = new ConstraintTeleportingBlender(opt_timestep * 1e-3, verbose ? &cout : 0, "ctb  ");
  }
  else if ("unconstrained" == opt_blender) {
    blender = new UnconstrainedBlender(opt_timestep * 1e-3);
  }
  else if ("bouncing" == opt_blender) {
    blender = new ConstraintBouncingBlender(opt_timestep * 1e-3, 1e-2);
  }
  else {
    errx(EXIT_FAILURE, "invalid blender '%s' (have: teleporting, unconstrained, bouncing)", opt_blender.c_str());
  }
  
  if ("eegoal" == opt_compound) {
    interactive_compound = new EEGoalCompoundObjective(planar_robot);
  }
  else if ("elastic" == opt_compound) {
    interactive_compound = new ElasticLinksCompoundObjective(planar_robot);
  }
  else {
    errx(EXIT_FAILURE, "invalid compound '%s' (use 'eegoal' or 'elastic')", opt_compound.c_str());
  }
  interactive_compound->init(dimx, dimy);
}


void cleanup()
{
  delete interactive_compound;
  delete blender;
}


int main(int argc, char ** argv)
{
  if (0 != atexit(cleanup)) {
    err(EXIT_FAILURE, "atexit");
  }
  parse_options(argc, argv);
  init_gui(&argc, &argv);
  gtk_main();
}
