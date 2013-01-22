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

#include <gtk/gtk.h>
#include <cmath>
#include <iostream>

using namespace std;


static double const dimx(10.);
static double const dimy(8.);
static double const deg(M_PI / 180.);

static GtkWidget * gw(0);
static gint gw_width(800), gw_height(640);
static gint gw_sx, gw_sy, gw_x0, gw_y0;
static int play(0);


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
  
  cairo_set_source_rgb (cr, 0.5, 0.5, 0.5);
  cairo_set_line_width (cr, 2.0);
  cairo_rectangle (cr, gw_x0 - 2, gw_y0 + 2, dimx * gw_sx + 4, dimy * gw_sy - 4);
  cairo_stroke (cr);
  
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


gint idle (gpointer data)
{
  if (play) {
    update ();
  }
  return TRUE;
}


int main (int argc, char ** argv)
{
  GtkWidget *window, *vbox, *hbox, *btn;
  
  gtk_init (&argc, &argv);
  //  init ();
  
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
  gtk_main ();
  
  return 0;
}
