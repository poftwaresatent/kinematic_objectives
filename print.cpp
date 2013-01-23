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

#include "print.hpp"
#include <iostream>
#include <sstream>
#include <stdio.h>


namespace pbmockup {  
  
  
  string pstring(Vector const & vv)
  {
    ostringstream os;
    print(vv, os, "", "", true);
    return os.str();
  }
  
    
  string pstring(Matrix const & mm, string const & prefix)
  {
    ostringstream os;
    print(mm, os, "", prefix);
    return os.str();
  }
    
    
  void print(Vector const & vv, ostream & os,
	     string const & title, string const & prefix,
	     bool nonl)
  {
    print((Matrix const &) vv, os, title, prefix, true, nonl);
  }
  
  
  string pstring(double vv)
  {
    static int const buflen(32);
    static char buf[buflen];
    memset(buf, 0, sizeof(buf));
    if (isinf(vv)) {
      snprintf(buf, buflen-1, " inf    ");
    }
    else if (isnan(vv)) {
      snprintf(buf, buflen-1, " nan    ");
    }
    else if (fabs(fmod(vv, 1)) < 1e-6) {
      snprintf(buf, buflen-1, "%- 7d  ", static_cast<int>(rint(vv)));
    }
    else {
      snprintf(buf, buflen-1, "% 6.4f  ", vv);
    }
    string str(buf);
    return str;
  }
  
  
  void print(Matrix const & mm, ostream & os,
	     string const & title, string const & prefix,
	     bool vecmode, bool nonl)
  {
    char const * nlornot("\n");
    if (nonl) {
      nlornot = "";
    }
    if ( ! title.empty()) {
      os << prefix << title << nlornot;
    }
    if ((mm.rows() <= 0) || (mm.cols() <= 0)) {
      os << prefix << "   (empty)" << nlornot;
    }
    else {
      
      if (vecmode) {
	if ( ! prefix.empty())
	  os << prefix;
	os << "  ";
	for (int ir(0); ir < mm.rows(); ++ir) {
	  os << pstring(mm.coeff(ir, 0));
	}
	os << nlornot;
	
      }
      else {

	for (int ir(0); ir < mm.rows(); ++ir) {
	  if ( ! prefix.empty())
	    os << prefix;
	  os << "  ";
	  for (int ic(0); ic < mm.cols(); ++ic) {
	    os << pstring(mm.coeff(ir, ic));
	  }
	  os << nlornot;
	}
	  
      }
    }
  }

}
