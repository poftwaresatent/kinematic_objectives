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

#include "baerlocher_algorithm.hpp"
#include <iostream>

using namespace kinematic_elastic;


int main (int argc, char ** argv)
{
  size_t const ndof(2);
  Model const model(ndof);
  
  for (double bm(1.0e-3); bm <= 1.0; bm *= 1.2) {
    
    cout << "# bm: " << bm << "\n";
    
    Vector state(ndof);
    state << 0.3, -0.2;
    
    tasklist_t tasklist;
    tasklist.push_back(task_s(ndof, 1, bm));
    tasklist.push_back(task_s(ndof, 1, bm * 0.5 * M_PI / 180.0));
    
    tasklist[0].desired << 1.2;
    tasklist[1].desired << 35.0 * M_PI / 180.0;
    
    tasklist[1].Jacobian << 0.0, 1.0; // constant in this case
    
    for (size_t ii(0); ii < 10000; ++ii) {
      double const q0(state.coeff(0));
      double const q1(state.coeff(1));
      
      tasklist[0].current <<
	cos(q0) + cos(q0 + q1);
      tasklist[0].Jacobian <<
	-sin(q0) - sin(q0 + q1),
	-sin(q0 + q1);
      
      tasklist[1].current <<
	q1;
      
      Vector dq = baerlocher_algorithm (model, state, tasklist);
      
      dump(state, tasklist, dq);
      
      state += dq;
      
    }
    
    cout << "\n\n";
  }
  
  return 0;
}
