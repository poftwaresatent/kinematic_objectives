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

#ifndef KINEMATIC_OBJECTIVES_PRIORITIZATION_HPP
#define KINEMATIC_OBJECTIVES_PRIORITIZATION_HPP

#include <kinematic_objectives/types.h>
#include <iosfwd>
#include <vector>

namespace kinematic_objectives {
  
  class Objective;
  
  
  /**
     This will become an abstract base class at some point. Hopefully,
     anyway. But for now it just implements the scheme from
     [Siciliano:1991]. Generally recognized as having good convergence
     properties. But Chiaverini:1997 clearly shows that it suffers
     from undesirable algoritmic singularities, as discussed in a
     practical context by Mistry:2007.
  */
  class Prioritization
  {
  public:
    Prioritization();
    
    void processObjective(Matrix const & N_in,
			  Vector const & bias_in,
			  Objective const * objective,
			  Matrix & N_updater,
			  Vector & bias_out);
    
    void processCompound(/** initial nullspace projector */
			 Matrix const & N_init,
			 /** hierarchy of objectives to fuse/blend */
			 vector<Objective*> const & objectives,
			 /** resulting (fused/blended) bias */
			 Vector & bias_res,
			 /** accumulated nullspace projector at the end of the fusion */
			 Matrix & N_res);
    
    /** stream for debug output (use 0 for silent operation) */
    ostream * dbgos_;
    
    /** prefix for debug output (prepended to each line) */
    string dbgpre_;
  };
  
}

#endif // KINEMATIC_OBJECTIVES_PRIORITIZATION_HPP
