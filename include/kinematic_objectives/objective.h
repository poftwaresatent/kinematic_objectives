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

#ifndef KINEMATIC_OBJECTIVES_OBJECTIVE_HPP
#define KINEMATIC_OBJECTIVES_OBJECTIVE_HPP

#include <kinematic_objectives/types.h>
#include <limits>


namespace kinematic_objectives {
  
  class KinematicModel;
  
  
  /**
     Working title, name and contents will depend on what turns out to
     actually work. Could also be termed a ProgressMonitor or
     BlenderFeedback or something along those lines.
     
     Ideas for things to store here:
     
     From blender.cpp perform_prioritization():
     - J_bar = J * N
     - N_up (nullspace updater: N -= N_up at each hierarchy level)
     - sigma (eigenvalues of J_bar, might be interesting to
       distinguish before / after regularization)
     - maybe eigenvalues of J
     - maybe eigenvalues of nullspace
     - probably input / output spaces of J_bar SVD (or some data
       derived from them)
     - maybe the same for the non-projected J (or something derived)
     - maybe compensated bias, or conversely the amount to which
       higher levels "help" achieving this one
     - biased and unbiased "objective update" (this is the gist behind
       Chiaverini's reconstruction error I think)
     - number of dimensions left after the objective
          
     From Blender::updateCompoundObjective():
     - current position (and velocity?) in objective space
     - 

     
     


  */
  struct Achievability {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    Vector placeholder;
  };
  
  
  class Objective
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    virtual ~Objective() {}
    
    virtual void init(KinematicModel const & model) { }
    
    /**
       \todo [low] consider moving this into a lower subclass, it
       seems to be relevant only for constraints.
    */
    virtual bool isActive() const { return true; }
    
    /**
       Subclasses have to set the bias_ and jacobian_ fields.
    */
    virtual void update(KinematicModel const & model) = 0;
    
    /**
       \todo [low] find a somewhat less hacky way for stacking multiple objectives on top of each other.
    */
    void stack(Objective const & t1, Objective const & t2);
    
    /**
       \todo [low] find a somewhat less hacky way for stacking multiple objectives on top of each other.
    */
    template<typename iterator_t>
    void stack(iterator_t begin, iterator_t end)
    {
      size_t ttnrows(0);
      for (iterator_t ii(begin); ii != end; ++ii) {
	ttnrows += (*ii)->jacobian_.rows();
      }
      size_t const ndof((*begin)->jacobian_.cols());
      bias_.resize(ttnrows);
      jacobian_.resize(ttnrows, ndof);
      for (size_t row(0); begin != end; row += (*begin++)->jacobian_.rows()) {
	bias_.block(   row, 0, (*begin)->jacobian_.rows(),    1) = (*begin)->bias_;
	jacobian_.block(row, 0, (*begin)->jacobian_.rows(), ndof) = (*begin)->jacobian_;
      }
    }
    
    inline Vector const & getBias() const
    { return bias_; }
    
    inline Matrix const & getJacobian() const
    { return jacobian_; }
    
    inline Achievability const & getAchievability() const
    { return achievability_; }
    
  protected:
    /**
       For target-based objectives, this is basically "desired -
       current." For potential-field-based (a.k.a. gradient-based)
       objectives, this is just the gradient. For constraints, this is
       directly the desired joint position change, whereas for hard
       and soft objectives, this is (usually) interpreted as an
       acceleration.
    */
    Vector bias_;
    
    Matrix jacobian_;
    
    Achievability achievability_;
  };
  
}

#endif // KINEMATIC_OBJECTIVES_OBJECTIVE_HPP
