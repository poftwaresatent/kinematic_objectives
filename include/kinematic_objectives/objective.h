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
#include <kinematic_objectives/pseudo_inverse.h>
#include <limits>


namespace kinematic_objectives {
  
  class KinematicModel;
  
  
  /**
     \note
     IDEAS for feedback info about objectives and their blending... so
     this is just a working title. Name and contents will depend on
     what turns out to actually work. Could also be termed a
     ProgressMonitor or BlenderFeedback or something along those
     lines.
     
     Potenital (other) things to store somewhere:
     - N_up (nullspace updater: N -= N_up at each hierarchy level)
     - biased and unbiased "objective update" (this is the gist behind
       Chiaverini's reconstruction error I think)
     - number of dimensions left after the objective
     - current position in objective space (not explicitly required
       from Objective API for now)
     - current velocity in objective space (trivial to compute from J
       and qdot)
       
     Another description of desirables: for example, which unilateral
     constraints were activated, the resulting constraint nullspace,
     whether that created any conflicts with hard objectives, whether
     there were any algorithmic singularities between constraints and
     hard objectives. For soft objectives, it should at least be
     possible to specify the remaining space of motion freedom they
     had available collectively.
  */
  class Objective
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    string const name_;
    
    
    explicit Objective(string const & name)
      : name_(name) { }
    
    virtual ~Objective()
    { /* nop */ }
    
    virtual void init(KinematicModel const & model)
    { /* nop */ }
    
    /**
       Default implementation: always returns true.
       
       \todo [low] consider moving this into a lower subclass, it
       seems to be relevant only for constraints.
    */
    virtual bool isActive() const
    { return true; }
    
    /**
       Subclasses have to set the bias_ and jacobian_ fields.
    */
    virtual void update(KinematicModel const & model) = 0;
    
    /**
       Subclasses have to implement this method so that
       Achievability::residual_error_magnitude_ can be set to a
       relevant number. The problem is that each objective lives in a
       different space, so there is no unifying metric. Even worse,
       some objectives may live in spaces that do not have a sensible
       metric, for instance when they mix rotational and translational
       degrees of freedom. So it is left to the implementer to decide
       how best to characterize a residual error.
       
       \note The definition of residual error comes from
       [Chiaverini:1997] and is \f$ e = \dot{x} - J\dot{q} \f$.
       
     */
    virtual double computeResidualErrorMagnitude(Vector const & ee) const = 0;
    
    /**
       Turn this objective into a stack of the two given
       objectives. This simply stacks the Jacobians and the bias
       vectors of t1 onto those of t2.
    */
    void stack(Objective const & t1, Objective const & t2);
    
    /**
       Turn this objective into a stack of a collection of
       objectives. The collection (or range) to use is simply given as
       iterators in a semi-generic way. Should work as-is for
       std::vector<Objective*> and std::list<Objective*>.
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
    
  public:
    // XXXX all mutable public fields should probably be collected
    // elsewhere separately
    
    void clearFeedback() const	// rfct
    {
      jbar_svd_.original_sigma.resize(0);
      jbar_svd_.regularized_sigma.resize(0);
      jbar_svd_.output_space.resize(0, 0);
      jbar_svd_.input_space.resize(0, 0);
      bias_comp_.resize(0);
      jbar_.resize(0, 0);
      jbar_inv_.resize(0, 0);
    }
    
    mutable PseudoInverseFeedback jbar_svd_;
    mutable Vector bias_comp_;
    mutable Matrix jbar_;
    mutable Matrix jbar_inv_;
  };
  
}

#endif // KINEMATIC_OBJECTIVES_OBJECTIVE_HPP
