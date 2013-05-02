/** \file hoOpticalFlowSolver.h
    \brief Abstract class for a CPU-based optical flow registration solver.

    hoOpticalFlowSolver is derived from class opticalFlowSolver 
    and implements the computation of the spatial and temporal gradients.
    A pure virtual function is expected to implement the specific algorithm (Horn-Schunck, Cornelius-Kanade).
*/

#pragma once

#include "hoNDArray.h"
#include "hoNDArray_operators.h"
#include "hoNDArray_elemwise.h"
#include "hoRegistration_utils.h"
#include "opticalFlowSolver.h"
#include "cpureg_export.h"

namespace Gadgetron{

  template<class T, unsigned int D> class EXPORTCPUREG hoOpticalFlowSolver 
    : public opticalFlowSolver< hoNDArray<T>,D >
  {  
  public:
  
    hoOpticalFlowSolver() : opticalFlowSolver< hoNDArray<T>,D >() {}   
    virtual ~hoOpticalFlowSolver() {}
    
  protected:

    // Inherited and still pure virtual...
    virtual boost::shared_ptr< hoNDArray<T> > core_solver( hoNDArray<T> *gradient_image, hoNDArray<T> *stencil_image ) = 0;      

    // CPU-based computation of the spatial and temporal image gradient
    //
    
    virtual void core_grad_spatial( T *fixed_image, T *moving_image, T *gradient_image, 
				    typename uintd<D>::Type matrix_size_moving, 
				    unsigned int number_of_batches_fixed, 
				    unsigned int number_of_batches_moving );
    
    virtual void core_grad_temporal( T *fixed_image, T *moving_image, T *gradient_image, 
				     typename uintd<D>::Type matrix_size_moving, 
				     unsigned int number_of_batches_fixed, 
				     unsigned int number_of_batches_moving );
  };  
}
