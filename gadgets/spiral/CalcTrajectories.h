#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_spiral_export.h"

#include <ismrmrd/ismrmrd.h>

namespace Gadgetron{

  class EXPORTGADGETS_SPIRAL CalcTrajectories : 
  public Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
    {
    public:
      GADGET_DECLARE(CalcTrajectories);

      CalcTrajectories();
      ~CalcTrajectories();

    protected:

      virtual int process_config(ACE_Message_Block* mb);

      virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
			  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
      
      //hoNDArray<float> SpiralGIRFPredictedTraj(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1);  
      
    private:
      
      std::vector<float> bufferX;
      std::vector<float> bufferY;
      std::vector<float> bufferZ;
      std::vector<float> GradWaveform;
    };

  //void EXPORTGADGETS_SPIRAL 
  hoNDArray<float> SpiralGIRFPredictedTraj(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, std::vector<float> bufferX,  std::vector<float> bufferY, std::vector<float> bufferZ, std::vector<float> GradWaveform); 
}
