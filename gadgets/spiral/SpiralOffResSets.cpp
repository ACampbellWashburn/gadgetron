#include "SpiralOffResSets.h"
#include "hoNDArray_fileio.h"
#include "hoNDArray_utils.h"

namespace Gadgetron{

  SpiralOffResSets::SpiralOffResSets() 
    : Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >()
    , buffer_(ACE_Message_Queue_Base::DEFAULT_HWM * 10000, ACE_Message_Queue_Base::DEFAULT_LWM * 10000)
  {
  }

  SpiralOffResSets::~SpiralOffResSets() {}

  int SpiralOffResSets::process_config(ACE_Message_Block* mb)
  {
    //filename_ = filename.value();
    return GADGET_OK;
  }

  int SpiralOffResSets::close(unsigned long flags) {
    
    GDEBUG("SpiralOffResSets::close...\n");
    GDEBUG("Number of items on Q: %d\n", buffer_.message_count());

    int ret = Gadget::close(flags);
    unsigned int readouts_buffered = buffer_.message_count();

    if( readouts_buffered == 0 )
      return GADGET_OK;
    
    // Get the array size from the dimensions of the first buffer entry
    //

    ACE_Message_Block* mbq;
    if (buffer_.dequeue_head(mbq) < 0) {
      GDEBUG("Message dequeue failed\n");
      return GADGET_FAIL;
    }

    GadgetContainerMessage< hoNDArray< std::complex<float> > > *daq = AsContainerMessage<hoNDArray< std::complex<float> > >(mbq);
    
    if (!daq) {
      GDEBUG("Unable to interpret data on message queue\n");
      return GADGET_FAIL;
    }

    hoNDArray< std::complex<float> > *entry = daq->getObjectPtr();
    std::vector<size_t> dims_profile = *entry->get_dimensions();
    std::vector<size_t> dims = dims_profile;
    dims.push_back(readouts_buffered);

    // Allocate array for result
    //

    hoNDArray< std::complex<float> > result( &dims );

    // And copy over the first profile
    //

    {
      hoNDArray< std::complex<float> > tmp( &dims_profile, result.get_data_ptr() );
      tmp = *entry;
    }

    mbq->release();
    
    // Copy the remaining profiles to the array
    //
    
    for (unsigned int i = 1; i < readouts_buffered; i++) {
      
      if (buffer_.dequeue_head(mbq) < 0) {
        GDEBUG("Message dequeue failed\n");
        return GADGET_FAIL;
      }
      
      daq = AsContainerMessage<hoNDArray< std::complex<float> > >(mbq);
      
      if (!daq) {
        GDEBUG("Unable to interpret data on message queue\n");
        return GADGET_FAIL;
      }
      
      entry = daq->getObjectPtr();
      hoNDArray< std::complex<float> > tmp( &dims_profile, result.get_data_ptr()+i*entry->get_number_of_elements() );
      tmp = *entry;
      mbq->release();
    }      
  
    // Reshape to get the coil dimension as the last
    //
  
    std::vector<size_t> order; order.push_back(0); order.push_back(2); order.push_back(1);
    result = *permute( &result, &order);

    // Write out the result
    //
  
    //write_nd_array< std::complex<float> >( &result, filename_.c_str() );
  
    return GADGET_OK;
  }
  
  int SpiralOffResSets::
  process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
          GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  {
    
    // Noise should have been consumed by the noise adjust, but just in case...
    //
    
    bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT).isSet(m1->getObjectPtr()->flags);
    if (is_noise) {
      m1->release();
      return GADGET_OK;
    }

    int gw; 
    gw = m1->getObjectPtr()->user_int[7];
    GDEBUG("ACW guidewire = %i\n",gw);
   
    //ACW test
    //boost::shared_ptr<hoNDArray<float_complext>> test=m2->getObjectPtr()->data_;  
       
    //hoNDArray< std::complex<float> >* pdata = m2->getObjectPtr();     
    //write_nd_array(pdata,"/write_test.txt"); 
 

    GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* ncm1= new GadgetContainerMessage<ISMRMRD::AcquisitionHeader>();
    *ncm1->getObjectPtr()= *m1->getObjectPtr();
    ncm1->getObjectPtr()->idx.set =1;

    GadgetContainerMessage< hoNDArray< std::complex<float> > >* ncm2 = new GadgetContainerMessage< hoNDArray< std::complex<float> > >;
    *ncm2->getObjectPtr() = *m2->getObjectPtr();

    GadgetContainerMessage <hoNDArray< float > >* m3 = AsContainerMessage<hoNDArray<float> >(m2->cont());

    GadgetContainerMessage< hoNDArray<float> >* ncm3 = new GadgetContainerMessage <hoNDArray< float > >;
    *ncm3->getObjectPtr() = *m3->getObjectPtr();
    
    std::complex<float> *p1 = ncm2->getObjectPtr()->get_data_ptr();
    int num_coils = m1->getObjectPtr()->active_channels;
    int num_samples = m1->getObjectPtr()->number_of_samples;
    int ii = 0; 
    float dt = m1->getObjectPtr()->sample_time_us/1000000;//3.9e-6;
    GDEBUG("ACW sample time = %f\n",dt); 
    float PI = std::atan(1.0)*4;
    for(int c=0;c<num_coils;c++){
      float t = 0;       
      for(int s=0;s<num_samples;s++){
	std::complex<float> test = p1[ii];
	if(gw==0){
	  p1[ii]=test;
	}else if (gw==1){
	  std::complex<float> mult = exp(std::complex<float>(0,2*PI*(100)*t));
	  p1[ii]=test*mult;
	}else{
	  std::complex<float> mult = exp(std::complex<float>(0,2*PI*(0)*t)); 
	  p1[ii] = test*mult; 
	}
	t = t+dt;
	ii++;
      }
    }
    memcpy(ncm2->getObjectPtr()->get_data_ptr(),p1,ncm2->getObjectPtr()->get_number_of_elements()*sizeof(float)*2);


    std::complex<float> *p2 = m2->getObjectPtr()->get_data_ptr();
    ii=0;
    for(int c =0; c<num_coils;c++){
      float t =0; 
      for(int s=0;s<num_samples;s++){
	std::complex<float>test = p2[ii];
	std::complex<float>mult = exp(std::complex<float>(0,2*PI*(-100)*t));
	if(gw==0){
	  p2[ii] = 0; 
	}else if(gw==1){
	  p2[ii]=test*mult;
	}else{
	  p2[ii]=test*mult; 
	}
        t=t+dt;
	ii++;
      }
    }
    memcpy(m2->getObjectPtr()->get_data_ptr(),p2,m2->getObjectPtr()->get_number_of_elements()*sizeof(float)*2);


    ncm1->cont(ncm2);
    ncm2->cont(ncm3);

    
       
    if (this->next()->putq(m1) < 0) {
     GDEBUG("Unable to put data on queue\n");
     return GADGET_FAIL;
    }

    if (this->next()->putq(ncm1)<0){
     GDEBUG("Unable to put data on queue\n");
     return GADGET_FAIL; 
    }

   
    return GADGET_OK;
  }
  
  GADGET_FACTORY_DECLARE(SpiralOffResSets)
}
