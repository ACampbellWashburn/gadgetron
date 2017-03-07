#include "SpiralOffResSubtract.h"
#include "ismrmrd/xml.h"

#ifdef USE_OMP
#include <omp.h>
#endif 

namespace Gadgetron{

  SpiralOffResSubtract::SpiralOffResSubtract() {}

  SpiralOffResSubtract::~SpiralOffResSubtract() {}

  int SpiralOffResSubtract::process_config(ACE_Message_Block* mb)
  {

    ISMRMRD::IsmrmrdHeader h;
    ISMRMRD::deserialize(mb->rd_ptr(),h);
    
    if (h.encoding.size() != 1) {
      GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
      GDEBUG("This Gadget only supports one encoding space\n");
      return GADGET_FAIL;
    }

  ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
  ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
  ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;
  
  sets_ = e_limits.set ? e_limits.set->maximum + 1 : 1;
  
  if (sets_ > 2) {
    GDEBUG("Phase subtraction only implemented for two sets for now\n");
    GDEBUG("Number of sets detected: %d, bailing out.\n", sets_);
    return GADGET_FAIL;
  }
  
  buffer_ = boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>[sets_]); 
  
  size_t bsize = sizeof(GadgetContainerMessage< GadgetContainerMessage<ISMRMRD::ImageHeader> >)*10000;
  
  for( size_t i=0; i<sets_; i++ ){
    buffer_[i].high_water_mark(bsize);
    buffer_[i].low_water_mark(bsize);
  }
  
  return GADGET_OK;
  }

  int SpiralOffResSubtract::
  process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
	  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  {
    GDEBUG("ACW sets = %i\n",sets_);
    // We need two sets to make a phase subtraction
    if (sets_ < 2) {
      GDEBUG("ACW sets<2\n");      
      if (this->next()->putq(m1) < 0) {
	return GADGET_FAIL;
      }
      return GADGET_OK;
    }

    int gw; 
    gw = m1->getObjectPtr()->user_int[7];
    GDEBUG("Subtract ACW guidewire = %i\n",gw);
    //if(gw==0){
    // if(this->next()->putq(m1)<0){
    //	return GADGET_FAIL; 
    //}
    //}else{

    size_t set = m1->getObjectPtr()->set;

    // Enqueue until we have images from both sets
    //

    if( buffer_[set].enqueue_tail(m1) < 0 ){
      GDEBUG("Message enqueue failed\n");
      return GADGET_FAIL;
    };

    // Phase subtract 
    //
    GDEBUG("ACW test Phase Subtract\n"); 
    while( buffer_[0].message_count()>0 && buffer_[1].message_count()>0 ) {

      ACE_Message_Block *mbq1, *mbq2;

      if( buffer_[0].dequeue_head(mbq1) < 0 || buffer_[1].dequeue_head(mbq2) < 0 ) {
	GDEBUG("Message dequeue failed\n");
	if( buffer_[set].message_count() > 0 ) 
	  buffer_[set].dequeue_tail(mbq1); // or m1 will be attempted deleted twice
	return GADGET_FAIL;
      }
	
      GadgetContainerMessage<ISMRMRD::ImageHeader> *pm1 = 
	AsContainerMessage<ISMRMRD::ImageHeader>(mbq1);

      GadgetContainerMessage< hoNDArray< std::complex<float> > > *cpm1 = 
	AsContainerMessage<hoNDArray< std::complex<float> > >(mbq1->cont());

      GadgetContainerMessage<ISMRMRD::ImageHeader> *pm2 = 
	AsContainerMessage<ISMRMRD::ImageHeader>(mbq2);

      GadgetContainerMessage< hoNDArray< std::complex<float> > > *cpm2 = 
	AsContainerMessage<hoNDArray< std::complex<float> > >(mbq2->cont());
	
      // Some validity checks
      //
     
      if( pm1->getObjectPtr()->image_index != pm2->getObjectPtr()->image_index ) {
	GDEBUG("Mismatch in image indices detected (%d, %d). Bailing out.\n", 
		      pm1->getObjectPtr()->image_index, pm2->getObjectPtr()->image_index);
	pm1->release();
	if( buffer_[set].message_count() > 0 ){
	  pm2->release();		
	  buffer_[set].dequeue_tail(mbq1); // or m1 will be attempted deleted twice
	}
	return GADGET_FAIL;
      }
      
      if (cpm1->getObjectPtr()->get_number_of_elements() != cpm2->getObjectPtr()->get_number_of_elements()) {
	GDEBUG("Mismatch in number of elements detected. Bailing out.\n");
	pm1->release();
	if( buffer_[set].message_count() > 0 ){
	  pm2->release();
	  buffer_[set].dequeue_tail(mbq1); // or m1 will be attempted deleted twice
	}
	return GADGET_FAIL;
      }

      std::complex<float> *p1 = cpm1->getObjectPtr()->get_data_ptr();
      std::complex<float> *p2 = cpm2->getObjectPtr()->get_data_ptr();
      float minimum_val = 0.0;

      // if(gw==1){

      #ifdef USE_OMP
      #pragma omp parallel for
      #endif

     
      for( long i = 0; i < (long)m2->getObjectPtr()->get_number_of_elements(); i++ ) {
	std::complex<float> tmp; 
	if(gw==2){
	  tmp = p1[i]-p2[i];//Complex subtraction
	}else if(gw==1){//Magnitude subtraction
	  tmp = (std::abs(p1[i])-std::abs(p2[i])); 
	}else{
	  tmp = p2[i];//Antomical imaging 
	}
	if(std::real( tmp) < minimum_val){
	  minimum_val=std::real(tmp);
	}
	p2[i] = tmp;
      }
      if(gw==1){
	for(long i = 0; i<(long)m2->getObjectPtr()->get_number_of_elements(); i++){
	  p2[i] = p2[i]+minimum_val; //Magnitude subtraction
	}
      
      }

      
      pm1->release();	
      pm2->getObjectPtr()->set = 0;

      if (this->next()->putq(pm2) < 0) {
	if( buffer_[set].message_count() > 0 ) {
	  pm2->release();
	  buffer_[set].dequeue_tail(mbq1); // or m1 will be attempted deleted twice
	}
	return GADGET_FAIL;
      }

    }//ifdeg USE_OMP
    
    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(SpiralOffResSubtract)
}
