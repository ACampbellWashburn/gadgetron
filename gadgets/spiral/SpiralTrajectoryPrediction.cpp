#include "SpiralTrajectoryPrediction.h"
#include "hoNDArray_fileio.h"
#include "hoNDArray_utils.h"
#include "hoNDFFT.h"
#include <iostream>
#include "gadgetron_paths.h"
#include "ismrmrd/xml.h"

namespace Gadgetron{

  SpiralTrajectoryPrediction::SpiralTrajectoryPrediction() 
    : Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >()
    , buffer_(ACE_Message_Queue_Base::DEFAULT_HWM * 10000, ACE_Message_Queue_Base::DEFAULT_LWM * 10000)
  {
  }

  SpiralTrajectoryPrediction::~SpiralTrajectoryPrediction() {}

  int SpiralTrajectoryPrediction::process_config(ACE_Message_Block* mb)
  {

    ISMRMRD::IsmrmrdHeader h;
    ISMRMRD::deserialize(mb->rd_ptr(),h);
    //h.userParameters->userParameterString[0].value; 
    //filename_ = filename.value();
    return GADGET_OK;
  }

  int SpiralTrajectoryPrediction::close(unsigned long flags) {
    
    GDEBUG("SpiralTrajectoryPrediction::close...\n");
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
  
  int SpiralTrajectoryPrediction::
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


   

    GadgetContainerMessage <hoNDArray< float > >* m3 = AsContainerMessage<hoNDArray<float> >(m2->cont());

   
    int interleave = m1->getObjectPtr()->idx.kspace_encode_step_1;
    int num_coils = m1->getObjectPtr()->active_channels;
    int num_samples = m1->getObjectPtr()->number_of_samples;
    int ii = 0; 
    float dt = m1->getObjectPtr()->sample_time_us/1000000;//3.9e-6;
  

    //ACE_Message_Block* mbq;
    //if (buffer_.dequeue_head(mbq) < 0) {
      //GDEBUG("Message dequeue failed\n");
      //return GADGET_FAIL;
      //}

    //std::string xml = m1->getObjectPtr()->
    //ISMRMRD::IsmrmrdHeader h; 
    //ISMRMRD::deserialize(mbq->rd_ptr(),h); 


   
    int j; 
    FILE *file; 
    //FILE OPTION
    
    std::string FileName = get_gadgetron_home(); 
    FileName.append("/../GIRF/"); 
    FileName.append("trajectories.txt"); 
    file = fopen(FileName.c_str(),"r"); 
    float k_in[200000];
   j = 0; 
    while(!feof(file)){
      fscanf(file,"%f",&k_in[j]); 
      j++;
    }
    fclose(file);
    
    
    std::string FileNameX = get_gadgetron_home();
    FileNameX.append("/../GIRF/");
    FileNameX.append("GIRFxx_65.txt");
    file  = fopen(FileNameX.c_str(),"r");
    float GIRFx_in[200000];
    j = 0; 
    while(!feof(file)){
      fscanf(file,"%f",&GIRFx_in[j]);
      j++;
    }
    fclose(file);

    std::string FileNameY = get_gadgetron_home();
    FileNameY.append("/../GIRF/");
    FileNameY.append("GIRFyy_65.txt");
    file  = fopen(FileNameY.c_str(),"r");
    float GIRFy_in[200000];
    j = 0;
    while(!feof(file)){
      fscanf(file,"%f",&GIRFy_in[j]);
      j++;
    }
    fclose(file);  

    std::string FileNameZ = get_gadgetron_home();
    FileNameZ.append("/../GIRF/");
    FileNameZ.append("GIRFzz_65.txt");
    file  = fopen(FileNameZ.c_str(),"r");
    float GIRFz_in[200000];
    j = 0;
    while(!feof(file)){
      fscanf(file,"%f",&GIRFz_in[j]);
      j++;
    }
    fclose(file);

    uint16_t N_GIRF = int(GIRFx_in[0]);
    float dt_GIRF =GIRFx_in[1]*1e-6;

    std::complex<float> *GIRFxx = new std::complex<float> [N_GIRF];
    std::complex<float> *GIRFyy = new std::complex<float> [N_GIRF];
    std::complex<float> *GIRFzz = new std::complex<float> [N_GIRF];
    j =2; 
    for(int s = 0; s<N_GIRF; s++){
      GIRFxx[s] = std::complex<float>(GIRFx_in[j],GIRFx_in[j+1]);
      GIRFyy[s] = std::complex<float>(GIRFy_in[j],GIRFy_in[j+1]);
      GIRFzz[s] = std::complex<float>(GIRFz_in[j],GIRFz_in[j+1]);
      j = j+2; 
    }

     //-----------------TRAJECTORY PREDICITON-------------//
   

    //1. ACCUMULATE KX AND KY
    uint16_t N = std::round(N_GIRF*dt_GIRF/dt); 
    if(N%2){N = N-1;}
    std::complex<float> *kx = new std::complex<float> [N];
    std::complex<float> *ky = new std::complex<float> [N];
    std::complex<float> *kz = new std::complex<float> [N];
    memset(kx,0,N*sizeof(std::complex<float>));
    memset(ky,0,N*sizeof(std::complex<float>));
    memset(kz,0,N*sizeof(std::complex<float>));

    float R11 = m1->getObjectPtr()->read_dir[0];//R2??                                                                                                                                    
    float R12 = m1->getObjectPtr()->read_dir[1];
    float R13 = m1->getObjectPtr()->read_dir[2];
    float R21 = m1->getObjectPtr()->phase_dir[0];//R1?? 
    float R22 = m1->getObjectPtr()->phase_dir[1];
    float R23= m1->getObjectPtr()->phase_dir[2];
    float R31= m1->getObjectPtr()->slice_dir[0];
    float R32= m1->getObjectPtr()->slice_dir[1];
    float R33= m1->getObjectPtr()->slice_dir[2];

    float detR = R11*(R22*R33-R23*R32)-R12*(R21*R33-R23*R31)+R13*(R21*R32-R22*R31);
    float invR11 = 1/detR*(R22*R33-R23*R32);
    float invR12 = -1/detR*(R12*R33-R13*R32);
    float invR13 = 1/detR*(R12*R23-R13*R22);
    float invR21 = -1/detR*(R21*R33-R23*R31);
    float invR22 = 1/detR*(R11*R33-R13*R31);
    float invR23 = -1/detR*(R11*R23-R13*R21);
    float invR31 = 1/detR*(R21*R32-R22*R31);
    float invR32 = -1/detR*(R11*R32-R12*R31);
    float invR33 = 1/detR*(R11*R22-R12*R21);

//ISMRMRD OPTION                         
                                                                                                                                                       
    float *p3 = m3->getObjectPtr()->get_data_ptr();
    /*
    for(int s=0;s<num_samples;s++){
      kx[s] = std::complex<float> (R11*p3[ii]+R12*p3[ii+1],0);
      ky[s] = std::complex<float>(R21*p3[ii]+R22*p3[ii+1],0);
      kz[s] = std::complex<float>(R31*p3[ii]+R32*p3[ii+1],0);

      ii = ii+3;

      kx[num_samples*2-s-1] = kx[s];
      ky[num_samples*2-s-1] = ky[s];
      kz[num_samples*2-s-1] = kz[s];
    }   
*/

//File option
    for(int s = 0; s<num_samples; s++){
      float kxx = k_in[interleave*num_samples*2+s]; 
      float kyy = k_in[interleave*num_samples*2+s+num_samples];

      kx[s] = std::complex<float>(R11*kxx+R12*kyy,0); 
      ky[s] = std::complex<float>(R21*kxx+R22*kyy,0);
      kz[s] = std::complex<float>(R31*kxx+R32*kyy,0); 
 
      kx[num_samples*2-s-1] = kx[s]; 
      ky[num_samples*2-s-1] = ky[s]; 
      kz[num_samples*2-s-1] = kz[s]; 
    }

    //

    std::complex<float> *GIRFx = new std::complex<float> [N];
    std::complex<float> *GIRFy = new std::complex<float> [N];
    std::complex<float> *GIRFz = new std::complex<float> [N];
    memset(GIRFx,0,N*sizeof(std::complex<float>)); 
    memset(GIRFy,0,N*sizeof(std::complex<float>));
    memset(GIRFz,0,N*sizeof(std::complex<float>));

    int fill =round( std::abs((N-N_GIRF)/2));
    GDEBUG("ACW fill = %i\n",fill);    
     if(N>N_GIRF){    
       for(int s = 0; s<N_GIRF; s++){
      	 GIRFx[fill+s] = GIRFxx[s]; 
	 GIRFy[fill+s] = GIRFyy[s];
	 GIRFz[fill+s] = GIRFzz[s];
       }
    }else{
       for(int s = 0; s<N; s++){
	 GIRFx[s] = GIRFxx[s+fill];
	 GIRFy[s] = GIRFyy[s+fill];
	 GIRFz[s] = GIRFzz[s+fill]; 
       }
    } 
    

     //4. FFT
                                                                                                             
    GadgetContainerMessage< hoNDArray< std::complex<float> > >* Ix =
      new GadgetContainerMessage<hoNDArray< std::complex<float> > >();
    Ix->getObjectPtr()->create(N);
     

    GadgetContainerMessage< hoNDArray< std::complex<float> > >* Iy =
      new GadgetContainerMessage<hoNDArray< std::complex<float> > >();
    Iy->getObjectPtr()->create(N);

    GadgetContainerMessage< hoNDArray< std::complex<float> > >* Iz =
      new GadgetContainerMessage<hoNDArray< std::complex<float> > >();
    Iz->getObjectPtr()->create(N);
       
    memcpy(Ix->getObjectPtr()->get_data_ptr(),kx,N*sizeof(std::complex<float>));  
    hoNDFFT<float>::instance()->fft(Ix->getObjectPtr(),0);

    memcpy(Iy->getObjectPtr()->get_data_ptr(),ky,N*sizeof(std::complex<float>));
    hoNDFFT<float>::instance()->fft(Iy->getObjectPtr(),0);

    memcpy(Iz->getObjectPtr()->get_data_ptr(),kz,N*sizeof(std::complex<float>));
    hoNDFFT<float>::instance()->fft(Iz->getObjectPtr(),0);

    hoNDFFT<float>::instance()->fftshift1D(*Ix->getObjectPtr());
    hoNDFFT<float>::instance()->fftshift1D(*Iy->getObjectPtr());
    hoNDFFT<float>::instance()->fftshift1D(*Iz->getObjectPtr());
    

    std::complex<float> *Ix1=Ix->getObjectPtr()->get_data_ptr();
    std::complex<float> *Iy1=Iy->getObjectPtr()->get_data_ptr();
    std::complex<float> *Iz1=Iz->getObjectPtr()->get_data_ptr();

    //5. APPLY ADC SHIFT
  float ADC_shift = 0.85e-6+0.5*dt; 
  float dw = float(1/(dt*N)); 
  GDEBUG("ACW ADC_shift = %f  and dw = %f and num_samples = %i \n",ADC_shift, dw,num_samples); 

    //6. CONVOLVE WITH GIRF IN EACH AXIS (CROSS TERMS?!)
  std::complex<float> *Ox = new std::complex<float> [N];
  std::complex<float> *Oy = new std::complex<float> [N];
  std::complex<float> *Oz = new std::complex<float> [N];
  float w =0;  
  float PI = std::atan(1.0)*4; 
  


  std::string FileNameT = get_gadgetron_home();
  FileNameT.append("/../");
  FileNameT.append("test.txt");
  FILE *file2; 
  file2  = fopen(FileNameT.c_str(),"wb");
 


 for(int s =0; s<N; s++){
   //Add GIRF to multiplcation
   //Sum cross terms
   Ox[s] = Ix1[s];//*GIRFx[s];//*std::exp(std::complex<float>(0,2*PI*w*ADC_shift)); 
   Oy[s] = Iy1[s];//*GIRFy[s];//*std::exp(std::complex<float>(0,2*PI*w*ADC_shift)); 
   Oz[s] = Iz1[s];//*GIRFz[s];//*std::exp(std::complex<float>(0,2*PI*w*ADC_shift)); 
    w = w+dw; 
    fprintf(file2,"%f \n %f \n %f \n",std::real(kx[s]),std::real(ky[s]),std::real(kz[s]));
    //GDEBUG("ACW w = %f\n", w); 
  }
fclose(file2);
  


    //7. CROP BANDWIDTH TO MATCH DWELL TIME

  if(N<N_GIRF){
  }else{
  }

    //8. IFFT
  memcpy(Ix->getObjectPtr()->get_data_ptr(),Ox,N*sizeof(std::complex<float>));
  hoNDFFT<float>::instance()->ifftshift1D(*Ix->getObjectPtr());
  hoNDFFT<float>::instance()->ifft(Ix->getObjectPtr(),0);

  memcpy(Iy->getObjectPtr()->get_data_ptr(),Oy,N*sizeof(std::complex<float>));
  hoNDFFT<float>::instance()->ifftshift1D(*Iy->getObjectPtr());  
  hoNDFFT<float>::instance()->ifft(Iy->getObjectPtr(),0);

  memcpy(Iz->getObjectPtr()->get_data_ptr(),Oz,N*sizeof(std::complex<float>));
  hoNDFFT<float>::instance()->ifftshift1D(*Iz->getObjectPtr());  
  hoNDFFT<float>::instance()->ifft(Iz->getObjectPtr(),0);

  

    //9. CORRECT POLARITY AND TAKE ONLY NUM_SAMPLES
  
  std::complex<float> *Ix2=Ix->getObjectPtr()->get_data_ptr();
  std::complex<float> *Iy2=Iy->getObjectPtr()->get_data_ptr();
  std::complex<float> *Iz2=Iz->getObjectPtr()->get_data_ptr();
  float Tx[num_samples];
  float Ty[num_samples]; 
  float Tz[num_samples];
  bool u;   
  for(int s = 0; s<num_samples; s++){
    u = (std::real(Ix2[s])>0);
    if(u){
      Tx[s] = std::abs(Ix2[s]); 
    }else{
      Tx[s] = -std::abs(Ix2[s]); 
    }

    u = (std::real(Iy2[s])>0);
    if(u){
      Ty[s] = std::abs(Iy2[s]);
    }else{
      Ty[s] =-std::abs(Iy2[s]);
    }

    u = (std::real(Iz2[s])>0);
    if(u){
      Tz[s] = std::abs(Iz2[s]);
    }else{
      Tz[s] =-std::abs(Iz2[s]);
    }
    
    //fprintf(file2,"%f \n %f \n %f \n",Tx[s],(Ty[s]),Tz[s]);   
  }
  //fclose(file2);


    //10. ROTATE BACK TO LOGICAL COORDINATES
  float max_kx; 
  float max_ky; 
  float kx_Pred[num_samples];
  float ky_Pred[num_samples];
  for (int s = 0; s<num_samples; s++){
    kx_Pred[s] = invR11*Tx[s]+invR12*Ty[s]+invR13*Tz[s]; 
    ky_Pred[s] = invR21*Tx[s]+invR22*Ty[s]+invR23*Tz[s];
    float kz = invR31*Tx[s]+invR32*Ty[s]+invR33*Tz[s];
    if(std::abs(kx_Pred[s])>max_kx){
      max_kx = std::abs(kx_Pred[s]); 
    }
    if(std::abs(ky_Pred[s])>max_ky){
      max_ky = std::abs(ky_Pred[s]); 
    }  
    //GDEBUG("ACW kz = %f \n",kz); 
  }

  float max_k = std::max(max_kx,max_ky);
    //11. CALCULATE DENSITY COMPENSTAION

  float weights[num_samples]; 
  float Gx; 
  float Gy; 
  for(int s = 0; s<num_samples; s++){
    if(s==0){
      Gx = (kx_Pred[s+1]-kx_Pred[s])/dt;
      Gy = (ky_Pred[s+1]-ky_Pred[s])/dt; 
    }else if(s==(num_samples-1)){
      Gx=(kx_Pred[s]-kx_Pred[s-1])/dt;
      Gy =(ky_Pred[s]-ky_Pred[s-1])/dt;
    }else{
      Gx = (kx_Pred[s+1]-kx_Pred[s-1])/(2*dt);
      Gy = (ky_Pred[s+1]-ky_Pred[s-1])/(2*dt); 
    }
    float w1 = std::sqrt(Gx*Gx+Gy*Gy); 
    float w2 = std::abs(std::sin(std::atan(Gx/Gy)-std::atan(kx_Pred[s]/ky_Pred[s])));
    weights[s] = w1*w2; 
  }
 

    //12. Scale to -0.5:0.5 and PASS ALONG


   ii=0;
    for(int s=0; s<num_samples; s++){
      //GDEBUG("ACW p3 = %f, kx = %f\n",p3[ii],kx_Pred[s]/max_k*0.5);
      p3[ii]=kx_Pred[s]/max_k*0.5; 
      p3[ii+1]=ky_Pred[s]/max_k*0.5;
      p3[ii+2]=weights[s]; //p3[ii+2];
      //fprintf(file2,"%f \n %f \n",p3[ii],p3[ii+1]);       
      ii = ii+3; 
    }
    //fclose(file2); 
  
    memcpy(m3->getObjectPtr()->get_data_ptr(),p3,m3->getObjectPtr()->get_number_of_elements()*sizeof(float));

    //-----------------END TRAJECTORY PREDICTION-----------------//

  //ncm1->cont(ncm2);
  //ncm2->cont(ncm3);

    
    m1->cont(m2); 
    m2->cont(m3);
    if (this->next()->putq(m1) < 0) {
     GDEBUG("Unable to put data on queue\n");
     return GADGET_FAIL;
    }

    //if (this->next()->putq(ncm1)<0){
    //GDEBUG("Unable to put data on queue\n");
    //return GADGET_FAIL; 
    //}

   
    return GADGET_OK;
  }
  
  GADGET_FACTORY_DECLARE(SpiralTrajectoryPrediction)
}
