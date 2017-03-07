#include "CalcTrajectories.h"
#include "hoNDArray_fileio.h"
#include "hoNDArray_utils.h"
#include "hoNDFFT.h"
#include <iostream>
#include "gadgetron_paths.h"
#include "ismrmrd/xml.h"

namespace Gadgetron{

  CalcTrajectories::CalcTrajectories() 
    : Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >()
  {    
    bufferX; 
    bufferY; 
    bufferZ;
    GradWaveform; 
  }

  CalcTrajectories::~CalcTrajectories() {}

  int CalcTrajectories::process_config(ACE_Message_Block* mb)
  {
   
    //GIRFs are stored as txt files
    FILE *file;    
    std::string FileNameX = get_gadgetron_home();
    FileNameX.append("/../GIRF/");
    FileNameX.append("GIRFx.txt");
    file  = fopen(FileNameX.c_str(),"r");
    if(!file){
      GDEBUG("Error: Unable to load file: %s \n",FileNameX.c_str());
      return GADGET_FAIL; 
    }else{
      while(!feof(file)){
	float tmp; 
	fscanf(file,"%f",&tmp);
	bufferX.push_back(tmp);
      }
    }
    fclose(file);

    std::string FileNameY = get_gadgetron_home();
    FileNameY.append("/../GIRF/");
    FileNameY.append("GIRFy.txt");
    file  = fopen(FileNameY.c_str(),"r");
    if(!file){
      GDEBUG("Error: Unable to load file: %s \n",FileNameY.c_str());
      return GADGET_FAIL;
    }else{
      while(!feof(file)){
	float tmp;
	fscanf(file,"%f",&tmp);
	bufferY.push_back(tmp);
      }
    }
    fclose(file);

    std::string FileNameZ = get_gadgetron_home();
    FileNameZ.append("/../GIRF/");
    FileNameZ.append("GIRFz.txt");
    file  = fopen(FileNameZ.c_str(),"r");
    if(!file){
      GDEBUG("Error: Unable to load file: %s \n",FileNameZ.c_str());
      return GADGET_FAIL;
    }else{
      while(!feof(file)){
	float tmp;
	fscanf(file,"%f",&tmp);
	bufferZ.push_back(tmp);
      }
    }
    fclose(file);

    //k-space trajectories are saved as txt file    
    std::string FileName = get_gadgetron_home();
    FileName.append("/../GIRF/");
    FileName.append("trajectories.txt");
    file = fopen(FileName.c_str(),"r");
    if(!file){
      GDEBUG("Error: Unable to load file: %s \n",FileName.c_str());
      return GADGET_FAIL;
    }else{
      while(!feof(file)){
	float tmp; 
	fscanf(file,"%f",&tmp); 
	GradWaveform.push_back(tmp);
      }
    }
    fclose(file);
  
    return GADGET_OK;
  }

  int CalcTrajectories::
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

    //Attach trajectories to data
    GadgetContainerMessage <hoNDArray< float > >* m3 = new GadgetContainerMessage<hoNDArray<float> >();  
    m1->getObjectPtr()->trajectory_dimensions = 3;     

    //Spiral GIRF Predicted Trajectories//
    *(m3->getObjectPtr()) = SpiralGIRFPredictedTraj(m1,bufferX, bufferY, bufferZ, GradWaveform);

    m1->cont(m2); 
    m2->cont(m3);
    if (this->next()->putq(m1) < 0) {
     GDEBUG("Unable to put data on queue\n");
     return GADGET_FAIL;
   }
    
    return GADGET_OK;
}


  //This would be modified for a specific application //
  hoNDArray<float> SpiralGIRFPredictedTraj(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, 
					   std::vector<float> bufferX, 
					   std::vector<float> bufferY, 
					   std::vector<float> bufferZ, 
					   std::vector<float> GradWaveform)
  {
  
/*--1. SET PARAMETERS--*/
   int interleave = m1->getObjectPtr()->idx.kspace_encode_step_1;
   int num_coils = m1->getObjectPtr()->active_channels;
   int num_samples = m1->getObjectPtr()->number_of_samples;
   int ii = 0;
   float dt = m1->getObjectPtr()->sample_time_us/1000000;   

   //Rotation matrices
   float R11 = m1->getObjectPtr()->read_dir[0];    
   float R12 = m1->getObjectPtr()->read_dir[1];
   float R13 = m1->getObjectPtr()->read_dir[2];
   float R21 = m1->getObjectPtr()->phase_dir[0];                 
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

   //GIRF parameters are stored in GIRF.txt files
   uint16_t N_GIRF = int(bufferX[0]);
   float dt_GIRF =bufferX[1]*1e-6; 

/*--2. NOMINAL TRAJECTORY--*/                           
   uint16_t N = std::round(N_GIRF*dt_GIRF/dt);
   if(N%2){N = N-1;} 
   
   hoNDArray< std::complex<float> > kx(N), ky(N), kz(N);
   std::fill(kx.begin(),kx.end(),0);
   std::fill(ky.begin(),ky.end(),0);                                                                                                                 
   std::fill(kz.begin(),kz.end(),0); 
   
   for(int s = 0; s<num_samples; s++){
     //Select current interleave
     float kxx = GradWaveform[interleave*num_samples*2+s];
     float kyy = GradWaveform[interleave*num_samples*2+s+num_samples];
     std::complex<float> test = std::complex<float>(R11*kxx+R12*kyy,0); 
     kx[s] = std::complex<float>(R11*kxx+R12*kyy,0);
     ky[s] = std::complex<float>(R21*kxx+R22*kyy,0);
     kz[s] = std::complex<float>(R31*kxx+R32*kyy,0);

     //make waveform periodic
     kx[num_samples*2-s-1] = kx[s];
     ky[num_samples*2-s-1] = ky[s];
     kz[num_samples*2-s-1] = kz[s];
   }

/*-- 3. GIRF: MATCH GIRF BANDWIDTH TO TRAJECTORIES--*/
   //This assumes trajectories are sampled the same as data. 
   //Otherwise, bandwith would be adjusted for both trajectories and GIRFs to match data sampling rate                                        
   
   hoNDArray< std::complex<float> > GIRFx(N), GIRFy(N), GIRFz(N);
   std::fill(GIRFx.begin(),GIRFx.end(),0);
   std::fill(GIRFy.begin(),GIRFy.end(),0);
   std::fill(GIRFz.begin(),GIRFz.end(),0);

   int fill =round( std::abs((N-N_GIRF)/2));
   if(N>N_GIRF){
     int j = 2; 
     for(int s = 0; s<N_GIRF; s++){
       GIRFx[fill+s] = std::complex<float>(bufferX[j],bufferX[j+1]); 
       GIRFy[fill+s] = std::complex<float>(bufferY[j],bufferY[j+1]); 
       GIRFz[fill+s] = std::complex<float>(bufferZ[j],bufferZ[j+1]);  
       j = j+2; 
     }
   }else{ 
     int j = 2 +fill; 
     for(int s = 0; s<N; s++){
       GIRFx[s] = std::complex<float>(bufferX[j],bufferX[j+1]); 
       GIRFy[s] = std::complex<float>(bufferY[j],bufferY[j+1]); 
       GIRFz[s] = std::complex<float>(bufferZ[j],bufferZ[j+1]); 
       j = j+2;      
     }
   }

/*-- 4. FFT --*/                                                                                                                                                                                                
   GadgetContainerMessage< hoNDArray< std::complex<float> > >* Ix =
     new GadgetContainerMessage<hoNDArray< std::complex<float> > >();
   *(Ix->getObjectPtr()) = kx;

   GadgetContainerMessage< hoNDArray< std::complex<float> > >* Iy =
     new GadgetContainerMessage<hoNDArray< std::complex<float> > >();
   *(Iy->getObjectPtr()) = ky;
   
   GadgetContainerMessage< hoNDArray< std::complex<float> > >* Iz =
     new GadgetContainerMessage<hoNDArray< std::complex<float> > >();
   *(Iz->getObjectPtr()) = kz;

   hoNDFFT<float>::instance()->fft(Ix->getObjectPtr(),0);
   hoNDFFT<float>::instance()->fft(Iy->getObjectPtr(),0);
   hoNDFFT<float>::instance()->fft(Iz->getObjectPtr(),0);

   hoNDFFT<float>::instance()->fftshift1D(*Ix->getObjectPtr());
   hoNDFFT<float>::instance()->fftshift1D(*Iy->getObjectPtr());
   hoNDFFT<float>::instance()->fftshift1D(*Iz->getObjectPtr());

   std::complex<float> *Ix1=Ix->getObjectPtr()->get_data_ptr();
   std::complex<float> *Iy1=Iy->getObjectPtr()->get_data_ptr();
   std::complex<float> *Iz1=Iz->getObjectPtr()->get_data_ptr();

/*-- 5. CONVOLVE NOMINAL TRAJECTORIES WITH GIRF --*/                                                                                                  
   
   //ADC_shift matches gradient clock to ADC clock
   float ADC_shift = 0.85e-6 + 0.5*dt;
   float dw = float(1/(dt*N));  
   float w = 0; 
   float PI = std::atan(1.0)*4; 

   hoNDArray< std::complex<float> > Ox(N),Oy(N),Oz(N);
 
   for(int s =0; s<N; s++){  
     Ox[s] = Ix1[s]*GIRFx[s]*std::exp(std::complex<float>(0,-2*PI*w*ADC_shift));  
     Oy[s] = Iy1[s]*GIRFy[s]*std::exp(std::complex<float>(0,-2*PI*w*ADC_shift));                                      
     Oz[s] = Iz1[s]*GIRFz[s]*std::exp(std::complex<float>(0,-2*PI*w*ADC_shift));
     w = w+dw;                                                                     
   }

/*-- 6. IFFT --*/                                                               

   *(Ix->getObjectPtr()) = Ox;
   hoNDFFT<float>::instance()->ifftshift1D(*Ix->getObjectPtr());
   hoNDFFT<float>::instance()->ifft(Ix->getObjectPtr(),0);
   
   *(Iy->getObjectPtr()) = Oy; 
   hoNDFFT<float>::instance()->ifftshift1D(*Iy->getObjectPtr());
   hoNDFFT<float>::instance()->ifft(Iy->getObjectPtr(),0);
   
   *(Iz->getObjectPtr()) = Oz;
   hoNDFFT<float>::instance()->ifftshift1D(*Iz->getObjectPtr());
   hoNDFFT<float>::instance()->ifft(Iz->getObjectPtr(),0);

/*-- 7. CORRECT WAVEFORM POLARITY AND TAKE ONLY NUM_SAMPLES --*/

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
   }

/*-- 8. ROTATE BACK TO LOGICAL COORDINATES --*/                                                                                                                                                                  
   float max_kx;
   float max_ky;
   hoNDArray< float > kx_Pred(num_samples), ky_Pred(num_samples); 
   for (int s = 0; s<num_samples; s++){
     kx_Pred[s] = invR11*Tx[s]+invR12*Ty[s]+invR13*Tz[s];
     ky_Pred[s] = invR21*Tx[s]+invR22*Ty[s]+invR23*Tz[s];
 
     if(std::abs(kx_Pred[s])>max_kx){
       max_kx = std::abs(kx_Pred[s]);
     }
     if(std::abs(ky_Pred[s])>max_ky){
       max_ky = std::abs(ky_Pred[s]);
     }                                                                                                                                                
   }

   float max_k = std::max(max_kx,max_ky);

/*-- 9. CALCULATE DENSITY COMPENSTAION --*/                                                                                                           

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

/*-- 10. PUT RESULS IN TRAJ FOR OUTPUT --*/	 
  
   hoNDArray<float> traj(3*num_samples);  
   ii=0;
   for(int s=0; s<num_samples; s++){
     traj[ii] =float (kx_Pred[s]/max_k*0.5); 
     traj[ii+1] = float(ky_Pred[s]/max_k*0.5); 
     traj[ii+2] = float(weights[s]);                                         
     ii = ii+3;
   }
 
   return traj; 
}


  GADGET_FACTORY_DECLARE(CalcTrajectories)                                                                                                                                                                    
}
