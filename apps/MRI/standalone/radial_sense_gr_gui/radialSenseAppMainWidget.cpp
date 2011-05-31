#include "radialSenseAppMainWidget.h"

#include "NFFT.h"
#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "vector_td_operators.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"
#include "radial_utilities.h"
#include "cgOperatorNonCartesianSense.h"
#include "cuCGIdentityOperator.h"
#include "cuCG.h"
#include "b1_map.h"

#include "UIconstants.h"
#include "GLReconWidget.h"

#include <QtGui/QFileDialog>
#include <QtGui/QProgressDialog>
#include <QtGui/QMessageBox>
#include <QtCore/QSignalMapper>

#include <assert.h>

using namespace std;

void radialSenseAppMainWindow::resetPrivateData()
{
  if( statusLabel ) delete statusLabel;
  statusLabel = 0x0;
}

radialSenseAppMainWindow::radialSenseAppMainWindow(QWidget *parent) : QMainWindow(parent)
{
  setupUi(this);
  retranslateUi(this);

  resetPrivateData();

  matrixSizeSpinBox->setValue(MATRIX_SIZE_INITIAL_VALUE);
  oversampledMatrixSizeSpinBox->setValue(MATRIX_SIZE_OS_INITIAL_VALUE);
  numIterationsSpinBox->setValue(NUM_ITERATIONS_INITIAL_VALUE);
  regularizationWeightSpinBox->setValue(REG_WEIGHT_INITIAL_VALUE);
  kernelSizeSpinBox->setValue(KERNEL_SIZE_INITIAL_VALUE);

  // Menu actions
  connect(actionOpen_cplx_file, SIGNAL(triggered()), this, SLOT(open()));
  connect(actionSave_image, SIGNAL(triggered()), this, SLOT(saveImage()));
  connect(actionClose, SIGNAL(triggered()), this, SLOT(close()));
  connect(actionExit, SIGNAL(triggered()), qApp, SLOT(quit()));

  // Originally, the thought was to put multiple ReconWidgets in the app. 
  // This is why the SignalMapper is used rather than the basic signals below.

  // Connect to the reconWidgets' frameChanged slots
  QSignalMapper *signalMapper1 = new QSignalMapper(this);
  connect(reconWidget->projectionSelectionScrollBar, SIGNAL(valueChanged(int)), signalMapper1, SLOT(map()));
  signalMapper1->setMapping(reconWidget->projectionSelectionScrollBar, 1 );
  connect(signalMapper1, SIGNAL(mapped(int)), this, SLOT(centralProjectionChanged(int)));

  // Connect to the reconWidgets' projectionsPerFrameChanged slots
  QSignalMapper *signalMapper2 = new QSignalMapper(this);
  connect(reconWidget->numProjectionsScrollBar, SIGNAL(valueChanged(int)), signalMapper2, SLOT(map()));
  signalMapper2->setMapping(reconWidget->numProjectionsScrollBar, 1 );
  connect(signalMapper2, SIGNAL(mapped(int)), this, SLOT(projectionsPerFrameChanged(int)));

  // Allocate encoding operator for non-Cartesian Sense
  E = boost::shared_ptr< cgOperatorNonCartesianSense<float,2> >( new cgOperatorNonCartesianSense<float,2>() );  

  // Allocate preconditioner
  D = boost::shared_ptr< cuCGPrecondWeight<float_complext::Type> >( new cuCGPrecondWeight<float_complext::Type>() );

  // Allocate regularization image operator and corresponding rhs operator
  rhs_buffer = boost::shared_ptr< cgOperatorSenseRHSBuffer<float,2> >( new cgOperatorSenseRHSBuffer<float,2>() );
  R = boost::shared_ptr< cuCGImageOperator<float,float_complext::Type> >( new cuCGImageOperator<float,float_complext::Type>() );  
  R->set_weight( 1.0f );
  R->set_encoding_operator( rhs_buffer );

  // Setup solver
  cg.add_matrix_operator( E );  // encoding matrix
  cg.add_matrix_operator( R );  // regularization matrix
  cg.set_preconditioner ( D );  // preconditioning matrix
  cg.set_iterations( get_num_iterations() );
  cg.set_limit( 1e-6 );
  cg.set_output_mode( cuCG<float, float_complext::Type>::OUTPUT_SILENT );
}

/*
  Slots
*/

void radialSenseAppMainWindow::open()
{
  // Open dialog box
  QString filename = QFileDialog::getOpenFileName( this, tr("Open File"), "./", tr("Raw data (*.cplx)"));

  if( filename.size() == 0 )
    return; // Cancel

  // Close current file
  close();

  // Update status bar
  statusLabel = new QLabel(filename);	
  statusBar()->addWidget(statusLabel);

  // Chose startup frame
  reconWidget->projectionSelectionScrollBar->setValue(get_matrix_size().vec[0]>>2);
  reconWidget->numProjectionsScrollBar->setValue(34);

  // Read samples from disk
  host_samples = read_nd_array<float_complext::Type>(filename.toLatin1().constData());

  cout << endl << "loaded dataset with " << host_samples.get_number_of_elements() << " samples." << endl;

  // This is to prevent the user changing the matrix sizes before any data is initially loaded
  matrixSizeSpinBox->setEnabled(true); 
  oversampledMatrixSizeSpinBox->setEnabled(true);
  
  replan();
}

void radialSenseAppMainWindow::saveImage()
{ /*
  // Open dialog box
  QString filename = QFileDialog::getSaveFileName( this, tr("Save image to file"), "./", tr("Raw float data (*.raw)"));

  if( filename.size() == 0 )
    return; // Cancel

  // This code is copied from 'reconstruct' and slightly modified...

  <cut..>

  LOOP:

      // Save file
      cudaMemcpy( tmp, devPtr, prod(get_matrix_size())*sizeof(float), cudaMemcpyDeviceToHost );
      fwrite( tmp, prod(get_matrix_size()), sizeof(float), fout );

      // Report any errors not already caught...
      err = cudaGetLastError();
      if( err != cudaSuccess ){
	QMessageBox::critical( this, tr("Cuda error"), tr(cudaGetErrorString(err)) );
	actionExit->trigger();
      }
	
      END LOOP:

      reconWidget->projectionNumberSpinBox->setValue(reconWidget->projectionNumberSpinBox->value()+20);
    }

  fclose(fout);
  cudaFree(devPtr);
  */
}

void radialSenseAppMainWindow::close()
{	
  resetPrivateData();
}

void radialSenseAppMainWindow::replan()
{
  QProgressDialog progress("Calibrating", "", 0, 4, this);
  progress.setWindowModality(Qt::WindowModal);
  progress.setValue(0);
  progress.show();

  // Set GUI elements before the plan is created to avoid triggering unneccessary reconstructions
  unsigned int maxProjections = min(get_matrix_size().vec[0]<<2, (get_num_points_per_array_coil()/get_num_samples_per_projection())>>1);
  reconWidget->numProjectionsScrollBar->setMaximum(maxProjections);
  reconWidget->numProjectionsSpinBox->setMaximum(maxProjections);
  unsigned int maxCentralProjection = get_maximum_central_projection();
  reconWidget->projectionSelectionScrollBar->setMaximum(maxCentralProjection);
  reconWidget->projectionNumberSpinBox->setMaximum(maxCentralProjection);
  unsigned int minCentralProjection = get_num_projections_per_frame()>>1;
  reconWidget->projectionSelectionScrollBar->setMinimum(minCentralProjection);
  reconWidget->projectionNumberSpinBox->setMinimum(minCentralProjection);
									    
  progress.setValue(1);

  // Pass matrix size to GLReconWidget::initializeGL
  //	reconWidget->openglCanvas->setMatrixSize( get_matrix_size().vec[0], get_matrix_size().vec[1] );

  progress.setValue(2);

  const unsigned int samples_per_profile = get_num_samples_per_projection();
  const unsigned int num_profiles = get_num_points_per_array_coil() / samples_per_profile;
  const unsigned int profiles_per_frame = get_num_projections_per_frame();
  const unsigned int frames_per_reconstruction = NUM_FRAMES_PER_CSM_RECON;
  const unsigned int profiles_per_reconstruction = get_num_projections_per_frame()*frames_per_reconstruction;
  const unsigned int samples_per_reconstruction = profiles_per_reconstruction*samples_per_profile;

  // Density compensation weights are constant throughout all reconstrutions
  dcw  = compute_radial_dcw_golden_ratio_2d
    ( samples_per_profile, profiles_per_frame, (float)get_matrix_size_os().vec[0]/(float)get_matrix_size().vec[0], 
      get_one<float>()/((float)samples_per_profile/(float)max(get_matrix_size().vec[0],get_matrix_size().vec[1])) );
  
  progress.setValue(3);

  // Setup plan for convolution
  plan.setup( get_matrix_size(), get_matrix_size_os(), get_kernel_width() );
  
  // Temporary oversampled image buffer
  vector<unsigned int> image_os_dims = uintd_to_vector<2>(get_matrix_size_os()); 
  image_os_dims.push_back(frames_per_reconstruction); image_os_dims.push_back(get_num_coils());    
  cuNDArray<float_complext::Type> *image_os = new cuNDArray<float_complext::Type>(); 
  if( !image_os->create(image_os_dims) ){
    cerr << "Unable to allocate device memory for oversampled image" << endl;
    exit(1);
  }
  
  // Extract coil sensitivity maps and training data using all the data
  for( unsigned int iteration = 0; iteration < num_profiles/profiles_per_reconstruction; iteration++ ) {
    
    // Define trajectories
    boost::shared_ptr< cuNDArray<floatd2::Type> > traj = compute_radial_trajectory_golden_ratio_2d<float>
      ( samples_per_profile, profiles_per_frame, frames_per_reconstruction, iteration*profiles_per_reconstruction );
    
    // Preprocess
    plan.preprocess( traj.get(), NFFT_plan<float,2>::NFFT_PREP_BACKWARDS );
    traj.reset();
    
    // Upload data
    boost::shared_ptr< cuNDArray<float_complext::Type> > csm_data = 
      upload_data( iteration*profiles_per_reconstruction, samples_per_profile, samples_per_reconstruction,
		   num_profiles*samples_per_profile, get_num_coils(), &host_samples );
    
    // Accumulate k-space for CSM estimation
    plan.convolve( csm_data.get(), image_os, dcw.get(), NFFT_plan<float,2>::NFFT_BACKWARDS, (iteration==0) ? false : true );
    csm_data.reset();
  }
  
  // We now have 'frames_per_reconstruction' k-space images of each coil. Add these up.
  boost::shared_ptr< cuNDArray<float_complext::Type> > acc_image_os = cuNDA_sum<float_complext::Type>( image_os, 2 );    
  delete image_os; image_os = 0x0;
  
  // Complete gridding of k-space CSM image
  plan.fft( acc_image_os.get(), NFFT_plan<float,2>::NFFT_BACKWARDS );
  plan.deapodize( acc_image_os.get() );
  
  // Remove oversampling
  vector<unsigned int> image_dims = uintd_to_vector<2>(get_matrix_size()); image_dims.push_back(get_num_coils());
  cuNDArray<float_complext::Type> *image = new cuNDArray<float_complext::Type>(); 
  if( !image->create(image_dims) ){
    cerr << "Unable to allocate device memory for image" << endl;
    exit(1);
  }
  cuNDA_crop<float_complext::Type,2>( (get_matrix_size_os()-get_matrix_size())>>1, acc_image_os.get(), image );
  acc_image_os.reset();
  
  // Estimate CSM
  csm = estimate_b1_map<float,2>( image );

  // Setup encoding and regularization operators
  rhs_buffer->set_csm(csm);
  E->setup( get_matrix_size(), get_matrix_size_os(), get_kernel_width() ); 
  R->compute( image, uintd_to_vector<2>(get_matrix_size()), cg.get_cublas_handle() );
  delete image; image = 0x0; 

  // Define preconditioning weights
  update_preconditioning_weights();
  
  progress.setValue(4);

  if( E->set_csm(csm) < 0 ) {
    cout << "Failed to set csm on encoding matrix" << endl;
  }
  
  progress.setValue(5);

  // Trigger the #projections slot
  reconWidget->numProjectionsScrollBar->setValue(reconWidget->numProjectionsScrollBar->value()+1);
    
  // Perform reconstruction
  reconstruct();
}

void radialSenseAppMainWindow::update_preconditioning_weights()
{
  boost::shared_ptr< cuNDArray<float> > _precon_weights = cuNDA_ss<float,float_complext::Type>( csm.get(), 2 );
  cuNDA_axpy<float>( get_kappa(), R->get(), _precon_weights.get(), cg.get_cublas_handle() );  
  cuNDA_reciprocal_sqrt<float>( _precon_weights.get() );
  boost::shared_ptr< cuNDArray<float_complext::Type> > precon_weights = cuNDA_real_to_complext<float>( _precon_weights.get() );
  D->set_weights( precon_weights );
}

void radialSenseAppMainWindow::projectionsPerFrameChanged(int)
{
  // the integer is an 'id' not the slider value!

  unsigned int value = get_num_projections_per_frame();

  // Enforce even values
  if( value%2 ){
    value--;
    reconWidget->numProjectionsScrollBar->setValue(value);
    return;
  }

  // Remove the Qt lag of the slider rendering
  QApplication::processEvents();

  // The range of the frames slider/spinbox has changed
  unsigned int maxCentralProjection = get_maximum_central_projection();
  reconWidget->projectionSelectionScrollBar->setMaximum(maxCentralProjection);
  reconWidget->projectionNumberSpinBox->setMaximum(maxCentralProjection);
  reconWidget->projectionSelectionScrollBar->setSingleStep(value>>2);
  reconWidget->projectionNumberSpinBox->setSingleStep(value>>2);

  const unsigned int samples_per_profile = get_num_samples_per_projection();
  const unsigned int profiles_per_frame = get_num_projections_per_frame();
  
  // Density compensation weights are constant throughout all reconstrutions
  dcw  = compute_radial_dcw_golden_ratio_2d
    ( samples_per_profile, profiles_per_frame, (float)get_matrix_size_os().vec[0]/(float)get_matrix_size().vec[0], 
      get_one<float>()/((float)samples_per_profile/(float)max(get_matrix_size().vec[0],get_matrix_size().vec[1])) );
  
  // Set density compensation weights
  if( E->set_dcw(dcw) < 0 ) {
    cout << "Failed to set density compensation weights on encoding matrix" << endl;
  }
  
  // Reconstruct
  reconstruct();
}

void radialSenseAppMainWindow::centralProjectionChanged(int id)
{
  // the integer is an 'id' not the slider value!

  // Enforce even values
  unsigned int value = get_central_projection();
  if( value%2 ){
    value--;
    reconWidget->projectionSelectionScrollBar->setValue(value);
    return;
  }

  // Remove the lag of the slider rendering
  QApplication::processEvents();

  // Perform reconstruction
  reconstruct();
}

void radialSenseAppMainWindow::matrixSizeChanged()
{
  static unsigned int lastValue = MATRIX_SIZE_INITIAL_VALUE;

  unsigned int value = matrixSizeSpinBox->value();
  unsigned int value_os = oversampledMatrixSizeSpinBox->value();

  if( value == lastValue )
    return;
  else 
    lastValue = value;

  // Pass matrix size to GLReconWidget
  reconWidget->openglCanvas->setMatrixSize( value, value );
	
  if( value_os < value ){
    oversampledMatrixSizeSpinBox->setValue(value);
  }
  
  // and encoding matrix
  E->setup( get_matrix_size(), get_matrix_size_os(), get_kernel_width() );  
  
  replan();
}

void radialSenseAppMainWindow::matrixSizeOSChanged()
{
  static unsigned int lastValue = MATRIX_SIZE_OS_INITIAL_VALUE;

  unsigned int value = matrixSizeSpinBox->value();
  unsigned int value_os = oversampledMatrixSizeSpinBox->value();

  if( value_os == lastValue )
    return;
  else 
    lastValue = value_os;

  if( value_os < value ){
    oversampledMatrixSizeSpinBox->setValue(value);
    return;
  }
	
  if( value_os%2 ){
    value_os++;
    oversampledMatrixSizeSpinBox->setValue(value_os);
    return;
  }

  E->setup( get_matrix_size(), get_matrix_size_os(), get_kernel_width() );  

  reconstruct();
}

void radialSenseAppMainWindow::kernelWidthChanged()
{
  static double lastValue = KERNEL_SIZE_INITIAL_VALUE;

  double value = kernelSizeSpinBox->value();

  if( value == lastValue )
    return;
  else 
    lastValue = value;

  E->setup( get_matrix_size(), get_matrix_size_os(), get_kernel_width() );  
  
  reconstruct();
}

void radialSenseAppMainWindow::numIterationsChanged()
{
  static unsigned int lastValue = NUM_ITERATIONS_INITIAL_VALUE;

  unsigned int value = numIterationsSpinBox->value();

  if( value == lastValue )
    return;
  else 
    lastValue = value;

  cg.set_iterations( get_num_iterations() );

  reconstruct();
}

void radialSenseAppMainWindow::regularizationWeightChanged()
{
  static double lastValue = REG_WEIGHT_INITIAL_VALUE;

  double value = regularizationWeightSpinBox->value();

  if( value == lastValue )
    return;
  else 
    lastValue = value;

  // Update D
  update_preconditioning_weights();
  
  // Update operator R 
  R->set_weight( get_kappa() );

  reconstruct();
}

void radialSenseAppMainWindow::windowScaleChanged(double)
{
  reconstruct();
}

/*
  Reconstruct frame
*/

void radialSenseAppMainWindow::reconstruct()
{
  
  // Check if any data has been loaded
  if( host_samples.get_number_of_elements() == 0 )
    return;
  
  // See if there is any uncaught errors before starting
  cudaError_t err;
  err = cudaGetLastError();
  if( err != cudaSuccess ){
    QMessageBox::critical( this, tr("Cuda error"), tr(cudaGetErrorString(err)) );
    actionExit->trigger();
  }

  // Map result to OpenGL
  reconWidget->openglCanvas->mapPBO();

  // Be optimistic...
  bool success = true;

  const unsigned int samples_per_profile = get_num_samples_per_projection();
  const unsigned int num_profiles = get_num_points_per_array_coil() / samples_per_profile;
  const unsigned int profiles_per_frame = get_num_projections_per_frame();
  const unsigned int frames_per_reconstruction = 1; 
  const unsigned int profiles_per_reconstruction = get_num_projections_per_frame()*frames_per_reconstruction;
  const uintd2::Type matrix_size = get_matrix_size();
  const uintd2::Type matrix_size_os = get_matrix_size_os();
  const unsigned int num_coils = get_num_coils();
  const unsigned int samples_per_reconstruction = profiles_per_reconstruction*samples_per_profile;

  // Determine trajectories
  boost::shared_ptr< cuNDArray<floatd2::Type> > traj = compute_radial_trajectory_golden_ratio_2d<float>
    ( samples_per_profile, profiles_per_frame, frames_per_reconstruction,  get_first_projection() );
  
  // Upload data
  boost::shared_ptr< cuNDArray<float_complext::Type> > data = 
    upload_data( get_first_projection(), samples_per_profile, samples_per_reconstruction,
		 num_profiles*samples_per_profile, num_coils, &host_samples );
    
    // Set current trajectory and trigger NFFT preprocessing
    if( E->preprocess(traj.get()) < 0 ) {
      cout << "Failed to set trajectory on encoding matrix" << endl;
    }
        
    // Form rhs (use result array to save memory)
    vector<unsigned int> rhs_dims = uintd_to_vector<2>(matrix_size); rhs_dims.push_back(frames_per_reconstruction);
    cuNDArray<float_complext::Type> rhs; rhs.create(rhs_dims);
    E->mult_MH( data.get(), &rhs );
    
    //
    // Conjugate gradient solver
    //

    boost::shared_ptr< cuNDArray<float_complext::Type> > cgresult = cg.solve(&rhs);

    // Magnitudes image for visualization
    boost::shared_ptr< cuNDArray<float> > tmp_res = cuNDA_norm<float,float_complext::Type>(cgresult.get());
    cuNDA_normalize( tmp_res.get(), get_window_scale(), cg.get_cublas_handle() );

    // Copy to OpenGL/pbo
    cudaMemcpy( reconWidget->openglCanvas->getDevPtr(),
		tmp_res->get_data_ptr(),
		prod(matrix_size)*sizeof(float), cudaMemcpyDeviceToDevice );
       
    // Report any errors not already caught...
    err = cudaGetLastError();
    if( err != cudaSuccess ){
      QMessageBox::critical( this, tr("Cuda error"), tr(cudaGetErrorString(err)) );
      actionExit->trigger();
    }
    
    reconWidget->openglCanvas->unmapPBO();
    
    if( !success ){
      QMessageBox::critical( this, tr("Reconstruction error"), tr("Check console. Quitting.") );
      actionExit->trigger();
      exit(EXIT_FAILURE);
    }
    
    reconWidget->openglCanvas->updateGL();
}

/*
  "Gets..."
*/

uintd2::Type radialSenseAppMainWindow::get_matrix_size()
{
  int value = matrixSizeSpinBox->value();
  return uintd2( value, value );
}

uintd2::Type radialSenseAppMainWindow::get_matrix_size_os()
{
  int value = oversampledMatrixSizeSpinBox->value();
  return uintd2( value, value );
}

float radialSenseAppMainWindow::get_kernel_width()
{
  double value = kernelSizeSpinBox->value();
  return (float) value;	
}

float radialSenseAppMainWindow::get_window_scale()
{
  double value = windowScaleSpinBox->value();
  return (float) value;	
}

unsigned int radialSenseAppMainWindow::get_num_samples_per_projection()
{
  if( host_samples.get_number_of_dimensions() > 0 )
    return host_samples.get_size(0);
  else return 0;
}

unsigned int radialSenseAppMainWindow::get_first_projection()
{
  int value = reconWidget->projectionNumberSpinBox->value();
  value -= get_num_projections_per_frame()>>1;
  if( value<0 )
    value = 0;
  return value;
}

unsigned int radialSenseAppMainWindow::get_central_projection()
{
  int value = reconWidget->projectionSelectionScrollBar->value();
  return value;
}

unsigned int radialSenseAppMainWindow::get_maximum_central_projection()
{
  if( get_num_samples_per_projection() == 0 )
    return 0;
	
  unsigned int maxCentralProjection = get_num_points_per_array_coil()/get_num_samples_per_projection()-get_num_projections_per_frame()/2-get_num_projections_per_frame()%2;
  return maxCentralProjection;
}

unsigned int radialSenseAppMainWindow::get_num_projections_per_frame()
{
  int value = reconWidget->numProjectionsSpinBox->value();
  return value;
}

unsigned int radialSenseAppMainWindow::get_num_coils()
{
  if( host_samples.get_number_of_dimensions() < 3 )
    return 0;

  unsigned int val;
  if( host_samples.get_number_of_dimensions() == 3 )
    val = host_samples.get_size(2);
  else{
    printf("\nUnknown number of dimensions in dataset. Quitting.\n");
    exit(1);
  }
  
  return val;
}

unsigned int radialSenseAppMainWindow::get_num_points_per_reconstruction()
{
  unsigned int val = get_num_samples_per_projection()*get_num_projections_per_frame();
  return val;
}

hoNDArray<floatd2::Type>* radialSenseAppMainWindow::get_sample_values_array()
{
  return &host_samples;
}

unsigned int radialSenseAppMainWindow::get_num_points_per_array_coil()
{
  if(host_samples.get_number_of_dimensions()<2)
    return 0;

  unsigned int val = host_samples.get_size(0)*host_samples.get_size(1);
  return val;
}

unsigned int radialSenseAppMainWindow::get_num_iterations()
{
  int value = numIterationsSpinBox->value();
  return value;
}

inline float radialSenseAppMainWindow::get_kappa()
{
  double value = regularizationWeightSpinBox->value();
  return (float)value;
}

// Upload samples for one reconstruction from host to device
boost::shared_ptr< cuNDArray<float_complext::Type> > 
radialSenseAppMainWindow::upload_data( unsigned int profile_offset, unsigned int samples_per_profile, unsigned int samples_per_reconstruction, 
				       unsigned int total_samples_per_coil, unsigned int num_coils,
				       hoNDArray<float_complext::Type> *host_data )
{
  vector<unsigned int> dims; dims.push_back(samples_per_reconstruction); dims.push_back(num_coils);
  cuNDArray<float_complext::Type> *data = new cuNDArray<float_complext::Type>(); 
  if( !data->create( dims ) ){
    cerr << "Unable to allocate device memory for samples" << endl;
    exit(1);
  }
  
  for( unsigned int i=0; i<num_coils; i++ )
    cudaMemcpy( data->get_data_ptr()+i*samples_per_reconstruction, 
		host_data->get_data_ptr()+i*total_samples_per_coil+profile_offset*samples_per_profile, 
		samples_per_reconstruction*sizeof(float_complext::Type), cudaMemcpyHostToDevice );
  
  return boost::shared_ptr< cuNDArray<float_complext::Type> >(data);
}
