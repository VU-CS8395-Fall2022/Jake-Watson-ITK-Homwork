#include <iostream>
//#include <itkSpatialObjectToImageFilter.h>
//for homework
#include <itkAffineTransform.h>
#include <itkStatisticsImageFilter.h>
//#include <itkBSplineTransform.h>

//libraries  here
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
//#include <itkImageSeriesReader.h>
//#include <itkImageSeriesWriter.h>
#include <itkImageRegistrationMethod.h>
//#include <itkMultiResolutionImageRegistrationMethod.h>
//#include <itkLinearInterpolateImageFunction.h>
#include <itkMeanSquaresImageToImageMetric.h>
#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkResampleImageFilter.h>
#include <itkCommand.h>

//example
#include <itkBSplineTransform.h>
#include <itkBSplineTransformInitializer.h>
//#include <itkLBFGS2Optimizerv4.h>
//#include <itkBSplineTransformParametersAdaptor.h>
//#include <itkSquaredDifferenceImageFilter.h>
#include <itkIdentityTransform.h>
//#include <itkBSplineTransformInitializer.h"
//#include <itkTransformToDisplacementFieldFilter.h>

//mesh experiment for volume
// #include <itkSimplexMesh.h>
// #include <itkRegularSphereMeshSource.h>
// #include <itkTriangleMeshToSimplexMeshFilter.h>
// #include <itkSimplexMeshVolumeCalculator.h>
// #include <itkTriangleMeshToSimplexMeshFilter.h>
#include <itkBinaryThresholdImageFilter.h>

// things we could/should do differently
// if we are inter-subject, A_T1 to B_T1, affine isn't good enough, should be deformable
// if we are intra-subject, A_T1 to A_T2, sum of square differences is not a good cross-modality metric, try correlation or mutual information

// declare our image type 
typedef itk::Image < int, 3 > myImageType ;

// todo: implement this class

//a class is a user defined data type with its own data members and member functions
// class MyObserverCommandClass : public itk::Command //the class is called 'MyObserverCommandClass'
// {
// public:
//   itkNewMacro ( MyObserverCommandClass ) ;

//   typedef itk::RegularStepGradientDescentOptimizer OptimizerType ;
//   //typedef OptimizerType * OptimizerPointerType ;
//   void Execute ( itk::Object *caller, const itk::EventObject & event ) override
//   {
//     // std::cout << "HelloWorld" << std::endl ;
//     OptimizerType::Pointer optimizer = dynamic_cast < OptimizerType * > ( caller ) ;
//     unsigned int iterationNumber = optimizer->GetCurrentIteration () ;
//     double currentMetricValue = optimizer->GetValue () ;
//     std::cout << iterationNumber << " " << currentMetricValue << std::endl ;
//     return ;
//   }

//   void Execute ( const itk::Object *caller, const itk::EventObject & event ) override
//   {
//     std::cout << "Hello world  -const" << std::endl ;
//     return ;
//   }
// } ;

//FUNCTION FOR THRESHOLDING IMAGE
myImageType::Pointer ThresholdImage ( myImageType::Pointer inputImage, int lowT, int highT, int inside, int outside )
{
  typedef itk::BinaryThresholdImageFilter < myImageType, myImageType > ThresholdFilterType ;
  ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New() ;
  thresholdFilter->SetInput ( inputImage ) ;
  if ( lowT != -1 )
    {
      thresholdFilter->SetLowerThreshold ( lowT ) ;
    }
  if ( highT != -1 )
    {
      thresholdFilter->SetUpperThreshold ( highT ) ;
    }
  thresholdFilter->SetInsideValue ( inside ) ;
  thresholdFilter->SetOutsideValue ( outside ) ;
  thresholdFilter->Update() ;

  return thresholdFilter->GetOutput() ;
};


// class MyRegistrationObserverCommandClass : public itk::Command 
// {
// public:
//   itkNewMacro ( MyRegistrationObserverCommandClass ) ;

//   typedef itk::RegularStepGradientDescentOptimizer OptimizerType ;
//   typedef itk::MultiResolutionImageRegistrationMethod < myImageType, myImageType > RegistrationMethod ;
//   //typedef OptimizerType * OptimizerPointerType ;
//   void Execute ( itk::Object *caller, const itk::EventObject & event ) override
//   {
//     RegistrationMethod::Pointer regMethod = dynamic_cast < RegistrationMethod * > ( caller ) ;
//     int currentLevel = regMethod->GetCurrentLevel() ;
//     std::cout << "Starting registration level: " << currentLevel << std::endl ;

//     // optimizer is no longer the caller but we can get to it thru the regmethod
//     // need the optimizer so we can change reg settings per level, e.g. step length or nIterations
//     OptimizerType::Pointer optimizer = dynamic_cast < OptimizerType * > ( regMethod->GetModifiableOptimizer () ) ;

//     if ( currentLevel == 0 )
//       {
// 	optimizer->SetNumberOfIterations ( 60 ) ;
// 	std::cout << "Min: " << optimizer->GetMinimumStepLength () << std::endl ;
// 	std::cout << "Max: " << optimizer->GetMaximumStepLength () << std::endl ;
// 	optimizer->SetMinimumStepLength ( 0 ) ;
// 	optimizer->SetMaximumStepLength ( 0.125 ) ; // looks ok enough
//       }
//     else if ( currentLevel == 1 )
//       {
// 	optimizer->SetNumberOfIterations ( 20 ) ;
// 	std::cout << "Min: " << optimizer->GetMinimumStepLength () << std::endl ;
// 	std::cout << "Max: " << optimizer->GetMaximumStepLength () << std::endl ;
// 	optimizer->SetMinimumStepLength ( 0 ) ;
// 	optimizer->SetMaximumStepLength ( 0.0625 ) ;
//       }
//     else
//       {
// 	optimizer->SetNumberOfIterations ( 10 ) ;
// 	std::cout << "Min: " << optimizer->GetMinimumStepLength () << std::endl ;
// 	std::cout << "Max: " << optimizer->GetMaximumStepLength () << std::endl ;
// 	optimizer->SetMinimumStepLength ( 0 ) ;
// 	optimizer->SetMaximumStepLength ( 0.125 ) ;
//       }   
//   }

//   void Execute ( const itk::Object *caller, const itk::EventObject & event ) override
//   {
//     // std::cout << "HelloWorld" << std::endl ;
//   }
// };

//MAIN
int main (int argc, char **argv)
{
  if ( argc <= 7 )
    {
      std::cout << "Usage: " << argv[0] << " " << "inputMovingFileName inputLabelsAtlas inputFixedFileName1 inputFixedFileName2 inputFixedFileName3 inputFixedFileName4 inputFixedFileName5  outputLaff1 outputLaff2 outputLaff3 outputLaff4 outputLaff5 outputLm1 outputLm2 outputLm3 outputLm4 outputLm5" << std::endl ;
      return 0 ;
    }
  
  // declare our image reader type
  typedef itk::ImageFileReader < myImageType > myFileReaderType ;

  // read the files in

  const int i = 3;
  //this is the image_atlas file
  myFileReaderType::Pointer movingReader = myFileReaderType::New() ;
  movingReader->SetFileName ( argv[1] ) ;
  movingReader->Update() ;
  
  //this is the labels_atlas file
  myFileReaderType::Pointer labelsReader = myFileReaderType::New() ;
  labelsReader->SetFileName (argv[2] ) ;
  labelsReader->Update() ;

  //loop through the code for the five fixed images Im (subject images)
  for(int i=3;i<8;i++)
  {
     
  myFileReaderType::Pointer fixedReader = myFileReaderType::New() ;
  fixedReader->SetFileName ( argv[i] ) ;
  fixedReader->Update() ;

  // the registration algorithm

  // all the typedefs, typedef basically allows you to rename data something easier for you to read
  typedef itk::ImageRegistrationMethod < myImageType, myImageType > RegistrationMethod ;
  //typedef itk::MultiResolutionImageRegistrationMethod < myImageType, myImageType > RegistrationMethod ;
  typedef itk::AffineTransform < double, 3 > TransformType ;
  typedef itk::LinearInterpolateImageFunction < myImageType, double > InterpolatorType ;
  typedef itk::MeanSquaresImageToImageMetric < myImageType, myImageType > MetricType ;
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType ;
  //
  // all the pointers, a pointer is a variable that stores a memory address
  RegistrationMethod::Pointer regMethod = RegistrationMethod::New() ;
  TransformType::Pointer transform = TransformType::New() ;
  InterpolatorType::Pointer interpolator = InterpolatorType::New() ;
  MetricType::Pointer metric = MetricType::New() ;
  OptimizerType::Pointer optimizer = OptimizerType::New() ;
  
  // parameter setup for registration method
  // (pointer name) -> (variable name), giving pointers variable names
  regMethod->SetMovingImage ( movingReader->GetOutput() ) ;
  regMethod->SetFixedImage ( fixedReader->GetOutput() ) ;
  regMethod->SetTransform ( transform ) ;
  regMethod->SetInterpolator ( interpolator ) ;
  regMethod->SetMetric ( metric ) ;
  regMethod->SetOptimizer ( optimizer ) ;
  
  transform->SetIdentity() ;
  regMethod->SetInitialTransformParameters ( transform->GetParameters() ) ;
  //regMethod->SetNumberOfLevels ( 3 ) ;
  //regMethod->SetFixedImageRegion ( fixedReader->GetOutput()->GetLargestPossibleRegion() ) ;
  // look at transform params that might need to be set up
  //transform->SetIdentity() ;
  //regMethod->SetInitialTransformParameters ( transform->GetParameters() ) ;
  metric->SetTransform ( transform ) ;
  //metric->Initialize() ;
  //metric->SetNumberOfSpatialSamples(10000L ) ;
  // look at interpolator params that might need to be set up - looks ok to leave alone
   // look at metric params that might need to be set up
  //metric->SetTransform ( transform ) ;
  //metric->Initialize() ;
  // leaving alone for now - the reg method is taking care of these for us when we call update
    // optimizer params that might need to be set up

  optimizer->MinimizeOn() ;
  optimizer->SetNumberOfIterations ( 30 ) ;
  //std::cout << "Min: " << optimizer->GetMinimumStepLength () << std::endl ;
  //std::cout << "Max: " << optimizer->GetMaximumStepLength () << std::endl ;
  optimizer->SetMinimumStepLength ( 0 ) ;
  optimizer->SetMaximumStepLength ( 0.2 ) ;

  // set up an observer for the optimizer
  //MyObserverCommandClass::Pointer customCommand = MyObserverCommandClass::New() ;
  //optimizer->AddObserver ( itk::IterationEvent (), customCommand ) ;

  // set up second observer for the registration method
  //MyRegistrationObserverCommandClass::Pointer regObserver = MyRegistrationObserverCommandClass::New() ;
  //regMethod->AddObserver ( itk::IterationEvent(), regObserver ) ;
  
  // run the registration method
  regMethod->Update() ; 
  //std::cout << "GradMagTolerance: " << optimizer->GetGradientMagnitudeTolerance() << std::endl ;
  //std::cout << "Stop condition: " << optimizer->GetStopCondition() << std::endl ;

  //I believe this technically completes the affine registration of image_atlas to the five images
  
  // get the transform from regmethod and apply it to the moving image //ipek
  // HW: apply the transform from the affine registration to the image Iatlas and labels Latlas, to create new files Iaff,m atlas and Laff,m atlas
  typedef itk::ResampleImageFilter < myImageType, myImageType > ResampleFilterType ;
  ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New() ;
  resampleFilter->SetInput ( movingReader->GetOutput() ) ;
  resampleFilter->SetReferenceImage ( fixedReader->GetOutput() ) ;
  resampleFilter->UseReferenceImageOn() ;
  resampleFilter->SetDefaultPixelValue ( 0 ) ;

  transform->SetParameters ( regMethod->GetLastTransformParameters () ) ;
  resampleFilter->SetTransform ( transform ) ;
  resampleFilter->Update () ;

  //binarize transformed Labels Atlas for writing image to display in ITK
  typedef itk::ImageFileWriter < myImageType > myFileWriterType ;
  myFileWriterType::Pointer myFileWriter = myFileWriterType::New() ;
  //myFileWriter->SetFileName ( argv[i+5] ) ;
  //myFileWriter->SetInput ( threshIm->GetOutput() ) ; 
  //myFileWriter->Write();
  
  resampleFilter->SetInput ( labelsReader->GetOutput() ) ;
  resampleFilter->SetReferenceImage ( fixedReader->GetOutput() ) ;
  resampleFilter->UseReferenceImageOn() ;
  resampleFilter->SetDefaultPixelValue ( 0 ) ;

  transform->SetParameters ( regMethod->GetLastTransformParameters () ) ;
  resampleFilter->SetTransform ( transform ) ;
  resampleFilter->Update () ;
  
  myImageType::Pointer laffm = resampleFilter->GetOutput() ;
  myImageType::Pointer thresholdedImage = ThresholdImage ( laffm, 1, -1, 1, 0 );
  
  myFileWriter->SetFileName ( argv[i+5] ) ;
  myFileWriter->SetInput ( thresholdedImage ) ;
  myFileWriter->Write() ;
  
  
//   //PART 3
//   //BSpline for Deformable registration
   typedef itk::BSplineTransform < double, 3, 3 > transformTypeTwo ;
   transformTypeTwo::Pointer transformTwo = transformTypeTwo::New() ;

   myImageType::Pointer fixedIm = (fixedReader->GetOutput());
   
   transformTypeTwo::ParametersType params;// = transformTypeTwo::New() ;
   const unsigned int numberOfParams = transformTwo->GetNumberOfParameters();
   transformTypeTwo::ParametersType parameters(numberOfParams);
   parameters.Fill(0.0);

   regMethod->SetMovingImage ( resampleFilter->GetOutput() );
   regMethod->SetTransform ( transformTwo ) ;

   transformTwo->SetIdentity();

   transformTypeTwo::PhysicalDimensionsType fixedPhysicalDimensions;
   transformTypeTwo::MeshSizeType meshSize;
   transformTypeTwo::OriginType fixedOrigin;

   unsigned int numberOfGridNodesInOneDimension = 15;
   for (unsigned int j = 0; j < 3; j++)
     {
       fixedOrigin[j] = fixedIm->GetOrigin()[j];
       fixedPhysicalDimensions[j] = fixedIm->GetSpacing()[j] * static_cast<double>(fixedIm->GetLargestPossibleRegion().GetSize()[j] - 1 );
     }

   meshSize.Fill( numberOfGridNodesInOneDimension - 3 ); 
//   //parameters
//   typedef itk::BSplineTransformInitializer < transformTypeTwo, myImageType > InitializerType ;
//   InitializerType::Pointer transformInitializer = InitializerType::New();
//   //transformInitializer->SetTransform(transformTwo);
//   //typedef itk::LBFGS2Optimizerv4 bOptimizerType ;
//   //bOptimizerType::Pointer bOptimizer = bOptimizerType::New();

//   //Jackie help on optimizer
//   //use regulargradientstepoptimizer
//   //addobserver optional, iterations 1-- max step .1 min step .0001
   optimizer->SetNumberOfIterations(100);
   optimizer->SetMinimumStepLength(0.001);
   optimizer->SetMaximumStepLength(0.1);
//   regMethod->Update();
   
    transformTwo->SetTransformDomainOrigin( fixedOrigin );
    transformTwo->SetTransformDomainPhysicalDimensions( fixedPhysicalDimensions );
    transformTwo->SetTransformDomainDirection( fixedIm->GetDirection()  );
    transformTwo->SetTransformDomainMeshSize( meshSize );

    transformTwo->SetParameters( parameters );
    regMethod->SetInitialTransformParameters ( transformTwo->GetParameters() );

    //set metric
    metric->SetTransform( transformTwo );
    metric->SetNumberOfSpatialSamples ( 10000L );
    regMethod->Update();
//registration method updated for bspline - ready to go


//   // // PART 4
//   // //take the created transform and and apply to the atlas Labels
//   // //make sure resampleFilter is setup for the Labels file,  need to change transform parameters for Bspline
//   // //reference image should still be Im
    ResampleFilterType::Pointer lm = ResampleFilterType::New() ;
    lm->SetTransform(transformTwo);
    lm->SetInput(laffm);
    lm->SetReferenceImage(fixedReader->GetOutput() );
    lm->UseReferenceImageOn() ;
    lm->SetDefaultPixelValue(0);
//   resampleFilter->SetInput ( labelsReader->GetOutput() );
    transformTwo->SetParameters (regMethod->GetLastTransformParameters () );
    lm->Update() ;
//   resampleFilter->SetTransform ( transformTwo ) ;
    
//   // //binarize before writing file for display in ITK
   myImageType::Pointer lmout = lm->GetOutput() ;
   myImageType::Pointer threshImTwo = ThresholdImage( lmout, 1, -1, 1, 0 );
   myFileWriter->SetInput ( threshImTwo ) ;
   myFileWriter->SetFileName (argv[i+10]) ;
   myFileWriter->Write() ;
   
   //PART 5
   //use the statistics image filter library (from tips and tricks in assignment document) to compute the volume of the segmentation

   /* 
   myImageType::Pointer LmOut = labelsReader->GetOutput() ;
   typedef itk::StatisticsImageFilter < myImageType > statsim; 
   statsim::Pointer statisticsImFilter = statsim::New() ;
   statisticsImFilter->SetInput(LmOut) ;
   statisticsImFilter->Update() ;

   myImageType::SpacingType space ;
   space->GetSpacing();
   int  volumeLM = space[0] * space[1] * space[2];
   std::cout << " The volume of  segmentation: "<< volumeLM  << std::endl;
   */
 }
   return 0;
 }
