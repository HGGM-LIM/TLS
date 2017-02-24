/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: GradientMagnitudeImageFilter.cxx,v $
  Language:  C++
  Date:      $Date: 2005/08/31 13:55:21 $
  Version:   $Revision: 1.28 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif


#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkCurvatureAnisotropicDiffusionImageFilter.h"

//#include "itkClosingByReconstructionImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include "itkSpeedCalculatorImageFilter.h"
#include "itkCastImageFilter.h"


int main( int argc, char * argv[] )
{
  if( argc < 3 ) 
    { 
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  inputImage  outputImage " << std::endl;
    return EXIT_FAILURE;
    }

  // typedefs  
  typedef    signed short    InputPixelType;
  typedef    float    OutputPixelType;
  
  typedef itk::Image< InputPixelType,  3 >   InputImageType;
  typedef itk::Image< OutputPixelType, 3 >   OutputImageType;

  typedef itk::ImageFileReader< InputImageType >  ReaderType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
  
  // reader and writer
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );


  // filter
  /*typedef itk::BinaryBallStructuringElement< InputPixelType, 3 > BallType_3D;
  typedef itk::ClosingByReconstructionImageFilter< InputImageType,
  InputImageType, BallType_3D> FilterType;

  FilterType::Pointer closing = FilterType::New();
  BallType_3D kernel;
  kernel.SetRadius(1);
  kernel.CreateStructuringElement();
  closing->SetKernel(kernel);*/
  
  // filter
  typedef itk::CurvatureAnisotropicDiffusionImageFilter<InputImageType, OutputImageType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetNumberOfIterations( 2 );
  filter->SetTimeStep( 0.01 );
  filter->SetConductanceParameter( 3 );
  filter->UseImageSpacingOn();

  // speed image calculator
  typedef itk::SpeedCalculatorImageFilter< OutputImageType,
  OutputImageType >  SpeedFilterType;
  SpeedFilterType::Pointer speedFilter = SpeedFilterType::New();
      
  // define pipeline: reader->filter->speedFilter->writer
  filter->SetInput(  reader->GetOutput());
  filter->Update();

  // the scale parameter should be proportional 
  // to the expected airway size and it is
  // given in image spacing units (it is the 
  // sigma of the hessian filter)
  speedFilter->SetScale(2);
  speedFilter->SetInput( filter->GetOutput() );
  writer->SetInput( speedFilter->GetOutput() );
  
  // execute pipeline
  writer->Update();
  
}
