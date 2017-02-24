/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAdaptiveBsplinesImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/07/17 12:26:43 $
  Version:   $Revision: 1.4 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkFindTracheaImageFilter_txx
#define _itkFindTracheaImageFilter_txx

#include "itkFindTracheaImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include "itkImageRegionConstIteratorWithIndex.h"
//#include "itkImageRegionConstIterator.h"

#include "itkHuThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"

#include "itkImageFileWriter.h"

#include "math.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_matrix_inverse.h"



namespace itk
{

template <typename TInputImage, typename TOutputImage>
void
FindTracheaImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}

template <typename TInputImage, typename TOutputImage>
FindTracheaImageFilter<TInputImage, TOutputImage>
::FindTracheaImageFilter() 
{
  m_tracheaLabel = 0;  
}

template< typename TInputImage, typename TOutputImage >
void
FindTracheaImageFilter< TInputImage, TOutputImage >
::GenerateInputRequestedRegion ()
{
   // call the superclass's implementation of this method
  Superclass::GenerateInputRequestedRegion();

  if ( !this->GetInput() )
    {
    return;
    }

  // get pointers to the input 
  InputImagePointer  inputPtr  =   const_cast< InputImageType *>( this->GetInput() );

  // Request the entire input image
  InputImageRegionType inputRegion;
  inputRegion = inputPtr->GetLargestPossibleRegion();
  
  inputPtr->SetRequestedRegion(inputRegion);

  return;
}

template< typename TInputImage, typename TOutputImage >
bool FindTracheaImageFilter< TInputImage, TOutputImage >
::IsNotComputed( InputPixelType value)
{
  typename std::vector<Label>::iterator it, end;
  it = m_labelVector.begin();
  end = m_labelVector.end();
  
  // loop to check if labelID is in label vector
  while(it!=end)
    {
    if(value==(*it).LabelID)
      {
      return false;
      }
    ++it;
    }
  
  return true;
}

template< typename TInputImage, typename TOutputImage >
void 
FindTracheaImageFilter< TInputImage, TOutputImage >
::ComputeAndAddLabelStatistics( InputPixelType label, InputImageConstPointer labelImage)
{
  float size = 0;
  float CI = 0;
  float centerX = 0;
  float centerY = 0;
  typedef ImageRegionConstIterator< TInputImage> InputIteratorType;

  // iterator to go through label image
  InputIteratorType labelIt = InputIteratorType( labelImage, labelImage->GetLargestPossibleRegion() );
   
  for(labelIt.GoToBegin(); !labelIt.IsAtEnd(); ++labelIt)
    {
    if(labelIt.Get()==label)
      {
      size++;
      centerX = centerX + labelIt.GetIndex()[0];
      centerY = centerY + labelIt.GetIndex()[1];
      }
       
    }
  Label newLabel;
  newLabel.LabelID = label;
  newLabel.size = size;
  newLabel.centerX = centerX/size;
  newLabel.centerY = centerY/size;
  m_labelVector.push_back(newLabel);
  
  std::cout << "Region " << label << " size: " << size << " , center: " << newLabel.centerX << ", " << newLabel.centerY << std::endl; 
  
}


template< typename TInputImage, typename TOutputImage >
void
FindTracheaImageFilter< TInputImage, TOutputImage >
::GenerateData()
{
  
  // get input and output images
  InputImageConstPointer input = this->GetInput();
  OutputImagePointer output = this->GetOutput();
  typename TInputImage::SizeType imageSize = input->GetLargestPossibleRegion().GetSize();
  // allocate memory for output
  this->AllocateOutputs(); 
  
  // define iterators
  typedef ImageRegionConstIterator< TInputImage> InputIteratorType;
  typedef ImageRegionIterator< TOutputImage> OutIteratorType;
  OutIteratorType outIt = OutIteratorType( output, input->GetLargestPossibleRegion() );
  
  // apply Hu threshold on image
  typedef HuThresholdImageFilter<TInputImage, TInputImage> HuThresholdFilterType;
  typename HuThresholdFilterType::Pointer huThreshold = HuThresholdFilterType::New();

   huThreshold->SetInput( input );
   
   // connected components on binarized image
   typedef ConnectedComponentImageFilter<TInputImage, TInputImage> ConnectedComponentFilterType;
   typename ConnectedComponentFilterType::Pointer connectedComponent = ConnectedComponentFilterType::New();
   connectedComponent->SetInput( huThreshold->GetOutput() );
   connectedComponent->Update();

   InputImageConstPointer labelImage = connectedComponent->GetOutput();
   
   // iterator to go through label image
   InputIteratorType labelIt = InputIteratorType( labelImage, input->GetLargestPossibleRegion() );
   
   // look in middle region for trachea
   

   typename TInputImage::PixelType value;   
   
   // loop through image computing sizes and centers
   for(labelIt.GoToBegin(); !labelIt.IsAtEnd(); ++labelIt)
     {
     value = labelIt.Get();
     if(this->IsNotComputed(value))
       {
       this->ComputeAndAddLabelStatistics(value, labelImage);
	}
     }
     
   // loop over regions to find trachea
   // it must have certain size and position
   typename std::vector<Label>::iterator it, end;
   it = m_labelVector.begin();
   end = m_labelVector.end();
   Label candidate;
   bool search = true;
   while(it!=end && search)
     {
     candidate = *it;
     if(100 < candidate.size && candidate.size< 5000)
       {
       if(imageSize[0]*0.25 < candidate.centerX < imageSize[1]*0.75)
         {
         if(imageSize[1]*0.25 < candidate.centerY < imageSize[1]*0.75)
	   {
	   m_tracheaLabel = candidate.LabelID;
	   search = false;
	   }
	 }
       }
     ++it;
     }
    // error message if trachea not found
    if(m_tracheaLabel==0)
      {
      std::cerr << "Error! Trachea not found in first slice" << std::endl;
      exit(0);
      }
      
   // write output image with only trachea
   for(labelIt.GoToBegin(), outIt.GoToBegin(); !labelIt.IsAtEnd(); ++labelIt, ++outIt)
     {
     if(labelIt.Get()==m_tracheaLabel)
       {
       outIt.Set(1);
       }
     else
       {
       outIt.Set(0);
       }
     }
   std::cout << "Trachea is label: " << m_tracheaLabel << std::endl;
   
   // to check
   typedef itk::ImageFileWriter < OutputImageType >  WriterType;
   typename WriterType::Pointer writer = WriterType::New();
   writer->SetInput(output);
   writer->SetFileName("trachea.tiff");
   writer->Update();

}
  
  


}// end namespace itk

#endif
