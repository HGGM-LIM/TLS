/*
 * itkNOVAOperator.h
 *
 *  Created on: May 25, 2012
 *      Author: mceresa
 */


#ifndef __itkNOVAOperator_txx
#define __itkNOVAOperator_txx

#include "itkNOVAOperator.h"
#include "itkObject.h"

namespace itk
{

template <class TPixel, unsigned int VDimension, class TAllocator>
void
NOVAOperator <TPixel, VDimension, TAllocator>
::Fill(const CoefficientVector &coeff)
{
  this->InitializeToZero();

  // Note that this code is only good for 2d and 3d operators.  Places the
  // coefficients in the exact center of the neighborhood
  unsigned int i;
  int x,y,z, pos;
  unsigned int center = this->GetCenterNeighborhoodIndex();

  if (VDimension == 2)
    {
    i = 0;
    for (y = -1; y <= 1; y++ )
      {
      for (x = -1; x <= 1; x++)
        {
        pos = center + y * this->GetStride(1) + x * this->GetStride(0);
        this->operator[](pos) = static_cast<TPixel>(coeff[i]);
        i++;
        }
      }
    }
  else
    {
    itkExceptionMacro( << "The ND version of the NOVA operator is not yet implemented.  Currently only the 2D and 3D versions are available." );
    }

}

template <class TPixel, unsigned int VDimension, class TAllocator>
typename NOVAOperator<TPixel, VDimension, TAllocator>
::CoefficientVector
 NOVAOperator<TPixel, VDimension, TAllocator>
::GenerateCoefficients()
{
  std::vector<double> coeff;
  if (VDimension == 2)
    {
    coeff.push_back(1.0);  coeff.push_back(1.0);  coeff.push_back(1.0);
    coeff.push_back(1.0);  coeff.push_back(0.0);  coeff.push_back(1.0);
    coeff.push_back(1.0);  coeff.push_back(1.0);  coeff.push_back(1.0);
    }
  else
    {
    itkExceptionMacro( << "The ND version of the Nova operator has not been implemented.  Currently only 2D and 3D versions are available." );
    }

  return coeff;
}

} // namespace itk

#endif





