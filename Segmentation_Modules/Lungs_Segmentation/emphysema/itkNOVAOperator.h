/*
 * itkNOVAOperator.h
 *
 *  Created on: May 25, 2012
 *      Author: mceresa
 */


#ifndef __itkNOVAOperator_h
#define __itkNOVAOperator_h

#include "itkExceptionObject.h"
#include "itkNeighborhoodOperator.h"

namespace itk {


template<class TPixel, unsigned int VDimension=2,
  class TAllocator = NeighborhoodAllocator<TPixel> >
class ITK_EXPORT NOVAOperator
  : public NeighborhoodOperator<TPixel, VDimension, TAllocator>
{
public:
  /** Standard typedefs */
  typedef NOVAOperator                                         Self;
  typedef NeighborhoodOperator<TPixel, VDimension, TAllocator>  Superclass;

  itkTypeMacro(NOVAOperator, NeighborhoodOperator);

  NOVAOperator() {}
  NOVAOperator(const Self& other)
    : NeighborhoodOperator<TPixel, VDimension, TAllocator>(other)
    {}


  /** Creates the operator with a specified radius ("square", same length
   * on each side). The spatial location of the coefficients within the
   * operator is defined by the subclass implementation of the Fill method.
   * \sa CreateDirectional \sa Fill */
  // virtual void CreateToRadius(const unsigned long);
  /**
   * Assignment operator
   */
  Self &operator=(const Self& other)
  {
    Superclass::operator=(other);
    return *this;
  }
  /**
   * Prints some debugging information
   */
  virtual void PrintSelf(std::ostream &os, Indent i) const
    {
    os << i << "SobelOperator { this=" << this  << "}" << std::endl;
    Superclass::PrintSelf(os, i.GetNextIndent());
    }

protected:
  /**
   * Typedef support for coefficient vector type.  Necessary to
   * work around compiler bug on VC++.
   */
  typedef typename Superclass::CoefficientVector CoefficientVector;
  typedef typename Superclass::PixelType         PixelType;

  /**
   * Calculates operator coefficients.
   */
  CoefficientVector GenerateCoefficients();

  /**
   * Arranges coefficients spatially in the memory buffer.
   */
  void Fill(const CoefficientVector &c);

};

} // namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_NOVAOperator(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT NOVAOperator< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef NOVAOperator< ITK_TEMPLATE_2 x > \
                                                  SobelOperator##y; } \
  }


#if ITK_TEMPLATE_TXX
# include "itkNOVAOperator.txx"
#endif

#endif




