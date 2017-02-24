/*=========================================================================
 *
 * Restartable FastMarching Image Filter utilized to segment the airway
 *
 * based on the original ITK implementation
 *=========================================================================*/

#ifndef __itkFastMarchingRestartableImageFilter_txx
#define __itkFastMarchingRestartableImageFilter_txx

#include "itkFastMarchingRestartableImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"
#include <algorithm>
#include "log.h"
#include <assert.h>


namespace itk
{
template< class TLevelSet, class TSpeedImage >
FastMarchingRestartableImageFilter< TLevelSet, TSpeedImage >
::FastMarchingRestartableImageFilter()
{
  this->ProcessObject::SetNumberOfRequiredInputs(0);

  OutputSizeType outputSize;
  outputSize.Fill(16);
  typename LevelSetImageType::IndexType outputIndex;
  outputIndex.Fill(0);

  m_OutputRegion.SetSize(outputSize);
  m_OutputRegion.SetIndex(outputIndex);

  m_OutputOrigin.Fill(0.0);
  m_OutputSpacing.Fill(1.0);
  m_OutputDirection.SetIdentity();
  m_OverrideOutputInformation = false;

  m_AlivePoints = NULL;
  m_OutsidePoints = NULL;
  m_TrialPoints = NULL;
  m_ProcessedPoints = NULL;

  m_SpeedConstant = 1.0;
  m_InverseSpeed = -1.0;
  m_LabelImage = LabelImageType::New();

  m_LargeValue = static_cast< PixelType >( NumericTraits< PixelType >::max() / 2.0 );
  m_StoppingValue = static_cast< double >( m_LargeValue );
  m_CollectPoints = false;

  m_NormalizationFactor = 1.0;
  m_Restarted = false;
  latest_time = 0.0;
}

template< class TLevelSet, class TSpeedImage >
void
FastMarchingRestartableImageFilter< TLevelSet, TSpeedImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Alive points: " << m_AlivePoints.GetPointer() << std::endl;
  os << indent << "Trial points: " << m_TrialPoints.GetPointer() << std::endl;
  os << indent << "Speed constant: " << m_SpeedConstant << std::endl;
  os << indent << "Stopping value: " << m_StoppingValue << std::endl;
  os << indent << "Large Value: "
     << static_cast< typename NumericTraits< PixelType >::PrintType >( m_LargeValue )
     << std::endl;
  os << indent << "Normalization Factor: " << m_NormalizationFactor << std::endl;
  os << indent << "Collect points: " << m_CollectPoints << std::endl;
  os << indent << "OverrideOutputInformation: ";
  os << m_OverrideOutputInformation << std::endl;
  os << indent << "OutputRegion: " << m_OutputRegion << std::endl;
  os << indent << "OutputOrigin: " << m_OutputOrigin << std::endl;
  os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;
  os << indent << "OutputDirection: " << m_OutputDirection << std::endl;
}

template <class TLevelSet, class TSpeedImage>
std::vector< itk::Point< unsigned int, itk::FastMarchingRestartableImageFilter<TLevelSet, TSpeedImage>::SetDimension > >
FastMarchingRestartableImageFilter<TLevelSet,TSpeedImage>
::GetCurrentTrials()
{
	/*
	 * Calculate the points which are trial in this moment and return them.
	 * Useful to now which is the current wavefront
	 */

  pLogType l = get_logger("FM::GetCurrentTrials");
  l->infoStream() << "start()";

  std::vector< itk::Point< unsigned int, SetDimension > >  points;

  /*
   * This implementation is very low because it loops over the whole image each time
   * now that we have a more efficient HashMap implementation let's use it!
   *

  typedef ImageRegionIterator< LabelImageType > LabelIterator;
  typedef itk::Point<unsigned int, SetDimension> PointType;

  LabelIterator typeIt( m_LabelImage,
                        m_LabelImage->GetBufferedRegion() );

  for ( typeIt.GoToBegin(); !typeIt.IsAtEnd(); ++typeIt )
    {
    if ( typeIt.Get() == TrialPoint )
      {
      PointType p;
      for (int i=0;i<SetDimension;i++)
        p.SetElement(i, typeIt.GetIndex().GetElement(i));
      points.push_back(p);
      }

    }
	*/

  l->infoStream() << "HashMap size: " << HashLabelMap.size();
  std::map<std::string, int>::iterator it = HashLabelMap.begin();
  while (it!=HashLabelMap.end()){

	  if (it->second == TrialPoint) {
		  itk::Point< unsigned int, SetDimension > idx = str2point(it->first);
		  //std::cout << idx << " ";
		  points.push_back(idx);
	  }
	  it++;
  }


  return points;

}

template <class TLevelSet, class TSpeedImage>
std::vector< itk::Point< unsigned int, itk::FastMarchingRestartableImageFilter<TLevelSet, TSpeedImage>::SetDimension > >
FastMarchingRestartableImageFilter<TLevelSet,TSpeedImage>
::GetCurrentAlives()
{
/*
 * Calculate the points which are alive in this moment and return them.
 * Useful to now which is the current segmentation
 * NOTE: a quicker way to get the segmentation is to threshold the output image
 */

	pLogType l = get_logger("FM::GetCurrentAlives");
	  l->infoStream() << "start()";
  std::vector< itk::Point< unsigned int, SetDimension > >  points;

  typedef ImageRegionIterator< LabelImageType > LabelIterator;
  typedef itk::Point<unsigned int, SetDimension> PointType;

  LabelIterator typeIt( m_LabelImage,
                        m_LabelImage->GetBufferedRegion() );

  for ( typeIt.GoToBegin(); !typeIt.IsAtEnd(); ++typeIt )
    {
    if ( typeIt.Get() == AlivePoint )
      {
      PointType p;
      for (int i=0;i<SetDimension;i++)
        p.SetElement(i, typeIt.GetIndex().GetElement(i));
      points.push_back(p);
      }

    }

  return points;

}


template< class TLevelSet, class TSpeedImage >
void
FastMarchingRestartableImageFilter< TLevelSet, TSpeedImage >
::GenerateOutputInformation()
{
  // copy output information from input image
  Superclass::GenerateOutputInformation();

  // use user-specified output information
  if ( this->GetInput() == NULL || m_OverrideOutputInformation )
    {
    LevelSetPointer output = this->GetOutput();
    output->SetLargestPossibleRegion(m_OutputRegion);
    output->SetOrigin(m_OutputOrigin);
    output->SetSpacing(m_OutputSpacing);
    output->SetDirection(m_OutputDirection);
    }
}

template< class TLevelSet, class TSpeedImage >
void
FastMarchingRestartableImageFilter< TLevelSet, TSpeedImage >
::EnlargeOutputRequestedRegion(
  DataObject *output)
{
  // enlarge the requested region of the output
  // to the whole data set
  TLevelSet *imgData;

  imgData = dynamic_cast< TLevelSet * >( output );
  if ( imgData )
    {
    imgData->SetRequestedRegionToLargestPossibleRegion();
    }
  else
    {
    // Pointer could not be cast to TLevelSet *
    itkWarningMacro( << "itk::FastMarchingImageFilter"
                     << "::EnlargeOutputRequestedRegion cannot cast "
                     << typeid( output ).name() << " to "
                     << typeid( TLevelSet * ).name() );
    }
}

template< class TLevelSet, class TSpeedImage >
void
FastMarchingRestartableImageFilter< TLevelSet, TSpeedImage >
::AllocateOutput(LevelSetImageType *output)
{
	// allocate memory for the output buffer
	  output->SetBufferedRegion( output->GetRequestedRegion() );
	  output->Allocate();

	  // cache some buffered region information
	  m_BufferedRegion = output->GetBufferedRegion();
	  m_StartIndex = m_BufferedRegion.GetIndex();
	  m_LastIndex = m_StartIndex + m_BufferedRegion.GetSize();
	  typename LevelSetImageType::OffsetType offset;
	  offset.Fill(1);
	  m_LastIndex -= offset;

	  // allocate memory for the PointTypeImage
	  m_LabelImage->CopyInformation(output);
	  m_LabelImage->SetBufferedRegion( output->GetBufferedRegion() );
	  m_LabelImage->Allocate();

	  // set all output value to infinity
	  typedef ImageRegionIterator<LevelSetImageType>  OutputIterator;

	  OutputIterator outIt ( output, output->GetBufferedRegion() );

	  PixelType outputPixel;
	  outputPixel = m_LargeValue;

	  outIt.GoToBegin();

	  while( !outIt.IsAtEnd() )
	    {
	    outIt.Set(outputPixel);
	    ++outIt;
	    }

	  /*
	   * Set all points type to FarPoint
	   * without this the propagation will not grow further because we propagate only if the points
	   * in the label map are not AlivePoint
	   *
	   */

	  typedef ImageRegionIterator< LabelImageType > LabelIterator;

	  LabelIterator typeIt( m_LabelImage,
	                        m_LabelImage->GetBufferedRegion() );


	  typeIt.GoToBegin();
	  while( !typeIt.IsAtEnd() )
	    {
	    typeIt.Set(FarPoint);
	    ++typeIt;
	    }

}

template< class TLevelSet, class TSpeedImage >
void
FastMarchingRestartableImageFilter< TLevelSet, TSpeedImage >
::SetupPropagation(LevelSetImageType *output)
{

	pLogType l = get_logger("FM::SetupPropagation");
	l->infoStream() << "start()";

	  PixelType outputPixel;

	  // process input alive points
	  AxisNodeType node;
	  NodeIndexType idx;

	  if ( m_AlivePoints )
	    {
		  l->debugStream() << "Processing alive points";
	    typename NodeContainer::ConstIterator pointsIter = m_AlivePoints->Begin();
	    typename NodeContainer::ConstIterator pointsEnd = m_AlivePoints->End();

	    while( pointsIter != pointsEnd )
	      {
	      // get node from alive points container
	      node = pointsIter.Value();
	      idx = node.GetIndex();

	      // check if node index is within the output level set
	      if ( m_BufferedRegion.IsInside( idx ) )
	        {
	        // make this an alive point
	        m_LabelImage->SetPixel(idx, AlivePoint);
	        HashLabelMap[idx2str(idx)] = AlivePoint;

	        outputPixel = node.GetValue();
	        output->SetPixel(idx, outputPixel);
	        l->debugStream() << "Set output point: " << idx << " = " << outputPixel;
	        }

	      ++pointsIter;
	      }
	    }

	  if( m_OutsidePoints )
	    {
		  l->debugStream() << "Processing outside points";
	    typename NodeContainer::ConstIterator pointsIter = m_OutsidePoints->Begin();
	    typename NodeContainer::ConstIterator pointsEnd = m_OutsidePoints->End();

	    while( pointsIter != pointsEnd )
	      {
	      // get node from alive points container
	      node = pointsIter.Value();
	      idx = node.GetIndex();

	      // check if node index is within the output level set
	      if ( m_BufferedRegion.IsInside( idx ) )
	        {
	        // make this an outside point
	        m_LabelImage->SetPixel(idx, OutsidePoint );
	        HashLabelMap[idx2str(idx)] = OutsidePoint;

	        outputPixel = node.GetValue();
	        output->SetPixel(idx, outputPixel);
	        l->debugStream() << "Set output point: " << idx << " = " << outputPixel;
	        }

	      ++pointsIter;
	      }
	    }



	  // process the input trial points
	  if ( m_TrialPoints )
	    {
		  l->debugStream() << "Clearing " << m_TrialHeap.size() << " heap points... ";
		  // make sure the heap is empty
		 	  while ( !m_TrialHeap.empty() )
		 	    {
		 		//std::cout << "Removed from trial: " << m_TrialHeap.top() << std::endl;
		 	    m_TrialHeap.pop();
		 	    }
		  l->debugStream() << "Processing input trial points: " << m_TrialPoints->Size();
	    typename NodeContainer::ConstIterator pointsIter = m_TrialPoints->Begin();
	    typename NodeContainer::ConstIterator pointsEnd = m_TrialPoints->End();

	    while( pointsIter != pointsEnd )
	      {
	      // get node from trial points container
	      node = pointsIter.Value();
	      idx = node.GetIndex();

	      // check if node index is within the output level set
	      if ( m_BufferedRegion.IsInside( idx ) )
	        {
	        // make this an initial trial point
	        m_LabelImage->SetPixel(idx, InitialTrialPoint);
	        HashLabelMap[idx2str(idx)] = InitialTrialPoint;

	        outputPixel = node.GetValue();
	        output->SetPixel(idx, outputPixel);
	        //std::cout << "Set output point: " << idx << " = " << outputPixel << std::endl;

	        m_TrialHeap.push(node);
	        }

	      ++pointsIter;
	      }
	    }


}

template< class TLevelSet, class TSpeedImage >
void
FastMarchingRestartableImageFilter< TLevelSet, TSpeedImage >
::Initialize(LevelSetImageType *output)
{

	pLogType l = get_logger("FM::Initialize");
	  l->infoStream() << "start()";
  this->AllocateOutput(output);
  this->SetupPropagation(output);

}

template< class TLevelSet, class TSpeedImage >
void
FastMarchingRestartableImageFilter< TLevelSet, TSpeedImage >
::ClearTrials(void) {
	/*
	 * Set all trials to Alive, effectively blocking the propagation.
	 * Used to switch propagation wavefront
	 */
	pLogType l = get_logger("FM::ClearTrials");
	  l->infoStream() << "start()";
	  typedef ImageRegionIterator< LabelImageType > LabelIterator;

	  LabelIterator typeIt( m_LabelImage,
							m_LabelImage->GetBufferedRegion() );

	  typeIt.GoToBegin();
	  while( !typeIt.IsAtEnd() )
		{
		  if (typeIt.Get() == TrialPoint)
			  typeIt.Set(AlivePoint);
		++typeIt;
		}
}



template <class TLevelSet, class TSpeedImage>
void
FastMarchingRestartableImageFilter<TLevelSet,TSpeedImage>
::RestartPropagation(  )
{
/*
	  typedef ImageRegionIterator< LabelImageType > LabelIterator;

	  LabelIterator typeIt( m_LabelImage,
	                        m_LabelImage->GetBufferedRegion() );


	  typeIt.GoToBegin();
	  while( !typeIt.IsAtEnd() )
	    {
	    typeIt.Set(FarPoint);
	    ++typeIt;
	    }*/


   this->SetupPropagation(this->GetOutput());

}


template< class TLevelSet, class TSpeedImage >
void
FastMarchingRestartableImageFilter< TLevelSet, TSpeedImage >
::GenerateData()
{

  //Clear hash map so it contains only last points
  HashLabelMap.clear();
  pLogType l = get_logger("FM::GenerateData");
  if( m_NormalizationFactor < vnl_math::eps )
    {
    ExceptionObject err(__FILE__, __LINE__);
    err.SetLocation(ITK_LOCATION);
    err.SetDescription("Normalization Factor is null or negative");
    throw err;
    }

  LevelSetPointer output = this->GetOutput();
  SpeedImageConstPointer speedImage = this->GetInput();


  if ( !m_Restarted ) {

	  l->infoStream() << "Initializing a new propagation...";
	  this->Initialize(output);
	  this->latest_time = 0.0;
  }

  if ( m_CollectPoints )
    {
    m_ProcessedPoints = NodeContainer::New();
    }

  // process points on the heap
  AxisNodeType node;
  double currentValue;
  double oldProgress = 0;

  this->UpdateProgress(0.0); // Send first progress event

  l->infoStream() << "Current time is " << this->get_last_arrival_time() << " and stopping time is: " << m_StoppingValue;
  while ( !m_TrialHeap.empty() )
    {
    // get the node with the smallest value
    node = m_TrialHeap.top();
    m_TrialHeap.pop();

    //l->debugStream() << "pop! " << node.GetIndex() << " (still " << m_TrialHeap.size() << " points )";
    // does this node contain the current value ?
    currentValue = static_cast< double >( output->GetPixel( node.GetIndex() ) );


    if ( node.GetValue() == currentValue )
      {
    	l->debugStream() << "node value (" << node.GetValue() << ") == currentValue  (" << currentValue << ")";
      // is this node already alive ?
      if ( m_LabelImage->GetPixel( node.GetIndex() ) != AlivePoint )
        {
    	  l->debugStream() << "Node is not currently Alive";

    	  if ((currentValue > this->latest_time) && (currentValue != m_LargeValue))
    		  this->latest_time = currentValue;


        if ( currentValue > m_StoppingValue )
          {
        	l->debugStream() << "Node value "<< currentValue << " is higher than stopping value " <<  m_StoppingValue;
          this->UpdateProgress(1.0);

          break;
          }



        if ( m_CollectPoints )
          {
          m_ProcessedPoints->InsertElement(m_ProcessedPoints->Size(), node);
          }

        // set this node as alive
        m_LabelImage->SetPixel(node.GetIndex(), AlivePoint);
        HashLabelMap[idx2str(node.GetIndex())] = AlivePoint;

        // update its neighbors
        this->UpdateNeighbors(node.GetIndex(), speedImage, output);

        // Send events every certain number of points.
        const double newProgress = currentValue / m_StoppingValue;
        if ( newProgress - oldProgress > 0.01 ) // update every 1%
          {

          this->UpdateProgress(newProgress);
          oldProgress = newProgress;
          if ( this->GetAbortGenerateData() )
            {
            this->InvokeEvent( AbortEvent() );
            this->ResetPipeline();
            ProcessAborted e(__FILE__, __LINE__);
            e.SetDescription("Process aborted.");
            e.SetLocation(ITK_LOCATION);
            throw e;
            }
          }
        } else {
        	l->debugStream() << "Node is already Alive! Not doing anything...";
        }
      } else {
    	  l->debugStream() << "node value (" << node.GetValue() << ") <> currentValue  (" << currentValue << ")";
      }
    }
}

template< class TLevelSet, class TSpeedImage >
void
FastMarchingRestartableImageFilter< TLevelSet, TSpeedImage >
::UpdateNeighbors(
  const IndexType & index,
  const SpeedImageType *speedImage,
  LevelSetImageType *output)
{
	pLogType l = get_logger("FM::UpdateNeighbors");
	l->debugStream() << "start()";


  IndexType neighIndex = index;
  unsigned char label;

  for ( unsigned int j = 0; j < SetDimension; j++ )
    {
    // update left neighbor
    if ( index[j] > m_StartIndex[j] )
      {
      neighIndex[j] = index[j] - 1;
      }

    label = m_LabelImage->GetPixel(neighIndex);

    if ( ( label != AlivePoint ) &&
         ( label != InitialTrialPoint ) &&
         ( label != OutsidePoint ) )
      {
      this->UpdateValue(neighIndex, speedImage, output);
      } else
    	  l->debugStream() << "Won't update a point with label " << (int) label;

    // update right neighbor
    if ( index[j] < m_LastIndex[j] )
      {
      neighIndex[j] = index[j] + 1;
      }

    label = m_LabelImage->GetPixel(neighIndex);

    if ( ( label != AlivePoint ) &&
         ( label != InitialTrialPoint ) &&
         ( label != OutsidePoint ) )
      {
      this->UpdateValue(neighIndex, speedImage, output);
      } else
    	  l->debugStream() << "Won't update a point with label " << (int) label;

    //reset neighIndex
    neighIndex[j] = index[j];
    }
}

template< class TLevelSet, class TSpeedImage >
double
FastMarchingRestartableImageFilter< TLevelSet, TSpeedImage >
::UpdateValue(
  const IndexType & index,
  const SpeedImageType *speedImage,
  LevelSetImageType *output)
{
	pLogType l = get_logger("FM::UpdateValue");
l->debugStream() << "start()";
  IndexType neighIndex = index;

  PixelType neighValue;

  // just to make sure the index is initialized (really cautious)
  AxisNodeType node;
  node.SetIndex( index );

  for ( unsigned int j = 0; j < SetDimension; j++ )
    {
    node.SetValue(m_LargeValue);

    // find smallest valued neighbor in this dimension
    for ( int s = -1; s < 2; s = s + 2 )
      {
      neighIndex[j] = index[j] + s;

      // make sure neighIndex is not outside from the image
      if ( ( neighIndex[j] > m_LastIndex[j] ) ||
           ( neighIndex[j] < m_StartIndex[j] ) )
        {
        continue;
        }

      if ( m_LabelImage->GetPixel( neighIndex ) == AlivePoint )
        {
        neighValue = static_cast< PixelType >( output->GetPixel(neighIndex) );

        // let's find the minimum value given a direction j
        if ( node.GetValue() > neighValue )
          {
          node.SetValue(neighValue);
          node.SetIndex(neighIndex);
          }
        }
      }

    // put the minimum neighbor onto the heap
    m_NodesUsed[j] = node;
    m_NodesUsed[j].SetAxis(j);

    // reset neighIndex
    neighIndex[j] = index[j];
    }

  // sort the local list
  std::sort(m_NodesUsed, m_NodesUsed + SetDimension);

  // solve quadratic equation
  double solution = static_cast< double >( m_LargeValue );

  double aa( 0.0 );
  double bb( 0.0 );
  double cc( m_InverseSpeed );

  if ( speedImage )
    {
    cc = static_cast< double >( speedImage->GetPixel(index) ) / m_NormalizationFactor;
    cc = -1.0 * vnl_math_sqr(1.0 / cc);
    }

  OutputSpacingType spacing = this->GetOutput()->GetSpacing();

  double discrim;

  for ( unsigned int j = 0; j < SetDimension; j++ )
    {
    node = m_NodesUsed[j];
    const double value = static_cast< double >( node.GetValue() );
    l->debugStream() << "Value of the current node: " << value;

    if ( solution >= value )
      {
      const int axis = node.GetAxis();
      // spaceFactor = \frac{1}{spacing[axis]^2}
      const double spaceFactor = vnl_math_sqr(1.0 / spacing[axis]);
      aa += spaceFactor;
      bb += value * spaceFactor;
      cc += vnl_math_sqr(value) * spaceFactor;

      discrim = vnl_math_sqr(bb) - aa * cc;
      if ( discrim < 0.0 )
        {
        // Discriminant of quadratic eqn. is negative
        ExceptionObject err(__FILE__, __LINE__);
        err.SetLocation(ITK_LOCATION);
        err.SetDescription("Discriminant of quadratic equation is negative");
        throw err;
        }

      solution = ( vcl_sqrt(discrim) + bb ) / aa;
      l->debugStream() << "solution = " << solution;
      }
    else
      {
      break;
      }
    }

  if ( solution < m_LargeValue )
    {

	  l->debugStream() << "creating new trials " ;
    // write solution to m_OutputLevelSet
    PixelType outputPixel = static_cast< PixelType >( solution );
    output->SetPixel(index, outputPixel);

    // insert point into trial heap
    m_LabelImage->SetPixel(index, TrialPoint);
    HashLabelMap[idx2str(index)] = TrialPoint;
    node.SetValue( outputPixel );
    node.SetIndex( index );
    m_TrialHeap.push(node);
    } else
    	l->debugStream() << "solution too large and no new trials are created! ";

  return solution;
}
} // namespace itk

#endif

