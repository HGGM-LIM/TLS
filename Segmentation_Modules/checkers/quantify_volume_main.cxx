
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkMedianImageFilter.h"
#include "itkImageRegion.h"
#include "itkImageRegionConstIterator.h"

int main( int argc, char * argv[] )
{ 

  if (argc < 3){
    std::cout << "Usage: " << argv[0] << " file_name smooth_image" << std::endl;
    return EXIT_SUCCESS;
  }
  
  if ( (strcmp(argv[1],"--help") == 0) || ( (strcmp(argv[1],"-h") == 0) ) ){
    std::cout << "Usage: " << argv[0] << " file_name smooth_image" << std::endl;
    return EXIT_SUCCESS;
  }
  
  char * infilepath  = argv[1];
  int smoothMe  = atoi(argv[2]); //Whether to smooth or not
  
  std::cout << "Reading input serie  " << infilepath << std::endl;
  
  typedef signed short PixelType;
  const unsigned int Dimension=3;

  typedef unsigned long LabelType;

  typedef itk::Image<PixelType, Dimension>                       		ImageType;
  
  typedef ImageType::ConstPointer                         InputImageConstPointer;
  typedef itk::ImageRegionConstIterator<ImageType>                    	ConstIteratorType;
  typedef ImageType::IndexType                                     	IndexType;
  typedef ImageType::SizeType                                      	SizeType;
  typedef ImageType::RegionType                                    	RegionType;
  
  typedef itk::ImageFileReader< ImageType >                      	 	ReaderType;
	ReaderType::Pointer reader =                                     	ReaderType::New();
	
  typedef itk::MedianImageFilter< ImageType, ImageType >  MedianFilterType;
 
  
  double spacing[Dimension];
	
  //Reading the lungs data
  reader->SetFileName(infilepath);
  reader->Update();
  
  spacing[0] = reader->GetOutput()->GetSpacing()[0];
  spacing[1] = reader->GetOutput()->GetSpacing()[1];
  spacing[2] = reader->GetOutput()->GetSpacing()[2];
  
  InputImageConstPointer imagePtr = reader->GetOutput();
  
  // How big is the input image?
  RegionType inputRegion = imagePtr->GetLargestPossibleRegion();
  SizeType size = inputRegion.GetSize();
  IndexType startIndex = inputRegion.GetIndex();
  
  if (smoothMe) {

    MedianFilterType::Pointer median = MedianFilterType::New();
    
    ImageType::SizeType indexRadius;
    indexRadius[0] = 1; // radius along x
    indexRadius[1] = 1; // radius along y
    indexRadius[2] = 1; // radius along z
    median->SetRadius( indexRadius );
    median->SetInput(imagePtr);
    median->Update();

    imagePtr = median->GetOutput();
    
  } 
  
  ConstIteratorType in(imagePtr, inputRegion);
  
  float m_VolLung = 0;
    
  for (in.GoToBegin();!in.IsAtEnd();++in)
  {
    if (in.Get() != 0)
    {
    	m_VolLung++;
    }
  }

  std::cout << "Spacing: " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << std::endl;
  std::cout << "Volume (cc): " << m_VolLung/1000*spacing[0]*spacing[1]*spacing[2] << std::endl;
  
    
  return EXIT_SUCCESS;
  
}

