

#include "itkGDCMImageIO.h"
#include "itkDICOMSeriesFileNames.h"

#include "itkImage.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char * argv[] )
{
 
  char * infilename  = argv[1];

  typedef signed short PixelType;
  const unsigned int Dimension=3;

  typedef unsigned long LabelType;

  typedef itk::Image< PixelType, Dimension >                       	ImageType;
	
  typedef ImageType::IndexType                                     	IndexType;
  typedef ImageType::SizeType                                      	SizeType;
  typedef ImageType::RegionType                                    	OutputImageRegionType;
 
  typedef itk::ImageSeriesReader< ImageType >                      	ReaderType;
	ReaderType::Pointer reader =                                     	ReaderType::New();
 
  typedef itk::ImageFileWriter< ImageType >                        	WriterType;
	WriterType::Pointer writer =                                     	WriterType::New();

//--------- Reading the serie --------------
  
  typedef itk::GDCMImageIO                                         ImageIOType;
  ImageIOType::Pointer dicomIO =                                   ImageIOType::New();

  typedef itk::DICOMSeriesFileNames NamesGeneratorType;
	
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
  nameGenerator->SetFileNameSortingOrderToSortByImageNumber();

  nameGenerator->SetDirectory( infilename );
  typedef std::vector<std::string> seriesIdContainer;
  const seriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
	
  seriesIdContainer::const_iterator seriesItr = seriesUID.begin();
  seriesIdContainer::const_iterator seriesEnd = seriesUID.end();
  
  std::string seriesIdentifier;
  seriesIdentifier = seriesUID.begin()->c_str();
	
  typedef std::vector<std::string> fileNamesContainer;
  fileNamesContainer fileNames;
  fileNames = nameGenerator->GetFileNames(seriesIdentifier);

  reader->SetFileNames( fileNames );
  reader->SetImageIO( dicomIO );

  try {
    reader->Update();
  } catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught: can't read serie" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }
  
  typedef itk::MetaDataDictionary   DictionaryType;

  ImageType::Pointer inputImage = reader->GetOutput();
  DictionaryType & dictionary = dicomIO->GetMetaDataDictionary();
  inputImage->SetMetaDataDictionary(dictionary);

  writer->SetInput( inputImage );
  writer->SetFileName("./0input.dcm" );

  
  try
    {
  	writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught: can't write back serie to a single file" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

	
   
  return EXIT_SUCCESS;
  
}
