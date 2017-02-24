

#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"

#include "itkShapeLabelObject.h"
#include "itkLabelMap.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkShapeLabelMapFilter.h"
#include "itkShapeOpeningLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"

#include "itkBinaryShapeOpeningImageFilter.h"


typedef unsigned short USPixelType;
typedef signed short SSPixelType;
const unsigned int Dimension=2;

typedef unsigned long LabelType;

typedef itk::Image<USPixelType, Dimension>                              USImageType;
typedef itk::Image<SSPixelType, Dimension>                              SSImageType;

typedef itk::ImageFileReader< USImageType >                             EmphMaskReaderType;
typedef itk::ImageFileReader< SSImageType >                             LungReaderType;
typedef itk::ImageRegionConstIterator<USImageType>                        USConstIteratorType;
typedef itk::ImageRegionConstIterator<SSImageType>                        SSConstIteratorType;


typedef itk::ShapeLabelObject<LabelType, Dimension>                                     LabelObjectType;
typedef itk::LabelMap<LabelObjectType>                                    LabelCollectionType; 
typedef itk::BinaryImageToLabelMapFilter<USImageType, LabelCollectionType> ConverterType;
typedef itk::ShapeLabelMapFilter<LabelCollectionType>                               ShapeFilterType;

typedef itk::ShapeOpeningLabelMapFilter<LabelCollectionType>                   ShapeOpeningFilterType;
typedef itk::LabelMapToLabelImageFilter< LabelCollectionType, USImageType > L2IType;
typedef itk::ImageFileWriter< USImageType >                             USWriterType;



int main( int argc, char * argv[] )
{ 

if (argc < 2){
    std::cout << "Usage: " << argv[0] << " enf_mask_file_name [size_th]" << std::endl;
    return EXIT_SUCCESS;
    }
    
    if ( (strcmp(argv[1],"--help") == 0) || ( (strcmp(argv[1],"-h") == 0) ) ){
        std::cout << "Usage: " << argv[0] << " enf_mask_file_name [size_th]" << std::endl;
        return EXIT_SUCCESS;
        }
        
        
        const char * emphfilepath  = argv[1];
        const char * lavfile = "lavdata.txt";
        const unsigned int enf_th = 255;
        unsigned int size_th = 0;
        
        if (argc == 3)
            size_th = atoi(argv[2]);
        
        
        EmphMaskReaderType::Pointer emph_mask_reader = EmphMaskReaderType::New();
        
        
        emph_mask_reader->SetFileName(emphfilepath);
        
        double spacing[Dimension];
                
        
            
        try {
            emph_mask_reader->Update();
        } catch (itk::ExceptionObject & excp) {
            std::cerr << "  Problems reading the emphysema mask image"   << std::endl;
            std::cerr << excp << std::endl;
            return EXIT_FAILURE;
        }
            
            
        
        
        //Reading the lungs data
        //std::cout << "Spacing is: " <<  spacing[0] << " " <<  spacing[1] << " " << spacing[2] << std::endl;
            
            
 
        
            
        ConverterType::Pointer converter =  ConverterType::New();
        converter->SetInput(emph_mask_reader->GetOutput());
        converter->SetInputForegroundValue(enf_th);
        //Any nonzero value should be read as True
        converter->SetFullyConnected(0);
                
        ShapeFilterType::Pointer shape = ShapeFilterType::New();
        shape->SetInput( converter->GetOutput() );
        shape->Update();
        LabelCollectionType::Pointer collection = shape->GetOutput();
                
        //std::cout << "Now saving lav sizes on file " << lavfile << "..."<< std::endl;
        
        std::ofstream lav_file_op(lavfile);
        unsigned long emph_cluster_count = 0;
        
        for(unsigned int label=1; label<collection->GetNumberOfLabelObjects(); label++ )
        {
            LabelObjectType::Pointer labelObject = collection->GetLabelObject( label );
            //lav_file_op << labelObject->GetSize() * spacing[0] * spacing[1] * spacing[2] << " ";
            lav_file_op << labelObject->GetNumberOfPixels() << " ";
            
            if (size_th > 0) {
                if (labelObject->GetNumberOfPixels() > size_th) {
                    //std::cout << "Label: " << label << " Size: " << labelObject->GetSize() * spacing[0] * spacing[1] * spacing[2] << " Centroid: " << labelObject->GetCentroid() << std::endl;
                    emph_cluster_count+=labelObject->GetNumberOfPixels();
                }
            }
           
        }
        lav_file_op.close();
        
        
        return EXIT_SUCCESS;
        
        }
                
