

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

#include "dicom_utilities.h"


typedef unsigned short USPixelType;
typedef unsigned short SSPixelType;
const unsigned int Dimension=3;

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

if (argc < 3){
    std::cout << "Usage: " << argv[0] << " lung_image enf_mask_file_name [size_th_min] [size_th_max]" << std::endl;
    return EXIT_SUCCESS;
    }
    
    if ( (strcmp(argv[1],"--help") == 0) || ( (strcmp(argv[1],"-h") == 0) ) ){
        std::cout << "Usage: " << argv[0] << " lung_image enf_mask_file_name [size_th_min] [size_th_max]" << std::endl;
        return EXIT_SUCCESS;
        }
        
        const char * lungfilepath  = argv[1];
        const char * emphfilepath  = argv[2];
        const char * lavfile = "lavdata.txt";
        const unsigned int enf_th = 1;
        unsigned int size_th_min = 0;
        unsigned int size_th_max = 0;
        
        if (argc == 4)
            size_th_min = atoi(argv[3]);
            
        if (argc == 5)
            size_th_max = atoi(argv[4]);
        
        LungReaderType::Pointer lung_reader = LungReaderType::New();
        EmphMaskReaderType::Pointer emph_mask_reader = EmphMaskReaderType::New();
        
        lung_reader->SetFileName(lungfilepath);
        emph_mask_reader->SetFileName(emphfilepath);
        
        double spacing[Dimension];
                
        try {
            lung_reader->Update();
        } catch (itk::ExceptionObject & excp) {
            std::cerr << "  Problems reading the lung image"   << std::endl;
            std::cerr << excp << std::endl;
            return EXIT_FAILURE;
        }
            
        try {
            emph_mask_reader->Update();
        } catch (itk::ExceptionObject & excp) {
            std::cerr << "  Problems reading the emphysema mask image"   << std::endl;
            std::cerr << excp << std::endl;
            return EXIT_FAILURE;
        }
            
            
        spacing[0] = lung_reader->GetOutput()->GetSpacing()[0];
        spacing[1] = lung_reader->GetOutput()->GetSpacing()[1];
        spacing[2] = lung_reader->GetOutput()->GetSpacing()[2];
		
		/* 
		 * There is currently a bug in gdcm that always returns 1 as z spacing of multiframe enhanced ct images
		 * Just check the internal dicom tag in this case
		 */
		if (spacing[2] == 1.0) {
			//std::cout << "Reading dicom tag... " <<  std::endl;
			spacing[2]=get_slice_thickness_tag<SSImageType::Pointer>(lung_reader->GetOutput());
		}
        
        //Reading the lungs data
        std::cout << "Spacing is: " <<  spacing[0] << " " <<  spacing[1] << " " << spacing[2] << std::endl;
            
            
        SSConstIteratorType lung_it(lung_reader->GetOutput(), lung_reader->GetOutput()->GetLargestPossibleRegion());
        USConstIteratorType emph_it(emph_mask_reader->GetOutput(), emph_mask_reader->GetOutput()->GetLargestPossibleRegion());
        
        unsigned long vol_count= 0;
        
        for (lung_it.GoToBegin(); !lung_it.IsAtEnd(); ++lung_it) 
            if (lung_it.Get() != 0)
                vol_count++;
                
        unsigned long emph_count = 0;
        
        for (emph_it.GoToBegin(); !emph_it.IsAtEnd(); ++emph_it)
            if (emph_it.Get() == enf_th) 
                emph_count++;
        
            
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
            
            if (size_th_min > 0) {
                if (labelObject->GetNumberOfPixels() <= size_th_min) {
                    continue;
                }
            }
            
            if (size_th_max > 0) {
                if (labelObject->GetNumberOfPixels() > size_th_max) {
                    continue;
                }
            }
//            lav_file_op << labelObject->GetSize() << " ";
            lav_file_op << labelObject->GetNumberOfPixels() * spacing[0] * spacing[1] * spacing[2] << " ";
            emph_cluster_count+=labelObject->GetNumberOfPixels();
           
        }
        lav_file_op.close();
/*        
        if (size_th > 0) {
            ShapeOpeningFilterType::Pointer opening = ShapeOpeningFilterType::New();
            opening->SetInput(shape->GetOutput());
            opening->SetLambda(size_th);
            opening->Update();
                        
            L2IType::Pointer l2i = L2IType::New();
            l2i->SetInput(opening->GetOutput());
            
            USWriterType:: Pointer opening_writer = USWriterType::New();
            opening_writer->SetInput(l2i->GetOutput());
            opening_writer->SetFileName("opening.tiff");
            opening_writer->Update();
        }
  */      
        //std::cout << "Vol counter:  " << vol_count << " Emph count: "<< emph_count << std::endl;
        std::cout << "Emphysema Index:  " << (double) emph_count*100/(double)vol_count << std::endl;
        if ((size_th_min > 0) || (size_th_max > 0)) {
            std::cout << "Emphysema Index (considering part of the clusters):  " << (double) emph_cluster_count*100/(double)vol_count << std::endl;
        }
        
        return EXIT_SUCCESS;
        
        }
                
