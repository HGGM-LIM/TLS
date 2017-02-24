#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkConnectivityImageFilter.h"

typedef unsigned short USPixelType;
const unsigned int Dimension=3;

typedef unsigned long LabelType;

typedef itk::Image<USPixelType, Dimension>                              USImageType;
typedef itk::ImageFileReader< USImageType >                             EmphMaskReaderType;
typedef itk::ImageFileWriter< USImageType >                             EmphMaskWriterType;
typedef itk::ConnectivityImageFilter<USImageType, USImageType>			ConnectivityFilterType;


typedef itk::ImageRegionIterator<USImageType>       ImageIterator;

int main(int argc, char *argv[]) {

    
    if (argc < 3){
        std::cout << "Usage: " << argv[0] << " input_image output_image [connectivity_threshold]" << std::endl;
        return EXIT_SUCCESS;
    }
    
    if ( (strcmp(argv[1],"--help") == 0) || ( (strcmp(argv[1],"-h") == 0) ) ){
        std::cout << "Usage: " << argv[0] << " input_image output_image [connectivity_threshold]" << std::endl;
        return EXIT_SUCCESS;
    }
    
    char * emphfilepath_in=argv[1];
    char * emphfilepath_mask="connectivity_map.tiff";
    char * emphfilepath_out=argv[2];
    int th = 2;
    
    if (argc == 4)
        th = atoi(argv[3]);
    
    std::cout << "The connectivity threshold is " << th << std::endl;
    
    EmphMaskReaderType::Pointer emph_mask_reader = EmphMaskReaderType::New();
    emph_mask_reader->SetFileName(emphfilepath_in);
    
    try {
        emph_mask_reader->Update();
    } catch (itk::ExceptionObject & excp) {
        std::cerr << "  Problems reading the emphysema mask image"   << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
    }
    
    USImageType::Pointer input_image  =  emph_mask_reader->GetOutput();

    ConnectivityFilterType::Pointer conn = ConnectivityFilterType::New();

    conn->SetInput(input_image);
    USImageType::Pointer output_image = conn->GetOutput();
    

    EmphMaskWriterType::Pointer emph_mask_writer = EmphMaskWriterType::New();
    emph_mask_writer->SetFileName(emphfilepath_mask);
    emph_mask_writer->SetInput(output_image);
    
    try {
        emph_mask_writer->Update();
    } catch (itk::ExceptionObject & excp) {
        std::cerr << "  Problems writing the emphysema mask image"   << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
    }
    
    // Threshold the connectivity map
    ImageIterator out(output_image, output_image->GetRequestedRegion());
    
    for (out = out.Begin(); ! out.IsAtEnd(); ++out ) 
    {
        
        out.Set( out.Get()>th );
        
    }
    
    emph_mask_writer->SetFileName(emphfilepath_out);
    emph_mask_writer->SetInput(output_image);
    
    try {
        emph_mask_writer->Update();
    } catch (itk::ExceptionObject & excp) {
        std::cerr << "  Problems writing the emphysema mask image"   << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
    }
    
    
    
    
                                    
}
