#ifndef __itkEmphysemaCalculationsImageFilter_txx
#define __itkEmphysemaCalculationsImageFilter_txx

#include "vnl/vnl_vector_fixed.h"
#include "itkProgressReporter.h"
#include "itkSize.h"
#include "itkImageRegion.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionReverseIterator.h"
#include "itkEmphysemaCalculationsImageFilter.h"




namespace itk
{
template< class TInputImage, class TOutputImage >
EmphysemaCalculationsImageFilter< TInputImage, TOutputImage >
::EmphysemaCalculationsImageFilter()
{
  itkDebugMacro(<< "EmphysemaCalculationsImageFilter::EmphysemaCalculationsImageFilter() called");

	m_LowEmphValue = -950;
  m_Backup = 0;
  m_AutomaticallyResize = 1;
  
  m_tEmphysema = 0;
  m_tConverter = 0;
  m_tShape     = 0;
  m_tLavFile   = 0;
  m_tHistogram = 0;
  
  m_NumberOfBins = 255;
  m_MarginScale  = 10;
  m_VolLung      = 0;
  m_VolEmph      = 0; 
  m_Mean = 0;
  m_Smooth = 0;
  m_SmoothRadius=1;
  
  m_Size.Fill(0);
  m_Index.Fill(0);
  m_ClusterFile = "lavdata.txt";
}

template< class TInputImage, class TOutputImage >
void
EmphysemaCalculationsImageFilter< TInputImage, TOutputImage >
::GenerateData()
{

  //this->AllocateOutputs();

  struct timeb initialTime;
  struct timeb actualTime;
  struct timeb time0;
	ftime(&initialTime);
	ftime(&time0);
	
  std::string histfile = "hist.csv";
  double min = -1024;
  double max = 0;
  	
  // Get the input and output pointers
  InputImageConstPointer  inputPtr  = this->GetInput(0);
  //OutputImagePointer      enfisemaImage = this->GetOutput(0);

  double spacing[3];
  spacing[0] = inputPtr->GetSpacing()[0];
  spacing[1] = inputPtr->GetSpacing()[1];
  spacing[2] = inputPtr->GetSpacing()[2];

  //outputPtr->SetSpacing(inputPtr->GetSpacing());

  // How big is the input image?
  RegionType inputRegion = inputPtr->GetRequestedRegion();
  BinaryImageType::Pointer enfisemaImage = BinaryImageType::New();
  enfisemaImage->SetRegions(inputRegion);
  enfisemaImage->Allocate();
  //SizeType size = inputRegion.GetSize();
  //IndexType startIndex = inputRegion.GetIndex();
  
  ConstIteratorType in(inputPtr, inputRegion);
  
  //We now write the emphysema image and calculate basic image properties
  BinaryIteratorType emphIt(enfisemaImage, enfisemaImage->GetRequestedRegion());
  m_VolLung = 0;
  m_Mean = 0;
  double mld_ei = 0;
  for (in.GoToBegin(),emphIt.GoToBegin();!in.IsAtEnd(),!emphIt.IsAtEnd();++in,++emphIt){
	  PixelType p = in.Get();
	  if (p != 0) {
		  m_VolLung++;
		  m_Mean = m_Mean + p;
		  if (p > max)
			  max = p;
	  }
	  if (p <= m_LowEmphValue){
		  emphIt.Set(255);
		  mld_ei += p;
	  }
  }

  StructuringElementType structuringElement;
  structuringElement.SetRadius(m_SmoothRadius);
  structuringElement.CreateStructuringElement();

  BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();
  BinaryImageType::Pointer imagePtr = enfisemaImage;
  if (m_Smooth) {

	  erodeFilter->SetInput(enfisemaImage);
	  erodeFilter->SetKernel(structuringElement);
	  erodeFilter->SetErodeValue(255);
	  erodeFilter->Update();

	  imagePtr = erodeFilter->GetOutput();

    }

  this->GraftOutput(imagePtr);

  BinaryIteratorType emphIt2(imagePtr, imagePtr->GetRequestedRegion());
  // Now we will quantify the emphysema
  m_VolEmph = 0;
  for (emphIt2.GoToBegin();!emphIt2.IsAtEnd();++emphIt2)
	if (emphIt2.Get() > 0)
		m_VolEmph++;

  m_Mean = (long int)((m_Mean*1.0)/m_VolLung);
  mld_ei = (long int)((mld_ei*1.0)/m_VolEmph);
  ftime(&actualTime);
	m_tEmphysema = 1000.0 * (actualTime.time - initialTime.time) + actualTime.millitm - initialTime.millitm;
	ftime(&initialTime);
	
	max = -10.0;	
	std::cout << "Setting hist max: " << max << " to ensure that 0 valued pixels are not included in HIST calculation" << std::endl;

	std::cout << "Image parameters:" << "Min: " << min << ", Max: " << max << " en " << m_tEmphysema <<"mseg.\n";




	ConverterType::Pointer converter =  ConverterType::New();
	converter->SetInput(imagePtr);
  converter->SetInputForegroundValue(255);
  //We don't want clusters that are separeted only by a few pixels to be considered one
  converter->SetFullyConnected(0);
  ftime(&actualTime);
	m_tConverter = 1000.0 * (actualTime.time - initialTime.time) + actualTime.millitm - initialTime.millitm;
	ftime(&initialTime);
	
	ShapeFilterType::Pointer shape = 																	ShapeFilterType::New();
  shape->SetInput( converter->GetOutput() );
  shape->Update();
  LabelCollectionType::Pointer collection = shape->GetOutput();
  ftime(&actualTime);
	m_tShape = 1000.0 * (actualTime.time - initialTime.time) + actualTime.millitm - initialTime.millitm;
	ftime(&initialTime);
  
   //Sizes in mmc for the BI calculation
   float c1 = 15.625;
   float c2 = 125.0;
   float c3 = 357.911;
   float a1 = 0.0; float a2=0.0; float a3=0.0; float a4=0.0;
   float g1 = 4.0; float g2=4.0; float g3=4.0; float g4=4.0;
  
  std::cout << "Now saving lav sizes on file " << m_ClusterFile << "..."<< std::endl;
  
  std::ofstream lav_file_op(m_ClusterFile.c_str());

  for(unsigned int label=1; label<collection->GetNumberOfLabelObjects(); label++ )
  {
   	LabelObjectType::Pointer labelObject = collection->GetLabelObject( label );
   	lav_file_op << labelObject->GetNumberOfPixels() << " ";
   	
   	float label_size = labelObject->GetNumberOfPixels();
   	float pARc = label_size*1.0/m_VolLung*100;
   	
   	if (label_size <= c1){
            // Class 1 bullae
            a1 += pARc;
    } else if ((label_size > c1) and (label_size <= c2)) {
            // Class 2 bullae
            a2 += pARc;
    } else if ((label_size > c2) and (label_size <= c3)) {
           // Class 3 bullae
            a3 += pARc;
    } else if (label_size > c3) {
           //  Class 4 bullae
            a4 += pARc;
    }
  	
  }
  

  if ((a1 >= 0) and (a1 <=4))
    g1 = a1;
  if ((a2 >= 0) and (a2 <=4))
    g2 = a2;
  if ((a3 >= 0) and (a3 <=4))
    g3 = a3;
  if ((a4 >= 0) and (a4 <=4)) 
    g4 = a4;
    
//  std::cout << "a1: " << a1 << " a2: " << a2 << " a3: " << a3 << " a4: " << a4 << std::endl;
//  std::cout << "g1: " << g1 << " g2: " << g2 << " g3: " << g3 << " g4: " << g4 << std::endl;
//  std::cout << "EI% = " << a1+a2+a3+a4 << std::endl;
  float bindex=(g2 + 2*g3 + 3*g4)/2.4;

   
  lav_file_op.close();
  ftime(&actualTime);
	m_tLavFile = 1000.0 * (actualTime.time - initialTime.time) + actualTime.millitm - initialTime.millitm;
	ftime(&initialTime);
	
  typedef itk::Statistics::ScalarImageToHistogramGenerator<ImageType > HistogramGeneratorType;
	typedef HistogramGeneratorType::HistogramType HistogramType;
	HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
	histogramGenerator->SetInput(inputPtr);
	histogramGenerator->SetNumberOfBins(m_NumberOfBins);
	histogramGenerator->SetMarginalScale(m_MarginScale);
	histogramGenerator->SetHistogramMin(min);
	histogramGenerator->SetHistogramMax(max);
	histogramGenerator->Compute();

	const HistogramType * histogram = histogramGenerator->GetOutput();

	unsigned int binNumber = 0;
	unsigned long freq=0;
	unsigned long cumfreq=0;
	
	HistogramType::ConstIterator itr = histogram->Begin();
	HistogramType::ConstIterator end = histogram->End();
	
	while( itr != end )
	{
		freq = (unsigned long)itr.GetFrequency();
		cumfreq = cumfreq + freq;
		++itr;
	}

	unsigned long totfreq=cumfreq;
	long bin = 0;
	double percent = 0;
	long HIST = 0;
	freq=0;
	cumfreq=0;

	std::cout << "Now saving hist data on file " << histfile << "..." << std::endl;
	std::ofstream hist_file_op(histfile.c_str());
	hist_file_op << "bin,val,frequency,cumfreq,percentile"<< std::endl;

	itr = histogram->Begin();
	end = histogram->End();
	
	while( itr != end )
	{
		freq = (unsigned long)itr.GetFrequency();
		cumfreq = cumfreq + freq;
		percent = (cumfreq * 100.0) /totfreq;
		bin = (long)(min + binNumber*((max-min)/m_NumberOfBins));
		hist_file_op << binNumber << "," << bin << "," << freq << "," <<  cumfreq << "," << percent << std::endl;
		++itr;
		++binNumber;
		
		if ((percent > 15.0) && (HIST == 0.0))
			HIST = bin; 
	}

  hist_file_op.close();

	std::cout << "Volume: " << m_VolLung/1000*spacing[0]*spacing[1]*spacing[2] << std::endl;

  std::cout << "Emphysema Index (" << m_LowEmphValue  << "HU): " << (double) m_VolEmph*100/(double)m_VolLung << std::endl;
  std::cout << "MLD: " << m_Mean << std::endl;
  std::cout << "HIST(15): " << HIST << std::endl;
  std::cout << "MLD (" << m_LowEmphValue  << "HU): " << mld_ei << std::endl;
  std::cout << "BI: " << bindex << std::endl;
  ftime(&actualTime);
	m_tHistogram = 1000.0 * (actualTime.time - initialTime.time) + actualTime.millitm - initialTime.millitm;
  std::cout << "Tiempo total: " << (float)(actualTime.time-time0.time+(actualTime.millitm-time0.millitm)/1000.0) << "s." << std::endl;
}

template< class TInputImage, class TOutputImage >
void
EmphysemaCalculationsImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
    
}

} // end namespace

#endif
