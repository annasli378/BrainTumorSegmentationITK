#include <iostream> // std::cout, std::cin
#include <string> // std::string
#include <itkImage.h> // itk::Image
#include <fstream>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkImageSeriesReader.h>
#include <itkImageSeriesWriter.h>
#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>
#include <itkNumericSeriesFileNames.h>

#include<itkRescaleIntensityImageFilter.h>
#include<itkThresholdImageFilter.h>
#include<itkBinaryThresholdImageFilter.h>
#include<itkConnectedThresholdImageFilter.h >
#include<itkRegionGrowImageFilter.h>

#include<itkBinaryBallStructuringElement.h>
#include<itkBinaryErodeImageFilter.h>
#include<itkBinaryDilateImageFilter.h>

#include<itkMeanImageFilter.h>
#include<itkBinaryMorphologicalClosingImageFilter.h>
#include<itkConnectedComponentImageFilter.h>
#include<itkScalarConnectedComponentImageFilter.h>
#include<itkNeighborhoodConnectedImageFilter.h>
#include<itkConfidenceConnectedImageFilter.h>
#include<itkWatershedImageFilter.h>
#include<itkGradientAnisotropicDiffusionImageFilter.h>
#include<itkGradientMagnitudeImageFilter.h>
#include<itkEllipseSpatialObject.h>
#include<itkSpatialObjectToImageStatisticsCalculator.h>
#include<itkExtractImageFilter.h>

using ImageType = itk::Image<short, 3>; // standard DICOM
using ImageType2D = itk::Image<short, 2>; 
using ImageTypeF = itk::Image<float, 3>;
using ReaderType = itk::ImageFileReader<ImageType>;
using WriterType = itk::ImageFileWriter<ImageType>;

using SeriesReaderType = itk::ImageSeriesReader<ImageType>;
using SeriesWriterType = itk::ImageSeriesWriter<ImageType, ImageType2D>;

template<class TImage>
void SaveImage(typename TImage::Pointer image, std::string fname)
{
	using WriterType = itk::ImageFileWriter<TImage>;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);
	writer->SetFileName(fname);
	writer->Update();
}

void SaveToVTK()
{
	// TODO: Lab.2
		// ODCZYT I ZAPIS OBRAZOW:
		//ReaderType::Pointer reader = ReaderType::New();
		//reader->SetFileName("../dane/p1/IM000070.dcm"); //!
		//reader->Update();
		//ImageType::Pointer image = reader->GetOutput();

		//WriterType::Pointer writer = WriterType::New();
		//writer->SetFileName("../wyniki/kopia2d.dcm");
		//writer->SetInput(reader->GetOutput());
		//writer->Update();

		//METADANE OBRAZU
		//ImageType::Pointer image = reader->GetOutput();
		//ImageType::DirectionType dirs = image->GetDirection();
		//ImageType::PointType orgin = image->GetOrigin();
		//ImageType::SpacingType spacing = image->GetSpacing();
		//ImageType::RegionType region = image->GetLargestPossibleRegion();

		//std::cout << dirs << std::endl;
		//std::cout << orgin << std::endl;
		//std::cout << spacing << std::endl;
		//std::cout << region << std::endl;

	itk::GDCMSeriesFileNames::Pointer gdcmNames
		= itk::GDCMSeriesFileNames::New();
	gdcmNames->SetDirectory("../dane/p4/MR-4");
	itk::SerieUIDContainer series = gdcmNames->GetSeriesUIDs();

	for (int i = 0; i < series.size(); i++)
	{
		std::cout << series[i] << std::endl;
	}

	itk::FilenamesContainer fileNames = gdcmNames->GetFileNames(series[0]);

	SeriesReaderType::Pointer seriesReader = SeriesReaderType::New();
	seriesReader->SetFileNames(fileNames);
	seriesReader->Update();

	ImageType::Pointer image3D = seriesReader->GetOutput();

	std::cout << image3D << std::endl;

	WriterType::Pointer writer2 = WriterType::New();
	writer2->SetFileName("../wyniki/kostki/p4/MR4.vtk");
	writer2->SetInput(image3D);
	writer2->Update();

}

void RescaleInten(ImageType::Pointer image, std::string resname)
{
	using IntensityFilter = itk::ThresholdImageFilter<ImageType>;
	IntensityFilter::Pointer threshold = IntensityFilter::New();
	threshold->SetInput(image);
	threshold->SetLower(-1024);

	using RescaleFilter = itk::RescaleIntensityImageFilter<ImageType>;
	RescaleFilter::Pointer filter = RescaleFilter::New();
	filter->SetInput(threshold->GetOutput());
	filter->SetOutputMinimum(0);
	filter->SetOutputMaximum(255);
	filter->Update();

	SaveImage<ImageType>(filter->GetOutput(), resname); 

}


int* Coordinates(std::string fname)
{
	static int coordinates[6];
	std::fstream inputFile;
	inputFile.open(fname);

	if (inputFile.is_open())
	{

		int x_1, y_1, z_1, x_2, y_2, z_2;
		inputFile >> x_1;
		inputFile >> y_1;
		inputFile >> z_1;
		inputFile >> x_2;
		inputFile >> y_2;
		inputFile >> z_2;

		std::cout << "  " << x_1 << "  " << y_1 << "  " << z_1 << "  ";
		std::cout << "  " << x_2 << "  " << y_2 << "  " << z_2 << "  " << std::endl;

		coordinates[0] = x_1;
		coordinates[1] = y_1;
		coordinates[2] = z_1;
		coordinates[3] = x_2;
		coordinates[4] = y_2;
		coordinates[5] = z_2;

		inputFile.close();   //close the file object.
		return coordinates;
	}
	else
	{
		std::cout << "The text file could not be opened!!" << std::endl;
	}
}

short PixValue(ImageType::Pointer image, int x, int y, int z)
{
	ImageType::IndexType index;
	index[0] = x; index[1] = y; index[2] = z;
	signed short pixel_value = image->GetPixel(index);
	return pixel_value;
}

float StatisticalThreshold(ImageType::Pointer image, int x1, int y1, int z1)
{
	static short tablica[5][5];
	tablica[0][0] = PixValue(image, x1 - 2, y1 - 2, z1);
	tablica[0][1] = PixValue(image, x1 - 2, y1 - 2, z1);
	tablica[0][2] = PixValue(image, x1 - 2, y1 - 2, z1);
	tablica[0][3] = PixValue(image, x1 - 2, y1 - 2, z1);
	tablica[0][4] = PixValue(image, x1 - 2, y1 - 2, z1);

	tablica[1][0] = PixValue(image, x1 - 1, y1 - 1, z1);
	tablica[1][1] = PixValue(image, x1 - 1, y1 - 1, z1);
	tablica[1][2] = PixValue(image, x1 - 1, y1 - 1, z1);
	tablica[1][3] = PixValue(image, x1 - 1, y1 - 1, z1);
	tablica[1][4] = PixValue(image, x1 - 1, y1 - 1, z1);

	tablica[2][0] = PixValue(image, x1 , y1, z1);
	tablica[2][1] = PixValue(image, x1 , y1 , z1);
	tablica[2][2] = PixValue(image, x1 , y1 , z1);
	tablica[2][3] = PixValue(image, x1 , y1 , z1);
	tablica[2][4] = PixValue(image, x1 , y1 , z1);

	tablica[3][0] = PixValue(image, x1 + 1, y1 + 1, z1);
	tablica[3][1] = PixValue(image, x1 + 1, y1 + 1, z1);
	tablica[3][2] = PixValue(image, x1 + 1, y1 + 1, z1);
	tablica[3][3] = PixValue(image, x1 + 1, y1 + 1, z1);
	tablica[3][4] = PixValue(image, x1 + 1, y1 + 1, z1);

	tablica[4][0] = PixValue(image, x1 + 2, y1 + 2, z1);
	tablica[4][1] = PixValue(image, x1 + 2, y1 + 2, z1);
	tablica[4][2] = PixValue(image, x1 + 2, y1 + 2, z1);
	tablica[4][3] = PixValue(image, x1 + 2, y1 + 2, z1);
	tablica[4][4] = PixValue(image, x1 + 2, y1 + 2, z1);
	
	float sum = 0.0, mean, std = 0.0;

	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			sum += tablica[i][j];
		}
	}
	mean = sum / 25;
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			std += pow(tablica[i][j]-mean,2);
		}
	}


	return sqrt(std/25);
}

void MorphOperations(ImageType::Pointer image, std::string resname)
{
	const unsigned int Dimension = 3;
	using BallStrelType = itk::BinaryBallStructuringElement<ImageType::PixelType, Dimension>;
	BallStrelType ball;
	ball.SetRadius(2);
	ball.CreateStructuringElement();

	//std::cout << ball.GetSize() << std::endl;

	using DilType = itk::BinaryDilateImageFilter<ImageType, ImageType, BallStrelType>;
	DilType::Pointer dyl = DilType::New();
	dyl->SetInput(image);
	dyl->SetKernel(ball);
	dyl->SetBackgroundValue(0);
	dyl->SetForegroundValue(250);
	dyl->Update();

	using ClosingType = itk::BinaryMorphologicalClosingImageFilter<ImageType, ImageType, BallStrelType>;
	ClosingType::Pointer close = ClosingType::New();
	close->SetInput(dyl->GetOutput());
	close->SetKernel(ball);
	close->SetForegroundValue(250);
	close->Update();

	SaveImage<ImageType>(close->GetOutput(), resname);
}

void FiltracjaConnected(ImageType::Pointer image, std::string resname, int x1, int y1, int z1, int x2, int y2, int z2)
{
	//Wyg³adzenie obrazu we
	using MeanType = itk::MeanImageFilter<ImageType, ImageType>;
	MeanType::Pointer meaner = MeanType::New();
	meaner->SetInput(image);
	meaner->SetRadius(1);
	meaner->Update();

	// znalezienie wartoœci w punktach startoych
	ImageType::IndexType index;
	index[0] = x1;
	index[1] = y1;
	index[2] = z1;
	signed short pixel_value = image->GetPixel(index);
	

	/*
	using FilterType = itk::ExtractImageFilter<ImageType, ImageType2D>;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput(meaner->GetOutput());
	filter->SetDirectionCollapseToSubmatrix();

	
	ImageType::RegionType inputRegion = meaner->GetOutput()->GetBufferedRegion();
	ImageType::SizeType size = inputRegion.GetSize();
	size[2] = 1;
	ImageType::IndexType start = inputRegion.GetIndex();
	int sliceNumber = z1;
	start[2] = sliceNumber;

	ImageType::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	filter->SetExtractionRegion(desiredRegion);
	filter->Update();

	std::cout << filter->GetOutput()->GetImageDimension() << std::endl;

	//Obliczenie progów -- statystyczne
	/*using EllipseType = itk::EllipseSpatialObject<2>;
	EllipseType::Pointer ellipse = EllipseType::New();
	ellipse->SetRadiusInObjectSpace(3); // promieñ ko³a
	EllipseType::PointType offset;
	offset[0] = x1;
	offset[1] = y1;
	// tu trzeba umeijscowiæ obszar na obazku
	ellipse->SetCenterInObjectSpace(offset);
	ellipse->Update();
	std::cout << ellipse->GetCenterInObjectSpace() << std::endl;

	//Wyciêcie 1 slajdu do policzenia statystyk z obrazu 3D po wyg³adzeniu

	std::cout << filter->GetOutput() << std::endl;
	
	std::cout << (filter->GetOutput())->GetImageDimension() << std::endl;
	std::cout << ellipse->GetLargestPossibleRegion() << std::endl;

	
	//policzenie statystyki
	using CalculatorType = itk::SpatialObjectToImageStatisticsCalculator<ImageType2D, EllipseType>;
	CalculatorType::Pointer calculator = CalculatorType::New();
	calculator->SetImage(filter->GetOutput());
	calculator->SetSpatialObject(ellipse);
	calculator->Update();
	std::cout << "Sample mean = " << calculator->GetMean() << std::endl;
	std::cout << "Sample covariance = " << calculator->GetCovarianceMatrix();
	*/

	// próg znaleziony statystycznie
	float prog = StatisticalThreshold(image, x1, y1, z1);
	short war_prog = round(prog);
	//std::cout << pixel_value << std::endl;
	//std::cout << round(prog) << std::endl;
	short upper, lower;
	// obliczenie progu
	short prog_l = war_prog; //5;
	short prog_u = war_prog;//6;

	if (pixel_value > prog_l)
	{
		 lower = pixel_value - prog_l;
		 upper = pixel_value + prog_u;
	}
	else
	{
		lower = 10; //pixel_value - 10
		upper = 100;//pixel_value + 10
	}


	std::cout << pixel_value << std::endl;
	std::cout << lower << std::endl;
	std::cout << upper << std::endl;

	using ConnThrFilterType = itk::ConnectedThresholdImageFilter<ImageType, ImageType>;
	ConnThrFilterType::Pointer connThr = ConnThrFilterType::New();
	connThr->SetInput(image);
	connThr->SetConnectivity(ConnThrFilterType::FullConnectivity);

	connThr->SetLower(lower);
	connThr->SetUpper(upper);
	connThr->SetReplaceValue(250);

	ConnThrFilterType::IndexType seed1, seed2;
	seed1[0] = x1; seed1[1] = y1; seed1[2] = z1;
	connThr->AddSeed(seed1);
	seed2[0] = x2; seed2[1] = y2; seed2[2] = z2;
	connThr->AddSeed(seed2);
	connThr->Update();

	SaveImage<ImageType>(connThr->GetOutput(), resname);
}

void FiltracjaScalar(ImageType::Pointer image, std::string resname)
{
	using ScalarFilterType = itk::ScalarConnectedComponentImageFilter<ImageType, ImageType>;
	ScalarFilterType::Pointer scalar = ScalarFilterType::New();
	scalar->SetInput(image);
	scalar->SetDistanceThreshold(5);
	scalar->SetBackgroundValue(0);
	scalar->Update();
	SaveImage<ImageType>(scalar->GetOutput(), resname);
}

void FiltracjaConfidence(ImageType::Pointer image, std::string resname, int x1, int y1, int x2, int y2)
{
	using ConfFilterType = itk::ConfidenceConnectedImageFilter<ImageType, ImageType>;
	ConfFilterType::Pointer conf = ConfFilterType::New();
	conf->SetInput(image);
	conf->SetInitialNeighborhoodRadius(50);
	conf->SetMultiplier(5);

	ConfFilterType::IndexType seed1, seed2;
	seed1[0] = x1; seed1[1] = y1;
	conf->AddSeed(seed1);
	seed2[0] = x2; seed2[1] = y2;
	conf->AddSeed(seed2);

	conf->Update();
	SaveImage<ImageType>(conf->GetOutput(), resname);
}

void Watershed(ImageType::Pointer image, std::string resname)
{

	std::cout << "Watershed:" << std::endl;
	//using DiffusionFilterType = itk::GradientAnisotropicDiffusionImageFilter<ImageType, ImageType> ;
	using GradientMagnitudeFilterType = itk::GradientMagnitudeImageFilter<ImageType, ImageType>;
	GradientMagnitudeFilterType::Pointer gradientMagnitudeFilter = GradientMagnitudeFilterType::New();
	gradientMagnitudeFilter->SetInput(image);
	gradientMagnitudeFilter->Update();

	using WatershedFilterType=itk::WatershedImageFilter<ImageType> ;
	WatershedFilterType::Pointer watershed = WatershedFilterType::New();

	//DiffusionFilterType::Pointer diffusion = DiffusionFilterType::New();
	//diffusion->SetNumberOfIterations(20);
	//diffusion->SetConductanceParameter(3.0);
	//diffusion->SetTimeStep(0.0625);

	float threshold = 0.5;
	float level = 0.1;

	watershed->SetThreshold(threshold);
	watershed->SetLevel(level);
	watershed->SetInput(gradientMagnitudeFilter->GetOutput());
	watershed->Update();



}

void RozdzielenieEtykiet(ImageType::Pointer image, std::string resname)
{
	using ConnComFilterType = itk::ConnectedComponentImageFilter<ImageType, ImageType, ImageType>;
	ConnComFilterType::Pointer connCom = ConnComFilterType::New();
	connCom->SetInput(image);
	connCom->SetBackgroundValue(0);
	connCom->SetMaskImage(image);
	connCom->Update();

	SaveImage<ImageType>(connCom->GetOutput(), resname);



}

//int main(char *argt[],char *argv[]) 
int main(int argc, char *argv[]) // glowna funkcja programu
{
	try {
		if (argc < 2)
		{
			std::cout << "Input arguments" << std::endl;
			std::cout << "Usage: " << argv[0] << "filenamePath" << std::endl;
			return EXIT_FAILURE;
		}

		
		std::cout << argv[1] << std::endl;
		std::cout << argv[2] << std::endl;
		std::cout << argv[3] << std::endl;

		// ODCZYTANIE Z PARAMETRÓ WEJŒCIOWYCH
		// ŒCIE¯KI DO PLIKU ZE WSPÓ£RZÊDNYMI
		std::string guzPath = argv[1];
		//"../wyniki/rescaled/p1/wspolrzedneMR4.txt";
		// ŒCIE¯KI DO PRZESKALOWANEGO OBRAZU
		std::string vtkPath = argv[2];
							  //"../wyniki/rescaled/p1/MR4.vtk";

		std::string tmpPath = "../wyniki/tmp/tmp.vtk";

		std::string savePath = argv[3];
		// ODCZYT PLIKU VTK
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(vtkPath);
		reader->Update();

		ImageType::Pointer image = reader->GetOutput();

		// ODCZYT WSPÓ£RZÊDNYCH GUZA
		int* wspolrzedneGuza;
		wspolrzedneGuza = Coordinates(guzPath);

		/*for (int i = 0; i < 6; i++)
		{
		std::cout << wspolrzedneGuza[i] << std::endl;
		}*/

		FiltracjaConnected(image, tmpPath, wspolrzedneGuza[0], wspolrzedneGuza[1], wspolrzedneGuza[2],  wspolrzedneGuza[3], wspolrzedneGuza[4], wspolrzedneGuza[5]);
		
		ReaderType::Pointer reader2 = ReaderType::New();
		reader2->SetFileName(tmpPath);
		reader2->Update();
		
		MorphOperations(reader2->GetOutput(), savePath);



	}
	catch (itk::ExceptionObject& ex) {
		ex.Print(std::cout);
	}
	catch (std::exception& ex) {
		std::cout << ex.what() << std::endl;
	}
	catch (...) {
		std::cout << "Unknown error!" << std::endl;
	}
	std::cout << "Hit [Enter]...";
	std::cin.get();
	return EXIT_SUCCESS; // albo EXIT_FAILURE
}
