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
	// saving vtk with path
	using WriterType = itk::ImageFileWriter<TImage>;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);
	writer->SetFileName(fname);
	writer->Update();
}


int* Coordinates(std::string fname)
{
	// get coordinates from txt file, return table with int values
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
	// get pixel value from image with coordinates x,y,z
	ImageType::IndexType index;
	index[0] = x; index[1] = y; index[2] = z;
	signed short pixel_value = image->GetPixel(index);
	return pixel_value;
}

float StatisticalThreshold(ImageType::Pointer image, int x1, int y1, int z1)
{
	// calculate threshold from 5x5 pixel square surrounding the starting pixel
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
	// morphology: 2 operations dilataion & closing -> better mask
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

void TheshConnected(ImageType::Pointer image, std::string resname, int x1, int y1, int z1, int x2, int y2, int z2)
{
	//Smoothing input image
	using MeanType = itk::MeanImageFilter<ImageType, ImageType>;
	MeanType::Pointer meaner = MeanType::New();
	meaner->SetInput(image);
	meaner->SetRadius(2);
	meaner->Update();

	// find values in start pixels
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

	// find statistical threshold
	float prog = StatisticalThreshold(image, x1, y1, z1);
	short war_prog = round(prog);
	//std::cout << pixel_value << std::endl;
	//std::cout << round(prog) << std::endl;
	short upper, lower;
	// calculate upper and lower thres
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

	// ConnectedThresholdImageFilter with given thres & starting points
	using ConnThrFilterType = //itk::NeighborhoodConnectedImageFilter<ImageType, ImageType>;
		itk::ConnectedThresholdImageFilter<ImageType, ImageType>;
	ConnThrFilterType::Pointer connThr = ConnThrFilterType::New();
	connThr->SetInput(image);
	//connThr->SetConnectivity(ConnThrFilterType::FullConnectivity);
	//connThr->SetRadius({ {25,25,25} });
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

		std::cout << "Wczytane œcie¿ki:" << std::endl;
		std::cout << argv[1] << std::endl;
		std::cout << argv[2] << std::endl;
		std::cout << argv[3] << std::endl;

		// READING FROM INPUT PARAMETERS
		// PATH TO THE COORDINATE FILE
		std::string guzPath = argv[1];
		// PATH TO THE IMAGE FILE
		std::string vtkPath = argv[2];
		// PATH TO THE RESULT IMAGE FILE
		std::string savePath = argv[3];
		// PATH TO THE TEMP IMAGE FILE
		std::string tmpPath = "../wyniki/tmp/tmp.vtk";

		// READ VTK FILE
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(vtkPath);
		reader->Update();

		ImageType::Pointer image = reader->GetOutput();

		// READ TUMOR COORDINATEs
		std::cout << "Odczytane punkty:" << std::endl;
		int* wspolrzedneGuza;
		wspolrzedneGuza = Coordinates(guzPath);


		// FILTERING
		TheshConnected(image, tmpPath, wspolrzedneGuza[0], wspolrzedneGuza[1], wspolrzedneGuza[2],  wspolrzedneGuza[3], wspolrzedneGuza[4], wspolrzedneGuza[5]);
		
		ReaderType::Pointer reader2 = ReaderType::New();
		reader2->SetFileName(tmpPath);
		reader2->Update();
		
		// GET BETTER RESULTS -> MORPHOLOGY ON MASK IMAGE
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
