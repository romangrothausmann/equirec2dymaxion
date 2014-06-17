////program to reproject an equirectangular map to a dymaxion map


#include <itkImageFileReader.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageFileWriter.h>

#include "itkFilterWatcher2.h"

#define V 1

//#define Pi 3.14159265358979323846


/********************************************************/
/* This program is copyrighted by Robert W. Gray and    */
/* not be used in ANY for-profit project without his    */
/* say so.                                              */
/********************************************************/
/* C program to call Dymaxion Map conversion procedures */
/* It is intended to only show you how to call the      */
/* conversion routine so that you may write your own    */
/* driver.                                              */
/********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/********************************************************/
/* Define global variables and procedures               */
/********************************************************/

// double v_x[13], v_y[13], v_z[13];
// double center_x[21], center_y[21], center_z[21];
// double garc, gt, gdve, gel;

extern void init_stuff(void);
extern void convert_s_t_p(double lng, double lat, double *x, double *y);



int main(int argc, char **argv) {

    if ( argc != 4 )
        {
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Output_Image"
		  << " X-res"
                  << std::endl;
        return EXIT_FAILURE;
        }
    
    const char Dimension= 2;

    typedef itk::RGBPixel< unsigned char > InputPixelType;

    typedef itk::Image<InputPixelType, Dimension>  ImageType;

    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
 
    reader->SetFileName(argv[1]);
    FilterWatcher watcherI(reader, "reading");
    try
        { 
        reader->Update();
        }
    catch (itk::ExceptionObject &ex)
        { 
        if (!strcmp(ex.GetDescription(), "Filter does not have progress.")){
            std::cerr << ex << std::endl;
            //std::cerr << ex.GetDescription() << std::endl;
            return EXIT_FAILURE;
            }
        // else
        //   std::cerr << "Reader has no progress!" << std::endl;
        }



    typename ImageType::IndexType start;
    start.Fill(0);
    
    typename ImageType::SizeType size;
    //size.Fill(20);
    unsigned int xres= atoi(argv[3]);

    size[0] = xres;
    size[1] = sqrt(3*3 - 1.5*1.5) * xres / 5.5; // ratios taken from the info of http://www.rwgrayprojects.com/rbfnotes/maps/graymap6.html

    typename ImageType::RegionType region(start, size);

    typename ImageType::RegionType iregion= reader->GetOutput()->GetLargestPossibleRegion();
    int ixres= iregion.GetSize()[0];
    int iyres= iregion.GetSize()[1];

    typename ImageType::PixelType pixelValue;
    // pixelValue.SetRed(0);
    // pixelValue.SetGreen(0);
    // pixelValue.SetBlue(0);
    //pixelValue.Fill(itk::NumericTraits<ImageType::PixelType>::ZeroValue());

    typename ImageType::Pointer image = ImageType::New();
    image->SetRegions(region);
    image->Allocate();
    image->FillBuffer(pixelValue);

    std::cout << "Output canvas created!" << std::endl;
 
    itk::ImageRegionConstIterator<ImageType> cit(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());

    typename ImageType::IndexType pixelIndex; //= {{27,29,37}}
    double lng, lat, ox, oy;
    init_stuff();//very important!!!

    for (cit.GoToBegin(); !cit.IsAtEnd(); ++cit) {

        pixelIndex= cit.GetIndex();

        //printf("ix: %f; iy: %f\n", cit.GetIndex()[0], cit.GetIndex()[1]);

	// lng= 2 * Pi * static_cast<double>(pixelIndex[0]) / ixres - Pi;
	// lat=     Pi * static_cast<double>(pixelIndex[1]) / iyres - Pi/2;
	lng= 2 * 180 * static_cast<double>(pixelIndex[0]) / ixres - 180;
	lat=     180 * static_cast<double>(pixelIndex[1]) / iyres -  90;
	//printf("lng: %f; lat: %f\n", lng, lat);        

        convert_s_t_p(lng, -lat, &ox, &oy); //dymaxion projection

	//printf("ox: %f; oy: %f\n", ox, oy);        

	ox= ox * size[0] / 5.5;
	oy= oy * size[1] / sqrt(3*3 - 1.5*1.5);
	oy= size[1] - oy;
        if((ox < 0) || (ox >= size[0]) || (oy < 0) || (oy >= size[1])){
            printf("Resulting coords out of bonds: x= %f [0; %f], y= %f [0;%f]\n", ox, static_cast<double>(size[0]), oy, static_cast<double>(size[1]));
            }
        else{

            // pixelIndex[0]= ox;
            // pixelIndex[1]= oy;

            pixelValue = cit.Get();

            //pixelValue = image->GetPixel(pixelIndex); //somehow average previous values with new one

	    //printf("RGB: %d,%d,%d\n", static_cast<InputComponentType>(pixelValue.GetRed()), static_cast<InputComponentType>(pixelValue.GetGreen()), static_cast<InputComponentType>(pixelValue.GetBlue()));
	    //std::cout << "RGB: " <<  itk::NumericTraits<typename InputPixelType::ValueType>::PrintType(pixelValue.GetRed()) << "; " << itk::NumericTraits<typename InputPixelType::ValueType>::PrintType(pixelValue.GetGreen()) << "; " << itk::NumericTraits<typename InputPixelType::ValueType>::PrintType(pixelValue.GetBlue()) << std::endl;

	    InputPixelType::ValueType red = pixelValue.GetRed();
	    InputPixelType::ValueType green = pixelValue.GetGreen();
	    InputPixelType::ValueType blue = pixelValue.GetBlue();
	    //std::cout << "RGB: " <<  itk::NumericTraits<typename InputPixelType::ValueType>::PrintType(red) << "; " << itk::NumericTraits<typename InputPixelType::ValueType>::PrintType(green) << "; " << itk::NumericTraits<typename InputPixelType::ValueType>::PrintType(blue) << std::endl;



            pixelIndex[0]= ox;
            pixelIndex[1]= oy;

            image->SetPixel(pixelIndex, pixelValue);

	    //printf("ox: %f; oy: %f\n", ox, oy);        
	    //printf("ix: %f; iy: %f\n", pixelIndex[0], pixelIndex[1]);

            }

        }
    std::cerr << "Projection done." << std::endl;

    typedef itk::ImageFileWriter<ImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer, "writing");
    writer->SetFileName(argv[2]);
    writer->SetInput(image);
    //writer->UseCompressionOn();
    //writer->SetUseCompression(atoi(argv[4]));
    try
        { 
        writer->Update();
        }
    catch (itk::ExceptionObject &ex)
        { 
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    

    return EXIT_SUCCESS;
    }


