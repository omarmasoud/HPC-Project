#include <iostream>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include<mpi.h>
#include <ctime>// include this header 
#pragma once

#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>
using namespace std;
using namespace msclr::interop;
#define TagStartindex  1
#define  TagProcessorRows  2
#define TagIsLast  3
#define TagKernelSize  4
#define  TagLocalImage 5
#define  TaglocalimageSize  6

int* inputImage(int* w, int* h, System::String^ imagePath) //put the size of image in w & h
{
	int* input;


	int OriginalImageWidth, OriginalImageHeight;

	//*********************************************************Read Image and save it to local arrayss*************************	
	//Read Image and save it to local arrayss

	System::Drawing::Bitmap BM(imagePath);

	OriginalImageWidth = BM.Width;
	OriginalImageHeight = BM.Height;
	*w = BM.Width;
	*h = BM.Height;
	int* Red = new int[BM.Height * BM.Width];
	int* Green = new int[BM.Height * BM.Width];
	int* Blue = new int[BM.Height * BM.Width];
	input = new int[BM.Height * BM.Width];
	for (int i = 0; i < BM.Height; i++)
	{
		for (int j = 0; j < BM.Width; j++)
		{
			System::Drawing::Color c = BM.GetPixel(j, i);

			Red[i * BM.Width + j] = c.R;
			Blue[i * BM.Width + j] = c.B;
			Green[i * BM.Width + j] = c.G;

			input[i * BM.Width + j] = ((c.R + c.B + c.G) / 3); //gray scale value equals the average of RGB values

		}

	}
	return input;
}


void createImage(int* image, int width, int height, int index)
{
	System::Drawing::Bitmap MyNewImage(width, height);


	for (int i = 0; i < MyNewImage.Height; i++)
	{
		for (int j = 0; j < MyNewImage.Width; j++)
		{
			//i * OriginalImageWidth + j
			if (image[i * width + j] < 0)
			{
				image[i * width + j] = 0;
			}
			if (image[i * width + j] > 255)
			{
				image[i * width + j] = 255;
			}
			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}
	MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");
	std::cout << "result Image Saved " << index << endl;
}
void computePadding(int* image, int height, int width, int kernelSize, int& paddedHeight, int& paddedWidth, int*& paddedImage)
{
	int padding = (kernelSize - 1) / 2;  // Compute the padding size

	// Compute the padded dimensions
	paddedHeight = height + 2 * padding;
	paddedWidth = width + 2 * padding;

	// Allocate memory for the padded image array
	paddedImage = new int[paddedHeight * paddedWidth];

	// Initialize the padded image array with zeros


	for (int i = 0; i < paddedHeight; i++) {
		for (int j = 0; j < paddedWidth; j++) {

			if (i < padding || j < padding || j > width - 1 + padding || i > height - 1 + padding)
			{

				paddedImage[i * (paddedWidth)+j] = 0;
			}
			else
			{

				paddedImage[i * (paddedWidth)+j] = image[(i - padding) * (width)+(j - padding)];
			}
		}
	}


}

int* createHighPassKernel(int size)
{
	// Check if the size is odd
	if (size % 2 == 0)
		size++; // Increment size to make it odd

	// Calculate the center index of the kernel
	int center = size / 2;

	// Calculate the total number of elements in the kernel
	int kernelSize = size * size;

	// Create a 1D kernel array
	int* kernel = new int[kernelSize];

	// Set the elements of the kernel
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i == center && j == center)
				kernel[i * size + j] = kernelSize - 1;
			else
				kernel[i * size + j] = -1;
		}
	}

	return kernel;
}
int* applyKernel(int* kernel, int kernelSize, const int* paddedImage, int startIndex, int rowsPerProcessor, int paddedWidth, int& subsize)
{
	int subImageWidth = paddedWidth - (kernelSize - 1);
	int subImageHeight = rowsPerProcessor - (kernelSize - 1);
	printf("conv width output is %d\nconv height output is %d\n", subImageWidth, subImageHeight);
	subsize = subImageWidth * subImageHeight;
	int* subImage = new int[subsize];
	printf("start index is %d rowsperprocess is %d\n", startIndex, rowsPerProcessor);

	// Iterate over each row of the sub-image
	for (int i = 0; i < subImageHeight; i++)
	{
		// Calculate the current row index in the padded image
		int paddedRowIndex = startIndex + i;

		// Iterate over each column of the sub-image
		for (int j = 0; j < subImageWidth; j++)
		{
			// Calculate the current column index in the padded image
			int paddedColIndex = startIndex + j;

			// Apply the kernel to the corresponding pixel in the sub-image
			int sum = 0;

			// Iterate over each row of the kernel
			for (int k = 0; k < kernelSize; k++)
			{
				// Calculate the current row index in the kernel
				int kernelRowIndex = k;

				// Iterate over each column of the kernel
				for (int l = 0; l < kernelSize; l++)
				{
					// Calculate the current column index in the kernel
					int kernelColIndex = l;

					// Calculate the corresponding pixel index in the padded image
					int paddedPixelIndex = (paddedRowIndex + kernelRowIndex) * paddedWidth + (paddedColIndex + kernelColIndex);

					// Multiply the kernel value with the corresponding pixel value in the padded image
					int pixelValue = paddedImage[paddedPixelIndex];
					int kernelValue = kernel[kernelRowIndex * kernelSize + kernelColIndex];
					sum += pixelValue * kernelValue;
				}
			}

			// Store the result of applying the kernel to the corresponding pixel in the sub-image
			int subImagePixelIndex = i * subImageWidth + j;
			subImage[subImagePixelIndex] = sum;
		}
	}

	return subImage;
}



int main()
{
	int ImageWidth = 4, ImageHeight = 4;

	int start_s, stop_s, TotalTime = 0;
	int rank, size;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	System::String^ imagePath;
	std::string img;
	img = "..//Data//Input//test.png";

	imagePath = marshal_as<System::String^>(img);
	int* imageData = inputImage(&ImageWidth, &ImageHeight, imagePath);

	start_s = clock();
	/*

	*/
	stop_s = clock();
	TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
	int kernelsize = 0;
	int* kernel = NULL;
	int paddedheight;
	int paddedwidth;
	int* paddedimage = NULL;
	int startindex=0;
	int rowssent;
	int islast;
	int localimagesize;
	int rowsperprocessor;
	int subsize;
	MPI_Status ReveiveStatus;
	
	if (rank == 0) {
		kernelsize = 3;
		kernel = createHighPassKernel(kernelsize);
		
		printf("process zero printing unpadded height as %d", ImageHeight);
		computePadding(imageData, ImageHeight, ImageWidth, kernelsize, paddedheight, paddedwidth, paddedimage);
		MPI_Bcast(&kernelsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(kernel, kernelsize * kernelsize, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&paddedheight, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&paddedwidth, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(paddedimage, paddedheight * paddedwidth, MPI_INT, 0, MPI_COMM_WORLD);
		printf("process zero printing padded height as %d", paddedheight);
		
		rowsperprocessor = paddedheight / (size);
		printf("\nlsuasdiuasduiauiscbidasuduias %d ",(sizeof(paddedimage)/ sizeof(paddedimage[0])));
		for (int i = 1; i < size; i++) {
			startindex = rowsperprocessor * paddedwidth * (i-1);
			if (i == size - 1) {
				rowsperprocessor = paddedheight - (rowsperprocessor * (size - 2));
				islast = 1;
			}
			else {
				rowsperprocessor = paddedheight / size;
				islast = 0;
			}
			localimagesize = rowsperprocessor * paddedwidth;
			if (i != 0) {

				MPI_Send(&startindex, 1, MPI_INT, i, TagStartindex, MPI_COMM_WORLD);
				MPI_Send(&rowsperprocessor, 1, MPI_INT, i, TagProcessorRows, MPI_COMM_WORLD);
				MPI_Send(&islast, 1, MPI_INT, i, TagIsLast, MPI_COMM_WORLD);
				MPI_Send(&localimagesize, 1, MPI_INT, i, TaglocalimageSize, MPI_COMM_WORLD);

			}

		}
	}
	else {

		MPI_Bcast(&kernelsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
		kernel = new int[kernelsize * kernelsize];
		MPI_Bcast(kernel, kernelsize * kernelsize, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&paddedheight, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&paddedwidth, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(paddedimage, paddedheight * paddedwidth, MPI_INT, 0, MPI_COMM_WORLD);
		printf("another process than zero printing padded height as %d", paddedheight);


		MPI_Recv(&startindex, 1, MPI_INT, 0, TagStartindex, MPI_COMM_WORLD, &ReveiveStatus);
		MPI_Recv(&rowsperprocessor, 1, MPI_INT, 0, TagProcessorRows, MPI_COMM_WORLD, &ReveiveStatus);
		MPI_Recv(&islast, 1, MPI_INT, 0, TagIsLast, MPI_COMM_WORLD, &ReveiveStatus);
		MPI_Recv(&localimagesize, 1, MPI_INT, 0, TaglocalimageSize, MPI_COMM_WORLD, &ReveiveStatus);
	}




	int* localImage = applyKernel(kernel, kernelsize, paddedimage, startindex, rowsperprocessor, paddedwidth, subsize);
	createImage(localImage, ImageWidth, ImageHeight, 1);
	if(rank!=0){

		MPI_Send(&subsize, 1, MPI_INT, 0, TaglocalimageSize, MPI_COMM_WORLD);
		MPI_Send(localImage, subsize, MPI_INT, 0, TagLocalImage, MPI_COMM_WORLD);
		MPI_Send(&startindex, 1, MPI_INT, 0, TagStartindex, MPI_COMM_WORLD);
		printf("process %d sent sub image \n", rank);
	
	}







	if (rank == 0) {
		//createImage(paddedimage, paddedwidth, paddedheight, 1);
		//std::cout << "time: " << TotalTime <<std:: endl;
	}


	std::free(imageData);
	MPI_Finalize();
	return 0;

}



