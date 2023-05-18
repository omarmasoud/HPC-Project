#include <iostream>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include<mpi.h>
#include <ctime>// include this header 
#include<vector>
#include<algorithm>

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
		size+=1; // Increment size to make it odd

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


int* applyKernel(int* kernel, int kernelSize, const int* paddedImage, int startIndex, int rowsPerProcessor, int paddedWidth, int paddedHeight, int& subsize) {




	/*int subImageWidth = paddedWidth - (kernelSize - 1);
	int subImageHeight = rowsPerProcessor - (kernelSize - 1);*/
	vector<int> result;
	//// intit result  = 0 ;
	//for (int i = 0; i < subImageWidth * subImageHeight; i++) {
	//	result[i] = 0;
	//}
	//printf("startindex is %d\n", startIndex);
	for (int i = startIndex; i < startIndex + (rowsPerProcessor * paddedWidth); i++) {

		int x = i % paddedWidth;
		int y = i / paddedWidth;

		int bottomRightX = x + kernelSize - 1;
		int bottomRightY = y + kernelSize - 1;
	
		if (bottomRightX > paddedWidth - 1 || bottomRightY > paddedHeight - 1) {
			
			/*if (bottomRightY > paddedHeight - 1) {
				printf("true bottom right y");
			}
			else{ printf("false bottom right y"); }*/
			//printf("x is %d y is %d\n", x, y);
			continue;
			
		}


		else {
			int temp = 0;
			for (int j = y; j < y + kernelSize; j++) {
				
				for (int k = x; k < x + kernelSize; k++) {
					 temp += paddedImage[j * paddedWidth + k] * kernel[(j - y) * kernelSize + (k - x)];
					
					
					
				}
			}
			result.push_back(temp);
		}
	}
	//printf("kernel applied vector size is  %d\n",result.size());
	subsize = result.size();
	int* returnresult=new int[result.size()];
	copy(result.begin(), result.end(), returnresult);
	return returnresult;
}


int main()
{
	int ImageWidth = 4, ImageHeight = 4;

	
	int rank, size;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	

	
	int kernelsize = 0;
	int* kernel = NULL;
	int paddedheight;
	int paddedwidth;
	int* paddedimage ;
	int startindex=0;
	int rowssent;
	int islast;
	int localimagesize;
	int rowsperprocessor;
	int subsize;
	MPI_Status ReveiveStatus;
	int* resultingimage;
	int remainingrows;

	int start_s, stop_s, TotalTime = 0;
	if (rank == 0) {
		


		System::String^ imagePath;
		std::string img;
		img = "..//Data//Input//lena.png";
		


		
		

		cout << "enter the kernel size you want to create" << endl;
		cin >> kernelsize;
		int choice;
		cout << "enter 1 to work on lena image by default and 2 to enter your image name" << endl;
		cin >> choice;
		if (choice != 1) {
			string imgname;
			cout << "enter your image name" << endl;
			cin >> imgname;
			img = "..//Data//Input//" + imgname + ".png";

		}
		imagePath = marshal_as<System::String^>(img);
		int* imageData = inputImage(&ImageWidth, &ImageHeight, imagePath);
		//kernelsize = 7;
		start_s = clock();

		
		kernel = createHighPassKernel(kernelsize);

		if (kernelsize % 2 == 0)
			kernelsize += 1;
		
		
		
		
		
		computePadding(imageData, ImageHeight, ImageWidth, kernelsize, paddedheight, paddedwidth, paddedimage);

		MPI_Bcast(&kernelsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(kernel, kernelsize * kernelsize, MPI_INT, 0, MPI_COMM_WORLD);

		MPI_Bcast(&paddedheight, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&paddedwidth, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(paddedimage, paddedheight * paddedwidth, MPI_INT, 0, MPI_COMM_WORLD);
		
		////// Distribution
		rowsperprocessor = paddedheight / (size);
		remainingrows = paddedheight % size;
		localimagesize = rowsperprocessor * paddedwidth;
		//printf("\nlsuasdiuasduiauiscbidasuduias %d ",(sizeof(paddedimage)/ sizeof(paddedimage[0])));
		for (int i = 1; i < size; i++) {
			startindex = rowsperprocessor * paddedwidth * (i-1);
			/*if (i == size - 1) {
				rowsperprocessor = paddedheight - (rowsperprocessor * (size - 1));
				islast = 1;
			}
			else {
				rowsperprocessor = paddedheight / size;
				islast = 0;
			}*/
			
			

			MPI_Send(&startindex, 1, MPI_INT, i, TagStartindex, MPI_COMM_WORLD);
			MPI_Send(&rowsperprocessor, 1, MPI_INT, i, TagProcessorRows, MPI_COMM_WORLD);
			//MPI_Send(&islast, 1, MPI_INT, i, TagIsLast, MPI_COMM_WORLD);
			MPI_Send(&localimagesize, 1, MPI_INT, i, TaglocalimageSize, MPI_COMM_WORLD);

			

		}
		if (size != 1) {
			startindex += localimagesize;
			rowsperprocessor += remainingrows;
			//islast = 0;
			localimagesize = rowsperprocessor * paddedwidth;
		}
	}
	else {
		
		MPI_Bcast(&kernelsize, 1, MPI_INT, 0, MPI_COMM_WORLD);

		kernel = new int[kernelsize * kernelsize];
		
		MPI_Bcast(kernel, kernelsize * kernelsize, MPI_INT, 0, MPI_COMM_WORLD);
	
		MPI_Bcast(&paddedheight, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		MPI_Bcast(&paddedwidth, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
		paddedimage = new int[paddedwidth * paddedheight];

		MPI_Bcast(paddedimage, paddedheight * paddedwidth, MPI_INT, 0, MPI_COMM_WORLD);
	//	printf("another process than zero printing padded height as %d", paddedheight);
	

		MPI_Recv(&startindex, 1, MPI_INT, 0, TagStartindex, MPI_COMM_WORLD, &ReveiveStatus);

		
		MPI_Recv(&rowsperprocessor, 1, MPI_INT, 0, TagProcessorRows, MPI_COMM_WORLD, &ReveiveStatus);
		//MPI_Recv(&islast, 1, MPI_INT, 0, TagIsLast, MPI_COMM_WORLD, &ReveiveStatus);
		MPI_Recv(&localimagesize, 1, MPI_INT, 0, TaglocalimageSize, MPI_COMM_WORLD, &ReveiveStatus);
		
	}
	

	int* localImage = applyKernel(kernel, kernelsize, paddedimage, startindex, rowsperprocessor, paddedwidth,paddedheight, subsize);
	if (rank != 0) {


		MPI_Send(&subsize, 1, MPI_INT, 0, TaglocalimageSize, MPI_COMM_WORLD);
		MPI_Send(localImage, subsize, MPI_INT, 0, TagLocalImage, MPI_COMM_WORLD);

	}
	if (rank == 0) {
		int receivedsubsize;
		int* recievedimg=new int[ImageHeight*ImageWidth];
		
		for (int i = 1; i < size; i++) {
			MPI_Recv(&receivedsubsize, 1, MPI_INT, i, TaglocalimageSize, MPI_COMM_WORLD, &ReveiveStatus);
			//recievedimg = new int[receivedsubsize];
			MPI_Recv(&recievedimg[(i-1)*receivedsubsize], receivedsubsize, MPI_INT, i, TagLocalImage, MPI_COMM_WORLD,&ReveiveStatus);
		}
		for (int i = (size - 1) * receivedsubsize; i < ImageHeight * ImageWidth; i++) {
			recievedimg[i] = localImage[i - (size - 1) * receivedsubsize];
		}
		createImage(recievedimg, ImageWidth, ImageHeight, 1);
		stop_s = clock();
		TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;

		printf("total time elapsed by %d processes to go through image of dimensions %d by %d is %d ms \n", size, ImageWidth, ImageHeight, TotalTime);
	}
	


	
	


	




	//std::free(localImage);
	MPI_Finalize();
	return 0;

}



