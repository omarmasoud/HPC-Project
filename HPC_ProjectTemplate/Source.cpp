#include <iostream>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include <ctime>// include this header 
#include <omp.h>
#pragma once

#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>
using namespace std;
using namespace msclr::interop;

void printKernel(double*** kernel, int width, int height);
void freeKernel(double*** kernel, int width);
void inputKernel(double*** kernel, int* kernelSize);

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
	int *Red = new int[BM.Height * BM.Width];
	int *Green = new int[BM.Height * BM.Width];
	int *Blue = new int[BM.Height * BM.Width];
	input = new int[BM.Height*BM.Width];
	for (int i = 0; i < BM.Height; i++)
	{
		for (int j = 0; j < BM.Width; j++)
		{
			System::Drawing::Color c = BM.GetPixel(j, i);

			Red[i * BM.Width + j] = c.R;
			Blue[i * BM.Width + j] = c.B;
			Green[i * BM.Width + j] = c.G;

			input[i*BM.Width + j] = ((c.R + c.B + c.G) / 3); //gray scale value equals the average of RGB values

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
			if (image[i*width + j] < 0)
			{
				image[i*width + j] = 0;
			}
			if (image[i*width + j] > 255)
			{
				image[i*width + j] = 255;
			}
			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i*MyNewImage.Width + j], image[i*MyNewImage.Width + j], image[i*MyNewImage.Width + j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}
	MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");
	cout << "result Image Saved " << index << endl;
}


int main()
{
	int ImageWidth = 4, ImageHeight = 4;

	int start_s, stop_s, TotalTime = 0;

	System::String^ imagePath;
	std::string img;
	img = "..//Data//Input//test.png";

	imagePath = marshal_as<System::String^>(img);
	int* imageData = inputImage(&ImageWidth, &ImageHeight, imagePath);
	
	//start_s = clock();
	
	int kernelSize = 0;
	double** kernel;

	inputKernel(&kernel, &kernelSize);
	printKernel(&kernel,kernelSize,kernelSize);
	freeKernel(&kernel, kernelSize);

	//stop_s = clock();
	//TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
	//createImage(imageData, ImageWidth, ImageHeight, 1);
	//cout << "time: " << TotalTime << endl;

	//free(imageData);
	system("pause");
	return 0;

}

void printKernel(double*** kernel, int width, int height) {
	cout << "kernel is" << endl;
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
		
			cout << (*kernel)[i][j] <<" ";
		}
		cout << endl;
	}
}

void freeKernel(double*** kernel, int width) {
	for (int i = 0; i < width; i++) {
		delete[] (*kernel)[i];
	}
	delete[] *kernel;
}

void inputKernel(double*** kernel, int* kernelSize) {
	cout << "Enter Kernel Size: ";
	cin >> *kernelSize;
	cout << *kernelSize << endl;
	*kernel = new double* [*kernelSize];

	for (int i = 0; i < *kernelSize; i++) {
		(*kernel)[i] = new double[*kernelSize];
		for (int j = 0; j < *kernelSize; j++) {
			cout << "Kernel [" << i << "][" << j << "]: ";
			cin >> (*kernel)[i][j];
		}
	}

}