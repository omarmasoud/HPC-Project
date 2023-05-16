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
int	calculatePaddingValue(int kernelSize);
double* inputImage(int* w, int* h, int* originalImageWidth, int* originalImageHeight, int padding, System::String^ imagePath);
void createImage(double* image, int width, int height, int index);
void applyFilter(double** image, double**result, int imageWidth, int imageHeight, double*** kernel, int kernelSize);

int main()
{
	int ImageWidth = 0, ImageHeight = 0;
	int originalImageWidth, originalImageHeight = 0;

	int start_s, stop_s, TotalTime = 0;

	System::String^ imagePath;
	std::string img;
	img = "..//Data//Input//test.png";

	imagePath = marshal_as<System::String^>(img);	
	int kernelSize = 0;
	double** kernel;

	inputKernel(&kernel, &kernelSize);
	printKernel(&kernel,kernelSize,kernelSize);
	int padding = calculatePaddingValue(kernelSize);
	double* imageData = inputImage(&ImageWidth, &ImageHeight,&originalImageWidth, &originalImageHeight, padding, imagePath);

	double* resultImage = new double[originalImageWidth * originalImageHeight];
	//applyFilter(&imageData,originalImageWidth, originalImageHeight,&kernel,kernelSize);
	//start_s = clock();
	//stop_s = clock();
	//TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
	//createImage(imageData, ImageWidth, ImageHeight, 1);
	//cout << "time: " << TotalTime << endl;
	createImage(imageData, ImageWidth, ImageHeight, 1);

	freeKernel(&kernel, kernelSize);
	free(imageData);
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

int calculatePaddingValue(int kernelSize) {
	return (kernelSize - 1) / 2;
}

double* inputImage(int* w, int* h, int * originalImageWidth, int* originalImageHeight, int padding, System::String^ imagePath){
	double* input;
	int OriginalImageWidth, OriginalImageHeight;

	System::Drawing::Bitmap BM(imagePath);

	OriginalImageWidth = BM.Width;
	OriginalImageHeight = BM.Height;
	cout << OriginalImageWidth << " " << OriginalImageHeight<<endl;

	*w = BM.Width + 2 * padding;
	*h = BM.Height + 2 * padding;
	*originalImageWidth = OriginalImageWidth;
	*originalImageHeight = OriginalImageHeight;

	double* Red = new double[(*h) * (*w)];
	double* Green = new double[(*h) * (*w)];
	double* Blue = new double[(*h) * (*w)];
	input = new double[(*h) * (*w)];

	for (int i = 0; i < *w; i++) { 
		for (int j = 0; j < *h; j++) { 

			if (i < padding || j < padding || i > BM.Width - 1 + padding || j > BM.Height - 1 + padding) {
				Red[i * (*h)+ j] = 0;
				Blue[i * (*h) + j] = 0;
				Green[i * (*h) + j] = 0;
				input[i * (*h) + j] = 0; //gray scale value equals the average of RGB values
			}
			else {
				System::Drawing::Color c = BM.GetPixel(i - padding, j - padding);
				Red[i * (*h) + j] = c.R;
				Blue[i * (*h) + j] = c.B;
				Green[i * (*h) + j] = c.G;
				input[i * (*h) + j] = ((c.R + c.B + c.G) / 3); //gray scale value equals the average of RGB values
			}
		}
	}

	return input;
}

void createImage(double* image, int width, int height, int index){
	System::Drawing::Bitmap MyNewImage(width, height);

	for (int i = 0; i < MyNewImage.Width; i++) {
		for (int j = 0; j < MyNewImage.Height; j++) {
			if (image[i * height + j] < 0)
			{
				image[i * height + j] = 0;
			}
			if (image[i * height + j] > 255)
			{
				image[i * height + j] = 255;
			}
			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i * MyNewImage.Height + j], image[i * MyNewImage.Height + j], image[i * MyNewImage.Height + j]);
			MyNewImage.SetPixel(i, j, c);
		}
	}

	MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");
	cout << "result Image Saved " << index << endl;
}

void applyFilter(double** image, int imageWidth, int imageHeight, double*** kernel, int kernelSize) {


}