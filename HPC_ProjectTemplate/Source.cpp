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
double*** inputImage(int* w, int* h, int* originalImageWidth, int* originalImageHeight, int padding, System::String^ imagePath);
void createImage(double**** image, int width, int height, int index);
void applyFilter(double**** image, int imageWidth, int imageHeight, double****result, int originalImageWidth, int orginalImageHeight, double*** kernel, int kernelSize);
void initResultArray(double**** image, int width, int height);
void freeImage(double**** image, int width, int height);

int main()
{
	int ImageWidth = 0, ImageHeight = 0;
	int originalImageWidth, originalImageHeight = 0;

	int start_s, stop_s, TotalTime = 0;

	System::String^ imagePath;
	std::string img;
	img = "..//Data//Input//lena.png";

	imagePath = marshal_as<System::String^>(img);
	int kernelSize = 0;
	double** kernel;

	inputKernel(&kernel, &kernelSize);
	printKernel(&kernel, kernelSize, kernelSize);

	int padding = calculatePaddingValue(kernelSize);
	double*** imageData = inputImage(&ImageWidth, &ImageHeight, &originalImageWidth, &originalImageHeight, padding, imagePath);
	
	double*** resultImage;
	initResultArray(&resultImage,originalImageWidth,originalImageHeight);

	applyFilter(&imageData,ImageWidth, ImageHeight,&resultImage, originalImageWidth, originalImageHeight,&kernel,kernelSize);


	//start_s = clock();
	//applyFilter(&imageData, ImageWidth * ImageHeight,ImageWidth,ImageHeight ,&resultImage,originalImageWidth,originalImageWidth,&kernel,kernelSize);
	//stop_s = clock();
	//TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
	//cout << "time: " << TotalTime << endl;
	createImage(&resultImage, originalImageWidth, originalImageHeight, 1);

	freeKernel(&kernel, kernelSize);
	freeImage(&resultImage, originalImageWidth, originalImageHeight);
	freeImage(&imageData, ImageWidth, ImageHeight);
	//free(imageData);
	//delete[] resultImage;
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

double*** inputImage(int* w, int* h, int * originalImageWidth, int* originalImageHeight, int padding, System::String^ imagePath){
	double*** input;
	int OriginalImageWidth, OriginalImageHeight;

	System::Drawing::Bitmap BM(imagePath);

	OriginalImageWidth = BM.Width;
	OriginalImageHeight = BM.Height;
	cout << OriginalImageWidth << " " << OriginalImageHeight<<endl;

	*w = BM.Width + 2 * padding;
	*h = BM.Height + 2 * padding;
	*originalImageWidth = OriginalImageWidth;
	*originalImageHeight = OriginalImageHeight;

	input = new double**[3];
	input[0] = new double* [*h];
	input[1] = new double* [*h];
	input[2] = new double* [*h];

	for (int i = 0; i < *h; i++) { 
		input[0][i] = new double[*w];
		input[1][i] = new double[*w];
		input[2][i] = new double[*w];

		for (int j = 0; j < *w; j++) { 

			if (i < padding || j < padding || j > BM.Width - 1 + padding || i > BM.Height - 1 + padding) {
				input[0][i][j] = 0;
				input[1][i][j] = 0;
				input[2][i][j] = 0;
			}
			else {
				System::Drawing::Color c = BM.GetPixel(j - padding, i - padding);
				input[0][i][j] = c.R;
				input[1][i][j] = c.G;
				input[2][i][j] = c.B;
			}
		}
	}
	return input;
}

void createImage(double**** image, int width, int height, int index){
	System::Drawing::Bitmap MyNewImage(width, height);

	for (int i = 0; i < MyNewImage.Height; i++) {
		for (int j = 0; j < MyNewImage.Width; j++) {
			for (int k = 0; k < 3; k++) {
			
				if ((*image)[k][i][j] < 0) {
					(*image)[k][i][j] = 0;

				}
				else if ((*image)[k][i][j] > 255) {
				
					(*image)[k][i][j] = 255;
				}
			}

			System::Drawing::Color c = System::Drawing::Color::FromArgb((*image)[0][i][j], (*image)[1][i][j], (*image)[2][i][j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}

	MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");
	cout << "result Image Saved " << index << endl;
}

void applyFilter(double**** image, int imageWidth, int imageHeight, double**** result, int originalImageWidth, int orginalImageHeight, double*** kernel, int kernelSize){

	int padding = (imageWidth - originalImageWidth) / 2;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < imageHeight - kernelSize + 1; j++) { //0 -> 7
			for (int k = 0; k < imageWidth - kernelSize + 1; k++) { // 0 -> 7
			

				for (int l = j; l < j + kernelSize; l++) {
					for (int m = k; m < k + kernelSize; m++) {
						(*result)[i][j][k] += (*image)[i][l][m] * (*kernel)[l - j][m - k];
					}
				}
			}
		}
	}
}

void initResultArray(double**** image, int width, int height) {
	(*image) = new double** [3];
	(*image)[0] = new double* [height];
	(*image)[1] = new double* [height];
	(*image)[2] = new double* [height];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < height; j++) {
			(*image)[i][j] = new double[width];
			for (int k = 0; k < width; k++) {
				(*image)[i][j][k] = 0;
			}
		}
	}
}

void freeImage(double**** image, int width, int height) {

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < height; j++) {
			delete[](*image)[i][j];
		}
		delete[](*image)[i];
	}
	delete[](*image);
}
