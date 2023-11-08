// pngtest.cpp : Defines the entry point for the console application.
//

#include "PngProc.h"
#include <string.h>
#include <stdio.h>
#include <random>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

void AWGN_png(unsigned char* pOut, unsigned char* pIn, size_t nWidth, size_t nHeight, double s) {
	random_device rd;
	mt19937 gen(rd());
	normal_distribution<> X(0, s);
	double noise;
	for (size_t y = 0; y < nHeight; ++y)
	{
		for (size_t x = 0; x < nWidth; ++x)
		{
			noise = X(gen);
			if ((*pIn + noise) <= 0) {
				*pOut++ = 0;
			}
			else if ((*pIn + noise) >= 256) {
				*pOut++ = 255;
			}
			else if (((*pIn + noise) > 0) && ((*pIn + noise) < 256)) {
				*pOut++ = noise + *pIn++;
			}
		}
	}
	return;
}

void Impulse_png(unsigned char* pOut, unsigned char* pIn, size_t nWidth, size_t nHeight, double P, int a = 0, int b = 255) {
	double tmp;
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> distrib(0, 100);
	for (size_t y = 0; y < nHeight; ++y)
	{
		for (size_t x = 0; x < nWidth; ++x)
		{
			tmp = distrib(gen);
			if (tmp <= (100*P / 2)) {
				*pOut++ = a;
				pIn++;
			}
			if ((tmp <= (100*P)) && (tmp > (100 * P / 2))) {
				*pOut++ = b;
				pIn++;
			}
			if (tmp > (100 * P)) {
				*pOut++ = *pIn++;
			}
		}
	}
	return;
}

void Preprocessing(unsigned char* pIn, size_t nWidth, size_t nHeight, vector<int> &Pbimg) {
	
	for (size_t y = 0; y < nHeight; ++y)
	{
		for (size_t x = 0; x < nWidth; ++x)
		{
			Pbimg[nWidth*y + x] = *pIn++;
		}
	}
}

void Histogram4Gray8bpp(vector<int> &pbImg, int iHeight, int iWidth, int iWidthBytes, unsigned long* pulHist){

	// начальное обнуление массива гистограммы
	memset(pulHist, 0, 256 * sizeof(*pulHist));
	// вычисление гистограммы
	for (int y = 0; y < iHeight; ++y)
	{
		for (int x = 0; x < iWidth; ++x)
		{
			pulHist[pbImg[iWidth * y + x]]++;
		}
	}
}

double mean_value(unsigned long* pulHist, int nWidth, int nHeight, int bins) {
	double mean=0.;
	for (int i = 0; i < bins; i++) {
		double temp = mean;
		mean = temp + i * ((double)pulHist[i] / nWidth / nHeight);
	}
	cout << "average value: " << mean << "\n";
	return mean;
}

double variance(unsigned long* pulHist, int nWidth, int nHeight, int bins, double avg) {
	double var=0.;
	for (int i = 0; i < bins; i++) {
		double temp = var;
		var = temp + pow((i-avg), 2) * ((double)pulHist[i] / nWidth / nHeight);
	}
	cout << "variance: " << var << "\n";
	return var;
}

int quartile(unsigned long* pulHist, int nWidth, int nHeight, int bins, double num) {
	int k = 0;
	int tmp;
	for (int i = 0; i < bins; i++) {
		tmp = k;
		k = tmp + pulHist[i];
		if (k >= nWidth*nHeight * num) {
			return i;
		}
	}
	return 0;
}

void ImageProcessingGray(unsigned char* pOut, unsigned char* pIn, size_t nWidth, size_t nHeight) {
	random_device rd;
	mt19937 gen(rd());
	normal_distribution<> X(0, 20);
	for (size_t y = 0; y < nHeight; ++y)
	{
		for (size_t x = 0; x < nWidth; ++x)
		{
			*pOut++ = X(gen) + *pIn++;
		}
	}
	return;
}

int main(int argc, char* argv[])
{
	class CBitsPtrGuard
	{
	public:
		CBitsPtrGuard(unsigned char** pB) : m_ppBits(pB) { }
		~CBitsPtrGuard() { if (*m_ppBits) delete *m_ppBits, *m_ppBits = 0; }
	protected:
		unsigned char** m_ppBits;
	};

	// parse input parameters
	char	szInputFileName[256];
	char    szOutputFileName[256];
	if (argc < 2)
		printf("\nformat: pngtest <input_file> [<output_file>]");
	else 
	{
		strcpy(szInputFileName, argv[1]);
		if (argc > 2)
			strcpy(szOutputFileName, argv[2]);
		else
		{
			strcpy(szOutputFileName, szInputFileName);
			strcat(szOutputFileName, "_out.png");
		}
	}


	size_t nReqSize = NPngProc::readPngFile(szInputFileName, 0, 0, 0, 0);
	if (nReqSize == NPngProc::PNG_ERROR)
	{
		printf("\nError ocured while pngfile was read");
		return -1;
	}
	

	unsigned char* pInputBits = new unsigned char[nReqSize];
	if (!pInputBits)
	{
		printf("\nCan't allocate memory for image, required size is %u", nReqSize);
		return -1;
	}
	CBitsPtrGuard InputBitsPtrGuard(&pInputBits);


	unsigned char* pOutputBits = new unsigned char[nReqSize];
	if (!pOutputBits)
	{
		printf("\nCan't allocate memory for image, required size is %u", nReqSize);
		return -1;
	}


	CBitsPtrGuard OutputBitsPtrGuard(&pOutputBits);

	size_t nWidth, nHeight;
	unsigned int nBPP;

	size_t nRetSize = NPngProc::readPngFileGray(szInputFileName, pInputBits, &nWidth, &nHeight/*, &nBPP*/);
	nBPP = 8;


	// ASSERT(nRetSize == nReqSize);

	// TODO: image processing 
	double s = 20.;
	double P = 0.1;
	vector<int> Pbimg(nWidth * nHeight);
	const int bins = 256;
	unsigned long pulHist[bins];

	Preprocessing(pInputBits, nWidth, nHeight, Pbimg);
	Histogram4Gray8bpp(Pbimg, nHeight, nWidth, nBPP, pulHist);

	for (int i = 0; i < bins; i++) {
		cout << pulHist[i] << " ";
	}
	cout << "\n";

	double avg = mean_value(pulHist, nWidth, nHeight, bins);
	variance(pulHist, nWidth, nHeight, bins, avg);
	cout << "quartile distance: " << quartile(pulHist, nWidth, nHeight, bins, 0.75) - quartile(pulHist, nWidth, nHeight, bins, 0.25) << "\n";


	AWGN_png(pOutputBits, pInputBits, nWidth, nHeight, s);
	//Impulse_png(pOutputBits, pInputBits, nWidth, nHeight, P);
	
	//ImageProcessingGray(pOutputBits, pInputBits, nWidth, nHeight); 

	if (NPngProc::writePngFile(szOutputFileName, pOutputBits, nWidth, nHeight, nBPP) == NPngProc::PNG_ERROR)
	{
		printf("\nError ocuured during png file was written");
		return -1;
	}

	return 0;
}