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

void Neighbour(const unsigned char* pbIn, int lWidthIn, int lHeightIn, unsigned char* pbOut, int lWidthOut, int lHeightOut) {
	for (int i = 0; i < lHeightOut; ++i)
	{
		double yy = (double)i * (double)lHeightIn / lHeightOut;
		int y = (int)yy;
		for (int j = 0; j < lWidthOut; ++j)
		{
			double xx = (double)j * (double)lWidthIn / lWidthOut;
			int x = (int)xx;
			int lP00 = pbIn[y * lWidthIn + x],
				lP01 = pbIn[y * lWidthIn + x + 1],
				lP10 = pbIn[(y + 1) * lWidthIn + x],
				lP11 = pbIn[(y + 1) * lWidthIn + x + 1];
			if (((int)(xx + 0.5) - x) > 0) {
				if (((int)(yy + 0.5) - y) > 0) {
					pbOut[i * lWidthOut + j] = lP11;
				}
				else {
					pbOut[i * lWidthOut + j] = lP10;
				}
			}
			else {
				if (((int)(yy + 0.5) - y) > 0) {
					pbOut[i * lWidthOut + j] = lP01;
				}
				else {
					pbOut[i * lWidthOut + j] = lP00;
				}
			}
		}
	}
}

void ResizeBilinear(const unsigned char* pbIn, int lWidthIn, int lHeightIn, unsigned char* pbOut, int lWidthOut, int lHeightOut) {
	for (int i = 0; i < lHeightOut; ++i)
	{
		double yy = (double)i * (double)lHeightIn / lHeightOut;
		int y = (int)yy;
		double u = yy - (double)y;
		for (int j = 0; j < lWidthOut; ++j)
		{
			double xx = (double)j * (double)lWidthIn / lWidthOut;
			int x = (int)xx;
			double v = xx - (double)x; 
			int lP00 = pbIn[y * lWidthIn + x],
				lP01 = pbIn[y * lWidthIn + x + 1],
				lP10 = pbIn[(y + 1) * lWidthIn + x],
				lP11 = pbIn[(y + 1) * lWidthIn + x + 1];
			pbOut[i * lWidthOut + j] = (unsigned char)((1. - u) * (1. - v) * lP00 + u * (1. - v) * lP10 + v * (1. - u) * lP01 + u * v * lP11);
		}
	}
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

	int outWidth = 4096;
	int outHeight = 2048;
	const int OutSize = 4096*2048;

	unsigned char* pOutputBits = new unsigned char[OutSize];
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
	vector<int> Pbimg(nWidth * nHeight);

	//ResizeBilinear(pInputBits, nWidth, nHeight, pOutputBits, outWidth, outHeight);
	Neighbour(pInputBits, nWidth, nHeight, pOutputBits, outWidth, outHeight);

	if (NPngProc::writePngFile(szOutputFileName, pOutputBits, outWidth, outHeight, nBPP) == NPngProc::PNG_ERROR)
	{
		printf("\nError ocuured during png file was written");
		return -1;
	}

	return 0;
}