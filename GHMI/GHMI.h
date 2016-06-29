#pragma once
#include <iostream>
#include <opencv2\opencv.hpp>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <afxwin.h>
#include <stdlib.h>

using namespace std;
using namespace cv;

#define N 5
#define imagenum 39
#define InvDim 18
#define _STR(s) #s
#define STR(s) _STR(s) 

typedef long double Invtype;

typedef enum method
{
	No_Method = -1,
	Eu_Norm = 0,
	Chi_Norm
}Discriminationmethod;

class GHMI
{
public:
	GHMI();
	~GHMI();
	void FindAllFile(CString path, CString* filenames, int& count);
	double H(int p, double x);
	vector<Invtype> Calculate18GHMIs4SingleImage(vector<Invtype> &Invariant, Mat &I, double sigma);
	vector<Invtype> CalImgGHMIs(String &filename);
};