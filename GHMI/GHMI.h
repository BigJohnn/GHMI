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
//Copyright(c) 2016, Hou Peihong
//All rights reserved.
//
//Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met :
//
//1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//
//2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and / or other materials provided with the distribution.
//
//3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
//BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
//GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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