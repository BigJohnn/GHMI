#include "GHMI.h"

GHMI::GHMI()
{

}

GHMI::~GHMI()
{

}

void GHMI::FindAllFile(CString path, CString* filenames, int& count)
{
	CFileFind finder;
	BOOL working = finder.FindFile(path + "\\*.*");
	while (working)
	{
		working = finder.FindNextFile();
		if (finder.IsDots())
			continue;
		if (finder.IsDirectory())
		{
			//FindAllFile(finder.GetFilePath(), filenames, count);
		}
		else
		{
			CString filename = finder.GetFileName();
			filenames[count++] = filename;
		}
	}
}

//Calculate Hermite Polynomials by its recursive relationship.
double GHMI::H(int p, double x)
{
	if (p == 0)
		return 1;
	else if (p == 1)
		return 2 * x;
	else
		return 2 * x*H(p - 1, x) - 2 * (p - 1)*H(p - 2, x);
}

vector<Invtype> GHMI::Calculate18GHMIs4SingleImage(vector<Invtype> &Invariant, Mat &I, double sigma)
{
	double M[6][6] = { 0.0 };

	//Calculate--GHMs--00~55
	for (int q = 0; q <= N; q++)
	{
		for (int p = 0; p <= N; p++)
		{
			for (int x = 0; x < I.rows; x++)
			{
				uchar* data = I.ptr<uchar>(x);

				for (int y = 0; y < I.cols; y++)
				{
					M[p][q] += data[y]*1.0/255 * H(p, x / sigma)*H(q, y / sigma)*exp(-(x*x + y*y) / 2 / (sigma*sigma));
				}
			}
		}
	}

	//int K = I.cols; //For K * K image I, I.cols = I.rows = K.
	//
	//double c = 2.0 / K - 1;
	//for (int n = 0; n <= N; n++)
	//{
	//	for (int i = 0; i < K; i++)
	//	{
	//		double s = 0.0;
	//		uchar* data = I.ptr<uchar>(i);
	//		for (int j = 0; j < K; j++)
	//		{
	//			s += c*data[j] * H(n, i);  //data[j] is the pixel value of I(i,j).
	//		}
	//		for (int m = 0; m <= N; m++)
	//		{
	//			M[m][n] += c*s*H(m, i *1.0 / sigma);
	//			//cout << "M[m][n] " << m << " " << n << " " << M[m][n] << endl;
	//		}
	//	}
	//}

	Invtype Inv[InvDim] = { 0.0 };

	Inv[0] = M[2][0] + M[0][2];

	Inv[1] = (M[3][0] + M[1][2])*(M[3][0] + M[1][2]) + (M[2][1] + M[0][3])*(M[2][1] + M[0][3]);

	Inv[2] = (M[2][0] - M[0][2])*((M[3][0] + M[1][2])*(M[3][0] + M[1][2]) - (M[2][1] + M[0][3])*
		(M[2][1] + M[0][3])) + 4 * M[1][1] * (M[3][0] + M[1][2])*(M[2][1] + M[0][3]);

	Inv[3] = M[1][1] * ((M[3][0] + M[1][2])*(M[3][0] + M[1][2]) - (M[2][1] + M[0][3])*(M[2][1] + M[0][3])) -
		(M[2][0] - M[0][2])*(M[3][0] + M[1][2])*(M[2][1] + M[0][3]);

	Inv[4] = (M[3][0] - 3 * M[1][2])*(M[3][0] + M[1][2])*((M[3][0] + M[1][2])*(M[3][0] + M[1][2]) - 3 * (M[2][1] + M[0][3])*(M[2][1] + M[0][3])) +
		(M[0][3] - 3 * M[2][1])*(M[2][1] + M[0][3])*((M[2][1] + M[0][3])*(M[2][1] + M[0][3]) - 3 * (M[3][0] + M[1][2])*(M[3][0] + M[1][2]));

	Inv[5] = (M[3][0] - 3 * M[1][2])*(M[2][1] + M[0][3])*((M[2][1] + M[0][3])*(M[2][1] + M[0][3]) - 3 * (M[3][0] + M[1][2])*(M[3][0] + M[1][2])) +
		(3 * M[2][1] - M[0][3])*(M[3][0] + M[1][2])*((M[3][0] + M[1][2])*(M[3][0] + M[1][2]) - 3 * (M[2][1] + M[0][3])*(M[2][1] + M[0][3]));

	Inv[6] = M[4][0] + 2 * M[2][2] + M[0][4];

	Inv[7] = (M[4][0] - M[0][4])*((M[3][0] + M[1][2])*(M[3][0] + M[1][2]) - (M[2][1] + M[0][3])*(M[2][1] + M[0][3])) +
		4 * (M[3][1] + M[1][3])*(M[3][0] + M[1][2])*(M[2][1] + M[0][3]);

	Inv[8] = (M[3][1] + M[1][3])*((M[3][0] + M[1][2])*(M[3][0] + M[1][2]) - (M[2][1] + M[0][3])*(M[2][1] + M[0][3])) -
		(M[4][0] - M[0][4])*(M[3][0] + M[1][2])*(M[2][1] + M[0][3]);

	Inv[9] = (M[4][0] - 6 * M[2][2] + M[0][4])*((M[3][0] + M[1][2])*(M[3][0] + M[1][2])*(M[3][0] + M[1][2])*(M[3][0] + M[1][2]) - 6 * (M[3][0] + M[1][2])*
		(M[3][0] + M[1][2])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3]) + (M[2][1] + M[0][3])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3])) +
		16 * (M[3][1] - M[1][3])*(M[3][0] + M[1][2])*(M[2][1] + M[0][3])*((M[3][0] + M[1][2])*(M[3][0] + M[1][2]) - (M[2][1] + M[0][3])*(M[2][1] + M[0][3]));

	Inv[10] = (M[4][0] - 6 * M[2][2] + M[0][4])*(M[3][0] + M[1][2])*(M[2][1] + M[0][3])*((M[2][1] + M[0][3])*(M[2][1] + M[0][3]) - (M[3][0] + M[1][2])*(M[3][0] + M[1][2])) +
		(M[3][1] - M[1][3])*((M[3][0] + M[1][2])*(M[3][0] + M[1][2])*(M[3][0] + M[1][2])*(M[3][0] + M[1][2]) - 6 * (M[3][0] + M[1][2])*(M[3][0] + M[1][2])*
		(M[2][1] + M[0][3])*(M[2][1] + M[0][3]) + (M[2][1] + M[0][3])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3]));

	Inv[11] = (M[5][0] + 2 * M[3][2] + M[1][4])*(M[5][0] + 2 * M[3][2] + M[1][4]) + (M[4][1] + 2 * M[2][3] + M[0][5])*(M[4][1] + 2 * M[2][3] + M[0][5]);

	Inv[12] = (M[5][0] + 2 * M[3][2] + M[1][4])*(M[3][0] + M[1][2]) + (M[4][1] + 2 * M[2][3] + M[0][5])*(M[2][1] + M[0][3]);

	Inv[13] = (M[4][1] + 2 * M[2][3] + M[0][5])*(M[3][0] + M[1][2]) - (M[5][0] + 2 * M[3][2] + M[1][4])*(M[2][1] + M[0][3]);

	Inv[14] = (M[5][0] - 2 * M[3][2] - 3 * M[1][4])*((M[3][0] + M[1][2])*(M[3][0] + M[1][2])*(M[3][0] + M[1][2]) - 3 * (M[3][0] + M[1][2])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3])) -
		(3 * M[4][1] + 2 * M[2][3] - M[0][5])*((M[2][1] + M[0][3])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3]) - 3 * (M[2][1] + M[0][3])*(M[3][0] + M[1][2])*(M[3][0] + M[1][2]));

	Inv[15] = (M[5][0] - 2 * M[3][2] - 3 * M[1][4])*((M[2][1] + M[0][3])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3]) - 3 * (M[2][1] + M[0][3])*(M[3][0] + M[1][2])*(M[3][0] + M[1][2])) +
		3 * (M[4][1] + 2 * M[2][3] - M[0][5])*((M[3][0] + M[1][2])*(M[3][0] + M[1][2])*(M[3][0] + M[1][2]) - 3 * (M[3][0] + M[1][2])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3]));

	Inv[16] = (M[5][0] - 10 * M[3][2] + 5 * M[4][1])*(pow(M[3][0] + M[1][2], 5) - 10 * (M[3][0] + M[1][2])*(M[3][0] + M[1][2])*(M[3][0] + M[1][2])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3]) +
		5 * (M[3][0] + M[1][2]) * (M[2][1] + M[0][3])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3])) + (5 * M[4][1] - 10 * M[2][3] + M[0][5])*(
		pow(M[2][1] + M[0][3], 5) - 10 * (M[3][0] + M[1][2])*(M[3][0] + M[1][2])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3]) + 5 * (M[2][1] + M[0][3])*
		(M[3][0] + M[1][2])*(M[3][0] + M[1][2])*(M[3][0] + M[1][2])*(M[3][0] + M[1][2]));

	Inv[17] = (M[0][5] - 10 * M[2][3] + 5 * M[4][1])*(pow(M[3][0] + M[1][2], 5) - 10 * (M[3][0] + M[1][2])*(M[3][0] + M[1][2])*(M[3][0] + M[1][2])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3]) +
		5 * (M[3][0] + M[1][2]) * (M[2][1] + M[0][3])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3])) - (5 * M[1][4] - 10 * M[3][2] + M[5][0])*(
		pow(M[2][1] + M[0][3], 5) - 10 * (M[3][0] + M[1][2])*(M[3][0] + M[1][2])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3])*(M[2][1] + M[0][3]) + 5 * (M[2][1] + M[0][3])*
		(M[3][0] + M[1][2])*(M[3][0] + M[1][2])*(M[3][0] + M[1][2])*(M[3][0] + M[1][2]));

	double S = 0.0;

	for (int i = 0; i <= InvDim; i++)
	{
		S += Inv[i] * Inv[i];
	}

	double D = sqrt(S);
	if (0 == D)
		D = 1;

	for (int i = 0; i < InvDim; i++)
	{
			Inv[i] /= D;
		Invariant.push_back(Inv[i]);
	}
	return Invariant;
}

vector<Invtype> GHMI::CalImgGHMIs(String &filename)
{
	Mat src = imread(filename, 0);

	//int a = 
	Mat I, Igrad, Ilaplace;

	/*Mat Idiffx, Idiffy;
	Mat Idiff2x, Idiff2y;*/

	resize(src, I, Size(100, 100));

	Sobel(I, Igrad, I.depth(), 1, 1);
	Laplacian(I, Ilaplace, I.depth());

	if (!Igrad.isContinuous())
		Igrad.reshape(1, Igrad.cols*Igrad.rows);


	if (!Ilaplace.isContinuous())
		Ilaplace.reshape(1, Ilaplace.cols*Ilaplace.rows);

	vector<Invtype> Invariant;

	//Get Invariant whose dimension is 4*InvDim(18)*3¡Ö216.
	double sigma = 0.8;
	for (int i = 0; i < 4; i++)
	{
		Invariant = GHMI::Calculate18GHMIs4SingleImage(Invariant, I, sigma);
		Invariant = GHMI::Calculate18GHMIs4SingleImage(Invariant, Igrad, sigma);
		Invariant = GHMI::Calculate18GHMIs4SingleImage(Invariant, Ilaplace, sigma);
		sigma += 0.2;
	}

	return Invariant;
}