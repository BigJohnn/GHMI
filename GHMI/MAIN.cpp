#include "GHMI.h"

int main(int argc, char** argv)
{
	GHMI G;
	//String Dir = "F:\\GaussianHermite\\testingimages_salient_max\\";
	//int filenum0,filenum1;
	double EuclideanDist[imagenum][imagenum] = { 0.0 };
	double ChiSquareDist[imagenum][imagenum] = { 0.0 };

	char* mat_size = STR(imagenum);
	char* suffix = ".csv";

	char* outputdata = "order.txt";

	/*******Select the methods to distinguish.*********/
	Discriminationmethod method = Chi_Norm;		//Chi_Norm | Eu_Norm

	String Dir = (String)"HPH-" + STR(Chi_Norm) + "-" + (String)mat_size + "x" + (String)mat_size + (String)suffix;

	ofstream fout(Dir.c_str());
	wofstream foutd(outputdata);

	CString filenames[1024];
	int count = 0;
	char path[MAX_PATH] = "F:\\GaussianHermite\\testingimages_salient_max\\";
	vector<char*> Cstr;
	char* cstr;

	G.FindAllFile(path, filenames, count);

	for (int i = 1; i <= count; i++)
	{
		wcout << filenames[i - 1].GetString() << "\t" << i << endl;
		foutd << filenames[i - 1].GetString() << "\t" << i << endl;
		size_t len = wcslen(filenames[i - 1].GetString()) + 1;
		size_t converted = 0;
		cstr = (char*)malloc(len*sizeof(char));
		wcstombs_s(&converted, cstr, len, filenames[i - 1].GetString(), _TRUNCATE);
		Cstr.push_back(cstr);
	}

	vector< vector<Invtype> > Invs;
	vector<Invtype> inv;
	for (int i = 0; i < imagenum; i++)
	{
		String filename = (String)path + (String)Cstr[i];
		inv = G.CalImgGHMIs(filename);
		Invs.push_back(inv);
	}

	///************************
	fout << ""	"";
	for (int k = 1; k <= imagenum; k++)
		fout << "," << k;
	fout << endl;

	for (size_t filenum1 = 0; filenum1 != Invs.size(); filenum1++)
	{
		fout << filenum1 + 1 << ",";

		for (size_t filenum0 = 0; filenum0 < filenum1; filenum0++)
		{
			vector<Invtype> Inv0, Inv1;
			Inv0 = Invs[filenum0];
			Inv1 = Invs[filenum1];

			switch (method)
			{
			case Eu_Norm:
			{
				for (int i = 0; i < InvDim * 12; i++)
				{
					EuclideanDist[filenum0][filenum1] += (Inv0[i] - Inv1[i])*(Inv0[i] - Inv1[i]);
				}

				EuclideanDist[filenum0][filenum1] = sqrt(EuclideanDist[filenum0][filenum1]);

				fout << setiosflags(ios::fixed) << setprecision(4) << EuclideanDist[filenum0][filenum1] << ",";
			}
				break;
			case Chi_Norm:
			{
				for (int i = 0; i < InvDim * 12; i++)
				{
					if (0 == abs(Inv0[i] + Inv1[i]))
						ChiSquareDist[filenum0][filenum1] = (Inv0[i] - Inv1[i])*(Inv0[i] - Inv1[i]);
					else
						ChiSquareDist[filenum0][filenum1] += (Inv0[i] - Inv1[i])*(Inv0[i] - Inv1[i]) / abs(Inv0[i] + Inv1[i]);
				}

				fout << setiosflags(ios::fixed) << setprecision(4) << ChiSquareDist[filenum0][filenum1]*pow(10,60) << ",";
			}
				break;
			default:
				break;
			}
		}
		fout << endl;
	}

	/*Mat Ilaplace;
	Mat I = imread(Dir + "1.jpg", 0);
	Sobel(I, Igrad, I.depth(), 1, 1);
	Laplacian(I, Ilaplace, I.depth());
	namedWindow("image1", 1);
	imshow("image1", Ilaplace);
	waitKey(0);*/
	//system("pause");
	return 0;
}

