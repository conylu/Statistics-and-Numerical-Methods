#include <iostream>
#include <fstream>
#include <math.h>
//#include <Eigen/Eigenvalues>

//using namespace Eigen;
using namespace std;
void r_calculate(double* r, double** a, double* x, double* b, int m);
void boundary(double* v, int n, int m);
void boundary_zero(double* v, int n, int m);
void jacobi(double** a, double* b, double* x0, double* x1, int m);
void mat_produce(double* v,double** A,int n);//把n*n行的向量的矩阵转换成n*n的矩阵
double norm2(double* v,int n);//计算r的2范数


void r_calculate(double* r, double** a, double* x, double* b, int m){
	for (int i = 0; i < m; i++) {	//初始化r[i]
		r[i] = b[i];
		for (int j = 0; j < m; j++) {
			r[i] -= a[i][j] * x[j];
		}
	}
}

void boundary(double* v, int n, int m)//考虑边界条件
{
	for (int i = 0; i < m; i++) {
		if (i / n == 0)
			v[i] = -v[i + n];
		else if (i / n == n - 1)
			v[i] = -v[i - n];
		else if (i % n == 0)
			v[i] = -v[i + 1];
		else if (i % n == n - 1)
			v[i] = -v[i - 1];
	}
}

void boundary_zero(double* v, int n, int m)//将边界上的r设为0，否则影响2范数的计算
{
	for (int i = 0; i < m; i++) {
		if (i / n == 0 || i / n == n - 1 || i % n == 0 || i % n == n - 1) {
			v[i] = 0;
		}
	}
}

void jacobi(double** a, double* b, double* x0, double* x1, int m){
	for (int i = 0; i < m; i++) {
		x1[i] = 0;
		for (int j = 0; j < i ; j++) {
			x1[i] -= a[i][j]*x0[j];
		}
		for (int j = i + 1; j < m; j++) {
			x1[i] -= a[i][j]*x0[j];
		}
		x1[i] =(x1[i]+b[i])/a[i][i];
	}
}

void mat_produce(double* v, double** A,int n) {
	for (int i = 1; i < n*n; i++) {
		A[i / n][i % n] = v[i];
	}
}

double norm2(double* v, int m)
{
	double sum = 0.0;
	for (int i = 0; i < m; i++) {
		sum += pow(v[i], 2.0);
	}
	return pow(sum, 0.5);
}

int main() {
	/*
	总体思路
	1. 设置任意的x0,计算r0=b-a*x0
	2. 调用不同的迭代函数，使之||r0||2收敛到所需精度。
	*/

	/*参数设置*/
	int n0 = 10;	//设置格点大小
	int ncount = 10000;//循环次数
	double precision = pow(10, -6);

	/*后面常用的量*/
	int n = n0 + 2;
	int m = pow(n0 + 2, 2);

	/*分配各个变量的储存空间*/
	double* phi0 = new double[m];//储存初始解or旧解
	double* b = new double[m];//储存A*phi=b中的b,在本例中为Rho
	double** a = new double*[m];//储存矩阵A
	for (int i = 0; i < m; i++) {
		a[i] = new double[m];
	}                
	double* r = new double[m];//储存剩余矢量
	double* phi1 = new double[m];//储存新解
	double* rnorm2 = new double [ncount];//循环后的2范数
	double** phimat = new double* [n];//解的矩阵
	for (int i = 0; i < n; i++) {
		phimat[i] = new double[n];
	}
	double** rmat = new double* [n];//剩余矢量矩阵
	for (int i = 0; i < n; i++) {
		rmat[i] = new double[n];
	}

	/*初始化各个变量*/

    //初始phi0，全猜1
	for (int i = 0; i < m; i++) {
		phi0[i] = 0;	
	}
	boundary(phi0, n, m);

    //初始化b[i]。边界上的b[i]并不会影响后续的计算，无需特别设置
	for (int i = 0; i < m; i++) {	
		//int xi = -10 + ((i / n )- 0.5);
		//int yi = -10 + ((i % n ) -0.5);
		b[i] = exp(-(pow(-(n-2)/2 + ((i / n ) - 0.5), 2) + pow(-(n-2)/2 + ((i % n) - 0.5), 2)));
	}

    //对a[][]的初始化为先把每个元素都设为0，再在下一个循环中设置非0元素
	for (int i = 0; i < m; i++) {	
		for (int j = 0; j < m; j++) {
			a[i][j] = 0;
		}
	}
	for (int i = n; i < m - n; i++) {//注意溢出
		a[i][i] = -4;
		a[i][i - 1] = 1;
		a[i][i + 1] = 1;
		a[i][i - n] = 1;
		a[i][i + n] = 1;
	}

	//计算最初的r,及r的二范数
	r_calculate(r, a, phi0, b, m);
	boundary_zero(r, n, m);
	rnorm2[0] = norm2(r, m);
	cout << "R=" << rnorm2[0] << '\n';
	//fout << "R=" << rnorm2[0] << '\n';


	//用不同的方法计算并输出到文件
	ofstream fout;
	fout.open("norm2.txt");
	for (int count = 0; count < ncount; count++) {
		jacobi(a,b,phi0, phi1, m);
		boundary(phi1, n, m);
		mat_produce(phi1, phimat, n);
		//for (int i = 0; i < n; i++) {
		//	for (int j = 0; j < n; j++) {
		//		cout << phimat[i][j] << '\t';
		//	}
		//	cout << endl;
		//}

		r_calculate(r, a, phi1, b, m);
		boundary_zero(r, n, m);
		//mat_produce(r, rmat, n);
		rnorm2[count+1] = norm2(r, m);
		cout << "R=" << rnorm2[count + 1] << '\n';
		fout << "R=" << rnorm2[count + 1] << '\n';
		for (int i = 0; i < m; i++) {
			phi0[i] = phi1[i];
		}
	}
	fout.close();
	for (int count = 0; count < ncount; count++) {
		if (rnorm2[count] < precision) {
		cout << "所需的收敛次数为" << count << endl;
		break;
		}
	}

	
	


}


