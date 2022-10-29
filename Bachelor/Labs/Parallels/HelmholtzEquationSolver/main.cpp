#define _USE_MATH_DEFINES
#include <math.h>

#include <string>
#include <vector>
#include "Matrix.h"
#include <iostream>
#include <fstream>
#include <functional>
#include <mpi.h>

typedef double(*Function)(double x, double y, double k);
typedef std::function<Matrix<double>(double, Matrix<double>, int, int, int, int, int, double, int&)> MPI_Type;

#define max_iter 500

void OutputFile(const std::string& filePath, const Matrix<double>& A, double h, size_t step = 1)
{
	std::ofstream file(filePath);

	if (!file.is_open())
	{
		std::cout << "OutputFile : File wasn't open" << std::endl;
		return;
	}

	for (size_t i = 0; i < A.GetSize().first; i += step)
	{
		for (size_t j = 0; j < A.GetSize().second; j += step)
		{
			file << i * h << "\t" << j * h << "\t" << A[i][j] << std::endl;
		}
	}

	file.close();
}

double Error(const Matrix<double>& y, const Matrix<double>& yPrev) 
{
	int N = y.GetSize().first, M = y.GetSize().second;

	int i = N / 2., j = M / 2.;

	return fabs(y[i][j] - yPrev[i][j]);
}
double Error(const Matrix<double>& y) 
{
	int N = y.GetSize().first, M = y.GetSize().second;
	double h = 1. / (N - 1);

	auto u = [](double x, double y) 
	{
		return (1. - x) * x * sin(M_PI * y);
	};

	int i = N / 2., j = M / 2.;

	return fabs(y[i][j] - u(j * h, i * h));
}

void HelmholtzSolver(MPI_Type method, double k, Function f, int N, int myid, int np, double eps, const std::string& methodName)
{
	std::vector<int> nPerProc, dataDisp;
	double h = 1. / (N - 1), t = 0, c = 1. / (4. + pow(k * h, 2));
	int NLocal, dataDispLocal, iter = 0;
	Matrix<double> y;

	if (myid == 0)
	{
 
		y.Resize(N);
		nPerProc.resize(np, N / np);
		dataDisp.resize(np);
    
		for (size_t i = 0; i < N % np; ++i)
		{
			++nPerProc[i];
		}
		for (size_t i = 1; i < np; ++i)
		{
			dataDisp[i] = dataDisp[i - 1] + nPerProc[i - 1];
		}
	}

	MPI_Scatter(nPerProc.data(), 1, MPI_INT, &NLocal, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(dataDisp.data(), 1, MPI_INT, &dataDispLocal, 1, MPI_INT, 0, MPI_COMM_WORLD);

	Matrix<double> fMtx(NLocal, N);

	for (size_t i = 0; i < NLocal; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			fMtx[i][j] = (h * h * c) * f((dataDispLocal + i) * h, j * h, k);
		}
	}
   //MPI_Barrier(MPI_COMM_WORLD);

	Matrix<double> localY = method(k, fMtx, N, NLocal, dataDispLocal, myid, np, eps, iter);
   //MPI_Barrier(MPI_COMM_WORLD);
  t = -MPI_Wtime();
	/*if (myid == 0)
	{
		for (size_t i = 0; i < NLocal; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				y[i][j] = localY[i][j];
			}
		}
		for (size_t i = 1; i < np; ++i)
		{
			//for (size_t j = 0; j < nPerProc[i]; ++j)
			{
				MPI_Recv(y[dataDisp[i]], nPerProc[i]*N, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			}
		}
	}
	else
	{
		//for (size_t i = 0; i < NLocal; ++i)
		{
			MPI_Send(localY[0], N*NLocal, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		}
	}
 */
 
 MPI_Gather(localY[0], N*NLocal, MPI_DOUBLE, y[0], N*NLocal, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   t += MPI_Wtime();
	if (myid == 0)
	{
				std::cout << methodName << " : time = " << t << std::endl;

		OutputFile("OP" + methodName + ".txt", y, h, 32);
	}
}

Matrix<double> JacobiHelmholtzSolver(double k, Matrix<double> fMtx, int N, int NLocal, int dataDispLocal, int myid, int np, double eps, int& iter)
{
	double h = 1. / (N - 1), c = 1. / (4. + pow(k * h, 2)), err = 100;

	Matrix<double> localY(NLocal, N), localPrevY(NLocal, N);
	std::vector<double> localYTop(N), localYLow(N);

	iter = 0;
	do
	{
		++iter;

		for (int id = 0; id < np - 1; id++)
		{
			if (myid == id)
			{
				MPI_Send(localPrevY[NLocal - 1], N, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD);
			}
			if (myid == id + 1)
			{
				MPI_Recv(localYLow.data(), N, MPI_DOUBLE, myid - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			}
		}
		for (int id = np - 1; id > 0; id--)
		{
			if (myid == id)
			{
				MPI_Send(localPrevY[0], N, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD);
			}
			if (myid == id - 1)
			{
				MPI_Recv(localYTop.data(), N, MPI_DOUBLE, myid + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			}
		}


		for (size_t i = 1; i < NLocal - 1; ++i)
		{
			for (size_t j = 1; j < N - 1; ++j)
			{
				localY[i][j] = (localPrevY[i + 1][j] + localPrevY[i - 1][j] + 
								localPrevY[i][j + 1] + localPrevY[i][j - 1]) * c + fMtx[i][j];
			}
		}
		if (myid != 0)
		{
			for (size_t j = 1; j < N - 1; ++j)
			{
				localY[0][j] = (localPrevY[1][j] + localYLow[j] + 
								localPrevY[0][j + 1] + localPrevY[0][j - 1]) * c + fMtx[0][j];
			}
		}
		if (myid != np - 1)
		{
			for (size_t j = 1; j < N - 1; ++j)
			{
				localY[NLocal - 1][j] = (localYTop[j] + localPrevY[NLocal - 1 - 1][j] + 
										 localPrevY[NLocal - 1][j + 1] + localPrevY[NLocal - 1][j - 1]) * c + fMtx[NLocal - 1][j];
			}
		}

		std::swap(localPrevY, localY);

		double localMax = Error(localPrevY, localY);							// |y_(k+1) - y_(k)|
		//double localMax = Error(localPrevY);									// |y_(k+1) - u(x_(k+1), y_(k+1))|

		MPI_Reduce(&localMax, &err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Bcast(&err, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	} 
	//while (iter < max_iter);
	while (err > eps);
 if (myid==0)
 {
   std::cout << "iter=" << iter << std::endl;
 }

	return localPrevY;
}

Matrix<double> SRJacobiHelmholtzSolver(double k, Matrix<double> fMtx, int N, int NLocal, int dataDispLocal, int myid, int np, double eps, int& iter)
{
	double h = 1. / (N - 1), c = 1. / (4. + pow(k * h, 2)), err = 100;
	int dest = 0, source = 0;

	Matrix<double> localY(NLocal, N), localPrevY(NLocal, N);
	std::vector<double> localYTop(N), localYLow(N);

	if (np > 1)
	{
		dest = myid + 1, source = myid - 1;
		if (myid == 0)
		{
			source = np - 1;
		}
		else if (myid == np - 1)
		{
			dest = 0;
		}
	}
	int sendCount = (myid - (np - 1)) ? N : 0, recvCount = myid ? N : 0;

	iter = 0;
	do 
	{
		++iter;
    if(np > 1)
    {
  		if (myid == 0)
  		{
  			MPI_Sendrecv(localPrevY[NLocal - 1], sendCount, MPI_DOUBLE, myid + 1, 0, localYLow.data(), recvCount, MPI_DOUBLE, np - 1, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
  			MPI_Sendrecv(localPrevY[0], recvCount, MPI_DOUBLE, np - 1, 1, localYTop.data(), sendCount, MPI_DOUBLE, myid + 1, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
  		}
  		else if (myid == np - 1)
  		{
  			MPI_Sendrecv(localPrevY[NLocal - 1], sendCount, MPI_DOUBLE, 0, 0, localYLow.data(), recvCount, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
  			MPI_Sendrecv(localPrevY[0], recvCount, MPI_DOUBLE, myid - 1, 1, localYTop.data(), sendCount, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
  		}
  		else
  		{
  			MPI_Sendrecv(localPrevY[NLocal - 1], sendCount, MPI_DOUBLE, myid + 1, 0, localYLow.data(), recvCount, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
  			MPI_Sendrecv(localPrevY[0], recvCount, MPI_DOUBLE, myid - 1, 1, localYTop.data(), sendCount, MPI_DOUBLE, myid + 1, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
  		}
   }

		for (size_t i = 1; i < NLocal - 1; ++i)
		{
			for (size_t j = 1; j < N - 1; ++j)
			{
				localY[i][j] = (localPrevY[i + 1][j] + localPrevY[i - 1][j] + 
								localPrevY[i][j + 1] + localPrevY[i][j - 1]) * c + fMtx[i][j];
			}
		}

		if (myid != 0)
		{
			for (size_t j = 1; j < N - 1; ++j)
			{
				localY[0][j] = (localPrevY[1][j] + localYLow[j] + 
								localPrevY[0][j + 1] + localPrevY[0][j - 1]) * c + fMtx[0][j];
			}
		}
		if (myid != np - 1) 
		{
			for (size_t j = 1; j < N - 1; ++j)
			{
				localY[NLocal - 1][j] = (localYTop[j] + localPrevY[NLocal - 2][j] + 
										 localPrevY[NLocal - 1][j + 1] + localPrevY[NLocal - 1][j - 1]) * c + fMtx[NLocal - 1][j];
			}
		}

		std::swap(localPrevY, localY);

		//double localMax = Error(localPrevY, localY);							// |y_(k+1) - y_(k)|
		//double localMax = Error(localPrevY);									// |y_(k+1) - u(x_(k+1), y_(k+1))|

		//MPI_Reduce(&localMax, &err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		//MPI_Bcast(&err, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	while (iter < max_iter);
	//while (err > eps);

	return localY;
}

Matrix<double> IJacobiHelmholtzSolver(double k, Matrix<double> fMtx, int N, int NLocal, int dataDispLocal, int myid, int np, double eps, int& iter)
{
	double h = 1. / (N - 1), c = 1. / (4. + pow(k * h, 2)), err = 100;
	Matrix<double> localY(NLocal, N), localPrevY(NLocal, N);
	std::vector<double> localYTop(N), localYLow(N);
	MPI_Request* reqSend, * reqRecv;
	int nReqSend = 0, nReqRecv = 0;

	if (myid == 0 || myid == np - 1)
	{
		reqSend = new MPI_Request[1], reqRecv = new MPI_Request[1];
	}
	else
	{
		reqSend = new MPI_Request[2], reqRecv = new MPI_Request[2];
	}

	iter = 0;
	do
	{
		iter++;
		nReqSend = 0, nReqRecv = 0;

		if (myid != 0)
		{
			MPI_Isend(localY[0], N, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD, reqSend + nReqSend);
			nReqSend++;

			MPI_Irecv(localYLow.data(), N, MPI_DOUBLE, myid - 1, MPI_ANY_TAG, MPI_COMM_WORLD, reqRecv + nReqRecv);
			nReqRecv++;
		}
		if (myid != np - 1)
		{
			MPI_Isend(localY[NLocal - 1], N, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD, reqSend + nReqSend);
			nReqSend++;

			MPI_Irecv(localYTop.data(), N, MPI_DOUBLE, myid + 1, MPI_ANY_TAG, MPI_COMM_WORLD, reqRecv + nReqRecv);
			nReqRecv++;
		}

		for (int i = 1; i < NLocal - 1; ++i)
		{
			for (int j = 1; j < N - 1; ++j)
			{
				localY[i][j] = (localPrevY[i + 1][j] + localPrevY[i - 1][j] +
								localPrevY[i][j + 1] + localPrevY[i][j - 1]) * c + fMtx[i][j];
			}
		}

		MPI_Waitall(nReqRecv, reqRecv, MPI_STATUSES_IGNORE);

		if (myid != 0)
		{
			for (int j = 1; j < N - 1; ++j)
			{
				localY[0][j] = (localPrevY[1][j] + localYLow[j] +
								localPrevY[0][j + 1] + localPrevY[0][j - 1]) * c + fMtx[0][j];
			}
		}
		if (myid != np - 1)
		{
			for (int j = 1; j < N - 1; ++j)
			{
				localY[NLocal - 1][j] = (localYTop[j] + localPrevY[NLocal - 2][j] +
										 localPrevY[NLocal - 1][j + 1] + localPrevY[NLocal - 1][j - 1]) * c + fMtx[NLocal - 1][j];
			}
		}

		MPI_Waitall(nReqSend, reqSend, MPI_STATUSES_IGNORE);

		std::swap(localPrevY, localY);

		//double localMax = Error(localPrevY, localY);							// |y_(k+1) - y_(k)|
		//double localMax = Error(localPrevY);									// |y_(k+1) - u(x_(k+1), y_(k+1))|

		//MPI_Reduce(&localMax, &err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		//MPI_Bcast(&err, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	} 
	while (iter < max_iter);
	//while (err > eps);

	delete[] reqSend;
	delete[] reqRecv;

	return localY;
}

Matrix<double> SeidelHelmholtzSolver(double k, Matrix<double> fMtx, int N, int NLocal, int dataDispLocal, int myid, int np, double eps, int& iter)
{
	double h = 1. / (N - 1), c = 1. / (4. + pow(k * h, 2)), localPrev, err = 100;

	Matrix<double> localY(NLocal, N);
	std::vector<double> localYTop(N), localYLow(N);

	iter = 0;
	do
	{
		++iter;
		localPrev = localY[int(NLocal / 2.)][int(N / 2.)];

		for (int id = 0; id < np - 1; id++)
		{
			if (myid == id)
			{
				MPI_Send(localY[NLocal - 1], N, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD);
			}
			if (myid == id + 1)
			{
				MPI_Recv(localYLow.data(), N, MPI_DOUBLE, myid - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			}
		}
		for (int id = np - 1; id > 0; id--)
		{
			if (myid == id)
			{
				MPI_Send(localY[0], N, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD);
			}
			if (myid == id - 1)
			{
				MPI_Recv(localYTop.data(), N, MPI_DOUBLE, myid + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			}
		}

		for (size_t i = 1; i < NLocal - 1; ++i)
		{
			for (size_t j = (i % 2 + 1); j < N - 1; j += 2)
			{
				localY[i][j] = (localY[i + 1][j] + localY[i - 1][j] +
								localY[i][j + 1] + localY[i][j - 1]) * c + fMtx[i][j];
			}
		}
		for (size_t i = 1; i < NLocal - 1; ++i)
		{
			for (size_t j = ((i + 1) % 2 + 1); j < N - 1; j += 2)
			{
				localY[i][j] = (localY[i + 1][j] + localY[i - 1][j] + 
								localY[i][j + 1] + localY[i][j - 1]) * c + fMtx[i][j];
			}
		}

		if (myid != 0)
		{
			for (size_t j = 1; j < N - 1; j += 2)
			{
				localY[0][j] = (localY[1][j] + localYLow[j] + 
								localY[0][j + 1] + localY[0][j - 1]) * c + fMtx[0][j];
			}
			for (size_t j = 2; j < N - 1; j += 2)
			{
				localY[0][j] = (localY[1][j] + localYLow[j] + 
								localY[0][j + 1] + localY[0][j - 1]) * c + fMtx[0][j];
			}
		}
		if (myid != np - 1) 
		{
			for (size_t j = ((NLocal - 1) % 2 + 1); j < N - 1; j += 2)
			{
				localY[NLocal - 1][j] = (localYTop[j] + localY[NLocal - 2][j] + 
										 localY[NLocal - 1][j + 1] + localY[NLocal - 1][j - 1]) * c + fMtx[NLocal - 1][j];
			}
			for (size_t j = (NLocal % 2 + 1); j < N - 1; j += 2)
			{
				localY[NLocal - 1][j] = (localYTop[j] + localY[NLocal - 2][j] + 
										 localY[NLocal - 1][j + 1] + localY[NLocal - 1][j - 1]) * c + fMtx[NLocal - 1][j];
			}
		}

		double localMax = fabs(localPrev - localY[int(NLocal / 2.)][int(N / 2.)]);		// |y_(k+1) - y_(k)|					
		//double localMax = Error(localY);													// |y_(k+1) - u(x_(k+1), y_(k+1))|

		MPI_Reduce(&localMax, &err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Bcast(&err, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	} 
	//while (iter < max_iter);
	while (err > eps);
   if (myid==0)
 {
   std::cout << "iter=" << iter << std::endl;
 }
	return localY;
}

Matrix<double> SRSeidelHelmholtzSolver(double k, Matrix<double> fMtx, int N, int NLocal, int dataDispLocal, int myid, int np, double eps, int& iter)
{
	double h = 1. / (N - 1), c = 1. / (4. + pow(k * h, 2)), localPrev, err = 100;
	int dest = 0, source = 0;
	Matrix<double> localY(NLocal, N);
	std::vector<double> localYTop(N), localYLow(N);

	if (np > 1)
	{
		dest = myid + 1, source = myid - 1;
		if (myid == 0)
		{
			source = np - 1;
		}
		else if (myid == np - 1)
		{
			dest = 0;
		}
	}
	int sendCount = (myid - (np - 1)) ? N : 0, recvCount = myid ? N : 0;

	iter = 0;
	do 
	{
		iter++;
		localPrev = localY[int(NLocal / 2.)][int(N / 2.)];
    if(np > 1)
    {
      //std::cout<<"write smth"<<std::endl;
  		if (myid == 0)
  		{
  			MPI_Sendrecv(localY[NLocal - 1], sendCount, MPI_DOUBLE, myid + 1, 0, localYLow.data(), recvCount, MPI_DOUBLE, np - 1, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
  			MPI_Sendrecv(localY[0], recvCount, MPI_DOUBLE, np - 1, 1, localYTop.data(), sendCount, MPI_DOUBLE, myid + 1, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
  		}
  		else if (myid == np - 1)
  		{
  			MPI_Sendrecv(localY[NLocal - 1], sendCount, MPI_DOUBLE, 0, 0, localYLow.data(), recvCount, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
  			MPI_Sendrecv(localY[0], recvCount, MPI_DOUBLE, myid - 1, 1, localYTop.data(), sendCount, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
  		}
  		else
  		{
  			MPI_Sendrecv(localY[NLocal - 1], sendCount, MPI_DOUBLE, myid + 1, 0, localYLow.data(), recvCount, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
  			MPI_Sendrecv(localY[0], recvCount, MPI_DOUBLE, myid - 1, 1, localYTop.data(), sendCount, MPI_DOUBLE, myid + 1, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
  		}
    }
		for (size_t i = 1; i < NLocal - 1; ++i)
		{
			for (size_t j = (i % 2 + 1); j < N - 1; j += 2)
			{
				localY[i][j] = (localY[i + 1][j] + localY[i - 1][j] +
								localY[i][j + 1] + localY[i][j - 1]) * c + fMtx[i][j];
			}
		}
		for (size_t i = 1; i < NLocal - 1; ++i)
		{
			for (size_t j = ((i + 1) % 2 + 1); j < N - 1; j += 2)
			{
				localY[i][j] = (localY[i + 1][j] + localY[i - 1][j] +
								localY[i][j + 1] + localY[i][j - 1]) * c + fMtx[i][j];
			}
		}

		if (myid != 0)
		{
			for (size_t j = 1; j < N - 1; j += 2)
			{
				localY[0][j] = (localY[1][j] + localYLow[j] + localY[0][j + 1] + 
								localY[0][j - 1]) * c + fMtx[0][j];
			}
			for (size_t j = 2; j < N - 1; j += 2)
			{
				localY[0][j] = (localY[1][j] + localYLow[j] + 
								localY[0][j + 1] + localY[0][j - 1]) * c + fMtx[0][j];
			}
		}
		if (myid != np - 1)
		{
			for (size_t j = ((NLocal - 1) % 2 + 1); j < N - 1; j += 2)
			{
				localY[NLocal - 1][j] = (localYTop[j] + localY[NLocal - 2][j] + 
										 localY[NLocal - 1][j + 1] + localY[NLocal - 1][j - 1]) * c + fMtx[NLocal - 1][j];
			}
			for (size_t j = (NLocal % 2 + 1); j < N - 1; j += 2)
			{
				localY[NLocal - 1][j] = (localYTop[j] + localY[NLocal - 2][j] + 
										 localY[NLocal - 1][j + 1] + localY[NLocal - 1][j - 1]) * c + fMtx[NLocal - 1][j];
			}
		}

		//double localMax = fabs(localPrev - localY[int(NLocal / 2.)][int(N / 2.)]);		// |y_(k+1) - y_(k)|					
		//double localMax = Error(localY);													// |y_(k+1) - u(x_(k+1), y_(k+1))|

		//MPI_Reduce(&localMax, &err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		//MPI_Bcast(&err, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	} 
	while (iter < max_iter);
	//while (err > eps);

	return localY;
}

Matrix<double> ISeidelHelmholtzSolver(double k, Matrix<double> fMtx, int N, int NLocal, int dataDispLocal, int myid, int np, double eps, int& iter)
{
	double h = 1. / (N - 1), c = 1. / (4. + pow(k * h, 2)), localPrev, err = 100;
	Matrix<double> localY(NLocal, N);
	std::vector<double> localYTop(N), localYLow(N);
	MPI_Request* reqSend, * reqRecv;
	int nReqSend = 0, nReqRecv = 0;

	if (myid == 0 || myid == np - 1)
	{
		reqSend = new MPI_Request[1], reqRecv = new MPI_Request[1];
	}
	else
	{
		reqSend = new MPI_Request[2], reqRecv = new MPI_Request[2];
	}

	iter = 0;
	do 
	{
		iter++;
		nReqSend = 0, nReqRecv = 0;
		localPrev = localY[int(NLocal / 2.)][int(N / 2.)];

		if (myid != 0)
		{
			MPI_Isend(localY[0], N, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD, reqSend + nReqSend);
			nReqSend++;

			MPI_Irecv(localYLow.data(), N, MPI_DOUBLE, myid - 1, MPI_ANY_TAG, MPI_COMM_WORLD, reqRecv + nReqRecv);
			nReqRecv++;
		}
		if (myid != np - 1)
		{
			MPI_Isend(localY[NLocal - 1], N, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD, reqSend + nReqSend);
			nReqSend++;

			MPI_Irecv(localYTop.data(), N, MPI_DOUBLE, myid + 1, MPI_ANY_TAG, MPI_COMM_WORLD, reqRecv + nReqRecv);
			nReqRecv++;
		}

		for (size_t i = 1; i < NLocal - 1; ++i)
		{
			for (size_t j = (i % 2 + 1); j < N - 1; j += 2)
			{
				localY[i][j] = (localY[i + 1][j] + localY[i - 1][j] +
								localY[i][j + 1] + localY[i][j - 1]) * c + fMtx[i][j];
			}
		}
		for (size_t i = 1; i < NLocal - 1; ++i)
		{
			for (size_t j = ((i + 1) % 2 + 1); j < N - 1; j += 2)
			{
				localY[i][j] = (localY[i + 1][j] + localY[i - 1][j] +
								localY[i][j + 1] + localY[i][j - 1]) * c + fMtx[i][j];
			}
		}

		MPI_Waitall(nReqRecv, reqRecv, MPI_STATUSES_IGNORE);

		if (myid != 0)
		{
			for (size_t j = 1; j < N - 1; j += 2)
			{
				localY[0][j] = (localY[1][j] + localYLow[j] + 
								localY[0][j + 1] + localY[0][j - 1]) * c + fMtx[0][j];
			}
			for (size_t j = 2; j < N - 1; j += 2)
			{
				localY[0][j] = (localY[1][j] + localYLow[j] + 
								localY[0][j + 1] + localY[0][j - 1]) * c + fMtx[0][j];
			}
		}
		if (myid != np - 1)
		{
			for (size_t j = ((NLocal - 1) % 2 + 1); j < N - 1; j += 2)
			{
				localY[NLocal - 1][j] = (localYTop[j] + localY[NLocal - 2][j] + 
										 localY[NLocal - 1][j + 1] + localY[NLocal - 1][j - 1]) * c + fMtx[NLocal - 1][j];
			}
			for (size_t j = (NLocal % 2 + 1); j < N - 1; j += 2)
			{
				localY[NLocal - 1][j] = (localYTop[j] + localY[NLocal - 2][j] + 
										 localY[NLocal - 1][j + 1] + localY[NLocal - 1][j - 1]) * c + fMtx[NLocal - 1][j];
			}
		}

		MPI_Waitall(nReqSend, reqSend, MPI_STATUSES_IGNORE);

		double localMax = fabs(localPrev - localY[int(NLocal / 2.)][int(N / 2.)]);		// |y_(k+1) - y_(k)|					
		//double localMax = Error(localY);													// |y_(k+1) - u(x_(k+1), y_(k+1))|

		MPI_Reduce(&localMax, &err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Bcast(&err, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	//while (iter < max_iter);
	while (err > eps);
  if(myid==0)
  {
    std::cout << "iter=" << iter << std::endl;
  }
	delete[] reqSend;
	delete[] reqRecv;

	return localY;
}

int main(int argc, char* argv[])
{
	int N = 10260, np, myid;
	double k = 150*100, eps = 1e-4;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	auto f = [](double x, double y, double k)
	{
		return 2. * sin(M_PI * y) + k * k * (1. - x) * x * sin(M_PI * y) + M_PI * M_PI * (1. - x) * x * sin(M_PI * y);
	};

	switch (argv[1][0] - '0')
	{
	case 1:
		HelmholtzSolver(JacobiHelmholtzSolver, k, f, N, myid, np, eps, "JacobiHelmholtzSolver");
		break;
	case 2:		
		HelmholtzSolver(SRJacobiHelmholtzSolver, k, f, N, myid, np, eps, "SRJacobiHelmholtzSolver");
		break;
	case 3:
		HelmholtzSolver(IJacobiHelmholtzSolver, k, f, N, myid, np, eps, "IJacobiHelmholtzSolver");
		break;
	case 4:
		HelmholtzSolver(SeidelHelmholtzSolver, k, f, N, myid, np, eps, "SeidelHelmholtzSolver");
		break;
	case 5:
		HelmholtzSolver(SRSeidelHelmholtzSolver, k, f, N, myid, np, eps, "SRSeidelHelmholtzSolver");
		break;
	case 6:
		HelmholtzSolver(ISeidelHelmholtzSolver, k, f, N, myid, np, eps, "ISeidelHelmholtzSolver");
		break;
	}

	MPI_Finalize();

	return 0;
}