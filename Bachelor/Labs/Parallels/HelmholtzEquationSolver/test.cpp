//решаем уравнение Гельмгольца вида -d^2(u)/dx^2-d^2(u)/dy^2+k^2*u=f с нулевыми ГУ

#define PI 3.14159265358979323846264338328

#define JACOBI
//#define SEIDEL


//#define SEND
#define SENDRECV
//#define ISEND

#ifdef ISEND
  #define INIT_ISEND
#endif

#include <math.h>

#include <string>
#include <vector>

#include<iostream>
#include <fstream>
#include<mpi.h>

using namespace std;

typedef double Type;  

Type f(const Type x, const Type y, const Type k)
{
	return 2. * sin(PI * y) + k * k * (1. - x) * x * sin(PI * y) + PI * PI * (1. - x) * x * sin(PI * y);
}
Type u(const Type x, const Type y)
{
	return (1. - x) * x * sin(PI * y);
}


struct Matrix {
	int size_i;
	std::vector<Type> data;

	Matrix(const int size=0, const Type val=0) {
		size_i = size;
		data.resize(size * size, val);
	}

	Matrix(const int size_x, const int size_y, const Type val = 0) {
		size_i = size_y;
		data.resize(size_x * size_y, val);
	}

	void Resize(const int size , const Type val = 0) {
		size_i = size;
		data.resize(size * size, val);
	}

	void Print(const int N, const int M) {
		for (size_t i = 0; i < N; i++)
		{
			for (size_t j = 0; j < M; j++)
			{
				std::cout << (*this)(i, j) << ' ';
			}
			cout << '\n';
		}
	}

	Type& operator()(const int i, const int j) {
		return data[size_i * i + j];
	}

	Type* operator()(const int i) {
		return data.data() + size_i * i;
	}

};

int WriteFile(const std::string& name_file, const int step_draw, const Type h, Matrix& A)
{
	std::ofstream ofile;
	ofile.open(name_file);

	const int N = A.size_i;

	for (int i = 0; i < N; i += step_draw)
		for (int j = 0; j < N; j += step_draw)
			ofile << i * h << ' ' << j * h << ' ' << A(i, j) << '\n';

	ofile.close();
	return 0;
}

Type NormError(const Type h, Matrix& U) {

	const int n = U.size_i;

	Type max = 0;
	Type buf = -1;

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j) {
			buf = fabs(U(i, j) - u(i * h, j * h));
			if (buf > max)
				max = buf;
		}
	return max;
}

int main(int argc, char** argv)
{
	//const Type eps = 1e-5;
	const int max_iter = 500;
	const int n = 2048;				   //кол-во узлов
		
	
	Type h = 1. / (n-1);		       // шаг
	Type k = 3 / (20 * h);
	Type c = 1. / (k * k * h * h + 4.);  

	int np, myid;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);  //номер исполняемового процесса
	MPI_Comm_size(MPI_COMM_WORLD, &np);    //кол-во исполняемых процессов

	MPI_Status stat;

	std::vector<int> nLoc(0);   //кол-во разбиений матрицы U для каждого процесса
	int n_loc;                  //Кол-во строк марицы U_loc для текущего процесса

	Matrix U;		        	//Численное решение, найденное с заданной точностью

	std::vector<int> shift(0);  // сдвиги начала области для i-го процессора
	int shift_loc;              // сдвиг для конкртеного процессора


	Type time = -MPI_Wtime();

	// Выделение памяти
	if (myid == 0)
	{
		U.Resize(n, 0);

		nLoc.resize(np, n / np);           //кол-во строк под каждый процесс	
		for (int i = 0; i < n % np; i++)  
			++nLoc[i];   // если число процессоров не кратно числу узлов (добавить 1 в первые n%np процессоров 

		shift.resize(np, 0);
		for (int i = 1; i < np; i++) 
			shift[i] = shift[i - 1] + nLoc[i - 1];

	}

	
	//отправляем каждому процессу свое значение nLoc, записывая в переменную n_loc
	MPI_Scatter(nLoc.data(), 1, MPI_INT, &n_loc, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(shift.data(), 1, MPI_INT, &shift_loc, 1, MPI_INT, 0, MPI_COMM_WORLD);



	Matrix yk(n_loc, n, 0);     // локальное решение на предыудщем шаге
	Matrix yk1(n_loc, n, 0);    // локальное решение 
	Matrix right(n_loc, n, 0);  // правая часть

	for (int i = 0; i < n_loc; i++)
		for (int j = 0; j < n; j++)
			right(i, j) = (h * h * c) * f((shift_loc + i) * h, j * h, k);


	std::vector<Type> y_up(n, 0);    // решение от процессора "сверху" 
	std::vector<Type> y_down(n, 0);  // решение от процессора "снизу"


#ifdef ISEND
	MPI_Request* rq_send;  //запрос на отправку
	MPI_Request* rq_recv;  //запрос на прием

	int n_rq_send = 0;    // номер запроса отпраки
	int n_rq_recv = 0;    // номер запроса приема
	
	MPI_Status* st;		 // статус перессылки

	if (myid == 0 || myid == np - 1) // первый и последний делает одно сообщение (остальные 2)
	{
		rq_send = new MPI_Request;
		rq_recv = new MPI_Request;
		st = new MPI_Status;
	}
	else {
		rq_send = new MPI_Request[2];
		rq_recv = new MPI_Request[2];
		st = new MPI_Status[2];
	}
#else
	MPI_Status st;  // для блокирующих пересылок запросы не нужны
#endif

#ifdef INIT_ISEND

	if (myid < np - 1)
	{
		
		MPI_Send_init(yk(n_loc-1), n, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD, rq_send + n_rq_send);
		MPI_Recv_init(y_up.data(), n, MPI_DOUBLE, myid + 1, MPI_ANY_TAG, MPI_COMM_WORLD, rq_recv + n_rq_recv);
		n_rq_send++;
		n_rq_recv++;
	}
	if (myid > 0)
	{
		MPI_Send_init(yk(0), n, MPI_DOUBLE, myid - 1, 1, MPI_COMM_WORLD, rq_send + n_rq_send);
		MPI_Recv_init(y_down.data(), n, MPI_DOUBLE, myid - 1, MPI_ANY_TAG, MPI_COMM_WORLD, rq_recv + n_rq_recv);
		n_rq_send++;
		n_rq_recv++;
	}
#endif // INIT_ISEND

#ifdef SENDRECV
	int dest_up;
	int dest_down;

	int sourse_up;
	int sourse_down;

	int dest_up_size;
	int dest_down_size;

	int sourse_up_size;
	int sourse_down_size;
	{
		if (np > 1) {
			if (myid == 0) {
				dest_up = 1;
				sourse_up = 1;

				dest_down = np - 1;
				sourse_down = np - 1;

				dest_up_size = n;
				dest_down_size = 0;

				sourse_up_size = n;
				sourse_down_size = 0;

			}
			else if (myid == np - 1) {
				dest_up = 0;
				sourse_up = 0;

				dest_down = myid - 1;
				sourse_down = myid - 1;

				dest_up_size = 0;
				dest_down_size = n;

				sourse_up_size = 0;
				sourse_down_size = n;
			}
			else {
				dest_up = myid + 1;
				sourse_up = myid + 1;

				dest_down = myid - 1;
				sourse_down = myid - 1;

				dest_up_size = n;
				dest_down_size = n;

				sourse_up_size = n;
				sourse_down_size = n;
			}
		}
	}

#endif // SENDRECV

#ifdef JACOBI
	int iter = 0;
	do {

#ifdef SEND
		{
			if (myid < np - 1)// кроме последнего
			{
				// отправка (принятие) правого конца	
				MPI_Send(yk(n_loc - 1), n, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD);
        MPI_Recv(y_up.data(), n, MPI_DOUBLE, myid + 1, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
				
			//	MPI_Send(yk(n_loc - 1), n, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD);     
				
				
			}
			if (myid > 0)  // кроме первого
			{
				// отправка принятие левого конца
				MPI_Recv(y_down.data(), n, MPI_DOUBLE, myid - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
				MPI_Send(yk(0), n, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD);
				

			}
		}
#endif // SEND

#ifdef SENDRECV
		{
			if (np > 1) {

				MPI_Sendrecv(yk(n_loc - 1), dest_up_size, MPI_DOUBLE, dest_up, 0,
					y_down.data(), sourse_down_size, MPI_DOUBLE, sourse_down, MPI_ANY_TAG, MPI_COMM_WORLD, &st);

				MPI_Sendrecv(yk(0), dest_down_size, MPI_DOUBLE, dest_down, 1,
					y_up.data(), sourse_up_size, MPI_DOUBLE, sourse_up, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
			}
		}
#endif // SENDRECV

#ifdef ISEND
#ifdef INIT_ISEND
		MPI_Startall(n_rq_send, rq_send);
		MPI_Startall(n_rq_recv, rq_recv);

#else
		n_rq_send = 0;
		n_rq_recv = 0;

		if (myid < np - 1)// кроме последнего
		{
			// создаем запрос rq_send под номером n_rq_send и увеличивыем счетчик отправки
			MPI_Isend(yk(n_loc - 1), n, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD, rq_send + n_rq_send);
			n_rq_send++;
		
			// создаем запрос rq_recv под номером n_rq_recv и увеличивыем счетчик принятия
			MPI_Irecv(y_up.data(), n, MPI_DOUBLE, myid + 1, MPI_ANY_TAG, MPI_COMM_WORLD, rq_recv + n_rq_recv);
			n_rq_recv++;
		}
		if (myid > 0)  // кроме первого
		{
			MPI_Isend(yk(0), n, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD, rq_send + n_rq_send);
			MPI_Irecv(y_down.data(), n, MPI_DOUBLE, myid - 1, MPI_ANY_TAG, MPI_COMM_WORLD, rq_recv + n_rq_recv);
			n_rq_send++;
			n_rq_recv++;
		}
#endif //INIT_ISEND
#endif // ISEND


		const int N = n_loc - 1;

		for (int i = 1; i < N; i++) {
			for (int j = 1; j < n - 1; j++) 
				yk1(i, j) = c * (yk(i + 1, j) + yk(i - 1, j) + yk(i, j + 1) + yk(i, j - 1)) + right(i, j);
		}

#ifdef ISEND
// для неблокирующей перессылки ждем, пока все сообщения будут приняты
		MPI_Waitall(n_rq_recv, rq_recv, st);
#endif 

		if (myid < np - 1) {
			for (int j = 1; j < n - 1; j++)
				yk1(N, j) = c * (y_up[j] + yk(N - 1, j) + yk(N, j + 1) + yk(N, j - 1)) + right(N, j);
		}
		if (myid > 0) {
			for (int j = 1; j < n - 1; j++)
				yk1(0, j) = c * (yk(1, j) + y_down[j] + yk(0, j + 1) + yk(0, j - 1)) + right(0, j);
		}

#ifdef ISEND
		// для неблокирующей перессылки ждем, пока все сообщения будут отправлены
		MPI_Waitall(n_rq_send, rq_send, st); 
#endif 

		std::swap(yk.data, yk1.data);
		iter++;

	} while (iter < max_iter);

#endif // JACOBI


#ifdef SEIDEL

	int iter = 0;
	do {


#ifdef SEND

		if (myid < np - 1)// кроме последнего
		{
			MPI_Send(yk(n_loc - 1), n, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD);            // отправить свой "верх"
			MPI_Recv(y_up.data(), n, MPI_DOUBLE, myid + 1, MPI_ANY_TAG, MPI_COMM_WORLD, &st);  // принять себе "верх"??
		}
		if (myid > 0)  // кроме первого
		{
			MPI_Recv(y_down.data(), n, MPI_DOUBLE, myid - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
			MPI_Send(yk(0), n, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD);
		}

#endif // SEND

#ifdef SENDRECV
		{
			if (myid == 0 && np > 1) {
				MPI_Sendrecv(yk(n_loc - 1), n, MPI_DOUBLE, myid + 1, 0,
					y_up.data(), n, MPI_DOUBLE, myid + 1, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
			}

			else if (myid == np - 1) {

				MPI_Sendrecv(yk(0), n, MPI_DOUBLE, myid - 1, 1,
					y_down.data(), n, MPI_DOUBLE, myid - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
			}
			else {

				MPI_Sendrecv(yk(n_loc - 1), n, MPI_DOUBLE, myid + 1, 0,
					y_down.data(), n, MPI_DOUBLE, myid - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &st);


				MPI_Sendrecv(yk(0), n, MPI_DOUBLE, myid - 1, 1,
					y_up.data(), n, MPI_DOUBLE, myid + 1, MPI_ANY_TAG, MPI_COMM_WORLD, &st);

	}
		}
#endif // SENDRECV

#ifdef ISEND
#ifdef INIT_ISEND
		MPI_Startall(n_rq_send, rq_send);
		MPI_Startall(n_rq_recv, rq_recv);

#else
		n_rq_send = 0;
		n_rq_recv = 0;

		if (myid < np - 1)// кроме последнего
		{
			// создаем запрос rq_send под номером n_rq_send и увеличивыем счетчик отправки
			MPI_Isend(yk(n_loc - 1), n, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD, rq_send + n_rq_send);
			n_rq_send++;

			// создаем запрос rq_recv под номером n_rq_recv и увеличивыем счетчик принятия
			MPI_Irecv(y_up.data(), n, MPI_DOUBLE, myid + 1, MPI_ANY_TAG, MPI_COMM_WORLD, rq_recv + n_rq_recv);
			n_rq_recv++;
		}
		if (myid > 0)  // кроме первого
		{
			MPI_Isend(yk(0), n, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD, rq_send + n_rq_send);
			MPI_Irecv(y_down.data(), n, MPI_DOUBLE, myid - 1, MPI_ANY_TAG, MPI_COMM_WORLD, rq_recv + n_rq_recv);
			n_rq_send++;
			n_rq_recv++;
		}
#endif // INIT_ISEND
#endif // ISEND

		int N = n_loc - 1;
		for (int i = 1; i < N; i++)
			for (int j = (i % 2 + 1); j < n - 1; j += 2)
				yk(i, j) = c * (yk(i + 1, j) + yk(i - 1, j) + yk(i, j + 1) + yk(i, j - 1)) + right(i, j);

		for (int i = 1; i < N; i++)
			for (int j = ((i + 1) % 2 + 1); j < n - 1; j += 2)
				yk(i, j) = c * (yk(i + 1, j) + yk(i - 1, j) + yk(i, j + 1) + yk(i, j - 1)) + right(i, j);


#ifdef ISEND
		// для неблокирующей перессылки ждем, пока все сообщения будут приняты
		MPI_Waitall(n_rq_recv, rq_recv, st);
#endif 

		if (myid < np - 1) {
			for (int j = (N % 2 + 1); j < n - 1; j += 2)
				yk(N, j) = c * (y_up[j] + yk(N - 1, j) + yk(N, j + 1) + yk(N, j - 1)) + right(N, j);

			for (int j = ((N + 1) % 2 + 1); j < n - 1; j += 2)
				yk(N, j) = c * (y_up[j] + yk(N - 1, j) + yk(N, j + 1) + yk(N, j - 1)) + right(N, j);
		}
		if (myid > 0) {
			for (int j = 1; j < n - 1; j += 2)
				yk(0, j) = c * (yk(1, j) + y_down[j] + yk(0, j + 1) + yk(0, j - 1)) + right(0, j);

			for (int j = 2; j < n - 1; j += 2)
				yk(0, j) = c * (yk(1, j) + y_down[j] + yk(0, j + 1) + yk(0, j - 1)) + right(0, j);
		}


#ifdef ISEND
		// для неблокирующей перессылки ждем, пока все сообщения будут отправлены
		MPI_Waitall(n_rq_send, rq_send, st);
#endif 

		iter++;

	} while (iter < max_iter);
#endif // SEIDEL


	//Сбор решения на 0-й процесс
	if (myid == 0)
	{
			
		for (int i = 0; i < n_loc; ++i)
			for (int j = 0; j < n; ++j)
				U(i, j) = yk(i, j);
				
		
		for (int p = 1; p < np; ++p)
			for (int i = 0; i < nLoc[p]; ++i)
				MPI_Recv(U(i + shift[p]), n, MPI_DOUBLE, p, 42, MPI_COMM_WORLD, &stat);
		
	}
	else
	{
		for (int i = 0; i < n_loc; ++i)
			MPI_Send(yk(i), n, MPI_DOUBLE, 0, 42, MPI_COMM_WORLD);
	}

	time += MPI_Wtime();

	if (myid == 0) {
		cout << "Number of iter = " << iter << '\n';
		cout << "Time = " << time << '\n';		
		cout << "Error: " << NormError(h, U) << '\n';

		WriteFile("Solve.dat", 32, h, U);
	}

	MPI_Finalize();
	return 0;


	// почему не работает?   Invalid MPI_Request, error stack.
#ifdef ISEND
	
	for (int i = 0; i < n_rq_send; i++)
		MPI_Request_free(rq_send + i);

	for (int i = 0; i < n_rq_recv; i++)
		MPI_Request_free(rq_recv + i);
#endif // ISEND

	MPI_Finalize();
	return 0;
}

int mainExampleMPI(int argc, char** argv)
{
	int np, myid;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);


	if (myid == 0)
	{
		std::cout << "This is main comp:\n";
	}
	else
	{
		std::cout << "This isn't main comp: " << myid << "\n";
	}

	MPI_Finalize();
}

