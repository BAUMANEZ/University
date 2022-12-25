int main(int argc, char** argv)
{
	const double eps = 1.e-6;
	const int N = 10000;
	const int n = N + 1;
	const double h = 1.0 / N;
	const double k = 16000;

	int rank;
	int n_p;
	int m;


	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &n_p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	std::vector<double> y;
	std::vector<double> u;
	std::vector<int> num_row_p(n_p);
	std::vector<int> displs(n_p);

	std::vector<double> y_local;
	std::vector<double> y0_local;

	if (rank == 0)
	{
	
		for (int i = 0; i < n_p; ++i) num_row_p[i] = (n / n_p);

		for (int i = 0; i < n % n_p; ++i) ++num_row_p[i];

		displs.resize(n_p);
		displs[0] = 0;
		for (int i = 1; i < n_p; ++i) displs[i] = displs[i - 1] + num_row_p[i - 1] * n;
		for (int i = 0; i < n_p; ++i) num_row_p[i] += 2;
		num_row_p[0] -= 1;
		num_row_p[n_p - 1] -= 1;
	}
	MPI_Bcast(num_row_p.data(), n_p, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displs.data(), n_p, MPI_INT, 0, MPI_COMM_WORLD);

	double start, finish;
	double time1;
	if (rank == 0)
	{
		cout << "N = " << N << endl;
		cout << "n_p = " << n_p << endl;
		cout << endl;
		y.resize(n * n, 0);
		u.resize(n * n);
		U_vect(u, h, n);

		start = MPI_Wtime();
		Jacobi(y, m, n, h, k, eps, f);
		finish = MPI_Wtime();
		time1 = finish - start;
		double err = error(y, u, n);

		cout << "Jacobi         " << "   iterations =  " << m << "  time =   " << time1
			<< "  err  " << err << std::endl;
		cout << std::endl;

	}
	MPI_Barrier(MPI_COMM_WORLD);

	double time2;
	int num_el = (num_row_p[rank] - (((rank == 0) || (rank == n_p - 1)) ? 1 : 2)) * n;
	std::vector<int> num_elements_p = num_row_p;

	for (int i = 0; i < num_elements_p.size(); ++i)
		num_elements_p[i] *= n;

	y_local.resize(num_row_p[rank] * n, 0.0);
	for (int i = 0; i < y_local.size(); ++i)
		y_local[i] = 0.0;

	y0_local.resize(num_row_p[rank] * n, 0.0);
	for (int i = 0; i < y_local.size(); ++i)
		y_local[i] = 0.0;


	int row_num_global = 0;
	for (int i = 0; i < rank; ++i)
		row_num_global += num_row_p[i];
	row_num_global -= (rank == 0) ? 0 : rank * 2;

	int num_row = num_row_p[rank];

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0)
		start = MPI_Wtime();
	Jacobi_Send_Recv(y_local, n, h, eps, k, f, num_row, rank, n_p, m, row_num_global);
	if (rank == 0)
		finish = MPI_Wtime();

	MPI_Gatherv((rank == 0) ? y_local.data() : y_local.data() + n, num_el, MPI_DOUBLE, y.data(), num_elements_p.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		time2 = finish - start;
		double err = error(y, u, n);
		std::cout << "Jacobi Send + Recv " << "  iterations = " << m << "  time =  " << time2
			<< "  err  " << err << "  Speed-up =  " << time1 / time2 << std::endl;
		std::cout << std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);


	y_local.resize(num_row_p[rank] * n, 0.0);
	for (int i = 0; i < y_local.size(); ++i)
		y_local[i] = 0.0;

		y0_local.resize(num_row_p[rank] * n, 0.0);
	for (int i = 0; i < y_local.size(); ++i)
		y_local[i] = 0.0;

	if (rank == 0)
		start = MPI_Wtime();

	MPI_Barrier(MPI_COMM_WORLD);

	Jacobi_SendRecv(y_local, n, h, eps, k, f, num_row_p[rank], rank, n_p, m, row_num_global);

	if (rank == 0)
		finish = MPI_Wtime();

	time2 = finish - start;

	MPI_Gatherv((rank == 0) ? y_local.data() : y_local.data() + n, num_el, MPI_DOUBLE, y.data(), num_elements_p.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		double err = error(y, u, n);
		std::cout << "Jacobi SendRecv " << "  iterations =  " << m << "  time =  " << time2
			<< "  err  " << err << "  Speed-up =  " << time1 / time2 << std::endl;
		std::cout << std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);


	y_local.resize(num_row_p[rank] * n, 0.0);
	for (int i = 0; i < y_local.size(); ++i)
		y_local[i] = 0.0;

	y0_local.resize(num_row_p[rank] * n, 0.0);
	for (int i = 0; i < y_local.size(); ++i)
		y_local[i] = 0.0;

	if (rank == 0)
		start = MPI_Wtime();

	Jacobi_Isend_Irecv(y_local, n, h, eps, k, f, num_row_p[rank], rank, n_p, m, row_num_global);

	if (rank == 0)
		finish = MPI_Wtime();

	time2 = finish - start;

	MPI_Gatherv((rank == 0) ? y_local.data() : y_local.data() + n, num_el, MPI_DOUBLE, y.data(), num_elements_p.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		double err = error(y, u, n);
		std::cout << "Jacobi Isend + Irecv " << "  iterations =  " << m << "  time =  " << time2
			<< "  err  " << err << "  Speed-up =  " << time1 / time2 << std::endl;
		std::cout << std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);


	if (rank == 0)
	{

		start = MPI_Wtime();
		red_black_it(y, m, n, h, k, eps, f);
		finish = MPI_Wtime();
		time1 = finish - start;

		double err;
		err = error(y, u, n);
		std::cout << std::endl << "red black it seq   " << "  iterations =  " << m << "  time =  " << time1
			<< "  err  " << err << std::endl;
		std::cout << std::endl;

	}
	MPI_Barrier(MPI_COMM_WORLD);


	y_local.resize(num_row_p[rank] * n, 0.0);
	for (int i = 0; i < y_local.size(); ++i)
		y_local[i] = 0.0;

	if (rank == 0)
		start = MPI_Wtime();

	red_black_it_Send_Recv(y_local, n, h, eps, k, f, num_row_p[rank], rank, n_p, m, row_num_global);

	if (rank == 0)
		finish = MPI_Wtime();

	time2 = finish - start;

	MPI_Gatherv((rank == 0) ? y_local.data() : y_local.data() + n, num_el, MPI_DOUBLE, y.data(), num_elements_p.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		double err = error(y, u, n);
		std::cout << "red black Send + Recv " << "  iterations =  " << m << "  time =  " << time2
			<< "  err  " << err << "  Speed-up =  " << time1 / time2 << std::endl;
		std::cout << std::endl;
	}


	MPI_Barrier(MPI_COMM_WORLD);

	y_local.resize(num_row_p[rank] * n, 0.0);
	for (int i = 0; i < y_local.size(); ++i)
		y_local[i] = 0.0;

	if (rank == 0)
		start = MPI_Wtime();

	red_black_it_SendRecv(y_local, n, h, eps, k, f, num_row_p[rank], rank, n_p, m, row_num_global);

	if (rank == 0)
		finish = MPI_Wtime();

	time2 = finish - start;

	MPI_Gatherv((rank == 0) ? y_local.data() : y_local.data() + n, num_el, MPI_DOUBLE, y.data(), num_elements_p.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		double err = error(y, u, n);
		std::cout << "red black SendRecv " << "  iterations =  " << m << "  time =  " << time2
			<< "  err  " << err << "  Speed-up =  " << time1 / time2 << std::endl;
		std::cout << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	y_local.resize(num_row_p[rank] * n, 0.0);
	for (int i = 0; i < y_local.size(); ++i)
		y_local[i] = 0.0;

	if (rank == 0)
		start = MPI_Wtime();

	red_black_it_Isend_Irecv(y_local, n, h, eps, k, f, num_row_p[rank], rank, n_p, m, row_num_global);

	if (rank == 0)
		finish = MPI_Wtime();

	time2 = finish - start;

	MPI_Gatherv((rank == 0) ? y_local.data() : y_local.data() + n, num_el, MPI_DOUBLE, y.data(), num_elements_p.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		double err = error(y, u, n);
		std::cout << "red black Isend_Irecv " << "  iterations =  " << m << "  time =  " << time2
			<< "  err  " << err << "  Speed-up =  " << time1 / time2 << std::endl;
		std::cout << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
}


double f(double x, double y)
{
	double k = 16000;
	return (2 * sin(M_PI * y) + k * k * (1 - x) * x * sin(M_PI * y) + M_PI * M_PI * (1 - x) * x * sin(M_PI * y));
}


double f_Real(double x, double y)
{
	return (1 - x) * x * sin(M_PI * y);
}

double error(const std::vector<double>& y, const std::vector<double>& u, int n)
{
	double norma = 0.0;
	for (int i = 0; i < n * n; ++i)
		norma += (y[i] - u[i]) * (y[i] - u[i]);
	return sqrt(norma);
}



void U_vect(std::vector<double>& y_real, double h, int n)
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			y_real[i * n + j] = u(j * h, i * h);
}


void Jacobi(std::vector<double>& y, int& m, int n, double h, double k, double eps, double (*f)(double, double))
{
	for (int i = 0; i < n * n; ++i)
		y[i] = 0;

	std::vector<double> y_ = y;

	double a = 1 / (4.0 + h * h * k * k);

	m = 0;
	double err = 0.0;
	do
	{
		y.swap(y_);
		err = 0.0;

		for (int i = 1; i < n - 1; ++i)
			for (int j = 1; j < n - 1; ++j)
				y[i * n + j] = a * (y_[i * n + (j - 1)] + y_[(i - 1) * n + j] + y_[(i + 1) * n + j] + y_[i * n + (j + 1)] + h * h * f(j * h, i * h));

		for (int i = 0; i < n * n; ++i) err += sqr(y[i] - y_[i]);
		err = sqrt(err);

		++m;
	} while (err > eps);

}

void red_black_it(std::vector<double>& y, int& m, int n, double h, double k, double eps, double (*f)(double, double))
{
	for (int i = 0; i < n * n; ++i)
		y[i] = 0;

	std::vector<double> y_ = y;

	double a = 1 / (4.0 + h * h * k * k);

	m = 0;
	double err = 0.0;
	do
	{
		y.swap(y_);

		err = 0.0;

		for (int i = 1; i < n - 1; ++i)
			for (int j = (i % 2) + 1; j < n - 1; j += 2)
				y[i * n + j] = a * (y_[i * n + (j - 1)] + y_[(i - 1) * n + j] + y_[(i + 1) * n + j] + y_[i * n + (j + 1)] + h * h * f(j * h, i * h));

		for (int i = 1; i < n - 1; ++i)
			for (int j = ((i + 1) % 2) + 1; j < n - 1; j += 2)
				y[i * n + j] = a * (y[i * n + (j - 1)] + y[(i - 1) * n + j] + y[(i + 1) * n + j] + y[i * n + (j + 1)] + h * h * f(j * h, i * h));

		for (int i = 0; i < n * n; ++i) err += sqr(y[i] - y_[i]);
		err = sqrt(err);

		++m;

	} while (err > eps);

}


void Jacobi_Send_Recv(std::vector<double>& y, int n, double h, double eps, double k, double (*f)(double, double),
	int num_row_p, int rank, int n_p, int& m, int row_global)
{
	std::vector<double> y_ = y;
	double a = 1 / (4.0 + h * h * k * k);


	m = 0;

	int num_up_el = (rank != n_p - 1) ? n : 0; // количество элементов, которые необходимо отправить(получить) к верхнему блоку, либо 0, либо n
	int num_down_el = (rank != 0) ? n : 0; // количество элементов, которые необходимо отправить(получить) к нижнему блоку, либо 0, либо n

	int num_p_up = (rank != n_p - 1) ? rank + 1 : 0; // процесс, которому передаем верхнюю строчку
	int num_p_down = (rank != 0) ? rank - 1 : n_p - 1; // процесс от которого получаем нижнюю строку

	double err;
	double local_err;
	do
	{
		y.swap(y_);

		double* up_row = y_.data() + (num_row_p - 1) * n; // указатель на начало верхней строки
		double* down_row = y_.data(); // указатель на начало нижней строки


		MPI_Send(up_row - n, num_up_el, MPI_DOUBLE, num_p_up, 0, MPI_COMM_WORLD);
		MPI_Recv(down_row, num_down_el, MPI_DOUBLE, num_p_down, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);


		MPI_Send(down_row + n, num_down_el, MPI_DOUBLE, num_p_down, 0, MPI_COMM_WORLD);
		MPI_Recv(up_row, num_up_el, MPI_DOUBLE, num_p_up, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);


		for (int i = 1; i < num_row_p - 1; ++i)
			for (int j = 1; j < n - 1; ++j)
				y[i * n + j] = a * (y_[i * n + (j - 1)] + y_[(i - 1) * n + j] + y_[(i + 1) * n + j] + y_[i * n + (j + 1)] + h * h * f(j * h, (i + row_global) * h));

		local_err = 0.0;

		int i_begin = (rank != 0) ? n : 0;
		int i_end = (rank != n_p - 1) ? ((num_row_p - 1) * n) : num_row_p * n;

		for (int i = i_begin; i < i_end; ++i) local_err += sqr(y[i] - y_[i]);

		MPI_Allreduce(&local_err, &err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		err = sqrt(err);
		++m;
	} while (err > eps);

}

void red_black_it_Send_Recv(std::vector<double>& y, int n, double h, double eps, double k, double (*f)(double, double),
	int num_row_p, int rank, int n_p, int& m, int row_global)
{
	std::vector<double> y_ = y;

	double a = 1 / (4.0 + h * h * k * k);



	int num_up_el = (rank != n_p - 1) ? n : 0; // количество элементов, которые необходимо отправить(получить) к верхнему блоку, либо 0, либо n
	int num_down_el = (rank != 0) ? n : 0; // количество элементов, которые необходимо отправить(получить) к нижнему блоку, либо 0, либо n

	int num_p_up = (rank != n_p - 1) ? rank + 1 : 0; // процесс, которому передаем верхнюю строчку
	int num_p_down = (rank != 0) ? rank - 1 : n_p - 1; // процесс от которого получаем нижнюю строку


	m = 0;
	double err;
	double local_err;
	do
	{
		y_.swap(y);
		double* up_row = y_.data() + (num_row_p - 1) * n; // указатель на начало верхней строки
		double* down_row = y_.data(); // указатель на начало нижней строки

		MPI_Send(up_row - n, num_up_el, MPI_DOUBLE, num_p_up, 0, MPI_COMM_WORLD);
		MPI_Recv(down_row, num_down_el, MPI_DOUBLE, num_p_down, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

		MPI_Send(down_row + n, num_down_el, MPI_DOUBLE, num_p_down, 0, MPI_COMM_WORLD);
		MPI_Recv(up_row, num_up_el, MPI_DOUBLE, num_p_up, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);


		for (int i = 1; i < num_row_p - 1; ++i)
			for (int j = ((i + row_global) % 2) + 1; j < n - 1; j += 2)
				y[i * n + j] = a * (y_[i * n + (j - 1)] + y_[(i - 1) * n + j] + y_[(i + 1) * n + j] + y_[i * n + (j + 1)] + h * h * f(j * h, (i + row_global) * h));

		double* up_row_y = y.data() + (num_row_p - 1) * n; // указатель на начало верхней строки (нижняя строка верхнего блока)
		double* down_row_y = y.data(); // указатель на начало нижней строки

		MPI_Send(up_row_y - n, num_up_el, MPI_DOUBLE, num_p_up, 0, MPI_COMM_WORLD);
		MPI_Recv(down_row_y, num_down_el, MPI_DOUBLE, num_p_down, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

		MPI_Send(down_row_y + n, num_down_el, MPI_DOUBLE, num_p_down, 0, MPI_COMM_WORLD);
		MPI_Recv(up_row_y, num_up_el, MPI_DOUBLE, num_p_up, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

		for (int i = 1; i < num_row_p - 1; ++i)
			for (int j = (((i + row_global) + 1) % 2) + 1; j < n - 1; j += 2)
				y[i * n + j] = a * (y[i * n + (j - 1)] + y[(i - 1) * n + j] + y[(i + 1) * n + j] + y[i * n + (j + 1)] + h * h * f(j * h, (i + row_global) * h));

		local_err = 0.0;

		int i_begin = (rank != 0) ? n : 0;
		int i_end = (rank != n_p - 1) ? ((num_row_p - 1) * n) : num_row_p * n;

		for (int i = i_begin; i < i_end; ++i) local_err += sqr(y[i] - y_[i]);

		MPI_Allreduce(&local_err, &err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		err = sqrt(err);
		++m;
	} while (err > eps);

}

void Jacobi_SendRecv(std::vector<double>& y, int n, double h, double eps, double k, double (*f)(double, double),
	int num_row_p, int rank, int n_p, int& m, int row_global)
{
	std::vector<double> y_ = y;

	double a = 1 / (4.0 + h * h * k * k);

	int num_up_el = (rank != n_p - 1) ? n : 0; // количество элементов, которые необходимо отправить(получить) к верхнему блоку, либо 0, либо n
	int num_down_el = (rank != 0) ? n : 0; // количество элементов, которые необходимо отправить(получить) к нижнему блоку, либо 0, либо n

	int num_p_up = (rank != n_p - 1) ? rank + 1 : 0; // процесс, которому передаем верхнюю строчку
	int num_p_down = (rank != 0) ? rank - 1 : n_p - 1; // процесс от которого получаем нижнюю строку

	double err;
	double local_err;
	m = 0;
	do
	{
		y_.swap(y);

		double* up_row = y_.data() + (num_row_p - 1) * n; // указатель на начало верхней строки (нижняя строка верхнего блока)
		double* down_row = y_.data(); // указатель на начало нижней строки


		MPI_Sendrecv(up_row - n, num_up_el, MPI_DOUBLE, num_p_up, 1, y_.data(), num_down_el, MPI_DOUBLE, num_p_down, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		MPI_Sendrecv(down_row + n, num_down_el, MPI_DOUBLE, num_p_down, 2, up_row, num_up_el, MPI_DOUBLE, num_p_up, 2, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

		for (int i = 1; i < num_row_p - 1; ++i)
			for (int j = 1; j < n - 1; ++j)
				y[i * n + j] = a * (y_[i * n + (j - 1)] + y_[(i - 1) * n + j] + y_[(i + 1) * n + j] + y_[i * n + (j + 1)] + h * h * f(j * h, (i + row_global) * h));

		local_err = 0.0;

		int i_begin = (rank != 0) ? n : 0;
		int i_end = (rank != n_p - 1) ? ((num_row_p - 1) * n) : num_row_p * n;

		for (int i = i_begin; i < i_end; ++i) local_err += sqr(y[i] - y_[i]);

		MPI_Allreduce(&local_err, &err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		err = sqrt(err);
		++m;
	} while (err > eps);

}


void red_black_it_SendRecv(std::vector<double>& y, int n, double h, double eps, double k, double (*f)(double, double),
	int num_row_p, int rank, int n_p, int& m, int row_global) {
	std::vector<double> y_ = y;

	double a = 1 / (4.0 + h * h * k * k);


	int num_up_el = (rank != n_p - 1) ? n : 0; // количество элементов, которые необходимо отправить(получить) к верхнему блоку, либо 0, либо n
	int num_down_el = (rank != 0) ? n : 0; // количество элементов, которые необходимо отправить(получить) к нижнему блоку, либо 0, либо n

	int num_p_up = (rank != n_p - 1) ? rank + 1 : 0; // процесс, которому передаем верхнюю строчку
	int num_p_down = (rank != 0) ? rank - 1 : n_p - 1; // процесс от которого получаем нижнюю строку

	m = 0;
	double err;
	double local_err;
	do
	{
		y.swap(y_);

		double* up_row = y_.data() + (num_row_p - 1) * n; // указатель на начало верхней строки (нижняя строка верхнего блока)
		double* down_row = y_.data(); // указатель на начало нижней строки


		MPI_Sendrecv(up_row - n, num_up_el, MPI_DOUBLE, num_p_up, 1, down_row, num_down_el, MPI_DOUBLE, num_p_down, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		MPI_Sendrecv(down_row + n, num_down_el, MPI_DOUBLE, num_p_down, 2, up_row, num_up_el, MPI_DOUBLE, num_p_up, 2, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);


		for (int i = 1; i < num_row_p - 1; ++i)
			for (int j = ((i + row_global) % 2) + 1; j < n - 1; j += 2)
				y[i * n + j] = a * (y_[i * n + (j - 1)] + y_[(i - 1) * n + j] + y_[(i + 1) * n + j] + y_[i * n + (j + 1)] + h * h * f(j * h, (i + row_global) * h));


		double* up_row_y = y.data() + (num_row_p - 1) * n; // указатель на начало верхней строки (нижняя строка верхнего блока)
		double* down_row_y = y.data(); // указатель на начало нижней строки


		MPI_Sendrecv(up_row_y - n, num_up_el, MPI_DOUBLE, num_p_up, 1, down_row_y, num_down_el, MPI_DOUBLE, num_p_down, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		MPI_Sendrecv(down_row_y + n, num_down_el, MPI_DOUBLE, num_p_down, 2, up_row_y, num_up_el, MPI_DOUBLE, num_p_up, 2, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

		for (int i = 1; i < num_row_p - 1; ++i)
			for (int j = (((i + row_global) + 1) % 2) + 1; j < n - 1; j += 2)
				y[i * n + j] = a * (y[i * n + (j - 1)] + y[(i - 1) * n + j] + y[(i + 1) * n + j] + y[i * n + (j + 1)] + h * h * f(j * h, (i + row_global) * h));

		local_err = 0.0;

		int i_begin = (rank != 0) ? n : 0;
		int i_end = (rank != n_p - 1) ? ((num_row_p - 1) * n) : num_row_p * n;

		for (int i = i_begin; i < i_end; ++i) local_err += sqr(y[i] - y_[i]);

		MPI_Allreduce(&local_err, &err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		err = sqrt(err);
		++m;
	} while (err > eps);

}

void Jacobi_Isend_Irecv(std::vector<double>& y, int n, double h, double eps, double k, double (*f)(double, double),
	int num_row_p, int rank, int n_p, int& m, int row_global)
{
	std::vector<double> y_ = y;

	double a = 1 / (4.0 + h * h * k * k);


	int num_up_el = (rank != n_p - 1) ? n : 0;
	int num_down_el = (rank != 0) ? n : 0;

	int num_p_up = (rank != n_p - 1) ? rank + 1 : 0;
	int num_p_down = (rank != 0) ? rank - 1 : n_p - 1;


	MPI_Request req_up_send_y, req_up_recv_y;
	MPI_Request req_down_send_y, req_down_recv_y;

	double* up_row_y = y.data() + (num_row_p - 1) * n;
	double* down_row_y = y.data();

	MPI_Send_init(up_row_y - n, num_up_el, MPI_DOUBLE, num_p_up, 1, MPI_COMM_WORLD, &req_up_send_y);
	MPI_Recv_init(down_row_y, num_down_el, MPI_DOUBLE, num_p_down, 1, MPI_COMM_WORLD, &req_down_recv_y);

	MPI_Send_init(down_row_y + n, num_down_el, MPI_DOUBLE, num_p_down, 2, MPI_COMM_WORLD, &req_down_send_y);
	MPI_Recv_init(up_row_y, num_up_el, MPI_DOUBLE, num_p_up, 2, MPI_COMM_WORLD, &req_up_recv_y);

	MPI_Request req_up_send_y0, req_up_recv_y0;
	MPI_Request req_down_send_y0, req_down_recv_y0;

	double* up_row_y0 = y_.data() + (num_row_p - 1) * n; // указатель на начало верхней строки (нижняя строка верхнего блока)
	double* down_row_y0 = y_.data();


	MPI_Send_init(up_row_y0 - n, num_up_el, MPI_DOUBLE, num_p_up, 3, MPI_COMM_WORLD, &req_up_send_y0);
	MPI_Recv_init(down_row_y0, num_down_el, MPI_DOUBLE, num_p_down, 3, MPI_COMM_WORLD, &req_down_recv_y0);

	MPI_Send_init(down_row_y0 + n, num_down_el, MPI_DOUBLE, num_p_down, 4, MPI_COMM_WORLD, &req_down_send_y0);
	MPI_Recv_init(up_row_y0, num_up_el, MPI_DOUBLE, num_p_up, 4, MPI_COMM_WORLD, &req_up_recv_y0);

	double err;
	double local_err;
	m = 0;
	do
	{
		y.swap(y_);

		if (m % 2 == 0)
		{
			MPI_Start(&req_down_recv_y);
			MPI_Start(&req_down_send_y);
			MPI_Start(&req_up_recv_y);
			MPI_Start(&req_up_send_y);
		}
		else
		{
			MPI_Start(&req_down_recv_y0);
			MPI_Start(&req_down_send_y0);
			MPI_Start(&req_up_recv_y0);
			MPI_Start(&req_up_send_y0);
		}

		for (int i = 2; i < num_row_p - 2; ++i)
			for (int j = 1; j < n - 1; ++j)
				y[i * n + j] = a * (y_[i * n + (j - 1)] + y_[(i - 1) * n + j] + y_[(i + 1) * n + j] + y_[i * n + (j + 1)] + h * h * f(j * h, (i + row_global) * h));

		if (m % 2 == 0)
		{
			MPI_Wait(&req_down_recv_y, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_down_send_y, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_up_recv_y, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_up_send_y, MPI_STATUSES_IGNORE);
		}
		else
		{
			MPI_Wait(&req_down_recv_y0, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_down_send_y0, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_up_recv_y0, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_up_send_y0, MPI_STATUSES_IGNORE);
		}

		int i = 1;
		for (int j = 1; j < n - 1; ++j)
			y[i * n + j] = a * (y_[i * n + (j - 1)] + y_[(i - 1) * n + j] + y_[(i + 1) * n + j] + y_[i * n + (j + 1)] + h * h * f(j * h, (i + row_global) * h));

		i = num_row_p - 2;
		for (int j = 1; j < n - 1; ++j)
			y[i * n + j] = a * (y_[i * n + (j - 1)] + y_[(i - 1) * n + j] + y_[(i + 1) * n + j] + y_[i * n + (j + 1)] + h * h * f(j * h, (i + row_global) * h));

		local_err = 0.0;

		int i_begin = (rank != 0) ? n : 0;
		int i_end = (rank != n_p - 1) ? ((num_row_p - 1) * n) : num_row_p * n;

		for (int i = i_begin; i < i_end; ++i) local_err += sqr(y[i] - y_[i]);

		MPI_Allreduce(&local_err, &err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		err = sqrt(err);

		++m;
	} while (err > eps);

}

void red_black_it_Isend_Irecv(std::vector<double>& y, int n, double h, double eps, double k, double (*f)(double, double),
	int num_row_p, int rank, int n_p, int& m, int row_global)
{
	std::vector<double> y_ = y;

	double a = 1 / (4.0 + h * h * k * k);


	int num_up_el = (rank != n_p - 1) ? n : 0;
	int num_down_el = (rank != 0) ? n : 0;

	int num_p_up = (rank != n_p - 1) ? rank + 1 : 0;
	int num_p_down = (rank != 0) ? rank - 1 : n_p - 1;


	MPI_Request req_up_send_y, req_up_recv_y;
	MPI_Request req_down_send_y, req_down_recv_y;

	double* up_row_y = y.data() + (num_row_p - 1) * n;
	double* down_row_y = y.data();

	MPI_Send_init(up_row_y - n, num_up_el, MPI_DOUBLE, num_p_up, 1, MPI_COMM_WORLD, &req_up_send_y);
	MPI_Recv_init(down_row_y, num_down_el, MPI_DOUBLE, num_p_down, 1, MPI_COMM_WORLD, &req_down_recv_y);

	MPI_Send_init(down_row_y + n, num_down_el, MPI_DOUBLE, num_p_down, 2, MPI_COMM_WORLD, &req_down_send_y);
	MPI_Recv_init(up_row_y, num_up_el, MPI_DOUBLE, num_p_up, 2, MPI_COMM_WORLD, &req_up_recv_y);

	MPI_Request req_up_send_y0, req_up_recv_y0;
	MPI_Request req_down_send_y0, req_down_recv_y0;

	double* up_row_y0 = y_.data() + (num_row_p - 1) * n; // указатель на начало верхней строки (нижняя строка верхнего блока)
	double* down_row_y0 = y_.data();


	MPI_Send_init(up_row_y0 - n, num_up_el, MPI_DOUBLE, num_p_up, 3, MPI_COMM_WORLD, &req_up_send_y0);
	MPI_Recv_init(down_row_y0, num_down_el, MPI_DOUBLE, num_p_down, 3, MPI_COMM_WORLD, &req_down_recv_y0);

	MPI_Send_init(down_row_y0 + n, num_down_el, MPI_DOUBLE, num_p_down, 4, MPI_COMM_WORLD, &req_down_send_y0);
	MPI_Recv_init(up_row_y0, num_up_el, MPI_DOUBLE, num_p_up, 4, MPI_COMM_WORLD, &req_up_recv_y0);


	m = 0;
	double err;
	double local_err;
	do
	{
		y.swap(y_);

		int i = 1;
		for (int j = ((i + row_global) % 2) + 1; j < n - 1; ++j)
			y[i * n + j] = a * (y_[i * n + j - 1] + y_[i * n + j + 1] + y_[(i - 1) * n + j] + y_[(i + 1) * n + j] + h * h * f(j * h, (i + row_global) * h));

		i = num_row_p - 2;
		for (int j = ((i + row_global) % 2) + 1; j < n - 1; ++j)
			y[i * n + j] = a * (y_[i * n + j - 1] + y_[i * n + j + 1] + y_[(i - 1) * n + j] + y_[(i + 1) * n + j] + h * h * f(j * h, (i + row_global) * h));

		if (m % 2 == 1)
		{
			MPI_Start(&req_down_recv_y);
			MPI_Start(&req_down_send_y);
			MPI_Start(&req_up_recv_y);
			MPI_Start(&req_up_send_y);
		}
		else
		{
			MPI_Start(&req_down_recv_y0);
			MPI_Start(&req_down_send_y0);
			MPI_Start(&req_up_recv_y0);
			MPI_Start(&req_up_send_y0);
		}


		for (int i = 2; i < num_row_p - 2; ++i)
			for (int j = ((i + row_global) % 2) + 1; j < n - 1; j += 2)
				y[i * n + j] = a * (y_[i * n + j - 1] + y_[i * n + j + 1] + y_[(i - 1) * n + j] + y_[(i + 1) * n + j] + h * h * f(j * h, (i + row_global) * h));

		if (m % 2 == 1)
		{
			MPI_Wait(&req_down_recv_y, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_down_send_y, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_up_recv_y, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_up_send_y, MPI_STATUSES_IGNORE);
		}
		else
		{
			MPI_Wait(&req_down_recv_y0, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_down_send_y0, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_up_recv_y0, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_up_send_y0, MPI_STATUSES_IGNORE);
		}



		i = 1;
		for (int j = (((i + row_global) + 1) % 2) + 1; j < n - 1; j += 2)
			y[i * n + j] = a * (y[i * n + j - 1] + y[i * n + j + 1] + y[(i - 1) * n + j] + y[(i + 1) * n + j] + h * h * f(j * h, (i + row_global) * h));

		i = num_row_p - 2;
		for (int j = (((i + row_global) + 1) % 2) + 1; j < n - 1; j += 2)
			y[i * n + j] = a * (y[i * n + j - 1] + y[i * n + j + 1] + y[(i - 1) * n + j] + y[(i + 1) * n + j] + h * h * f(j * h, (i + row_global) * h));


		if (m % 2 == 1)
		{
			MPI_Start(&req_down_recv_y);
			MPI_Start(&req_down_send_y);
			MPI_Start(&req_up_recv_y);
			MPI_Start(&req_up_send_y);
		}
		else
		{
			MPI_Start(&req_down_recv_y0);
			MPI_Start(&req_down_send_y0);
			MPI_Start(&req_up_recv_y0);
			MPI_Start(&req_up_send_y0);
		}

		for (int i = 2; i < num_row_p - 2; ++i)
			for (int j = (((i + row_global) + 1) % 2) + 1; j < n - 1; j += 2)
				y[i * n + j] = a * (y[i * n + j - 1] + y[i * n + j + 1] + y[(i - 1) * n + j] + y[(i + 1) * n + j] + h * h * f(j * h, (i + row_global) * h));

		if (m % 2 == 1)
		{
			MPI_Wait(&req_down_recv_y, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_down_send_y, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_up_recv_y, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_up_send_y, MPI_STATUSES_IGNORE);
		}
		else
		{
			MPI_Wait(&req_down_recv_y0, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_down_send_y0, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_up_recv_y0, MPI_STATUSES_IGNORE);
			MPI_Wait(&req_up_send_y0, MPI_STATUSES_IGNORE);
		}


		local_err = 0.0;

		int i_begin = (rank != 0) ? n : 0;
		int i_end = (rank != n_p - 1) ? ((num_row_p - 1) * n) : num_row_p * n;

		for (int i = i_begin; i < i_end; ++i) local_err += sqr(y[i] - y_[i]);

		MPI_Allreduce(&local_err, &err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		err = sqrt(err);

		++m;
	} while (err > eps);
}


