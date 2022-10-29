#pragma once

template<typename T>
class Matrix
{
	std::vector<T> data;
	int N, M;

public:
	Matrix() : N(0), M(0)
	{
		data.resize(0);
	}
	Matrix(int pN, T pVal = 0) : N(pN), M(pN)
	{
		data.resize(N * M, pVal);
	}
	Matrix(int pN, int pM, T pVal = 0) : N(pN), M(pM)
	{
		data.resize(N * M, pVal);
	}

	Matrix(const Matrix<T>& other) = default;
	Matrix(Matrix<T>&& other) = default;

	Matrix& operator=(const Matrix<T>& other) = default;
	Matrix& operator=(Matrix<T>&& other) = default;

	~Matrix() = default;

	void Resize(int pN, T pVal = 0) 
	{
		N = pN; M = pN;
		data.resize(N * M, pVal);
	}

	T* operator[](const int i)
	{
		return data.data() + M * i;
	}
	T* operator[](const size_t i)
	{
		return data.data() + M * i;
	}
	const T* operator[](const int i) const
	{
		return data.data() + M * i;
	}
	const T* operator[](const size_t i) const
	{
		return data.data() + M * i;
	}

	typename std::vector<T>::iterator begin()
	{
		return data.begin();
	}
	typename std::vector<T>::iterator end()
	{
		return data.end();
	}

	std::pair<int, int> GetSize() const
	{
		return std::make_pair(N, M);
	}

	std::vector<T>& GetData()
	{
		return data;
	}
};