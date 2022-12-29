#pragma once

#include <vector>

#include "containers_tools.hpp"
#include "pair.hpp"
#include "triple.hpp"

template<typename T>
class double_indexed_vector
{
private:
	std::vector<T> data;
	std::size_t xc, yc;

public:
	double_indexed_vector(
		const std::size_t __xc,
		const std::size_t __yc
	) : data(__xc * __yc), xc(__xc), yc(__yc) {}

	double_indexed_vector(
		const double_indexed_vector<T> &other
	) : data(other.data), xc(other.xc), yc(other.yc) {}

	double_indexed_vector(
		double_indexed_vector<T> &&other
	) : data(std::move(other.data)), xc(other.xc), yc(other.yc) {}

	inline std::vector<T> &std_vec() { return data; }
	inline const std::vector<T> &std_vec() const { return data; }

	inline auto get_xc() const { return xc; }
	inline auto get_yc() const { return yc; }

	inline T &node(const std::size_t i, const std::size_t j)
	{
		return data[i * yc + j];
	}

	inline const T &node(const std::size_t i, const std::size_t j) const
	{
		return data[i * yc + j];
	}

	inline triple<T> lower_triangle(const std::size_t i, const std::size_t j) const
	{
		return { node(i, j), node(i + 1, j), node(i, j + 1) };
	}

	inline triple<T> upper_triangle(const std::size_t i, const std::size_t j) const
	{
		return { node(i + 1, j + 1), node(i, j + 1), node(i + 1, j) };
	}

	inline auto begin() { return data.begin(); }
	inline auto begin() const { return data.cbegin(); }
	inline auto сbegin() const { return data.cbegin(); }

	inline auto end() { return data.end(); }
	inline auto end() const { return data.cend(); }
	inline auto сend() const { return data.cend(); }
};

template<typename T>
inline dims_info<2> dims(
	const double_indexed_vector<T> &arg
) {
	return { arg.get_xc(), arg.get_yc() };
}

template<typename double_indexed_container_t>
inline double_indexed_container_t init(
	const dims_info<2> d
) {
	return double_indexed_container_t(d[0], d[1]);
}

enum form { LOWER, UPPER };

struct position { std::size_t i, j; form f; };

template<typename T>
struct triangles
{
	const double_indexed_vector<T> &data;

	triangles(const double_indexed_vector<T> &d) : data(d) {}

	inline triple<T> operator[](const position p) const
	{
		switch (p.f)
		{
		case LOWER: return data.lower_triangle(p.i, p.j);
		case UPPER: return data.upper_triangle(p.i, p.j);
		default: return triple<T>();
		}
	}

	struct iterator
	{
	private:
		const double_indexed_vector<T> &d;
		position p;

	public:
		iterator(
			const double_indexed_vector<T> &__d,
			const position __p
		) : d(__d), p(__p) {}

		inline position pos() const { return p; }

		inline triple<T> operator*() const
		{
			switch (p.f)
			{
			case LOWER: return d.lower_triangle(p.i, p.j);
			case UPPER: return d.upper_triangle(p.i, p.j);
			default: return triple<T>();
			}
		}

		inline iterator &operator++()
		{
			if (p.f == LOWER) p.f = UPPER;
			else if (p.j < d.get_yc() - 2) { p.f = LOWER; ++p.j; }
			else { p.f = LOWER; ++p.i; p.j = 0; }
			return *this;
		}

		inline iterator operator++(int)
		{
			auto it = *this;
			if (p.f == LOWER) p.f = UPPER;
			else if (p.j < d.get_yc() - 2) { p.f = LOWER; ++p.j; }
			else { p.f = LOWER; ++p.i; p.j = 0; }
			return it;
		}

		inline bool operator==(const iterator other)
		{
			return &d == &other.d && p.i == other.p.i && p.j == other.p.j && p.f == other.p.f;
		}

		inline bool operator!=(const iterator other)
		{
			return &d != &other.d || p.i != other.p.i || p.j != other.p.j || p.f != other.p.f;
		}
	};

	inline auto begin() const { return iterator(data, { 0, 0, LOWER }); }

	inline auto end() const { return iterator(data, { data.get_xc() - 1, 0, LOWER }); }

	inline std::size_t size() const { return 2 * (data.get_xc() - 1) * (data.get_yc() - 1); }
};

template<typename T>
inline dims_info<2> dims(
	const triangles<T> &arg
) {
	return { arg.data.get_xc() - 1, 2 * (arg.data.get_yc() - 1) };
}

template<typename T>
using plane_mesh = double_indexed_vector<pair<T>>;

template<typename T>
auto make_mesh(
	const pair<T> a,
	const pair<T> b,
	const pair<T> c,
	const pair<T> d,
	const std::size_t nx,
	const std::size_t ny
) {
	const auto xc = nx + 1, yc = ny + 1;

	plane_mesh<T> result(xc, yc);

	for (std::size_t i = 0; i < xc; ++i)
	{
		for (std::size_t j = 0; j < yc; ++j)
		{
			const T xi = (T)i / nx;
			const T eta = (T) j / ny;

			const T wa = (1 - xi) * (1 - eta);
			const T wb = xi * (1 - eta);
			const T wc = xi * eta;
			const T wd = (1 - xi) * eta;

			result.node(i, j) = wa * a + wb * b + wc * c + wd * d;
		}
	}
	
	return result;
}