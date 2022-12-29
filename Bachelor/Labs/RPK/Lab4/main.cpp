#include <iostream>
#include <fstream>
#include <cmath>
#include <filesystem>

#include "containers_tools.hpp"
#include "pair.hpp"
#include "triple.hpp"
#include "math_tools.hpp"
#include "finite_element.hpp"
#include "double_indexed_vector.hpp"
#include "io_tools.hpp"

template<typename T>
auto mean_by_fe(const double_indexed_vector<pair<T>> &arg)
{
	const std::size_t nx = arg.get_xc(), ny = arg.get_yc();
	const std::size_t lx = nx - 1, ly = ny - 1;

	double_indexed_vector<T> result(nx + 1, ny + 1);

	result.node(0, 0) = arg.node(0, 0).first;
	result.node(nx, ny) = arg.node(lx, ly).second;
	result.node(0, ny) = (arg.node(0, ly).first + arg.node(0, ly).second) / 2.;
	result.node(nx, 0) = (arg.node(lx, 0).first + arg.node(lx, 0).second) / 2.;

	for (std::size_t i = 1; i < nx; ++i)
	{
		result.node(i, 0) = (
			arg.node(i - 1, 0).first +
			arg.node(i - 1, 0).second +
			arg.node(i, 0).first
		) / 3.;
		result.node(i, ny) = (
			arg.node(i - 1, ly).second +
			arg.node(i, ly).first +
			arg.node(i, ly).second
		) / 3.;
	}

	for (std::size_t j = 1; j < ny; ++j)
	{
		result.node(0, j) = (
			arg.node(0, j - 1).first +
			arg.node(0, j - 1).second +
			arg.node(0, j).first
		) / 3.;
		result.node(nx, j) = (
			arg.node(lx, j - 1).second +
			arg.node(lx, j).first +
			arg.node(lx, j).second
		) / 3.;
	}

	for (std::size_t i = 1; i < nx; ++i)
	{
		for (std::size_t j = 1; j < ny; ++j)
		{
			result.node(i, j) = (
				arg.node(i - 1, j - 1).second +
				arg.node(i, j - 1).second +
				arg.node(i - 1, j).second +
				arg.node(i, j).first +
				arg.node(i, j - 1).first +
				arg.node(i - 1, j).first
			) / 6.;
		}
	}

	return result;
}

template<typename T>
inline T f(const pair<T> p)
{
	return 10 + 5 * p.second * std::sin(2 * PI_template<T> *p.first);
}

using number = double;
using numpair = pair<number>;
using numtriple = triple<number>;

int main(int argc, char **argv)
{
	io_tools::set_inline_formats();

	const std::filesystem::path
		input = "input", output = "output";

	std::ofstream fout;

	const numpair A{0, 0}, B{10, 0}, C{11, 10}, D{0, 9};
	const std::size_t nx = 10, ny = 10;
	const auto mesh = make_mesh(A, B, C, D, nx, ny);
	fout.open(output / "mesh.txt");
	for (const auto e : mesh)
		fout << e << std::endl;
	fout.close();

	const auto values = map<double_indexed_vector<number>>(
		mesh, f<number>
	);
	fout.open(output / "values.txt");
	for (const auto e : values)
		fout << e << std::endl;
	fout.close();

	const auto centers = map<double_indexed_vector<numpair>>(
		triangles(mesh),
		static_cast<numpair(*)(triple<numpair>)>(mean<numpair>)
	);
	fout.open(output / "centers.txt");
	for (const auto c : centers)
		fout << c << std::endl;
	fout.close();

	const auto data_fe = pairs(
		map_pair<double_indexed_vector<numtriple>>(
			triangles(mesh), triangles(values), derx_dery_int_fe<number>
		)
	);

	fout.open(output / "derx_fe.txt");
	for (const auto e : data_fe)
		fout << e.first.first << std::endl << e.second.first << std::endl;
	fout.close();

	fout.open(output / "dery_fe.txt");
	for (const auto e : data_fe)
		fout << e.first.second << std::endl << e.second.second << std::endl;
	fout.close();

	fout.open(output / "int_fe.txt");
	for (const auto e : data_fe)
		fout << e.first.third << std::endl << e.second.third << std::endl;
	fout.close();

	const auto data_nodes = mean_by_fe(data_fe);

	fout.open(output / "derx_nodes.txt");
	for (const auto e : data_nodes)
		fout << e.first << std::endl;
	fout.close();

	fout.open(output / "dery_nodes.txt");
	for (const auto e : data_nodes)
		fout << e.second << std::endl;
	fout.close();

	fout.open(output / "int_nodes.txt");
	for (const auto e : data_nodes)
		fout << e.third << std::endl;
	fout.close();

	std::ifstream fin;

	numpair p;
	fin.open(input / "point.txt");
	fin >> p;
	fin.close();

	const auto it_fin_elem = find_first_iterator(
		triangles(mesh),
		[p](const triple<numpair> fe) { return in(fe, p); }
	);
	const auto fin_elem = *it_fin_elem;
	const auto vals = triangles(values)[it_fin_elem.pos()];

	fout.open(output / "point_fin_elem.txt");
	fout
		<< fin_elem.first  << std::endl
		<< fin_elem.second << std::endl
		<< fin_elem.third  << std::endl;
	fout.close();

	const auto [derx_p, dery_p, int_p] = derx_dery_int_fe(fin_elem, vals);
	fout.open(output / "point_val_derx_dery_int.txt");
	fout
		<< p      << std::endl
		<< f(p)   << std::endl
		<< derx_p << std::endl
		<< dery_p << std::endl
		<< int_p  << std::endl;
	fout.close();
	
	return 0;
}