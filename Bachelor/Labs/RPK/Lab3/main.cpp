#include <iostream>
#include <string>
#include <fstream>
#include <vector>

std::vector<double> sum_with_scale(const std::vector<double>& vec, const std::vector<double>& vec_right, double scale);
double det2(double a11, double a12, double a21, double a22);
std::vector<double> nodes(const std::vector<double>& vec1, const std::vector<double>& vec2, \
	const std::vector<double>& vec3, const std::vector<double>& vec4);

enum DomainType {
	DT_UNKNOWN = 0,
	DT_SEGMENT = 2,
	DT_QUADRANGLE = 4
};

enum SegmentElemType {
	SET_UNKNOWN = 0,
	SET_LINEAR = 1,
	SET_QUADRATIC = 2
};

enum QuadrangleElemType {
	QET_UNKNOWN = 0,
	QET_QUADRANGLE = 1,
	QET_TRIANGLE_SIDE_DIAGONAL = 2,
	QET_TRIANGLE_MAIN_DIAGONAL = 3,
	QET_TRIANGLE_BOTH_DIAGONALS = 4
};

struct DataMesh {
	DomainType NP = DT_UNKNOWN;
	std::vector<double> x, y, z;
	int NE1 = 0, NE2 = 0;
	SegmentElemType segm_type = SET_UNKNOWN;
	QuadrangleElemType quadr_type = QET_UNKNOWN;

	void reserve_mem();
	void read(std::string file_name);
	void print();
	void write(std::string file_name);

	void write_segment_linear(std::ofstream& fout);
	void write_segment_quadratic(std::ofstream& fout);
	void write_quadrangle_quadrangle(std::ofstream& fout);
	void write_quadrangle_side_diag(std::ofstream& fout);
	void write_quadrangle_main_diag(std::ofstream& fout);
	void write_quadrangle_both_diag(std::ofstream& fout);
};

int main() {
	DataMesh data;
	std::string file_name_segm = "Segm";
	std::string file_name_quad = "QuadrangX";
	//data.read("import\\In" + file_name_segm + "1.txt");
	//data.print();
	//data.write("export\\Out" + file_name_segm + "1.txt");
	//std::cout << "=================\n";
	//data.read("import\\In" + file_name_segm + "2.txt");
	//data.print();
	//data.write("export\\Out" + file_name_segm + "2.txt");
	//std::cout << "=================\n";
	data.read("In" + file_name_quad + ".txt");
	data.print();
	data.write("Out" + file_name_quad + "1.txt");
	//std::cout << "=================\n";
	//data.read("import\\In" + file_name_quad + "2.txt");
	//data.print();
	//data.write("export\\Out" + file_name_quad + "2.txt");
	//std::cout << "=================\n";
	//data.read("import\\In" + file_name_quad + "3.txt");
	//data.print();
	//data.write("export\\Out" + file_name_quad + "3.txt");
	//std::cout << "=================\n";
	//data.read("import\\In" + file_name_quad + "4.txt");
	//data.print();
	//data.write("export\\Out" + file_name_quad + "4.txt");
	return 0;
}

std::vector<double> sum_with_scale(const std::vector<double>& vec, const std::vector<double>& vec_right, double scale) {
	std::vector<double> vec_res = vec;
	for (int i = 0; i < vec.size(); ++i) {
		vec_res[i] += scale * vec_right[i];
	}
	return vec_res;
}

double det2(double a11, double a12, double a21, double a22) {
	return a11 * a22 - a12 * a21;
}

double get_param(double a11, double a12, double a21, double a22, double b1, double b2) {
	double denom = a11 * a22 - a12 * a21;
	//return { (a22 * b1 - a12 * b2) / denom, (a11 * b2 - a21 * b1) / denom };
	return (a22 * b1 - a12 * b2) / denom;
}

std::vector<double> nodes(const std::vector<double>& vec1, const std::vector<double>& vec2, \
	const std::vector<double>& vec3, const std::vector<double>& vec4) {

	std::vector<double> res(3);
	double t;

	if (det2(vec2[0] - vec1[0], vec3[0] - vec4[0], vec2[1] - vec1[1], vec3[1] - vec4[1]) != 0.0) {
		// без 3
		//std::cout << "-3\n";
		t = get_param(vec2[0] - vec1[0], vec3[0] - vec4[0], vec2[1] - vec1[1], vec3[1] - vec4[1], \
			vec3[0] - vec1[0], vec3[1] - vec1[1]);
	}
	else if (det2(vec2[0] - vec1[0], vec3[0] - vec4[0], vec2[2] - vec1[2], vec3[2] - vec4[2]) != 0.0) {
		// без 2
		//std::cout << "-2\n";
		t = get_param(vec2[0] - vec1[0], vec3[0] - vec4[0], vec2[2] - vec1[2], vec3[2] - vec4[2], \
			vec3[0] - vec1[0], vec3[2] - vec1[2]);
	}
	else if (det2(vec2[1] - vec1[1], vec3[1] - vec4[1], vec2[2] - vec1[2], vec3[2] - vec4[2]) != 0.0) {
		// без 1
		//std::cout << "-1\n";
		t = get_param(vec2[1] - vec1[1], vec3[1] - vec4[1], vec2[2] - vec1[2], vec3[2] - vec4[2], \
			vec3[1] - vec1[1], vec3[2] - vec1[2]);
	}
	else {
		std::cout << "Miss!\n";
		t = 0;
	}

	for (int i = 0; i < 3; ++i) {
		res[i] = vec1[i] + t * (vec2[i] - vec1[i]);
	}

	return res;
}

void DataMesh::reserve_mem() {
	x.resize(NP); y.resize(NP); z.resize(NP);
}

void DataMesh::read(std::string file_name) {
	int temp;
	std::ifstream fin;
	fin.open(file_name.c_str(), std::ios::in);
	if (!fin.is_open()) {
		std::cout << "File error!\n";
	}
	else {
		fin >> temp;
		NP = (DomainType)temp;
		reserve_mem();
		for (int i = 0; i < NP; ++i) {
			fin >> x[i] >> y[i] >> z[i];
		}
		fin >> NE1;
		if (NP == DT_QUADRANGLE) {
			fin >> NE2;
		}
		fin >> temp;
		if (NP == DT_SEGMENT) {
			segm_type = (SegmentElemType)temp;
		}
		else if (NP == DT_QUADRANGLE) {
			quadr_type = (QuadrangleElemType)temp;
		}
	}
	fin.close();
}

void DataMesh::print() {
	std::cout << "NP = " << NP << "\n";
	std::cout << "x\ty\tz\n";
	for (int i = 0; i < NP; ++i) {
		std::cout << x[i] << "\t" << y[i] << "\t" << z[i] << "\n";
	}
	if (NP == DT_SEGMENT) {
		std::cout << "NE1 = " << NE1 << "\ntype = " << segm_type << "\n";
	}
	else if (NP == DT_QUADRANGLE) {
		std::cout << "NE1 = " << NE1 << ", NE2 = " << NE2 << "\ntype = " << quadr_type << "\n";
	}
}

void DataMesh::write(std::string file_name) {
	std::ofstream fout(file_name.c_str());
	if (!fout.is_open()) {
		std::cout << "File error!\n";
	}
	else {
		switch (NP) {
		case DT_SEGMENT:
			switch (segm_type) {
			case(SET_LINEAR):
				write_segment_linear(fout);
				break;
			case(SET_QUADRATIC):
				write_segment_quadratic(fout);
				break;
			}
			break;

		case DT_QUADRANGLE:
			switch (quadr_type) {
			case(QET_QUADRANGLE):
				write_quadrangle_quadrangle(fout);
				break;
			case(QET_TRIANGLE_SIDE_DIAGONAL):
				write_quadrangle_side_diag(fout);
				break;
			case(QET_TRIANGLE_MAIN_DIAGONAL):
				write_quadrangle_main_diag(fout);
				break;
			case(QET_TRIANGLE_BOTH_DIAGONALS):
				write_quadrangle_both_diag(fout);
				break;
			}
			break;
		}
	}
	fout.close();
	std::cout << "Write complete :)\n";
}

void DataMesh::write_segment_linear(std::ofstream& fout) {
	int NC = 2; // количество контуров
	fout << NE1 << " " << NE1 + 1 << " " << NC << "\n";

	int ENP = 2; // количество узлов в элементе
	for (int i = 1; i <= NE1; ++i) {
		fout << i << " " << ENP << " " << i << " " << i + 1 << "\n";
	}

	std::vector<double> coord(3);
	coord[0] = x[0];
	coord[1] = y[0];
	coord[2] = z[0];
	std::vector<double> base_vec(3);
	base_vec[0] = (x[NP - 1] - x[0]) / NE1;
	base_vec[1] = (y[NP - 1] - y[0]) / NE1;
	base_vec[2] = (z[NP - 1] - z[0]) / NE1;

	for (int i = 0; i < NE1 + 1; ++i) {
		fout << i + 1 << " " << coord[0] << " " << coord[1] << " " << coord[2] << "\n";
		for (int j = 0; j < 3; ++j) {
			coord[j] += base_vec[j];
		}
	}

	int CPN = 1; // количество узлов на контуре
	fout << CPN << " " << CPN << "\n";
	fout << 1 << "\n" << NE1 + 1;
}

void DataMesh::write_segment_quadratic(std::ofstream& fout) {
	int NC = 2; // количество контуров
	fout << NE1 << " " << 2 * NE1 + 1 << " " << NC << "\n";

	int ENP = 3; // количество узлов в элементе
	std::vector<int> EP(3);
	EP[0] = 1;
	EP[1] = 3;
	EP[2] = 2;
	for (int i = 1; i <= NE1; ++i) {
		fout << i << " " << ENP << " " << EP[0] << " " << EP[1] << " " << EP[2] << "\n";
		EP[0] = EP[2];
		EP[1] += 2;
		EP[2] += 2;
	}

	std::vector<double> coord(3);
	coord[0] = x[0];
	coord[1] = y[0];
	coord[2] = z[0];
	std::vector<double> base_vec(3);
	base_vec[0] = (x[NP - 1] - x[0]) / NE1;
	base_vec[1] = (y[NP - 1] - y[0]) / NE1;
	base_vec[2] = (z[NP - 1] - z[0]) / NE1;

	fout << 1 << " " << coord[0] << " " << coord[1] << " " << coord[2] << "\n";
	for (int j = 0; j < 3; ++j) {
		coord[j] += base_vec[j];
	}

	for (int i = 1; i < 2 * NE1 + 1; i += 2) {
		fout << i + 1 << " " << coord[0] << " " << coord[1] << " " << coord[2] << "\n";
		for (int j = 0; j < 3; ++j) {
			coord[j] -= 0.5 * base_vec[j];
		}
		fout << i + 2 << " " << coord[0] << " " << coord[1] << " " << coord[2] << "\n";
		for (int j = 0; j < 3; ++j) {
			coord[j] += 1.5 * base_vec[j];
		}
	}

	int CPN = 1; // количество узлов на контуре
	fout << CPN << " " << CPN << "\n";
	fout << 1 << "\n" << 2 * NE1;
}

void print_vector(std::vector<std::vector<double>> vec) {
	for (int i = 0; i < vec.size(); ++i) {
		for (int j = 0; j < 3; ++j) {
			std::cout << vec[i][j] << " ";
		}
		std::cout << "\n";
	}
}

void DataMesh::write_quadrangle_quadrangle(std::ofstream& fout) {
	int NC = 1; // количество контуров
	fout << NE1 * NE2 << " " << (NE1 + 1) * (NE2 + 1) << " " << NC << "\n";

	int ENP = 4; // количество узлов в элементе

	for (int j = 1; j <= NE2; ++j) {
		for (int i = 1; i <= NE1; ++i) {
			fout << (j - 1) * NE1 + i << " " << ENP << " " << \
				(j - 1) * (NE1 + 1) + i << " " << (j - 1) * (NE1 + 1) + i + 1 << " " << \
				j * (NE1 + 1) + i + 1 << " " << j * (NE1 + 1) + i << "\n";
		}
	}

	std::vector<double> coord1(3);
	coord1[0] = x[0];
	coord1[1] = y[0];
	coord1[2] = z[0];
	std::vector<double> coord2(3);
	coord2[0] = x[1];
	coord2[1] = y[1];
	coord2[2] = z[1];
	std::vector<double> coord3(3);
	coord3[0] = x[2];
	coord3[1] = y[2];
	coord3[2] = z[2];
	std::vector<double> coord4(3);
	coord4[0] = x[3];
	coord4[1] = y[3];
	coord4[2] = z[3];

	std::vector<double> base_vec1(3);
	base_vec1 = sum_with_scale(sum_with_scale(base_vec1, coord2, 1.0 / NE1), coord1, -1.0 / NE1);
	std::vector<double> base_vec2(3);
	base_vec2 = sum_with_scale(sum_with_scale(base_vec2, coord3, 1.0 / NE2), coord2, -1.0 / NE2);
	std::vector<double> base_vec3(3);
	base_vec3 = sum_with_scale(sum_with_scale(base_vec3, coord3, 1.0 / NE1), coord4, -1.0 / NE1);
	std::vector<double> base_vec4(3);
	base_vec4 = sum_with_scale(sum_with_scale(base_vec4, coord4, 1.0 / NE2), coord1, -1.0 / NE2);

	std::vector<std::vector<double>> XYZ12(NE1 + 1);
	for (int i = 0; i < XYZ12.size(); ++i) {
		XYZ12[i] = sum_with_scale(coord1, base_vec1, i);
	}
	std::vector<std::vector<double>> XYZ23(NE2 + 1);
	for (int i = 0; i < XYZ23.size(); ++i) {
		XYZ23[i] = sum_with_scale(coord2, base_vec2, i);
	}
	std::vector<std::vector<double>> XYZ43(NE1 + 1);
	for (int i = 0; i < XYZ43.size(); ++i) {
		XYZ43[i] = sum_with_scale(coord4, base_vec3, i);
	}
	std::vector<std::vector<double>> XYZ14(NE2 + 1);
	for (int i = 0; i < XYZ14.size(); ++i) {
		XYZ14[i] = sum_with_scale(coord1, base_vec4, i);
	}

	/*std::cout << "XYZ12\n"; print_vector(XYZ12);
	std::cout << "XYZ23\n"; print_vector(XYZ23);
	std::cout << "XYZ43\n"; print_vector(XYZ43);
	std::cout << "XYZ14\n"; print_vector(XYZ14);*/

	for (int j = 1; j <= NE2 + 1; ++j) {
		for (int i = 1; i <= NE1 + 1; ++i) {
			std::vector<double> nodes_point = nodes(XYZ14[j - 1], XYZ23[j - 1], XYZ12[i - 1], XYZ43[i - 1]);
			fout << (j - 1) * (NE1 + 1) + i << " " << nodes_point[0] << " " << nodes_point[1] \
				<< " " << nodes_point[2] << "\n";
		}
	}

	int CPN = 2 * (NE1 + NE2); // количество узлов на контуре
	fout << CPN << "\n";
	for (int i = 1; i <= NE1 + 1; ++i) {
		fout << i << "\n";
	}
	for (int j = 2; j <= NE2 + 1; ++j) {
		fout << j * (NE1 + 1) << "\n";
	}
	for (int i = 2; i <= NE1 + 1; ++i) {
		fout << (NE2 + 1) * (NE1 + 1) - i + 1 << "\n";
	}
	for (int j = 1; j < NE2; ++j) {
		fout << (NE2 - j) * (NE1 + 1) + 1 << "\n";
	}
}

void DataMesh::write_quadrangle_side_diag(std::ofstream& fout) {
	int NC = 1; // количество контуров
	fout << 2 * NE1 * NE2 << " " << (NE1 + 1) * (NE2 + 1) << " " << NC << "\n";

	int ENP = 3; // количество узлов в элементе

	for (int j = 1; j <= NE2; ++j) {
		for (int i = 1; i <= NE1; ++i) {
			fout << 2 * (j - 1) * NE1 + i << " " << ENP << " " << \
				(j - 1) * (NE1 + 1) + i << " " << (j - 1) * (NE1 + 1) + i + 1 << " " << \
				j * (NE1 + 1) + i + 1 << "\n";
		}
		for (int i = NE1 + 1; i <= 2 * NE1; ++i) {
			fout << 2 * (j - 1) * NE1 + i << " " << ENP << " " << \
				(j - 1) * (NE1 + 1) + i - NE1 << " " << j * (NE1 + 1) + i + 1 - NE1 << " " << \
				j * (NE1 + 1) + i - NE1 << "\n";
		}
	}

	std::vector<double> coord1(3);
	coord1[0] = x[0];
	coord1[1] = y[0];
	coord1[2] = z[0];
	std::vector<double> coord2(3);
	coord2[0] = x[1];
	coord2[1] = y[1];
	coord2[2] = z[1];
	std::vector<double> coord3(3);
	coord3[0] = x[2];
	coord3[1] = y[2];
	coord3[2] = z[2];
	std::vector<double> coord4(3);
	coord4[0] = x[3];
	coord4[1] = y[3];
	coord4[2] = z[3];

	std::vector<double> base_vec1(3);
	base_vec1 = sum_with_scale(sum_with_scale(base_vec1, coord2, 1.0 / NE1), coord1, -1.0 / NE1);
	std::vector<double> base_vec2(3);
	base_vec2 = sum_with_scale(sum_with_scale(base_vec2, coord3, 1.0 / NE2), coord2, -1.0 / NE2);
	std::vector<double> base_vec3(3);
	base_vec3 = sum_with_scale(sum_with_scale(base_vec3, coord3, 1.0 / NE1), coord4, -1.0 / NE1);
	std::vector<double> base_vec4(3);
	base_vec4 = sum_with_scale(sum_with_scale(base_vec4, coord4, 1.0 / NE2), coord1, -1.0 / NE2);

	std::vector<std::vector<double>> XYZ12(NE1 + 1);
	for (int i = 0; i < XYZ12.size(); ++i) {
		XYZ12[i] = sum_with_scale(coord1, base_vec1, i);
	}
	std::vector<std::vector<double>> XYZ23(NE2 + 1);
	for (int i = 0; i < XYZ23.size(); ++i) {
		XYZ23[i] = sum_with_scale(coord2, base_vec2, i);
	}
	std::vector<std::vector<double>> XYZ43(NE1 + 1);
	for (int i = 0; i < XYZ43.size(); ++i) {
		XYZ43[i] = sum_with_scale(coord4, base_vec3, i);
	}
	std::vector<std::vector<double>> XYZ14(NE2 + 1);
	for (int i = 0; i < XYZ14.size(); ++i) {
		XYZ14[i] = sum_with_scale(coord1, base_vec4, i);
	}

	for (int j = 1; j <= NE2 + 1; ++j) {
		for (int i = 1; i <= NE1 + 1; ++i) {
			std::vector<double> nodes_point = nodes(XYZ14[j - 1], XYZ23[j - 1], XYZ12[i - 1], XYZ43[i - 1]);
			fout << (j - 1) * (NE1 + 1) + i << " " << nodes_point[0] << " " << nodes_point[1] \
				<< " " << nodes_point[2] << "\n";
		}
	}

	int CPN = 2 * (NE1 + NE2); // количество узлов на контуре
	fout << CPN << "\n";
	for (int i = 1; i <= NE1 + 1; ++i) {
		fout << i << "\n";
	}
	for (int j = 2; j <= NE2 + 1; ++j) {
		fout << j * (NE1 + 1) << "\n";
	}
	for (int i = 2; i <= NE1 + 1; ++i) {
		fout << (NE2 + 1) * (NE1 + 1) - i + 1 << "\n";
	}
	for (int j = 1; j < NE2; ++j) {
		fout << (NE2 - j) * (NE1 + 1) + 1 << "\n";
	}
}

void DataMesh::write_quadrangle_main_diag(std::ofstream& fout) {
	int NC = 1; // количество контуров
	fout << 2 * NE1 * NE2 << " " << (NE1 + 1) * (NE2 + 1) << " " << NC << "\n";

	int ENP = 3; // количество узлов в элементе

	for (int j = 1; j <= NE2; ++j) {
		for (int i = 1; i <= NE1; ++i) {
			fout << 2 * (j - 1) * NE1 + i << " " << ENP << " " << \
				(j - 1) * (NE1 + 1) + i << " " << (j - 1) * (NE1 + 1) + i + 1 << " " << \
				j * (NE1 + 1) + i << "\n";
		}
		for (int i = NE1 + 1; i <= 2 * NE1; ++i) {
			fout << 2 * (j - 1) * NE1 + i << " " << ENP << " " << \
				(j - 1) * (NE1 + 1) + i + 1 - NE1 << " " << j * (NE1 + 1) + i + 1 - NE1 << " " << \
				j * (NE1 + 1) + i - NE1 << "\n";
		}
	}

	std::vector<double> coord1(3);
	coord1[0] = x[0];
	coord1[1] = y[0];
	coord1[2] = z[0];
	std::vector<double> coord2(3);
	coord2[0] = x[1];
	coord2[1] = y[1];
	coord2[2] = z[1];
	std::vector<double> coord3(3);
	coord3[0] = x[2];
	coord3[1] = y[2];
	coord3[2] = z[2];
	std::vector<double> coord4(3);
	coord4[0] = x[3];
	coord4[1] = y[3];
	coord4[2] = z[3];

	std::vector<double> base_vec1(3);
	base_vec1 = sum_with_scale(sum_with_scale(base_vec1, coord2, 1.0 / NE1), coord1, -1.0 / NE1);
	std::vector<double> base_vec2(3);
	base_vec2 = sum_with_scale(sum_with_scale(base_vec2, coord3, 1.0 / NE2), coord2, -1.0 / NE2);
	std::vector<double> base_vec3(3);
	base_vec3 = sum_with_scale(sum_with_scale(base_vec3, coord3, 1.0 / NE1), coord4, -1.0 / NE1);
	std::vector<double> base_vec4(3);
	base_vec4 = sum_with_scale(sum_with_scale(base_vec4, coord4, 1.0 / NE2), coord1, -1.0 / NE2);

	std::vector<std::vector<double>> XYZ12(NE1 + 1);
	for (int i = 0; i < XYZ12.size(); ++i) {
		XYZ12[i] = sum_with_scale(coord1, base_vec1, i);
	}
	std::vector<std::vector<double>> XYZ23(NE2 + 1);
	for (int i = 0; i < XYZ23.size(); ++i) {
		XYZ23[i] = sum_with_scale(coord2, base_vec2, i);
	}
	std::vector<std::vector<double>> XYZ43(NE1 + 1);
	for (int i = 0; i < XYZ43.size(); ++i) {
		XYZ43[i] = sum_with_scale(coord4, base_vec3, i);
	}
	std::vector<std::vector<double>> XYZ14(NE2 + 1);
	for (int i = 0; i < XYZ14.size(); ++i) {
		XYZ14[i] = sum_with_scale(coord1, base_vec4, i);
	}

	for (int j = 1; j <= NE2 + 1; ++j) {
		for (int i = 1; i <= NE1 + 1; ++i) {
			std::vector<double> nodes_point = nodes(XYZ14[j - 1], XYZ23[j - 1], XYZ12[i - 1], XYZ43[i - 1]);
			fout << (j - 1) * (NE1 + 1) + i << " " << nodes_point[0] << " " << nodes_point[1] \
				<< " " << nodes_point[2] << "\n";
		}
	}

	int CPN = 2 * (NE1 + NE2); // количество узлов на контуре
	fout << CPN << "\n";
	for (int i = 1; i <= NE1 + 1; ++i) {
		fout << i << "\n";
	}
	for (int j = 2; j <= NE2 + 1; ++j) {
		fout << j * (NE1 + 1) << "\n";
	}
	for (int i = 2; i <= NE1 + 1; ++i) {
		fout << (NE2 + 1) * (NE1 + 1) - i + 1 << "\n";
	}
	for (int j = 1; j < NE2; ++j) {
		fout << (NE2 - j) * (NE1 + 1) + 1 << "\n";
	}
}

void DataMesh::write_quadrangle_both_diag(std::ofstream& fout) {
	int NC = 1; // количество контуров
	fout << 4 * NE1 * NE2 << " " << (NE1 + 1) * (NE2 + 1) + NE1 * NE2 << " " << NC << "\n";

	int ENP = 3; // количество узлов в элементе

	for (int j = 1; j <= NE2; ++j) {
		for (int i = 1; i <= NE1; ++i) {
			fout << 4 * (j - 1) * NE1 + 4 * (i - 1) + 1 << " " << ENP << " " << \
				(j - 1) * (NE1 + 1) + i << " " << (j - 1) * (NE1 + 1) + i + 1 << " " << \
				(NE2 + 1) * (NE1 + 1) + (j - 1) * NE1 + i << "\n";
			fout << 4 * (j - 1) * NE1 + 4 * (i - 1) + 2 << " " << ENP << " " << \
				(j - 1) * (NE1 + 1) + i + 1 << " " << j * (NE1 + 1) + i + 1 << " " << \
				(NE2 + 1) * (NE1 + 1) + (j - 1) * NE1 + i << "\n";
			fout << 4 * (j - 1) * NE1 + 4 * (i - 1) + 3 << " " << ENP << " " << \
				j * (NE1 + 1) + i + 1 << " " << j * (NE1 + 1) + i << " " << \
				(NE2 + 1) * (NE1 + 1) + (j - 1) * NE1 + i << "\n";
			fout << 4 * (j - 1) * NE1 + 4 * (i - 1) + 4 << " " << ENP << " " << \
				j * (NE1 + 1) + i << " " << (j - 1) * (NE1 + 1) + i << " " << \
				(NE2 + 1) * (NE1 + 1) + (j - 1) * NE1 + i << "\n";
		}
	}

	std::vector<double> coord1(3);
	coord1[0] = x[0];
	coord1[1] = y[0];
	coord1[2] = z[0];
	std::vector<double> coord2(3);
	coord2[0] = x[1];
	coord2[1] = y[1];
	coord2[2] = z[1];
	std::vector<double> coord3(3);
	coord3[0] = x[2];
	coord3[1] = y[2];
	coord3[2] = z[2];
	std::vector<double> coord4(3);
	coord4[0] = x[3];
	coord4[1] = y[3];
	coord4[2] = z[3];

	std::vector<double> base_vec1(3);
	base_vec1 = sum_with_scale(sum_with_scale(base_vec1, coord2, 1.0 / NE1), coord1, -1.0 / NE1);
	std::vector<double> base_vec2(3);
	base_vec2 = sum_with_scale(sum_with_scale(base_vec2, coord3, 1.0 / NE2), coord2, -1.0 / NE2);
	std::vector<double> base_vec3(3);
	base_vec3 = sum_with_scale(sum_with_scale(base_vec3, coord3, 1.0 / NE1), coord4, -1.0 / NE1);
	std::vector<double> base_vec4(3);
	base_vec4 = sum_with_scale(sum_with_scale(base_vec4, coord4, 1.0 / NE2), coord1, -1.0 / NE2);

	std::vector<std::vector<double>> XYZ12(NE1 + 1);
	for (int i = 0; i < XYZ12.size(); ++i) {
		XYZ12[i] = sum_with_scale(coord1, base_vec1, i);
	}
	std::vector<std::vector<double>> XYZ23(NE2 + 1);
	for (int i = 0; i < XYZ23.size(); ++i) {
		XYZ23[i] = sum_with_scale(coord2, base_vec2, i);
	}
	std::vector<std::vector<double>> XYZ43(NE1 + 1);
	for (int i = 0; i < XYZ43.size(); ++i) {
		XYZ43[i] = sum_with_scale(coord4, base_vec3, i);
	}
	std::vector<std::vector<double>> XYZ14(NE2 + 1);
	for (int i = 0; i < XYZ14.size(); ++i) {
		XYZ14[i] = sum_with_scale(coord1, base_vec4, i);
	}

	std::vector<std::vector<double>> XYZ;
	for (int j = 1; j <= NE2 + 1; ++j) {
		for (int i = 1; i <= NE1 + 1; ++i) {
			std::vector<double> nodes_point = nodes(XYZ14[j - 1], XYZ23[j - 1], XYZ12[i - 1], XYZ43[i - 1]);
			fout << (j - 1) * (NE1 + 1) + i << " " << nodes_point[0] << " " << nodes_point[1] \
				<< " " << nodes_point[2] << "\n";
			XYZ.push_back(nodes_point);
		}
	}
	for (int j = 1; j <= NE2; ++j) {
		for (int i = 1; i <= NE1; ++i) {
			std::vector<double> nodes_point = nodes(XYZ[(j - 1) * (NE1 + 1) + i - 1], XYZ[j * (NE1 + 1) + i + 1 - 1], \
				XYZ[(j - 1) * (NE1 + 1) + i + 1 - 1], XYZ[j * (NE1 + 1) + i - 1]);

			fout << (NE2 + 1) * (NE1 + 1) + (j - 1) * NE1 + i << " " << nodes_point[0] << " " << nodes_point[1] \
				<< " " << nodes_point[2] << "\n";
		}
	}

	int CPN = 2 * (NE1 + NE2); // количество узлов на контуре
	fout << CPN << "\n";
	for (int i = 1; i <= NE1 + 1; ++i) {
		fout << i << "\n";
	}
	for (int j = 2; j <= NE2 + 1; ++j) {
		fout << j * (NE1 + 1) << "\n";
	}
	for (int i = 2; i <= NE1 + 1; ++i) {
		fout << (NE2 + 1) * (NE1 + 1) - i + 1 << "\n";
	}
	for (int j = 1; j < NE2; ++j) {
		fout << (NE2 - j) * (NE1 + 1) + 1 << "\n";
	}
}