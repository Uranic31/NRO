#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <Eigen/Dense>

double preveri(double x_trenutni, double y_trenutni, double x_sosed, double y_sosed) {

	double pozicija = -1;

	if ((x_trenutni - x_sosed) < 1e-9 && (x_trenutni - x_sosed) > -1e-9) {

		if ((y_trenutni - y_sosed) > 0) {
			pozicija = 1;
		}
		else {
			pozicija = 3;
		}
	}

	else if ((y_trenutni - y_sosed) < 1e-9 && (y_trenutni - y_sosed) > -1e-9) {

		if ((x_trenutni - x_sosed) > 0) {
			pozicija = 0;
		}
		else {
			pozicija = 2;
		}
	}
	return pozicija;
}

Eigen::MatrixXd vectorOfVectorsToEigenMatrix(const std::vector<std::vector<double>>& A) {

	if (A.size() == 0 || A[0].size() == 0)
	{
		std::cerr << "Error: Empty input vector." << std::endl;
		return Eigen::MatrixXd();
	}

	int vrstice = A.size();
	int stolpci = A[0].size();

	Eigen::MatrixXd eigenMatrix(vrstice, stolpci);

	for (int i = 0; i < vrstice; ++i)
	{
		for (int j = 0; j < stolpci; ++j)
		{
			eigenMatrix(i, j) = A[i][j];
		}
	}

	return eigenMatrix;
}

Eigen::VectorXd vectorToEigenVector(const std::vector<double>& b) {
	if (b.size() == 0)
	{
		std::cerr << "Error: Empty input vector." << std::endl;
		return Eigen::VectorXd();
	}

	int elementi = b.size();

	Eigen::VectorXd eigenVector(elementi);

	for (int i = 0; i < elementi; ++i)
	{
		eigenVector(i) = b[i];
	}

	return eigenVector;
}


int main() {

	// Definiranje osnovnih velièin
	double deltaX = 1;
	double deltaY = 1;
	double k = 24;

	// Definiranje vektorjev in matrik
	std::vector<double> X;
	std::vector<double> Y;
	std::vector<std::vector<double>> celice;
	std::vector<std::vector<double>> vozlisca_robnih_pogojev;
	std::vector<double> tipi_robnih_pogojev;
	std::vector<double> vrednosti_robnih_pogojev;
	std::vector<double> vrednosti_prestopa_toplote;
	std::vector<double> vrednosti_toplotnega_toka;

	// BRANJE DATOTEKE

	std::string imedat = "primer1mreza.txt";
	std::ifstream datoteka;
	datoteka.open(imedat);

	if (!datoteka.is_open())
	{
		std::cerr << "Napaka pri odpiranju datoteke: " << imedat << std::endl;
		return 1;
	}
	std::string vrstica;
	std::getline(datoteka, vrstica);

	std::istringstream iss(vrstica);

	std::string prviElement;
	int n_nodes;

	iss >> prviElement >> n_nodes;

	for (int i = 0; i < n_nodes; i++)
	{
		std::string s;
		std::getline(datoteka, s); // preberemo vrstico
		std::replace(s.begin(), s.end(), ';', ' ');
		std::replace(s.begin(), s.end(), ',', ' ');
		std::istringstream iss(s);
		int node_id;
		double x;
		double y;
		iss >> node_id >> x >> y;

		X.push_back(x);
		Y.push_back(y);

	}

	std::string prazna;
	std::getline(datoteka, prazna);

	std::string vrstica2;
	std::getline(datoteka, vrstica2);

	std::istringstream iss1(vrstica2);

	std::string nekej;
	int n_celic;

	iss1 >> nekej >> n_celic;


	for (int i = 0; i < n_celic; i++)
	{
		std::string s;
		std::getline(datoteka, s); // preberemo vrstico
		std::replace(s.begin(), s.end(), ';', ' ');
		std::replace(s.begin(), s.end(), ',', ' ');
		std::istringstream iss(s);

		int nekej2;
		int node1_id;
		int node2_id;
		int node3_id;
		int node4_id;
		iss >> nekej2 >> node1_id >> node2_id >> node3_id >> node4_id;

		std::vector<double> row;
		row.push_back(node1_id);
		row.push_back(node2_id);
		row.push_back(node3_id);
		row.push_back(node4_id);

		celice.push_back(row);

	}

	// 1. ROBNI POGOJ

	std::string prazna2;
	std::getline(datoteka, prazna2);

	std::string st_pog;
	std::getline(datoteka, st_pog);

	std::istringstream iss_st_pog(st_pog);

	std::string robni;
	std::string pogoji;
	int n_pog;

	iss_st_pog >> robni >> pogoji >> n_pog;

	std::string tip_pog1;
	std::getline(datoteka, tip_pog1);

	tipi_robnih_pogojev.push_back(0);

	std::string pog1;
	std::getline(datoteka, pog1);

	std::istringstream iss_pog1(pog1);

	std::string nekej3;
	int temp1;

	iss_pog1 >> nekej3 >> temp1;

	vrednosti_robnih_pogojev.push_back(temp1);

	vrednosti_prestopa_toplote.push_back(-1);
	vrednosti_toplotnega_toka.push_back(-1);

	std::string st_pog1;
	std::getline(datoteka, st_pog1);

	int n_pog1 = std::stoi(st_pog1);


	std::vector<double> pogoj1;
	for (int i = 0; i < n_pog1; i++)
	{
		std::string s;
		std::getline(datoteka, s); // preberemo vrstico
		double pog1 = std::stoi(s);

		pogoj1.push_back(pog1);

		//std::cout << pogoj1[i] << std::endl;
	}
	vozlisca_robnih_pogojev.push_back(pogoj1);

	// 2. ROBNI POGOJ

	std::string prazna3;
	std::getline(datoteka, prazna3);

	std::string tip_pog2;
	std::getline(datoteka, tip_pog2);

	tipi_robnih_pogojev.push_back(2);

	std::string pog2;
	std::getline(datoteka, pog2);

	std::istringstream iss_pog2(pog2);

	std::string nekej4;
	int temp2;

	iss_pog2 >> nekej4 >> temp2;
	vrednosti_robnih_pogojev.push_back(temp2);

	std::string pog22;
	std::getline(datoteka, pog22);

	std::istringstream iss_pog22(pog22);

	std::string nekej5;
	std::string nekej6;
	int pr_top1;

	iss_pog22 >> nekej5 >> nekej6 >> pr_top1;
	vrednosti_prestopa_toplote.push_back(pr_top1);

	vrednosti_toplotnega_toka.push_back(-1);

	std::string st_pog2;
	std::getline(datoteka, st_pog2);

	int n_pog2 = std::stoi(st_pog2);


	std::vector<double> pogoj2;
	for (int i = 0; i < n_pog2; i++)
	{
		std::string s;
		std::getline(datoteka, s);
		double pog2 = std::stoi(s);

		pogoj2.push_back(pog2);

		//std::cout << pogoj2[i] << std::endl;
	}
	vozlisca_robnih_pogojev.push_back(pogoj2);

	// 3. ROBNI POGOJ

	std::string prazna4;
	std::getline(datoteka, prazna4);

	std::string tip_pog3;
	std::getline(datoteka, tip_pog3);

	tipi_robnih_pogojev.push_back(0);

	std::string pog3;
	std::getline(datoteka, pog3);

	std::istringstream iss_pog3(pog3);

	std::string nekej7;
	int temp3;

	iss_pog3 >> nekej7 >> temp3;
	vrednosti_robnih_pogojev.push_back(temp3);

	vrednosti_prestopa_toplote.push_back(-1);
	vrednosti_toplotnega_toka.push_back(-1);

	std::string st_pog3;
	std::getline(datoteka, st_pog3);

	int n_pog3 = std::stoi(st_pog3);


	std::vector<double> pogoj3;
	for (int i = 0; i < n_pog3; i++)
	{
		std::string s;
		std::getline(datoteka, s);
		double pog3 = std::stoi(s);

		pogoj3.push_back(pog3);
	}
	vozlisca_robnih_pogojev.push_back(pogoj3);

	// 4. ROBNI POGOJ

	std::string prazna5;
	std::getline(datoteka, prazna5);

	std::string tip_pog4;
	std::getline(datoteka, tip_pog4);

	tipi_robnih_pogojev.push_back(1);

	std::string pog4;
	std::getline(datoteka, pog4);

	std::istringstream iss_pog4(pog4);

	std::string nekej8;
	int t_tok;

	iss_pog4 >> nekej8 >> t_tok;
	vrednosti_toplotnega_toka.push_back(t_tok);

	vrednosti_robnih_pogojev.push_back(-1);
	vrednosti_prestopa_toplote.push_back(-1);

	std::string st_pog4;
	std::getline(datoteka, st_pog4);

	int n_pog4 = std::stoi(st_pog4);


	std::vector<double> pogoj4;
	for (int i = 0; i < n_pog4; i++)
	{
		std::string s;
		std::getline(datoteka, s);
		double pog4 = std::stoi(s);

		pogoj4.push_back(pog4);

	}
	vozlisca_robnih_pogojev.push_back(pogoj4);

	std::cout << "Datoteka uspesno prebrana!" << std::endl;


	// ISKANJE SOSEDNJIH VOZLIŠÈ

	std::vector<std::vector<double>>sosednja_vozlisca;

	for (double node_id = 0; node_id < n_nodes; node_id++) {

		std::vector<double>vozlisce_i_sosedi = { -1, -1, -1, -1 };

		for (int nd = 0; nd < n_celic; nd++) {

			std::vector<double> trenutna_celica = celice[nd];
			double vozlisce1 = trenutna_celica[0];
			double vozlisce2 = trenutna_celica[1];
			double vozlisce3 = trenutna_celica[2];
			double vozlisce4 = trenutna_celica[3];

			if (node_id == vozlisce1 || node_id == vozlisce2 || node_id == vozlisce3 || node_id == vozlisce4) {

				for (double vozl = 0; vozl < 4; vozl++) {

					double sosednje_vozlisce = trenutna_celica[vozl];

					if (sosednje_vozlisce != node_id) {

						double x_vozl = X[node_id];
						double y_vozl = Y[node_id];
						double x_sosed = X[sosednje_vozlisce];
						double y_sosed = Y[sosednje_vozlisce];

						double pozicija = preveri(x_vozl, y_vozl, x_sosed, y_sosed);

						if (pozicija != -1) {
							vozlisce_i_sosedi[pozicija] = sosednje_vozlisce;
						}
					}
				}
			}
		}

		sosednja_vozlisca.push_back(vozlisce_i_sosedi);
	}
	std::cout << "Iskanje sosedov vozlisc zakljuceno!" << std::endl;

	// GENERIRANJE MATRIKE A IN VEKTORJA b

	std::vector<std::vector<double>> A;
	for (int r = 0; r < n_nodes; r++)
	{
		std::vector<double> stlopci;
		for (int c = 0; c < n_nodes; c++)
		{
			stlopci.push_back(0);
		}
		A.push_back(stlopci);
	}

	std::vector<double> b;
	for (int r = 0; r < n_nodes; r++)
	{
		b.push_back(0);
	}

	std::vector<double> sosedi;
	double levi_sosed;
	double spodnji_sosed;
	double desni_sosed;
	double zgornji_sosed;


	for (int node_id = 0; node_id < n_nodes; node_id++)
	{
		sosedi = sosednja_vozlisca[node_id];

		levi_sosed = sosedi[0];
		spodnji_sosed = sosedi[1];
		desni_sosed = sosedi[2];
		zgornji_sosed = sosedi[3];

		if (levi_sosed != -1 && spodnji_sosed != -1 && desni_sosed != -1 && zgornji_sosed != -1)
		{
			A[node_id][levi_sosed] = 1;
			A[node_id][spodnji_sosed] = 1;
			A[node_id][desni_sosed] = 1;
			A[node_id][zgornji_sosed] = 1;
			A[node_id][node_id] = -4;
		}

		else
		{
			double tip_robnega_pogoja;
			double vrednost;
			double vrednost_prestopa;
			double vrednost_toplotnega_toka;

			for (int robni_pogoj_id = 0; robni_pogoj_id < n_pog; robni_pogoj_id++)
			{
				std::vector<double> vozlisca_v_trenutnem_rp = vozlisca_robnih_pogojev[robni_pogoj_id];

				for (int id_vozlisce_rp = 0; id_vozlisce_rp < vozlisca_v_trenutnem_rp.size(); id_vozlisce_rp++)
				{
					double vozlisce_v_trenutnem_rp = vozlisca_v_trenutnem_rp[id_vozlisce_rp];

					if (node_id == vozlisce_v_trenutnem_rp)
					{
						tip_robnega_pogoja = tipi_robnih_pogojev[robni_pogoj_id];
						vrednost = vrednosti_robnih_pogojev[robni_pogoj_id];
						vrednost_prestopa = vrednosti_prestopa_toplote[robni_pogoj_id];
						vrednost_toplotnega_toka = vrednosti_toplotnega_toka[robni_pogoj_id];
					}
				}
			}

			if (tip_robnega_pogoja == 0)
			{
				A[node_id][node_id] = 1;
				b[node_id] = vrednost;
			}
			else if (tip_robnega_pogoja == 1)
			{
				int stevilo_sosedov = 0;
				for (int st = 0; st < 4; st++)
				{
					if (sosedi[st] != -1)
					{
						stevilo_sosedov += 1;
					}
				}
				if (stevilo_sosedov == 3)
				{
					if (levi_sosed == -1)
					{
						A[node_id][node_id] -= 4;
						A[node_id][desni_sosed] += 2;
						A[node_id][spodnji_sosed] += 1;
						A[node_id][zgornji_sosed] += 1;
						b[node_id] = -2 * (vrednost * deltaX / k);
					}
					if (desni_sosed == -1)
					{
						A[node_id][node_id] -= 4;
						A[node_id][levi_sosed] += 2;
						A[node_id][spodnji_sosed] += 1;
						A[node_id][zgornji_sosed] += 1;
						b[node_id] = -2 * (vrednost * deltaX / k);
					}
					if (spodnji_sosed == -1)
					{
						A[node_id][node_id] -= 4;
						A[node_id][zgornji_sosed] += 2;
						A[node_id][levi_sosed] += 1;
						A[node_id][desni_sosed] += 1;
						b[node_id] = -2 * (vrednost * deltaX / k);
					}
					if (zgornji_sosed == -1)
					{
						A[node_id][node_id] -= 4;
						A[node_id][spodnji_sosed] += 2;
						A[node_id][levi_sosed] += 1;
						A[node_id][desni_sosed] += 1;
						b[node_id] = -2 * (vrednost * deltaX / k);
					}
				}
			}
			else if (tip_robnega_pogoja == 2)
			{
				int stevilo_sosedov = 0;
				for (int st = 0; st < 4; st++)
				{
					if (sosedi[st] != -1)
					{
						stevilo_sosedov += 1;
					}
				}
				if (stevilo_sosedov == 3)
				{
					if (levi_sosed == -1)
					{
						A[node_id][node_id] = A[node_id][node_id] - 2 * (vrednost_prestopa * deltaX / k + 2);
						A[node_id][desni_sosed] += 2;
						A[node_id][spodnji_sosed] += 1;
						A[node_id][zgornji_sosed] += 1;
						b[node_id] = b[node_id] - 2 * vrednost_prestopa * deltaX * vrednost / k;
					}
					if (desni_sosed == -1)
					{
						A[node_id][node_id] = A[node_id][node_id] - 2 * (vrednost_prestopa * deltaX / k + 2);
						A[node_id][levi_sosed] += 2;
						A[node_id][spodnji_sosed] += 1;
						A[node_id][zgornji_sosed] += 1;
						b[node_id] = b[node_id] - 2 * vrednost_prestopa * deltaX * vrednost / k;
					}
					if (spodnji_sosed == -1)
					{
						A[node_id][node_id] = A[node_id][node_id] - 2 * (vrednost_prestopa * deltaX / k + 2);
						A[node_id][zgornji_sosed] += 2;
						A[node_id][levi_sosed] += 1;
						A[node_id][desni_sosed] += 1;
						b[node_id] = b[node_id] - 2 * vrednost_prestopa * deltaX * vrednost / k;
					}
					if (zgornji_sosed == -1)
					{
						A[node_id][node_id] = A[node_id][node_id] - 2 * (vrednost_prestopa * deltaX / k + 2);
						A[node_id][spodnji_sosed] += 2;
						A[node_id][levi_sosed] += 1;
						A[node_id][desni_sosed] += 1;
						b[node_id] = b[node_id] - 2 * vrednost_prestopa * deltaX * vrednost / k;
					}
				}
			}
		}

	}
	std::cout << "Matrika A in vektor b sta generirana!" << std::endl;

	// REŠEVANJE Z GAUSS-SEIDELOVO METODO - 1000 ITERACIJ

	std::vector<double> T_Gauss;
	for (int T_zacetna = 0; T_zacetna < n_nodes; T_zacetna++)
	{
		T_Gauss.push_back(100);
	}

	int st_it = 1000;

#pragma omp parallel for
	for (int ii = 0; ii < st_it; ii++)
	{

#pragma omp critical
		for (int i = 0; i < n_nodes; i++)
		{
			double d = b[i];

			for (int j = 0; j < n_nodes; j++)
			{
				if (i != j)
				{
					d = d - A[i][j] * T_Gauss[j];
				}
			}

			T_Gauss[i] = d / A[i][i];

		}
	}

	std::cout << "Vrednosti temperature po Gauss-Seidlu: " << std::endl;

	for (int i = 0; i < T_Gauss.size(); i++)
	{
		std::cout << T_Gauss[i] << std::endl;
	}


	//TOLE NWM ÈE NE BI KR VN DAU K POMOJM NI TREBA
	// Rešim s knjižnico eigen3 in uporabim funkcijo "partialPivLu" - delna LU dekompozija

	Eigen::MatrixXd MatrikaA_eigen = vectorOfVectorsToEigenMatrix(A);
	std::cout << "Matrika A prevedena v eigen!" << std::endl;
	Eigen::VectorXd Vektorb_eigen = vectorToEigenVector(b);
	std::cout << "Vector b preveden v eigen!" << std::endl;

	Eigen::MatrixXd T_eigen;

	T_eigen = MatrikaA_eigen.partialPivLu().solve(Vektorb_eigen);
	std::cout << "PartialPivLU je zracunu!" << std::endl;

	std::cout << "Resitev temperature PartialPivLU:" << std::endl;
	std::cout << T_eigen << std::endl;


	// VIZUALIZACIJA PODATKOV - VTK

	// Roèni zapis rezultata po VTK - Gauss_Seidel metoda
	std::ofstream rezultat_dat("rezultat_vtk_Gauss.vtk");

	if (!rezultat_dat.is_open())
	{
		std::cerr << "Napaka pri ustvarjanju datoteke!" << std::endl;
		return 1;
	}

	// Glava
	rezultat_dat << "# vtk DataFile Version 3.0\n";
	rezultat_dat << "Mesh_1\n";
	rezultat_dat << "ASCII\n";
	rezultat_dat << "DATASET UNSTRUCTURED_GRID\n";

	// Toèke
	rezultat_dat << "POINTS " << n_nodes << " float\n";
	for (int koordinata_id = 0; koordinata_id < n_nodes; koordinata_id++)
	{
		rezultat_dat << X[koordinata_id] << " " << Y[koordinata_id] << " 0\n";
	}

	rezultat_dat << "\n";

	// Celice
	rezultat_dat << "CELLS " << n_celic << " " << n_celic * 5 << "\n";
	for (int celica_id = 0; celica_id < n_celic; celica_id++)
	{
		int vozl_id1 = celice[celica_id][0];
		int vozl_id2 = celice[celica_id][1];
		int vozl_id3 = celice[celica_id][2];
		int vozl_id4 = celice[celica_id][3];
		rezultat_dat << "4 " << vozl_id1 << " " << vozl_id2 << " " << vozl_id3 << " " << vozl_id4 << "\n";
	}

	rezultat_dat << "\n";

	// Tip celic
	rezultat_dat << "CELL_TYPES " << n_celic << "\n";
	for (int celica_id = 0; celica_id < n_celic; ++celica_id)
	{
		rezultat_dat << "9\n";
	}

	rezultat_dat << "\n";

	// Info toèk
	rezultat_dat << "POINT_DATA " << n_nodes << "\n";
	rezultat_dat << "SCALARS Temperature float 1\n";
	rezultat_dat << "LOOKUP_TABLE default\n";

	for (int koordinata_id = 0; koordinata_id < n_nodes; koordinata_id++)
	{
		rezultat_dat << T_Gauss[koordinata_id] << "\n";
	}

	// Roèni zapis rezultata po VTK - PartialPivLU (eigen3)
	std::ofstream rezultat_dat_eigen("rezultat_vtk_eigen.vtk");

	if (!rezultat_dat_eigen.is_open())
	{
		std::cerr << "Napaka pri ustvarjanju datoteke!" << std::endl;
		return 1;
	}

	// Glava
	rezultat_dat_eigen << "# vtk DataFile Version 3.0\n";
	rezultat_dat_eigen << "Mesh_1\n";
	rezultat_dat_eigen << "ASCII\n";
	rezultat_dat_eigen << "DATASET UNSTRUCTURED_GRID\n";

	// Toèke
	rezultat_dat_eigen << "POINTS " << n_nodes << " float\n";
	for (int koordinata_id = 0; koordinata_id < n_nodes; koordinata_id++)
	{
		rezultat_dat_eigen << X[koordinata_id] << " " << Y[koordinata_id] << " 0\n";
	}

	rezultat_dat_eigen << "\n";

	// Celice
	rezultat_dat_eigen << "CELLS " << n_celic << " " << n_celic * 5 << "\n";
	for (int celica_id = 0; celica_id < n_celic; celica_id++)
	{
		int vozl_id1 = celice[celica_id][0];
		int vozl_id2 = celice[celica_id][1];
		int vozl_id3 = celice[celica_id][2];
		int vozl_id4 = celice[celica_id][3];
		rezultat_dat_eigen << "4 " << vozl_id1 << " " << vozl_id2 << " " << vozl_id3 << " " << vozl_id4 << "\n";
	}

	rezultat_dat_eigen << "\n";

	// Tip celic
	rezultat_dat_eigen << "CELL_TYPES " << n_celic << "\n";
	for (int celica_id = 0; celica_id < n_celic; ++celica_id)
	{
		rezultat_dat_eigen << "9\n";
	}

	rezultat_dat_eigen << "\n";

	// Info toèk
	rezultat_dat_eigen << "POINT_DATA " << n_nodes << "\n";
	rezultat_dat_eigen << "SCALARS Temperature float 1\n";
	rezultat_dat_eigen << "LOOKUP_TABLE default\n";

	rezultat_dat_eigen << T_eigen << "\n";
}
