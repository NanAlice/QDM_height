#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <ctime>
#include <random>
#define ULL unsigned long long

using namespace std;

const int default_L = 4;
const double default_V = 0.9;
const double default_T = 1.0;
const int default_layer = 20;
const double default_J = 1.0;
const int default_loop_n = 250000;
const int default_equil_n = 10000;

const int dx[4] = {1, -1, 0, 0};
const int dy[4] = {0, 0, 1, -1};

unsigned random_seed = (unsigned)time(NULL);
//unsigned random_seed = 812;
std::default_random_engine e(random_seed);
std::uniform_int_distribution<int> random_4(0,3);
std::uniform_real_distribution<double> random_Double(0,1.0);

double random_double()
{
	return random_Double(e);
}
int random_int(int n)
{
	return rand() % n;
}

class Configuration {
public:
	int ***h;
	int ekn;
	int H, W, layer;
	double T, V;
	double J;
	map<ULL, int> state_map;
	int equil_n;
	double dtau, beta;
	double flip_try, flip_succ;
	int loop_n;
	double **M;
	FILE *energy_file, *M_file;
	Configuration(int hh = default_L, int ww = default_L, int Layer = default_layer, double TT = default_T, double VV = default_V, int Equil = default_equil_n): \
	H(hh), W(ww), layer(Layer * 2), T(TT), V(VV), equil_n(Equil) {
		loop_n = default_loop_n;
		J = default_J;
		h = new int **[layer];
		beta = 1.0 / T;
		dtau = beta / (0.5 * layer);
		for (int i = 0; i < layer; ++i) {
			h[i] = new int *[H];
			for (int j = 0; j < H; ++j) {
				h[i][j] = new int[W];
				for (int k = 0; k < W; ++k) {
					//h[i][j][k] = (5 - (k % 2) - (j % 4)) % 4;
					// take the left-down site to be the coordinate of plaquette
					// make a columnar state here
					// dimers along x-axis(between y and y + 1)
					// dimer between x & x + 1: if active, h(x + 1) = h(x) + 1
					// else h(x + 1) = h(x) - 1

					// dimer between y & y + 1: if active:
					// if (x + y + 1) is odd, h(y + 1) = h(y) - 1
					// else h(y + 1) = h(y) + 1
					// else, if not active,
					// if (x + y + 1) is odd, h(y + 1) = h(y) + 1
					// else h(y + 1) = h(y) - 1

					// -----------------------> y increase
					// 0 3 0 3 0 3 0 3
					// 3 2 3 2 3 2 3 2
					// 2 1 2 1 2 1 2 1
					// 1 0 1 0 1 0 1 0

					// x % 4 = 0,1: h[odd] = 3
					// else h[odd] = 1

					// x % 4 = 1,2: h[even] = 2
					// else h[even] = 0
					if ((j + k) % 2) { //odd
						if (j % 4 < 2) {
							h[i][j][k] = 3;
						}
						else {
							h[i][j][k] = 1;
						}
					}
					else {
						if (((j % 4) == 1) || ((j % 4) == 2)) {
							h[i][j][k] = 2;
						}
						else {
							h[i][j][k] = 0;
						}
					}
				}
			}
		}
		ekn = 0;
		state_map.clear();

		M = new double*[2];
		for (int i = 0; i < 2; ++i) {
			M[i] = new double[2];
			for (int j = 0; j < 2; ++j) {
				M[i][j] = 0.0;
			}
		}
	}
	void decode_layer(int z) {
		z = (z + layer) % layer;
		bool **dv, **dh;
		dv = new bool*[H];
		dh = new bool*[H];
		for (int x = 0; x < H; ++x) {
			dv[x] = new bool[W];
			dh[x] = new bool[W];
			for (int y = 0; y < W; ++y) {
				dv[x][y] = (((get(z, x, y - 1) - get(z, x, y) + 4) % 4) == ((((x + H) % H) + ((y + W) % W) + 1) % 2) * 2 + 1);
				dh[x][y] = (((get(z, x, y) - get(z, x - 1, y) + 4) % 4) == 1);
			}
		}
		for (int x = 0; x < H; ++x) {
			for (int y = 0; y < W; ++y) {
				cout << dv[x][y] << ' ';
			}
			cout << endl;
		}
		for (int x = 0; x < H; ++x) {
			for (int y = 0; y < W; ++y) {
				cout << dh[x][y] << ' ';
			}
			cout << endl;
		}
	}
	int get_bit(ULL hashval, int k) {
		return ((int)(hashval >> (2 * k))) & 3;
	}
	void print_state(bool **v) {
		for (int x = 0; x < H; ++x) {
			for (int y = 0; y < W; ++y) {
				cout << v[x][y] << ' ';
			}
			cout << endl;
		}
	}
	void decode_hash(ULL hashval, FILE *f) {
		bool **dv, **dh;
		
		dv = new bool*[H];
		dh = new bool*[H];
		for (int x = 0; x < H; ++x) {
			dv[x] = new bool[W];
			dh[x] = new bool[W];
			for (int y = 0; y < W; ++y) {
				dv[x][y] = (((get_bit(hashval, getn(x, y - 1)) - get_bit(hashval, getn(x, y)) + 4) % 4) == ((((x + H) % H) + ((y + W) % W) + 1) % 2) * 2 + 1);
				dh[x][y] = (((get_bit(hashval, getn(x, y)) - get_bit(hashval, getn(x - 1, y)) + 4) % 4) == 1);
			}
		}
		if (!check_config(dv, dh)) {
			cout << "error!" << endl;
			print_state(dv);
			print_state(dh);
		}
		for (int x = 0; x < H; ++x) {
			for (int y = 0; y < W; ++y) {
				fprintf(f, "%d ", dv[x][y]);
			}
			fprintf(f, "\n");
		}
		for (int x = 0; x < H; ++x) {
			for (int y = 0; y < W; ++y) {
				fprintf(f, "%d ", dh[x][y]);
			}
			fprintf(f, "\n");
		}

	}
	void print_layer(int l) {
		for (int i = 0; i < H; ++i) {
			for (int j = 0; j < W; ++j) {
				cout << h[l][i][j] << ' ';
			}
			cout << endl;
		}
		cout << endl;
		decode_layer(l);
	}
	void print_fp(int l) {
		for (int x = 0; x < H; ++x) {
			for (int y = 0; y < W; ++y) {
				if (FP(l, x, y)) {
					cout << x << ' ' << y << endl;
				}
			}
		}
	}
	int & get(int z, int x, int y) {
		return h[(z + layer) % layer][(x + H) % H][(y + W) % W];
	}
	int FP(int z, int x, int y) {
		bool left, right, up, down;
		// left: between x & x - 1
		// right: x & x + 1
		left = (((get(z, x, y) - get(z, x - 1, y) + 4) % 4) == 1);
		right = (((get(z, x + 1, y) - get(z, x, y) + 4) % 4) == 1);
		if (left && right) {
			return 1;
		}
		// up: between y & y + 1
		// down: y & y - 1
		up = (((get(z, x, y + 1) - get(z, x, y) + 4) % 4) == ((((x + H) % H) + ((y + W) % W) + 1) % 2) * 2 + 1);
		down = (((get(z, x, y - 1) - get(z, x, y) + 4) % 4) == ((((x + H) % H) + ((y + W) % W) + 1) % 2) * 2 + 1);
		if (up && down) {
			return 2;
		}
		return 0;
	}
	bool flippable(int z, int x, int y) {
		int fp = FP(z, x, y);
		flip_try += 1.0;
		if (!fp) {
			return false;
		}
		int zprev = (z - 1 + layer) % layer;
		int znext = (z + 1 + layer) % layer;
		int tx, ty;
		for (int i = 0; i < 4; ++i) {
			tx = (x + dx[i] + H) % H;
			ty = (y + dy[i] + W) % W;
			if (get(z, tx, ty) != get(z + 1, tx, ty)) {
				return false;
			}
			//if (get(z, tx, ty) != get(znext, tx, ty)) {
			//	return false;
			//}
		}
		flip_succ += 1.0;
		return true;
	}
	void flip_only(int z, int x, int y) {
		get(z, x, y) = (get(z, x, y) + 2) % 4;
	}
	double flip_prob(int z, int x, int y) {
		x = (x + H) % H;
		y = (y + W) % W;
		z = (z + layer) % layer;
		double dE = 0;
		int tx, ty;
		for (int i = 0; i < 4; ++i) {
			tx = (x + dx[i] + H) % H;
			ty = (y + dy[i] + W) % W;
			if (FP(z, tx, ty)) {
				dE -= V;
			}
		//	if (FP(z + 1, tx, ty)) {
		//		dE -= V;
		//	}
		}
		flip_only(z, x, y);
		//flip_only(z + 1, x, y);
		for (int i = 0; i < 4; ++i) {
			tx = (x + dx[i] + H) % H;
			ty = (y + dy[i] + W) % W;
			if (FP(z, tx, ty)) {
				dE += V;
			}
		//	if (FP(z + 1, tx, ty)) {
		//		dE += V;
		//	}
		}
		flip_only(z, x, y);
		//flip_only(z + 1, x, y);
		double res = exp(-dtau * dE);
		int zprev = (z - 1 + layer) % layer;
		int znext = (z + 2 + layer) % layer;
		if (get(z, x, y) == get(zprev, x, y)) {
			res *= tanh(dtau * J);
		}
		else {
			res /= tanh(dtau * J);
		}
		if (get(z + 1, x, y) == get(znext, x, y)) {
			res *= tanh(dtau * J);
		}
		else {
			res /= tanh(dtau * J);
		}
		return res;
	}
	void update() {
		int flipx, flipy;
		int can_flip;
		double dE, prob;
		for (int i = 0; i < layer; ++i) {
			for (int x = 0; x < H; ++x) {
				for (int y = 0; y < W; ++y) {
					if (((i & 1) == ((x + y) & 1)) && flippable(i, x, y)) {
						//dE = delta_E(i, x, y);
						if (random_double() < flip_prob(i, x, y)) {
							flip_only(i, x, y);
							flip_only(i + 1, x, y);
						}
					}
				}
			}
		}
	}
	void update_work() {
		ULL tmp;
		double energy = 0.0;
		double this_energy = 0.0;
		for (int i = 0; i < loop_n; ++i) {
			if (i % 1000 == 0) {
				cout << "loop number: " << i << endl;
			}
			update();
			if (i < equil_n) {
				continue;
			}
			for (int z = 0; z < layer; ++z) {
				tmp = hash_hv(z);
				if (state_map.find(tmp) != state_map.end()) {
					++state_map[tmp];
				}
				else {
					state_map[tmp] = 1;
				}
			}
			this_energy = Energy();
			Measure_M();
			fprintf(energy_file, "%lf\n", this_energy / (H * W));
			energy += this_energy;
			//if (!check_node()) {
			//	return;
			//}
		}
		cout << "average energy = " << energy / ((loop_n - equil_n) * (H * W)) << endl;
	}
	int getn(int x, int y) {
		x = (x + H) % H;
		y = (y + W) % W;
		return x * W + y;
	}
	ULL hash_config(int l) {
		ULL res = 0;
		for (int x = 0; x < H; ++x) {
			for (int y = 0; y < W; ++y) {
				res |= ((ULL)(get(l, x, y)) << (getn(x, y) * 2));
			}
		}
		return res;
	}
	ULL hash_hv(int z) {
		z = (z + layer) % layer;
		bool **dv, **dh;
		dv = new bool*[H];
		dh = new bool*[H];
		for (int x = 0; x < H; ++x) {
			dv[x] = new bool[W];
			dh[x] = new bool[W];
			for (int y = 0; y < W; ++y) {
				dv[x][y] = (((get(z, x, y - 1) - get(z, x, y) + 4) % 4) == ((((x + H) % H) + ((y + W) % W) + 1) % 2) * 2 + 1);
				dh[x][y] = (((get(z, x, y) - get(z, x - 1, y) + 4) % 4) == 1);
			}
		}
		ULL res = 0;
		for (int x = 0; x < H; ++x) {
			for (int y = 0; y < W; ++y) {
				res |= ((ULL)(dv[x][y]) << getn(x, y));
			}
		}
		for (int x = 0; x < H; ++x) {
			for (int y = 0; y < W; ++y) {
				res |= ((ULL)(dh[x][y]) << (getn(x, y) + H * W));
			}
		}
		return res;
	}
	void decode_hash_hv(ULL hashval, FILE *f) {
		for (int x = 0; x < H; ++x) {
			for (int y = 0; y < W; ++y) {
				fprintf(f, "%d ", int(hashval >> getn(x, y)) & 1);
			}
			fprintf(f, "\n");
		}
		for (int x = 0; x < H; ++x) {
			for (int y = 0; y < W; ++y) {
				fprintf(f, "%d ", int(hashval >> (getn(x, y) + H * W)) & 1);
			}
			fprintf(f, "\n");
		}
	}
	void output_states(FILE *f) {
		map<ULL, int>::iterator smp;
		for (smp = state_map.begin(); smp != state_map.end(); ++smp) {
			decode_hash_hv(smp->first, f);
			fprintf(f, "%d\n", smp->second);
		}
		cout << state_map.size() << endl;
	}

	bool check_config(bool **dv, bool **dh) {
		int **deg = new int*[H];
		for (int i = 0; i < H; ++i) {
			deg[i] = new int[W];
			for (int j = 0; j < W; ++j) {
				deg[i][j] = 0;
			}
		}
		for (int x = 0; x < H; ++x) {
			for (int y = 0; y < W; ++y) {
				if (dv[x][y]) {
					++deg[x][y];
					++deg[(x + 1) % H][y];
				}
				if (dh[x][y]) {
					++deg[x][y];
					++deg[x][(y + 1) % W];
				}
			}
		}
		for (int x = 0; x < H; ++x) {
			for (int y = 0; y < W; ++y) {
				if (deg[x][y] != 1) {
					return false;
				}
			}
		}
		return true;
	}
	int count_FP() {
		int res = 0;
		for (int z = 0; z < layer; ++z) {
			for (int x = 0; x < H; ++x) {
				for (int y = 0; y < W; ++y) {
					if (FP(z, x, y)) {
						++res;
					}
				}
			}
		}
		return res;
	}
	int count_diff() {
		int res = 0;
		for (int z = 0; z < layer; ++z) {
			for (int x = 0; x < H; ++x) {
				for (int y = 0; y < W; ++y) {
					if (get(z, x, y) != get(z + 1, x, y)) {
						++res;
					}
				}
			}
		}
		return res;
	}
	double Energy() {
		return (count_FP() * 0.5 * V - count_diff() * J) / (0.5 * layer);
		//return (count_FP() * V - ekn * J) / (H * W * layer * 1.0);
		//return -(ekn * J) / (H * W * layer * 1.0);
	}
	double modify_sum(double s)
	{
		return s;
		while (abs(s) > 1.01) {
			if (s > 0.99) {
				s -= 2.0;
			}
			else {
				s += 2.0;
			}
		}
		return s;
	}
	void getM(int l, double **M) {
		double **hn = new double*[H];
		for (int i = 0; i < H; ++i) {
			hn[i] = new double[W];
			for (int j = 0; j < W; ++j) {
				hn[i][j] = h[l][i][j] * 0.5;
				if (hn[i][j] > 1.01) {
					hn[i][j] -= 2.0;
				}
			}
		}
		double A = 0.0, B = 0.0, C = 0.0, D = 0.0;
		double del = -1.0;
		for (int i = 0; i < H / 2; ++i) {
			del *= (- 1.0);
			for (int j = 0; j < W / 2; ++j) {
				B += (hn[2 * i][2 * j] * del);
				D += (hn[2 * i][2 * j + 1] * del);
				A += (hn[2 * i + 1][2 * j] * del * (-1.0));
				C += (hn[2 * i + 1][2 * j + 1] * del * (-1.0));
			}
		}
		A = modify_sum(A); B = modify_sum(B);
		C = modify_sum(C); D = modify_sum(D);
		//cout << "A = " << A << " B = " << B << " C = " << C << " D = " << D << endl;
		M[0][0] = A - B - C + D;
		M[1][1] = A + B - C - D;
		M[0][1] = A - B - C - D;
		M[1][0] = - A + B - C - D;
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				M[i][j] = modify_sum(M[i][j]);
			}
		}
	}
	void Measure_M() {
		if (!M_file) {
			return;
		}
		double **tM, **mM;
		tM = new double*[2];
		mM = new double*[2];
		for (int i = 0; i < 2; ++i) {
			tM[i] = new double[2];
			mM[i] = new double[2];
			for (int j = 0; j < 2; ++j) {
				tM[i][j] = 0;
				mM[i][j] = 0;
			}
		}
		for (int l = 0; l < layer; ++l) {
			getM(l, tM);
			for (int i = 0; i < 2; ++i) {
				for (int j = 0; j < 2; ++j) {
					//mM[i][j] += tM[i][j];
					fprintf(M_file, "%lf ", tM[i][j]);
				}
			}
			fprintf(M_file, "\n");
		}
		/*for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				mM[i][j] /= layer;
				fprintf(M_file, "%lf ", mM[i][j]);
			}
		}
		fprintf(M_file, "\n");*/
	}
	void print_M() {
		double **tM, **mM;
		tM = new double*[2];
		mM = new double*[2];
		for (int i = 0; i < 2; ++i) {
			tM[i] = new double[2];
			mM[i] = new double[2];
			for (int j = 0; j < 2; ++j) {
				tM[i][j] = 0;
				mM[i][j] = 0;
			}
		}
		for (int l = 0; l < layer; ++l) {
			getM(l, tM);
			for (int i = 0; i < 2; ++i) {
				for (int j = 0; j < 2; ++j) {
					mM[i][j] += tM[i][j];
				}
			}
		}
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				cout << M[0][0] << ' ' << M[1][1] << endl;
				cout << M[1][0] << ' ' << M[0][1] << endl;
			}
		}
	}

};
void file_name(char *contents, char *time_string, char *ret) {
	//cout << contents << endl;
	//cout << time_string << endl;
	strcpy(ret, contents);
	strcat(ret, time_string);
	strcat(ret, ".dat");
	//cout << ret << endl;
	//cout << ret << endl;
	//cout <<strcat(strcat(contents, time_string), ".dat") << endl;
	//return strcat(strcat(contents, time_string), ".dat");
}

void print_log(FILE *f, Configuration &C) {
	fprintf(f, "T = %lf\n", C.T);
	fprintf(f, "V = %lf\n", C.V);
	fprintf(f, "loop = %d\n", default_loop_n);
	fprintf(f, "layer = %d\n", C.layer);
	fprintf(f, "H = %d\n", C.H);
	fprintf(f, "W = %d\n", C.W);
}
int main(int argc, char **argv) {
	srand(random_seed);
	Configuration C;
	if (argc > 1) {
		C.loop_n = atoi(argv[1]);
	}
	//C.print_layer(0);
	//C.decode_layer(0);
	char *time_string = new char[20];
	//cout << (int)time(NULL) << endl;
	sprintf(time_string, "%d", (int)time(NULL));
	char *M_file = new char[100];
	char *energy_file = new char[100];
	char *state_file = new char[100];
	file_name("Mij", time_string, M_file);
	file_name("energy", time_string, energy_file);
	file_name("state", time_string, state_file);
	FILE *energy_out = fopen(energy_file, "w");
	print_log(energy_out, C);
	FILE *M_out = fopen(M_file, "w");
	print_log(M_out, C);
	C.energy_file = energy_out;
	C.M_file = M_out;

	C.update_work();
	fclose(energy_out);
	//C.print_layer(0);
	//C.print_layer(0);
	//C.print_fp(0);

	FILE *state_out = fopen(state_file, "w");
	print_log(state_out, C);
	C.output_states(state_out);
	fclose(state_out);

	cout << "flip rate = " << C.flip_succ / C.flip_try << endl;

	return 0;
}