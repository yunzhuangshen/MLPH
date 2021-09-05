#ifndef LSCC_H
#define LSCC_H

#include<limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include<sys/times.h>
#include<unistd.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <string>
#include<math.h>
#include<assert.h>

#include <vector>

namespace LSCC_SOLVER {
#define DEBUG 0
#define	MAXV	1000000
#define MAXE	200000000
struct Edge1{
	long long v1;
	long long v2;
};


class LSCC 
{
long long nv,ne;
long long BMS=100;
tms start, finish;
double time_limit;
double real;
long long lbest;
long long M_iter = 0;

// heap memory
Edge1* edge;
long long* v_degree_tmp;
long long* adjaclen;
long long* neighbor_len;
long long* conf_change;
long long* time_stamp;
long long* temp_array;
long long* vectex;
long long* funch;
long long* address;
long long* tabuin;
long long* cruset;
long long* C0;
long long* C1;
long long* We;
long long* BC;
long long* TC1;
long long* FC1;
long long* Tbest;
long long* TTbest;
long long **neighbor;

double real_solve1=-1;
double real_solve2=-1;
char * File_Name;
const char *resPath = "results/";
const char *dataPath;
long long Max_Vtx,  Max_Iter;
long long f;
long long fbest;
long long len;
long long tm1;
long long tm2;
long long len0; // the length of C0
long long len1; // the length of C1
long long Iter; // the number of iterations taken
long long TABUL = 7;
long long Wf;
long long Wbest;
long long Waim;
long long Titer;
long long len_best = 0;
long long len_W;
long long Iteration[ 100 ];
double time_used[ 100 ];
long long len_used[ 100 ];
long long W_used[ 100 ];
char outfilename[100];
char filename[100];
long long len_improve;
long long len_time;
long long Wmode;
long long kkk;

long long edge_is(long long m, long long n);
void Initializing();
void dump_conf_change();
void neighbor_add(long long node);
void neighbor_drop(long long node);
bool is_forbiden_cc(long long node);
void dump_neighborhood();
void dump_cur_clique();
void Output(long long r);
long long randomInt( long long n );
void clearGamma();
long long selectC0( );
long long WselectC0( );
long long expand(long long SelN);
long long selectC1( );
long long WselectC1( );
long long plateau( long long SelN );
long long Mumi_Weigt();
long long backtract();
long long tabu( long long Max_Iter );
void verify();
void validate();
long long Max_Tabu();

public:
LSCC(){};
~LSCC();
std::vector<double> sol_objs;
std::vector<double> sol_times;
std::vector<std::vector<int>> sols;

long long lscc(long long nb_vtx, long long nb_edge, double cutoff, long long** AdjacentList, long long* Node_Degree, long long* Node_Weight);
};
}
#endif