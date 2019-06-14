#include "main.h"
#include "mrf.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <cstring>
#include <cstdio>
#define MAX(x,y) ((x<y)?y:x)
#define MIN(x,y) ((x<y)?x:y)
#include <unistd.h>

using namespace std;

//RES alignment data
int ** forbidden_res = NULL;
int *res_alignment = NULL;
double res_ub(0);
double res_lb(0);

//Optional parameters
//int n_limit_param(10000000);
int n_limit_param(INFINITY);
//int iter_limit_param(1000);
int iter_limit_param(INFINITY);
//int t_limit(1800);
int t_limit(36000);
int disp_level(1);
float epsilon=1;
//float epsilon= 0.00000001;

//various time check
double load_map_times(0.);
double res_alloc_time(0.);
double res_solve_time(0.);
double total_time(0.);

//
int LA;
int LB;


// scores
double* v_scores;
double* w_scores;

// score function
int use_w = 1;

// gap_open and gap_extend
double gap_open = 0.0;
double gap_extend = 0.0;

// output files
char* output_fname = NULL;
char* info_fname = NULL;


double f_vertex_mrf(int k, int i)
{
	return v_scores[i*LA+k];
}


// score pour l'alignement de deux arcs (A.(i,j),B.(k,l))
double f_edge_mrf(int k, int i, int l, int j)
{
	return w_scores[i*LA*LB*LB+j*LB*LB+k*LB+l];
}



void display_results(int LA, int LB, int nb_edges_A, int nb_edges_B, double self1, double self2, double res_lb, double res_ub, double total_time, double res_alloc_time, double res_solve_time, int nb_bb_nodes, res_alignment)
{
	cout << endl << "RESULT:" << endl;
        if(disp_level >= 1)
        {
                cout << "      |N1| = " << LA <<"\n";
                cout << "      |N2| = " << LB <<"\n";
                cout << "      |E1| = " << nb_edges_A <<"\n";
                cout << "      |E2| = " << nb_edges_B <<"\n";
                cout << "      UB = " << -res_lb <<"\n";
                cout << "      LB = " << -res_ub <<"\n";
                cout << "      Similarity_global = " << 2.0*-res_ub/(self1+self2) <<"\n";
                cout << "      Similarity_global_ub = " << 2.0*-res_lb/(self1+self2) <<"\n";
                cout << "      Similarity_local = " << -res_ub/MIN(self1,self2) <<"\n";
                cout << "      Test_Similarity = " << 2.0*-res_ub/(self1+self2) <<"\n";
                cout << "      T(sec) = " << total_time <<"\n";
                cout << "      Allocation (sec) = " << res_alloc_time << endl;
                cout << "      Solving (sec) = " << res_solve_time << endl;
                cout << "      Number visited nodes = " << nb_bb_nodes << endl;
        }
        cout << "\n";

        display_alignment(LB, res_alignment);
        display_aligned_nodes(LB, res_alignment);

	ofstream output_file;
        output_file.open(info_fname);
        output_file << "similarity_global,time,UB,LB,nb_visited_nodes,selfcomp1,selfcomp2" << endl;
        output_file << 2.0*-res_ub/(self1+self2) << "," << total_time << "," << -res_lb << "," << -res_ub << "," << nb_bb_nodes << "," << self1 << "," << self2 << endl;
        output_file.close();
}



int solve_prb(int ** forbidden, int * sol, double &alloc_time, double &solve_time, double &ub, double &lb, int& nb_bb_nodes, double dalih_bound = 0.0)
{
    std::cout << "started solving problem" << std::endl;
    int status(0);

    double (*score_vertex)(int,int);
    score_vertex = &f_vertex_mrf;
    double (*score_edge)(int,int,int,int);
    score_edge = &f_edge_mrf;

    int nb_row = LA;
    int nb_col = LB;

    int ** row_map = A.get_edges_map(); // TODO réécrire
    int ** col_map = B.get_edges_map();

    //start time
    long tic_per_sec = sysconf(_SC_CLK_TCK);
    struct tms start, end;
    int tic1(times(&start));

    parameters p    = parameters();
    p.max_iteration = iter_limit_param;
    p.max_node      = n_limit_param;
    p.time_limit    = t_limit;
    p.epsilon = epsilon;


    //create problem
    graph_apurva        g(nb_row, row_map, nb_col, col_map, score_vertex, score_edge, forbidden);
    std::cout << "graph apurva ok" << std::endl;
    lambda_mat_apurva   lm(g);
    std::cout << "lambda mat apurva ok" << std::endl;
    dp_mat_apurva       dp(g, gap_open, gap_extend);
    //dp_mat_apurva       dp(g);
    std::cout << "dp mat apurva ok" << std::endl;
    problem_apurva      prb(g,dp,lm,self1,self2);
    std::cout << "problem apurva ok" << std::endl;
    branch_and_bound    bandb;


    //compute allocation time
    int tic2 = times(&end);
    alloc_time = ((double)tic2 - (double)tic1) / (double)tic_per_sec;

    //branch and bound solve
    bandb.set_ub(-dalih_bound);
    bandb.solve(prb, p);
    status = bandb.get_solution_status();
    ub = bandb.get_ub();
    lb = bandb.get_lb();
    bandb.get_solution(sol);
    //get solving time
    solve_time = bandb.get_solve_time();
    nb_bb_nodes = bandb.get_nb_visited_nodes();


    return(status);

}


//display the alignment
void display_alignment(int nb_node2, int *res_alignment)
{
	int col(0);
	int nseg(0);
	while(col < nb_node2)
	{
		//find begin of segment
		while(col < nb_node2 && res_alignment[col] == -1)
		{
			col++;
		}
		if(col < nb_node2)
		{
			nseg++;
			col++;
			
			//find end of segment
			while( col < nb_node2 && res_alignment[col] == (res_alignment[col-1] + 1))
			{
				col++;
			}
		}
	}
	
	cout << "ALIGNMENT:\n" << nseg << " segments\n";
	
	int rowb, rowe,colb,cole;
	int *srowb = new int[nseg];
	int *srowe = new int[nseg];
	int *scolb = new int[nseg];
	int *scole = new int[nseg];
	
	/* Filling segments area */
	col=0;
	nseg=0;

	while(col < nb_node2)
	{
		//find begin of segment
		while(col < nb_node2 && res_alignment[col] == -1)
		{
			col++;
		}
		if(col < nb_node2)
		{
			colb = col;
			rowb = res_alignment[col];
			col++;
			
			//find end of segment
			while(col < nb_node2 && res_alignment[col] == (res_alignment[col-1] + 1))
			{
				col++;
			}
			cole = col-1;
			rowe = res_alignment[col-1];
			
			srowb[nseg] = rowb;
			srowe[nseg] = rowe;
			scolb[nseg] = colb;
			scole[nseg] = cole;
			nseg++;
		}
	}
	
	for(int i(0); i<nseg; ++i)
	{		
		cout << "[" << srowb[i] << " - " << srowe[i] << "] <-> [" << scolb[i] << " - " << scole[i] << "]\n";
	}
	cout << "\n";

	delete[] srowb;
	delete[] srowe;
	delete[] scolb;
	delete[] scole;
}



//display the aligned nodes
void display_aligned_nodes(int nb_node2, int *res_alignment)
{

	// Get the profit for each position
	int** edges_mapA = A.get_edges_map();
	int** edges_mapB = B.get_edges_map();
	vector<double> profits;
	for(int col1=0; col1<nb_node2; ++col1)
	{
		if(res_alignment[col1] != -1)
		{
			double profit = f_vertex_mrf(col1,res_alignment[col1]);
			for(int col2=0; col2<nb_node2; ++col2)
			{
				if(res_alignment[col2] != -1 && col2 != col1)
				{
					if(edges_mapB[col1][col2] == 1 && edges_mapA[res_alignment[col1]][res_alignment[col2]] == 1)
					{
						profit += 0.5*f_edge_mrf(col1,res_alignment[col1],col2,res_alignment[col2]);
					}
				}
			}
			profits.push_back(profit);
		}
	}

	cout << "ALIGNED_NODES:" << endl;

	ofstream output_file;
  	output_file.open(output_fname);
	output_file << "pos_ref,pos_2" << endl;

	vector<double>::iterator it = profits.begin();
	for(int col=0; col<nb_node2; ++col)
	{
		if(res_alignment[col] != -1)
		{
			cout << res_alignment[col] << " " << col << " " << *it << endl;
			output_file << res_alignment[col] << "," << col << endl;
			++it;
		}
	}
	cout << endl;
	output_file.close();

}


int call_from_python(double* v_scores_, double* w_scores_, int LA_, int LB_, double selfcompA_, double selfcompB_, int use_w_, double gap_open_, double gap_extend_, char*output_fname_, char* info_fname_, int n_limit_param_, int iter_limit_param_, int t_limit_, int disp_level_, float epsilon_)
{
	int status(0);
	v_scores = v_scores_;
	w_scores = w_scores_;
	LA = LA_;
	LB = LB_;
	selfcompA = selfcompA_;
	selfcompB = selfcompB_;
	use_w = use_w_;
	gap_open = gap_open_;
	gap_extend = gap_extend_;
	output_fname = output_fname_;
	info_fname = info_fname_;
	n_limit_param = n_limit_param_;
	iter_limit_param_ = iter_limit_param_;
	t_limit = t_limit_;
	disp_level = disp_level_;
	epsilon = epsilon_;

	// computation time
	long tic_per_sec = sysconf(_SC_CLK_TCK);
	struct tms start, end;
	int tic1(times(&start));


	// alloc residues contact map filter (TODO : check if deprecated)
	forbidden_res = new int* [nb_node2];
    	for(int col(0); col != nb_node2; ++col)
    	{
        	forbidden_res[col] = new int [nb_node1];
        	//by default all nodes are allowed
        	fill_n(forbidden_res[col], nb_node1, static_cast<int>(0));
    	}


	// alloc solution arrays
	res_alignment = new int[LB];
	for(int col(0); col != LB; ++col) // TODO pourquoi ?
    	{
                res_alignment[col] = -1;
    	}


	// solve
	int nb_bb_nodes = 0;
	solve_prb(forbidden_res, res_alignment, res_alloc_time, res_solve_time, res_ub, res_lb, nb_bb_nodes, -INFINITY);


	// computation time
	int tic2 = times(&end);
	total_time = ((double)tic2 - (double)tic1) / (double)tic_per_sec;

	display_results_and_print_to_files();

	for(int col(0); col != nb_node2; ++col)
	{
		delete[] forbidden_res[col];
    	}
    	delete[] forbidden_res;
    	delete[] res_alignment;

	return 0;
}


