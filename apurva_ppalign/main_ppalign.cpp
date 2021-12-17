#include "main.h"
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cstring>
#include <cstdio>
#define MAX(x,y) ((x<y)?y:x)
#define MIN(x,y) ((x<y)?x:y)
#include <unistd.h>

using namespace std;

//various time check
double load_map_times(0.);
double total_time(0.);


// scores and lengths are global variables because of f_vertex_mrf and f_edge_mrf
int LA;
int LB;
int q=21;

double alpha_w;
double offset_v;

float* v_A;
float* v_B;
float* w_A;
float* w_B;

// gap_open and gap_extend
//double gap_open = 0.0;
//double gap_extend = 0.0;

// position-specific insertion costs
float* insert_open_A;
float* insert_open_B;
float* insert_extend_A;
float* insert_extend_B;


// similarity function for vA_i and vB_k
float f_vertex_mrf(int k, int i)
{
	float score_vivk=0;
	for(int a=0; a<q; a++)
	{
		score_vivk+=v_A[i*q+a]*v_B[k*q+a];
	}
	return score_vivk-offset_v;
}




// similarity function for wA_ij and wB_kl
float f_edge_mrf(int k, int i, int l, int j)
{

	float score_wijwkl=0;
	for (int a=0; a<q; a++)
	{
		for (int b=0; b<q; b++)
		{
			if ((k<=l) && (i<=j))
			{
				score_wijwkl+=w_A[b+a*q+j*q*q+i*q*q*LA]*w_B[b+a*q+l*q*q+k*q*q*LB];
			}
			else if ((k<=l) && (i>j))
			{
				score_wijwkl+=w_A[a+b*q+i*q*q+j*q*q*LA]*w_B[b+a*q+l*q*q+k*q*q*LB];
			}
			else if ((k>l) && (i<=j))
			{
				score_wijwkl+=w_A[b+a*q+j*q*q+i*q*q*LA]*w_B[a+b*q+k*q*q+l*q*q*LB];
			}
			else
			{
				score_wijwkl+=w_A[a+b*q+i*q*q+j*q*q*LA]*w_B[a+b*q+k*q*q+l*q*q*LB];
			}
		}
	}
	return alpha_w*score_wijwkl;

}


/*
float f_edge_mrf(int k, int i, int l, int j)
{
	if (l<=k)
	{
		if (j<=i)
		{
			float score_wijwkl=0;
			for (int a=0; a<q; a++)
			{
				for (int b=0; b<q; b++)
				{
					score_wijwkl+=w_A[b+a*q+j*q*q+i*q*q*LA]*w_B[b+a*q+l*q*q+k*q*q*LB];
				}
			}
			return alpha_w*score_wijwkl;
		}
		else
		{
			return f_edge_mrf(k,j,l,i);
		}
	}
	else
	{
		return f_edge_mrf(l,i,k,j);
	}
}
*/

//display the alignment
void display_alignment(int *res_alignment)
{
	int col(0);
	int nseg(0);
	while(col < LB)
	{
		//find begin of segment
		while(col < LB && res_alignment[col] == -1)
		{
			col++;
		}
		if(col < LB)
		{
			nseg++;
			col++;
			
			//find end of segment
			while( col < LB && res_alignment[col] == (res_alignment[col-1] + 1))
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

	while(col < LB)
	{
		//find begin of segment
		while(col < LB && res_alignment[col] == -1)
		{
			col++;
		}
		if(col < LB)
		{
			colb = col;
			rowb = res_alignment[col];
			col++;
			
			//find end of segment
			while(col < LB && res_alignment[col] == (res_alignment[col-1] + 1))
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
void display_aligned_nodes(int *res_alignment, int** row_map, int** col_map, char* output_fname)
{

	// Get the profit for each position
	vector<double> profits;
	for(int col1=0; col1<LB; ++col1)
	{
		if(res_alignment[col1] != -1)
		{
			double profit = f_vertex_mrf(col1,res_alignment[col1]);
			for(int col2=0; col2<LB; ++col2)
			{
				if(res_alignment[col2] != -1 && col2 != col1)
				{
					if(col_map[col1][col2] == 1 && row_map[res_alignment[col1]][res_alignment[col2]] == 1)
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
	for(int col=0; col<LB; ++col)
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



void display_alignment_with_gaps(int* res_alignment, int* res_alignment_insert_before, char* output_fname)
{
	cout << "DISPLAY ALIGNMENT WITH GAPS" << endl;
	ofstream output_file;
  	output_file.open(output_fname);
	output_file << "pos_ref,pos_2" << endl;
	int row=0;
	for (int col=0; col<LB+1; col++)
	{
		for (int insertion=0; insertion<res_alignment_insert_before[col]; insertion++)
		{
			cout << row << " " << "-" << endl;
			output_file << row << "," << "-" << endl;
			row+=1;
		}
		if (col<LB)
		{
			if (res_alignment[col]!=-1)
			{
				row=res_alignment[col];
				cout << row << " " << col << endl;
				output_file << row << "," << col << endl;
				row+=1;
			}
			else
			{
				cout << "-" << " " << col << endl;
				output_file << '-' << "," << col << endl;
			}
		}

	}

	output_file.close();
}


int count_edges(int** edges_map, int L)
{
	int nb_edges=0;
	for (int i=0; i<L; i++)
	{
		for (int j=0; j<L; j++)
		{
			nb_edges+=edges_map[i][j];
		}
	}
	return nb_edges;
}


void display_results_and_print_to_files(int** row_map, int** col_map, double self1, double self2, double res_lb, double res_ub, double total_time, double res_alloc_time, double res_solve_time, int nb_bb_nodes, int* res_alignment, int* res_alignment_insert_before, int disp_level, char* aln_fname, char* info_fname, char* aln_with_gaps_fname, int status)
{
	if (status!=NOT_SIMILAR)
	{
		cout << endl << "RESULT:" << endl;
		if(disp_level >= 1)
		{
			cout << "      |N1| = " << LA <<"\n";
			cout << "      |N2| = " << LB <<"\n";
			cout << "      |E1| = " << count_edges(row_map, LA) <<"\n";
			cout << "      |E2| = " << count_edges(col_map, LB) <<"\n";
			cout << "      UB = " << -res_lb <<"\n";
			cout << "      LB = " << -res_ub <<"\n";
			cout << "      Similarity_global = " << 2.0*-res_ub/(self1+self2) <<"\n";
			cout << "      Similarity_global_ub = " << 2.0*-res_lb/(self1+self2) <<"\n";
			cout << "      Similarity_local = " << -res_ub/MIN(self1,self2) <<"\n";
			cout << "      T(sec) = " << total_time <<"\n";
			cout << "      Allocation (sec) = " << res_alloc_time << endl;
			cout << "      Solving (sec) = " << res_solve_time << endl;
			cout << "      Number visited nodes = " << nb_bb_nodes << endl;
		}
		cout << "\n";

		display_alignment(res_alignment);
		display_aligned_nodes(res_alignment, row_map, col_map, aln_fname);
		display_alignment_with_gaps(res_alignment, res_alignment_insert_before, aln_with_gaps_fname);

		ofstream output_file;
		output_file.open(info_fname);
		output_file << "similarity_global,solver_time,UB,LB,nb_visited_nodes,selfcomp1,selfcomp2" << endl;
		output_file << 2.0*-res_ub/(self1+self2) << "," << total_time << "," << -res_lb << "," << -res_ub << "," << nb_bb_nodes << "," << self1 << "," << self2 << endl;
		output_file.close();
	}
	else
	{
		cout << "Potts models are not similar enough." << endl;
		ofstream output_file;
		output_file.open(info_fname);
		output_file << "similarity_global,solver_time,UB,LB,nb_visited_nodes,selfcomp1,selfcomp2" << endl;
		output_file << "nan" << "," << total_time << "," << "nan" << "," << "nan" << "," << nb_bb_nodes << "," << self1 << "," << self2 << endl;
		output_file.close();
	}

}



int solve_prb(int ** forbidden, int * sol, int * sol_insert_before, double &alloc_time, double &solve_time, double &ub, double &lb, int& nb_bb_nodes, int** row_map, int** col_map, double self1, double self2, int iter_limit_param, int n_limit_param, double t_limit, double epsilon, double gamma, double theta, double stepsize_min, int nb_non_increasing_steps_max, double score_min, float* insert_open_A, float* insert_open_B, float* insert_extend_A, float* insert_extend_B, double dalih_bound=0.0)
{
    int status(0);

    float (*score_vertex)(int,int);
    score_vertex = &f_vertex_mrf;
    float (*score_edge)(int,int,int,int);
    score_edge = &f_edge_mrf;

    int nb_row = LA;
    int nb_col = LB;


    //start time
    long tic_per_sec = sysconf(_SC_CLK_TCK);
    struct tms start, end;
    int tic1(times(&start));

    parameters p    = parameters();
    p.max_iteration = iter_limit_param;
    p.max_node      = n_limit_param;
    p.time_limit    = t_limit;
    p.epsilon = epsilon;
    p.gamma = gamma;
    p.theta = theta;
    p.stepsize_min = stepsize_min;
    p.nb_non_increasing_steps_max = nb_non_increasing_steps_max;
    p.score_min = score_min;


    //create problem
    graph_apurva        g(nb_row, row_map, nb_col, col_map, score_vertex, score_edge, forbidden);
    lambda_mat_apurva   lm(g);
    dp_mat_apurva       dp(g, insert_open_A, insert_open_B, insert_extend_A, insert_extend_B);
    problem_apurva      prb(g,dp,lm,self1,self2,2*score_min/(self1+self2));
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
    bandb.get_solution(sol, sol_insert_before);
    //get solving time
    solve_time = bandb.get_solve_time();
    nb_bb_nodes = bandb.get_nb_visited_nodes();

    return(status);
}


// 1D array given to cython -> 2D array
int** unflatten(int* flat_array, int length)
{
	int** uf = new int*[length];
	for (int i=0; i<length; i++)
	{
		uf[i] = new int[length];
		for (int j=0; j<length; j++)
		{
			uf[i][j] = flat_array[i*length+j];
		}
	}
	return uf;
}


// free 2D array memory
void free_2d_int_array(int** array_2d, int length)
{
	for (int i=0; i<length; i++)
	{
		delete[] array_2d[i];
	}
    	delete[] array_2d;

}


extern "C" int call_from_python(float* v_A_, float* v_B_, float* w_A_, float* w_B_, int LA_, int LB_, int* edges_mapA, int* edges_mapB, double self1, double self2, float* insert_open_A_, float* insert_open_B_, float*insert_extend_A_, float* insert_extend_B_, char* aln_fname, char* info_fname, char* aln_with_gaps_fname, int n_limit_param, int iter_limit_param, double t_limit, int disp_level, double epsilon, double gamma, double theta, double stepsize_min, int nb_non_increasing_steps_max, double score_min, double alpha_w_, double offset_v_)
{
	int status(0);

	v_A = v_A_;
	v_B = v_B_;
	w_A = w_A_;
	w_B = w_B_;
	LA = LA_;
	LB = LB_;
	//gap_open = gap_open_;
	//gap_extend = gap_extend_;
	alpha_w = alpha_w_;
	offset_v = offset_v_;

	insert_open_A = insert_open_A_;
	insert_open_B = insert_open_B_;
	insert_extend_A = insert_extend_A_;
	insert_extend_B = insert_extend_B_;




	// computation time
	long tic_per_sec = sysconf(_SC_CLK_TCK);
	struct tms start, end;
	int tic1(times(&start));


	// alloc residues contact map filter
	int** forbidden_res = new int* [LB];
    	for(int col(0); col != LB; ++col)
    	{
        	forbidden_res[col] = new int [LA];
        	//by default all nodes are allowed
        	fill_n(forbidden_res[col], LA, static_cast<int>(0));
    	}


	// alloc solution arrays
	int* res_alignment = new int[LB];
	for(int col(0); col != LB; ++col)
    	{
                res_alignment[col] = -1;
    	}

	int* res_alignment_insert_before = new int[LB+1];
	for(int col(0); col != LB+1; ++col)
    	{
                res_alignment_insert_before[col] = 0;
    	}

	// solve
	int nb_bb_nodes = 0;
	double res_ub = 0;
	double res_lb = 0;
	double res_alloc_time(0.);
	double res_solve_time(0.);

	int** row_map = unflatten(edges_mapA, LA);
	int** col_map = unflatten(edges_mapB, LB);

	status = solve_prb(forbidden_res, res_alignment, res_alignment_insert_before, res_alloc_time, res_solve_time, res_ub, res_lb, nb_bb_nodes, row_map, col_map, self1, self2, iter_limit_param, n_limit_param, t_limit, epsilon, gamma, theta, stepsize_min, nb_non_increasing_steps_max, score_min, insert_open_A, insert_open_B, insert_extend_A, insert_extend_B, -INFINITY);


	// computation time
	int tic2 = times(&end);
	total_time = ((double)tic2 - (double)tic1) / (double)tic_per_sec;


	display_results_and_print_to_files(row_map, col_map, self1, self2, res_lb, res_ub, total_time, res_alloc_time, res_solve_time, nb_bb_nodes, res_alignment, res_alignment_insert_before, disp_level, aln_fname, info_fname, aln_with_gaps_fname, status);


	for(int col(0); col != LB; ++col)
	{
		delete[] forbidden_res[col];
    	}
    	delete[] forbidden_res;
    	delete[] res_alignment;

	free_2d_int_array(row_map, LA);
	free_2d_int_array(col_map, LB);
	
	return 0;
}


