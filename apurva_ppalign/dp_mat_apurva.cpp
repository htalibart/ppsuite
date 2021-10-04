#include "main.h" 
#include <omp.h>
#include <string.h>

using namespace std;

/* Compute the best sum of outgoing edge of node col1.row1
* -> deadly optimized
*/
double dp_mat_apurva :: value_arcs_out(int nb_next_col, int nb_next_row, int * next_col, int * next_row, double * coef_lambda_edge)
{
    int ind_col2, ind_row2;
    double value, vert, diag;

    //double dp_arc_tmp[nb_next_col+1][nb_next_row+1];
    //memset(dp_arc_tmp, 0, (nb_next_col+1)*(nb_next_row+1)*sizeof(double));

    //Main dynamic programming loop
    for(ind_col2 = 0; ind_col2 != nb_next_col; ++ind_col2)
    {
        //col2 = next_col[ind_col2];
        //For each node col2,row2, col2 in next_col(col1), row2 in next_row(row1)
        for (ind_row2 = 0; ind_row2 != nb_next_row; ++ind_row2)
        {
            // -> horizontal move
            value = dp_arc[ind_col2][ind_row2+1];
            // -> vertical move
            vert = dp_arc[ind_col2+1][ind_row2];
            if (vert < value)
                value = vert;
            // -> diagonal move, if col1,row1,col2,row2 is a valid edge
            diag = dp_arc[ind_col2][ind_row2] + coef_lambda_edge[ind_col2*nb_next_row + ind_row2];
            if(diag <= value)//< before 04/01/2012
                value = diag;
            dp_arc[ind_col2+1][ind_row2+1] = value;
        }
    }
    return(dp_arc[nb_next_col][nb_next_row]);
}


void dp_mat_apurva :: calc_arcs_out(int col1, int row1, graph_apurva & g, lambda_mat_apurva & lb_mat, int *lo, int *up)
{
    double value, vert, diag;
    int col2, row2, from, move;

    double ** coef_lambda_edge = lb_mat.get_coef_lambda_edge();
    int ind_node1 = g.get_ind_from_node(col1,row1);
    int nb_next_col = g.get_nb_next_col(col1, row1);
    int nb_next_row = g.get_nb_next_row(col1, row1);

    //Main dynamic programming loop same as in the upper function
    for(int ind_col2(0); ind_col2 != nb_next_col; ++ind_col2)
    {
        col2 = g.get_next_col(col1,row1,ind_col2);
        for (int ind_row2(0); ind_row2 != nb_next_row; ++ind_row2)
        {
            row2 = g.get_next_row(col1,row1,ind_row2);
            value = dp_arc[ind_col2][ind_row2+1];
            move = 0;
            vert = dp_arc[ind_col2+1][ind_row2];
            if (vert < value)
            {
                value = vert; move = -1;
            }
            if(row2 <= up[col2] && row2 >= lo[col2] && g.is_a_node(col2,row2))
            {
                diag = dp_arc[ind_col2][ind_row2] + coef_lambda_edge[ind_node1][ind_col2*nb_next_row + ind_row2];
                if(diag <= value) //< before 04/01/2012
                {
                    value = diag; move = 1;
                }
            }
            dp_arc[ind_col2+1][ind_row2+1] = value;
            arc_from[ind_col2+1][ind_row2+1] = move;
        }
    }

    /*
    * Retrieve solution in dp_arc_out
    * Backward dynamic programming
    */
    int j(nb_next_row);
    row2 = g.get_next_row(col1, row1, j-1);

    for(int i(nb_next_col); i != 0; --i)
    {
        from = arc_from[i][j];
        while(from == -1)
        {
            j--;
            from = arc_from[i][j];
        }

        if (from == 0)
        {
            dp_arc_out_col[col1][row1][i-1] = -1;
        }
        else
        {
            row2 = g.get_next_row(col1, row1, j-1);
            dp_arc_out_col[col1][row1][i-1] = row2;
            j--;
        }
    }

}

/**************************************************
 * FILL                                           *
 **************************************************/

void dp_mat_apurva :: fill(graph_apurva & g, lambda_mat_apurva & lb_mat, int * lo, int * up)
{
    /*
    long tic_per_sec = sysconf(_SC_CLK_TCK);
    struct tms start, end;
    int tic1(times(&start));
    */

    //get data once for all (faster)
    int nb_col = g.get_nb_col();
    //int nb_row(g.get_nb_row());
    int * nb_next_col(g.get_nb_next_col());
    int * nb_next_row(g.get_nb_next_row());
    int ** next_col(g.get_next_col());
    int ** next_row(g.get_next_row());
    int ** node_to_ind(g.get_node_to_ind());
    double * coef_lambda_node(lb_mat.get_coef_lambda_node());
    double ** coef_lambda_edge(lb_mat.get_coef_lambda_edge());

    int ** lambda_modified_nodes(lb_mat.get_lambda_modified_nodes());

    int col1;
    
    #pragma omp parallel for schedule(dynamic), default(none), num_threads(NUM_THREADS), private(col1), shared(nb_col, lo, up, node_to_ind, \
 next_row, next_col, lambda_modified_nodes, nb_next_col, nb_next_row, coef_lambda_edge, coef_lambda_node)
    for(col1 = 0; col1 < nb_col; ++col1)
    {
    	double score;
    	int ind_node1;
    	int i_nb_next_col, i_nb_next_row;

      // optimization : dont loop over every cell
      // thus cells out of range (lo - up) are not set to 0
      // but it was redundant as dp_mat_apurva :: solve already checks for out of range cells
      //for (int row1(0); row1 != nb_row; ++row1)
        for (int row1(lo[col1]); row1 <= up[col1]; ++row1)
        {
            //if col1.row1 is a valid node
            ind_node1 = node_to_ind[col1][row1];


            if( ind_node1 >= 0)
            {
                if(lambda_modified_nodes[col1][row1] != 0)
                {
                    //We need to recompute the score of node col1.row1

                    //Score of node col1,row1   (in apurva, S_col1.row1 = 0.)
                    score = 0;
                    i_nb_next_col = nb_next_col[col1];
                    i_nb_next_row = nb_next_row[row1];
                    //adding best sum of outgoing lambda-weighted edges
                    if(i_nb_next_col >= 1 && i_nb_next_row >=1)
                    {
                        score += value_arcs_out(i_nb_next_col, i_nb_next_row, next_col[col1], next_row[row1], coef_lambda_edge[ind_node1]);
                    }
                    //printf("%f \n",score);
                    //adding sum of lambda_node_row and sum of lambda_node_col, for node col1.row1, and all pred col, and all pred row
                    score -= coef_lambda_node[ind_node1];

                    dp_score[col1][row1] = score;

                    //Since score has been recomputed, setting node in order not to recompute it
                    lambda_modified_nodes[col1][row1] = 0;
                }
            }
            else  //node should not be taken
                dp_score[col1][row1] = 0.;
        }
    }

    /*
	//print dp_score
	for (int col_print=0; col_print<nb_col; col_print++)
	{
		for (int row_print=0; row_print<nb_row; row_print++)
		{
			cout << col_print << row_print << " " << dp_score[col_print][row_print] << endl;
		}
	}
*/


    /*
    int tic2 = times(&end);
    double solve_time = ((double)tic2 - (double)tic1) / (double)tic_per_sec;
    cout << "DP_fill done in " << solve_time << "(s)\n";
    */
}




double dp_mat_apurva :: solve_w_gapcosts(graph_apurva & g, int * sol, int * sol_insert_before, lambda_mat_apurva & lb_mat, int * lo, int * up)
{

    int nb_col(g.get_nb_col());
    int nb_row(g.get_nb_row());
    int ** node_to_ind(g.get_node_to_ind());

    double prev_M;
    double prev_GA;
    double prev_GB;

    for (int i=1; i<=nb_row; ++i)
	{
		for (int k=1; k<=nb_col; ++k)
		{
			// compute M
			prev_M = dp_M[k-1][i-1]+dp_score[k-1][i-1];
			prev_GA = dp_GA[k-1][i-1]+dp_score[k-1][i-1];
			prev_GB = dp_GB[k-1][i-1]+dp_score[k-1][i-1];

			if ( (prev_M<=prev_GA) && (prev_M<=prev_GB) && (lo[k-1] <= i-1 && up[k-1] >= i-1 && node_to_ind[k-1][i-1] >= 0) )
			{
				dp_M[k][i] = prev_M;
				dp_M_from[k][i] = 1;
			}
			else if ( (prev_GA<=prev_M) && (prev_GA<=prev_GB) )
			{
				dp_M[k][i] = prev_GA;
				dp_M_from[k][i] = -1;
			}
			else
			{
				dp_M[k][i] = prev_GB;
				dp_M_from[k][i] = 0;
			}


			// compute GA
			prev_M = dp_M[k-1][i]+get_insertion_open_after_row(i-1);
			prev_GA = dp_GA[k-1][i]+get_insertion_extend_after_row(i-1);
			prev_GB = dp_GB[k-1][i]+get_insertion_open_after_row(i-1);

			if ( (prev_M<=prev_GA) && (prev_M<=prev_GB) && (lo[k-1] <= i-1 && up[k-1] >= i-1 && node_to_ind[k-1][i-1] >= 0) )
			{
				dp_GA[k][i] = prev_M;
				dp_GA_from[k][i] = 1;
			}
			else if ( (prev_GA<=prev_M) && (prev_GA<=prev_GB) )
			{
				dp_GA[k][i] = prev_GA;
				dp_GA_from[k][i] = -1;
			}
			else
			{
				dp_GA[k][i] = prev_GB;
				dp_GA_from[k][i] = 0;
			}


			// compute GB
			prev_M = dp_M[k][i-1]+get_insertion_open_after_col(k-1);
			prev_GA = dp_GA[k][i-1]+get_insertion_open_after_col(k-1);
			prev_GB = dp_GB[k][i-1]+get_insertion_extend_after_col(k-1);

			if ( (prev_M<=prev_GA) && (prev_M<=prev_GB) && (lo[k-1] <= i-1 && up[k-1] >= i-1 && node_to_ind[k-1][i-1] >= 0) )
			{
				dp_GB[k][i] = prev_M;
				dp_GB_from[k][i] = 1;
			}
			else if ( (prev_GB<=prev_M) && (prev_GB <= prev_GA) )
			{
				dp_GB[k][i] = prev_GB;
				dp_GB_from[k][i] = 0;
			}
			else
			{
				dp_GB[k][i] = prev_GA;
				dp_GB_from[k][i] = -1;
			}

			/*
			cout << "dp_M["<<k<<"]["<<i<<"]="<<dp_M[k][i]<<endl;
			cout << "dp_GA["<<k<<"]["<<i<<"]="<<dp_GA[k][i]<<endl;
			cout << "dp_GB["<<k<<"]["<<i<<"]="<<dp_GB[k][i]<<endl;
			*/
		}
	}

    	double dp_max;
	int current_tr;
	if ((dp_M[nb_col][nb_row]<=dp_GA[nb_col][nb_row]) && (dp_M[nb_col][nb_row]<=dp_GB[nb_col][nb_row]))
	{
		dp_max=dp_M[nb_col][nb_row];
		current_tr=1;
	}
	else if ((dp_GA[nb_col][nb_row]<=dp_M[nb_col][nb_row]) && (dp_GA[nb_col][nb_row]<=dp_GB[nb_col][nb_row]))
	{
		dp_max=dp_GA[nb_col][nb_row];
		current_tr=-1;
	}
	else
	{
		dp_max=dp_GB[nb_col][nb_row];
		current_tr=0;
	}

	//cout << "dp_max=" << dp_max << endl;


	// reset sol_insert_before
	for (int col=0; col<nb_col+1; col++)
	{
		sol_insert_before[col]=0;
	}

	//Traceback
	int i=nb_row;
	int k=nb_col;

	while ( (i>0) && (k>0) )
	{
		//cout << "i="<<i<<",k="<<k<<endl;
		//cout << "current_tr=" << current_tr << endl;
		if (current_tr==1)
		{
			sol[k-1]=i-1;
			current_tr=dp_M_from[k][i];
			i--;
			k--;
		}
		else if (current_tr==0)
		{
			sol_insert_before[k]+=1;
			current_tr=dp_GB_from[k][i];
			i--;
		}
		else
		{
			sol[k-1]=-1;
			current_tr=dp_GA_from[k][i];
			k--;
		}
	}


	    // gaps in the beginning
	    while(i>0)
	    {
		sol_insert_before[k]+=1;
		--i;
	    }
	    while(k>0)
	    {
		sol[k-1] = -1;
		--k;
	    }

	    //Show solution
	    /*
	    cout << "Solution UB:" << endl;
	    for(int i=0; i<nb_col; ++i)
	    {
		cout << i << " " << sol[i] << endl;
	    }
	    cout << "Solution insert before UB:" << endl;
	    for(int i=0; i<nb_col+1; ++i)
	    {
		cout << i << " " << sol_insert_before[i] << endl;
	    }
	    */


	//cout << "( dp_max - lb_mat.get_sum_lambda_act() ) = " << ( dp_max - lb_mat.get_sum_lambda_act() ) <<endl;
	return( dp_max - lb_mat.get_sum_lambda_act() );
}

