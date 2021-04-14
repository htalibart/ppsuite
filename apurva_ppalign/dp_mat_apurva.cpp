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
    int tic2 = times(&end);
    double solve_time = ((double)tic2 - (double)tic1) / (double)tic_per_sec;
    cout << "DP_fill done in " << solve_time << "(s)\n";
    */
}

/**************************************************
 * SOLVE                                          *
 **************************************************/

double dp_mat_apurva :: solve(graph_apurva & g, int * sol, lambda_mat_apurva & lb_mat, int * lo, int * up)
{
    int nb_col(g.get_nb_col());
    int nb_row(g.get_nb_row());
    int ** node_to_ind(g.get_node_to_ind());

    double horiz, diag, vert;

    // First node, col0.row0
    if(dp_score[0][0] <= 0. && lo[0] <= 0 && up[0] >= 0 && node_to_ind[0][0] >= 0) //before dp_score[0][0] < 0.
    {
        dp[0][0] = dp_score[0][0];
        dp_from[0][0] = 1;
    }
    else
    {
        dp[0][0] = 0.;
        dp_from[0][0] = 0;
    }

    // First row
    for(int col(1); col != nb_col; ++col)
    {
        horiz = dp[col-1][0];
        diag = dp_score[col][0];
        if(diag <= horiz && lo[col] <= 0 && up[col] >= 0 && node_to_ind[col][0] >= 0) //was diag < horiz before 04/01/2012
        {
            dp[col][0] = diag;
            dp_from[col][0] = 1;
        }
        else
        {
            dp[col][0] = horiz;
            dp_from[col][0] = 0;
        }
    }

    // First col
    for(int row(1); row != nb_row; ++row)
    {
        vert = dp[0][row-1];
        diag = dp_score[0][row];
        if(diag <= vert && lo[0] <= row && up[0] >= row && node_to_ind[0][row] >= 0) //was diag < vert before 04/01/2012
        {
            dp[0][row] = diag;
            dp_from[0][row] = 1;
        }
        else
        {
            dp[0][row] = vert;
            dp_from[0][row] = -1;
        }
    }

    // Dynamic programming
    for(int col(1); col != nb_col; ++col)
    {
        for(int row(1); row != nb_row; ++row)
        {
            horiz = dp[col-1][row];
            vert = dp[col][row-1];
            diag = dp[col-1][row-1] + dp_score[col][row];

            if (vert < horiz)
            {
                dp[col][row] = vert;
                dp_from[col][row] = -1;
            }
            else
            {
                dp[col][row] = horiz;
                dp_from[col][row] = 0;
            }

            if(diag <= dp[col][row] && lo[col] <= row && up[col] >= row && node_to_ind[col][row] >= 0) //was diag < dp[col][row]
            {
                dp[col][row] = diag;
                dp_from[col][row] = 1;
            }
        }
    }

    //clear sol
    for(int i(0); i != nb_col; ++i)
        sol[i] = -1;

    //save sol
    int row(nb_row - 1), col(nb_col - 1), from;
    while(row >= 0 && col >= 0)
    {
        from = dp_from[col][row];
        if(from == -1)
            --row;
        else if(from == 0)
        {
            sol[col] = -1;
            --col;
        }
        else
        {
        	sol[col] = row;
            --col;
            --row;
        }
    }

    /*
    cout << "Relaxed solution : " << dp[nb_col - 1][nb_row - 1] << "\n";
    */
    //add the constant factor that stems from the activation constraints
    return( dp[nb_col - 1][nb_row - 1] - lb_mat.get_sum_lambda_act() );
}




double dp_mat_apurva :: solve_w_gapcosts(graph_apurva & g, int * sol, lambda_mat_apurva & lb_mat, int * lo, int * up)
{

    int nb_col(g.get_nb_col());
    int nb_row(g.get_nb_row());
    int ** node_to_ind(g.get_node_to_ind());

    // Initialization

    //Traceback symbols
    //diagonal: 1 (symbols aligned)
    //horizontal: 0 (gap in first sequence), Q (=dp_h)
    //vertical: -1 (gap in second sequence), P (=dp_v)

    //Dynamic programming
    double prev_aligned;
    double prev_gap_in_second;
    double prev_gap_in_first;
    for(int i=1; i<=nb_col; ++i) // second sequence
    {
    	for(int j=1; j<=nb_row; ++j) // first sequence
    	{
    		// compute P (alignment ends with gap in second sequence) 
    		prev_aligned = dp[i-1][j] + get_insertion_open_after_row(j-1);
    		prev_gap_in_second = dp_v[i-1][j] + get_insertion_extend_after_row(j-1);
    		prev_gap_in_first = dp_h[i-1][j] + get_insertion_open_after_row(j-1);
    		if(prev_aligned <= prev_gap_in_first && prev_aligned <= prev_gap_in_second && lo[i-1] <= j-1 && up[i-1] >= j-1 && node_to_ind[i-1][j-1] >= 0) //if values are identical, priority is to align
    		{
    			dp_v[i][j] = prev_aligned;
    			dp_v_from[i][j] = 1;
    		}
    		else if(prev_gap_in_second <= prev_gap_in_first) //if values are identical, priority is to extend a gap
    		{
    			dp_v[i][j] = prev_gap_in_second; //extend gap in second
    			dp_v_from[i][j] = -1;

    		}
    		else
    		{
    			dp_v[i][j] = prev_gap_in_first; //make new gap in first
    			dp_v_from[i][j] = 0;
    		}

    		// compute Q (alignment ends with gap in first sequence)
    		prev_aligned = dp[i][j-1] + get_insertion_open_after_col(i-1);
    		prev_gap_in_second = dp_v[i][j-1] + get_insertion_open_after_col(i-1);
    		prev_gap_in_first = dp_h[i][j-1] + get_insertion_extend_after_col(i-1);
    		if(prev_aligned <= prev_gap_in_first && prev_aligned <= prev_gap_in_second && lo[i-1] <= j-1 && up[i-1] >= j-1 && node_to_ind[i-1][j-1] >= 0) //if values are identical, priority is to align
    		{
    			dp_h[i][j] = prev_aligned;
    			dp_h_from[i][j] = 1;
    		}
    		else if(prev_gap_in_first < prev_gap_in_second) //if values are identical, priority is to put a gap in the second sequence
    		{
    			dp_h[i][j] = prev_gap_in_first; //extend gap in first
    			dp_h_from[i][j] = 0;
    		}
    		else
    		{
    			dp_h[i][j] = prev_gap_in_second; //make new gap in second
    			dp_h_from[i][j] = -1;
    		}


    		// compute D (alignment ends with two aligned characters)
    		prev_aligned = dp[i-1][j-1] + dp_score[i-1][j-1];
    		prev_gap_in_second = dp_v[i-1][j-1] + dp_score[i-1][j-1];
    		prev_gap_in_first = dp_h[i-1][j-1] + dp_score[i-1][j-1];
    		if(prev_aligned <= prev_gap_in_first && prev_aligned <= prev_gap_in_second && lo[i-1] <= j-1 && up[i-1] >= j-1 && node_to_ind[i-1][j-1] >= 0) //if values are identical, priority is to align
    		{
    			dp[i][j] = prev_aligned;
    			dp_from[i][j] = 1;
    		}
    		else if(prev_gap_in_first < prev_gap_in_second) //close gap in first seq; if values are identical, priority is to put gap in second
    		{
    			dp[i][j] = prev_gap_in_first;
    			dp_from[i][j] = 0;
    		}
    		else //close gap in second seq
    		{
    			dp[i][j] = prev_gap_in_second;
    			dp_from[i][j] = -1;
    		}
		
    	}
    }


    string current_tr;
    double dp_max;
    if ( (dp[nb_col][nb_row]<=dp_v[nb_col][nb_row]) && (dp[nb_col][nb_row]<=dp_h[nb_col][nb_row]) )
	{
		current_tr = "D";
		dp_max = dp[nb_col][nb_row];
	}
    else if ( (dp_v[nb_col][nb_row]<=dp[nb_col][nb_row]) && (dp_v[nb_col][nb_row]<=dp_h[nb_col][nb_row]) )
	{
		current_tr = "P";
		dp_max = dp_v[nb_col][nb_row];
	}
    else if ( (dp_h[nb_col][nb_row]<=dp[nb_col][nb_row]) && (dp_h[nb_col][nb_row]<=dp_v[nb_col][nb_row]) )
	{
		current_tr = "Q";
		dp_max = dp_h[nb_col][nb_row];
	}
    else
	{
		current_tr="IMPOSSIBLE";
	}



    //Traceback
    int i = nb_col;
    int j = nb_row;
    while(i>0 && j>0)
    {

    	if(current_tr == "D") //the two characters are aligned
    	{
    		sol[i-1] = j-1;

    		if(dp_from[i][j] == -1) //previous there was a gap in the second sequence
    		{
    			current_tr = "P";
    		}
    		else //previous there was a gap in the first sequence
    		{
    			current_tr = "Q";
    		}
    			--i;
    			--j;

    	}
    	else  if(current_tr == "P") //there is a gap in the second sequence
    	{

    		sol[i-1] = -1;

    		if(dp_v_from[i][j] == 1) //previous characters were aligned
    		{
    			current_tr = "D";
    		}
    		else if(dp_v_from[i][j] == -1) //previous there was a gap in the second sequence that was extended
    		{
    		}
    		else //previous there was a gap in the first sequence (an insertion followed by a deletion)
    		{
    			current_tr = "Q";
    		}
    		--i;
    	}
    	else if(current_tr == "Q") //there is a gap in the first sequence
    	{

    		if(dp_h_from[i][j] == 1) //previous two characters were aligned
    		{
    			current_tr = "D";
    		}
    		else if(dp_h_from[i][j] == 0) //previous there was also a gap in the first sequence (gap extension)
    		{
    		}
    		else //prevous there was a gap in the second sequence (an insertion followed by a deletion)
    		{
    			current_tr = "P";
    		}
    		--j;
    	}

    }
    //gaps in the beginning in first sequence
    while(j>0)
    {
    	--j;
    }
    //gaps in the beginning in second sequence
    while(i>0)
    {
    	sol[i-1] = -1;
    	--i;
    }

    //Show solution
    
    /*cout << "Solution:" << endl;
    for(int i=0; i<nb_col; ++i)
    {
    	cout << i << " " << sol[i] << endl;
    }*/


    //add the constant factor that stems from the activation constraints
    return( dp_max - lb_mat.get_sum_lambda_act() );
}
