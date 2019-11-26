#include "main.h"
using namespace std;
#include <unistd.h>

/**
* Check if constraint col1, col2, row2 is broken.
* The important point is that sum(Y) is either 0 or 1,
* so if one Y is activated, there is no need to check the other ones
**/
inline int problem_apurva :: is_brk_cst_col_node(int col1, int col2, int row2, int * sol)
{
    int nb_pred_col(g.get_nb_pred_col(col2,row2));
    int nb_pred_row(g.get_nb_pred_row(col2,row2));
    int row0(g.get_pred_row(col2,row2,0));
    int rowl(g.get_pred_row(col2,row2,nb_pred_row - 1));
    int ind_col1 = g.get_ind_col_col(col2,col1);

    int Sum_Y_edge(0);
    int X_node(0);

    //X_col_row   is the node taken in the solution ?
    if(sol[col2] == row2)
    {
        X_node = 1;
    }


    //Sum_Y_edge
    //->upper_part
    for(int ind_coln(nb_pred_col - 1); ind_coln > ind_col1; --ind_coln)
    {
        int coln(g.get_pred_col(col2,row2,ind_coln));
        if(sol[coln] == row0)
        {
            int ind_col2(g.get_ind_col_col(coln,col2));
            if( dp.get_arc_out(coln,row0,ind_col2) == row2)
            {
                Sum_Y_edge=1;
                ind_coln = ind_col1;
            }
        }
    }

    //->col_part
    if(Sum_Y_edge == 0)
    {
        for(int ind_rown(0); ind_rown < nb_pred_row; ++ind_rown)
        {
            int rown(g.get_pred_row(col2,row2,ind_rown));
            if(sol[col1] == rown)
            {
                int ind_col2(g.get_ind_col_col(col1,col2));
                if( dp.get_arc_out(col1,rown,ind_col2) == row2)
                {
                    Sum_Y_edge=1;
                    ind_rown = nb_pred_row;
                }
            }
        }
    }
    //->lower_part
    if(Sum_Y_edge == 0)
    {
        for(int ind_coln(ind_col1 - 1); ind_coln >= 0; --ind_coln)
        {
            int coln(g.get_pred_col(col2,row2,ind_coln));
            if(sol[coln] == rowl)
            {
                int ind_col2(g.get_ind_col_col(coln,col2));
                if( dp.get_arc_out(coln,rowl,ind_col2) == row2)
                {
                    Sum_Y_edge=1;
                    ind_coln = -1;
                }
            }
        }
    }

    return(X_node - Sum_Y_edge);
}




/**
* Check if constraint row1, col2, row2 is broken.
* The important point is that sum(Y) is either 0 or 1,
* so if one Y is activated, there is no need to check the other ones
**/
inline int problem_apurva :: is_brk_cst_row_node(int row1, int col2, int row2, int * sol)
{
    int nb_pred_col(g.get_nb_pred_col(col2,row2));
    int nb_pred_row(g.get_nb_pred_row(col2,row2));
    int col0(g.get_pred_col(col2,row2,0));
    int coll(g.get_pred_col(col2,row2,nb_pred_col - 1));
    int ind_row1 = g.get_ind_row_row(row2,row1);

    int Sum_Y_edge(0);
    int X_node(0);

    //X_col_row   is the node taken in the solution ?
    if(sol[col2] == row2)
    {
        X_node = 1;
    }
    //Sum_Y_edge
    //-> upper_part
    for(int ind_rown(0); ind_rown < ind_row1; ++ind_rown)
    {
        int rown(g.get_pred_row(col2,row2,ind_rown));
        if(sol[coll] == rown)
        {
            int ind_col2(g.get_ind_col_col(coll,col2));
            if( dp.get_arc_out(coll,rown,ind_col2) == row2)
            {
                Sum_Y_edge=1;
                ind_rown = ind_row1;
            }
        }
    }
    //-> row_part
    if(Sum_Y_edge == 0)
    {
        for(int ind_coln(0); ind_coln < nb_pred_col; ++ind_coln)
        {
            int coln(g.get_pred_col(col2,row2,ind_coln));
            if(sol[coln] == row1)
            {
                int ind_col2(g.get_ind_col_col(coln,col2));
                if( dp.get_arc_out(coln,row1,ind_col2) == row2)
                {
                    Sum_Y_edge=1;
                    ind_coln = nb_pred_col;
                }
            }
        }
    }
    //->lower_part
    if(Sum_Y_edge == 0)
    {
        for(int ind_rown(ind_row1 + 1); ind_rown < nb_pred_row; ++ind_rown)
        {
            int rown(g.get_pred_row(col2,row2,ind_rown));
            if(sol[col0] == rown)
            {
                int ind_col2(g.get_ind_col_col(col0,col2));
                if( dp.get_arc_out(col0,rown,ind_col2) == row2)
                {
                    Sum_Y_edge=1;
                    ind_rown = nb_pred_row;
                }
            }
        }
    }

    return(X_node - Sum_Y_edge);
}

/**
* Check if activation constraint row1, col2, row2 is broken.
* The important point is that sum(Y) is either 0 or 1,
* so if one Y is activated, there is no need to check the other ones
* This function returns 0 if for given row1, col2, row2, x_ik=0, sum (y_rsik-x_rs)=0
**/
inline int problem_apurva :: is_brk_cst_row_node_act(int row1, int col2, int row2, int * sol)
{
    //X_col_row is the node taken in the solution?
	//If yes then this is no broken constraint of this type
    if(sol[col2] == row2)
    {
        return 0;
    }

    int nb_pred_col(g.get_nb_pred_col(col2,row2));
    int nb_pred_row(g.get_nb_pred_row(col2,row2));
    int col0(g.get_pred_col(col2,row2,0));
    int coll(g.get_pred_col(col2,row2,nb_pred_col - 1));
    int ind_row1 = g.get_ind_row_row(row2,row1);

    int Sum_XY_edge(0);

    //Sum_Y_edge
    //-> upper_part
    for(int ind_rown(0); ind_rown < ind_row1; ++ind_rown)
    {
        int rown(g.get_pred_row(col2,row2,ind_rown));
        if(sol[coll] == rown && g.get_edge_coef(coll,rown,col2,row2) >= 0.0)
        {
        	Sum_XY_edge -= 1;
            int ind_col2(g.get_ind_col_col(coll,col2));
            if( dp.get_arc_out(coll,rown,ind_col2) == row2)
            {
                Sum_XY_edge += 1;
                ind_rown = ind_row1;
            }
        }
    }
    //-> row_part
    if(Sum_XY_edge == 0)
    {
        for(int ind_coln(0); ind_coln < nb_pred_col; ++ind_coln)
        {
            int coln(g.get_pred_col(col2,row2,ind_coln));
            if(sol[coln] == row1 && g.get_edge_coef(coln,row1,col2,row2) >= 0.0)
            {
            	Sum_XY_edge -= 1;
                int ind_col2(g.get_ind_col_col(coln,col2));
                if( dp.get_arc_out(coln,row1,ind_col2) == row2)
                {
                    Sum_XY_edge += 1;
                    ind_coln = nb_pred_col;
                }
            }
        }
    }
    //->lower_part
    if(Sum_XY_edge == 0)
    {
        for(int ind_rown(ind_row1 + 1); ind_rown < nb_pred_row; ++ind_rown)
        {
            int rown(g.get_pred_row(col2,row2,ind_rown));
            if(sol[col0] == rown && g.get_edge_coef(col0,rown,col2,row2) >= 0.0)
            {
            	Sum_XY_edge -= 1;
                int ind_col2(g.get_ind_col_col(col0,col2));
                if( dp.get_arc_out(col0,rown,ind_col2) == row2)
                {
                    Sum_XY_edge += 1;
                    ind_rown = nb_pred_row;
                }
            }
        }
    }

    return(Sum_XY_edge+1);
}













/*
* check if a col constraint is broken
* -> return X - Sum(Y)
*/
inline int problem_apurva :: check_cst_col_node(int col1, int col2, int row2, int * sol)
{
    int nb_pred_col(g.get_nb_pred_col(col2,row2));
    int nb_pred_row(g.get_nb_pred_row(col2,row2));
    int row0(g.get_pred_row(col2,row2,0));
    int rowl(g.get_pred_row(col2,row2,nb_pred_row - 1));
    int ind_col1 = g.get_ind_col_col(col2,col1);

    int Sum_Y_edge(0);
    int X_node(0);

    //X_col_row   is the node taken in the solution ?
    if(sol[col2] == row2)
    {
        X_node = 1;
    }
    //Sum_Y_edge
    //->upper_part
    for(int ind_coln(nb_pred_col - 1); ind_coln != ind_col1; --ind_coln)
    {
        int coln(g.get_pred_col(col2,row2,ind_coln));
        if(sol[coln] == row0)
        {
            int ind_col2(g.get_ind_col_col(coln,col2));
            if( dp.get_arc_out(coln,row0,ind_col2) == row2)
                Sum_Y_edge++;
        }
    }
    //->col_part
    for(int ind_rown(0); ind_rown != nb_pred_row; ++ind_rown)
    {
        int rown(g.get_pred_row(col2,row2,ind_rown));
        if(sol[col1] == rown)
        {
            int ind_col2(g.get_ind_col_col(col1,col2));
            if( dp.get_arc_out(col1,rown,ind_col2) == row2)
                Sum_Y_edge++;
        }
    }
    //->lower_part
    for(int ind_coln(ind_col1 - 1); ind_coln >= 0; --ind_coln)
    {
        int coln(g.get_pred_col(col2,row2,ind_coln));
        if(sol[coln] == rowl)
        {
            int ind_col2(g.get_ind_col_col(coln,col2));
            if( dp.get_arc_out(coln,rowl,ind_col2) == row2)
                Sum_Y_edge++;
        }
    }

    if(Sum_Y_edge > 1)
    {
        cout << "Fatal error in check_cst_col_node\n";
        exit(0);
    }

    return(X_node - Sum_Y_edge);
}

/*
* check if a row constraint is broken
* -> return X - Sum(Y)
*/
inline int problem_apurva :: check_cst_row_node(int row1, int col2, int row2, int * sol)
{
    int nb_pred_col(g.get_nb_pred_col(col2,row2));
    int nb_pred_row(g.get_nb_pred_row(col2,row2));
    int col0(g.get_pred_col(col2,row2,0));
    int coll(g.get_pred_col(col2,row2,nb_pred_col - 1));
    int ind_row1 = g.get_ind_row_row(row2,row1);

    int Sum_Y_edge(0);
    int X_node(0);

    //X_col_row   is the node taken in the solution ?
    if(sol[col2] == row2)
    {
        X_node = 1;
    }
    //Sum_Y_edge
    //-> upper_part
    for(int ind_rown(0); ind_rown != ind_row1; ++ind_rown)
    {
        int rown(g.get_pred_row(col2,row2,ind_rown));
        if(sol[coll] == rown)
        {
            int ind_col2(g.get_ind_col_col(coll,col2));
            if( dp.get_arc_out(coll,rown,ind_col2) == row2)
                Sum_Y_edge++;
        }
    }
    //-> row_part
    for(int ind_coln(0); ind_coln != nb_pred_col; ++ind_coln)
    {
        int coln(g.get_pred_col(col2,row2,ind_coln));
        if(sol[coln] == row1)
        {
            int ind_col2(g.get_ind_col_col(coln,col2));
            if( dp.get_arc_out(coln,row1,ind_col2) == row2)
                Sum_Y_edge++;
        }
    }
    //->lower_part
    for(int ind_rown(ind_row1 + 1); ind_rown != nb_pred_row; ++ind_rown)
    {
        int rown(g.get_pred_row(col2,row2,ind_rown));
        if(sol[col0] == rown)
        {
            int ind_col2(g.get_ind_col_col(col0,col2));
            if( dp.get_arc_out(col0,rown,ind_col2) == row2)
                Sum_Y_edge++;
        }
    }

    if(Sum_Y_edge > 1)
    {
        cout << "Fatal error in check_cst_col_node\n";
        exit(0);
    }

    return(X_node - Sum_Y_edge);
}


/**************************************************
 * GET_CORRECT_VALUE                              *
 **************************************************/
double problem_apurva :: get_correct_value(int * sol)
{

    double score(0.);
    int nb_col(g.get_nb_col());
    int nb_row(g.get_nb_row());
    double vertex_score = 0.0;
    double edge_score = 0.0;
    double gap_cost = 0.0;
    bool extend_gap = false;
    int last_row = -1;;

    for (int col(0); col != nb_col; ++col)
    {
        int row(sol[col]);
        if(row >= 0)
        {
            int nb_next_col(g.get_nb_next_col(col, row));
            score -= g.get_node_coef(col,row);
            vertex_score -= g.get_node_coef(col,row);

            for (int ind_ncol(0); ind_ncol != nb_next_col; ++ind_ncol)
            {
                int ncol(g.get_next_col(col, row, ind_ncol));
                int nrow(sol[ncol]);
                if(nrow >=0)
                {
                    if(g.is_an_edge(col,row,ncol,nrow))
                    {
                    	score += g.get_edge_coef(col,row,ncol,nrow);
                    	edge_score += g.get_edge_coef(col,row,ncol,nrow);
                    }
                }
            }

            if(row - last_row > 1) // count the row symbols aligned to a gap
            {
        		score += dp.get_gap_open();
        		gap_cost += dp.get_gap_open();
        		score += (row-last_row-2)*dp.get_gap_extend();
        		gap_cost +=  (row-last_row-2)*dp.get_gap_extend();
            }

            extend_gap = false;
            last_row = row;
        }
        else // count the col symbols aligned to a gap
        {
        	if(extend_gap)
        	{
        		score += dp.get_gap_extend();
        		gap_cost += dp.get_gap_extend();
        	}
        	else
        	{
        		score += dp.get_gap_open();
        		gap_cost += dp.get_gap_open();
        		extend_gap = true;
        	}
        }
    }

    // gaps at the end
    if(nb_row - last_row > 1) // count the row symbols aligned to a gap
    {
		score += dp.get_gap_open();
		gap_cost += dp.get_gap_open();
		score += (nb_row-last_row-2)*dp.get_gap_extend();
		gap_cost +=  (nb_row-last_row-2)*dp.get_gap_extend();
    }

    return(score);
}

/**************************************************
 * GET_SUB_GR_NORM                                *
 **************************************************/
int problem_apurva :: get_nz_sub_gr(int * sol)
{
    int nz_c(0);
    int nz_r(0);
    int nz_r_act(0);

    int ** nz_sgv_c = lb_mat.get_nz_sgv_c();
    int ** nz_sgv_r = lb_mat.get_nz_sgv_r();
    int ** nz_sgv_r_act = lb_mat.get_nz_sgv_r_act();

    int nb_col(g.get_nb_col());

    //Step 0 : find selected outgoing edges
    for(int col1(0); col1 != nb_col; ++col1)
    {
        int row1(sol[col1]);
        if(row1 >= 0)
        {
            int nb_next_row(g.get_nb_next_row(col1,row1));
            int nb_next_col(g.get_nb_next_col(col1,row1));
            if(nb_next_col >= 1 && nb_next_row >= 1)
                dp.calc_arcs_out(col1, row1, g, lb_mat, lo, up);
        }
    }

    //Step 1 : Find equations : - lambda (X - Sum Y) < 0  (lambda should be decreased : maximization)
    //That is x_ik=1 sum y_rsik=0
    for(int col1(0); col1 != nb_col; ++col1)
    {
        int row1(sol[col1]);
        if(row1 >= 0)
        {
            //Node col1,row1 is activated (x = 1)
            int nb_pred_col(g.get_nb_pred_col(col1,row1));
            int nb_pred_row(g.get_nb_pred_row(col1,row1));

            //find which equations has a Sum Y = 0
            if(nb_pred_col >=1 && nb_pred_row >=1)
            {
                //col checker
                for(int ind_col(0); ind_col != nb_pred_col; ind_col++)
                {
                    int col0(g.get_pred_col(col1,row1,ind_col));

                    if(lb_mat.get_node_col(col0,col1,row1) > 0.0)
                    {
                    	if(is_brk_cst_col_node(col0, col1, row1, sol) != 0)
                    	{
                            nz_sgv_c[nz_c][0] = col0;
                            nz_sgv_c[nz_c][1] = col1;
                            nz_sgv_c[nz_c][2] = row1;
                            nz_sgv_c[nz_c][3] = -1;
                            ++nz_c;
                        }
                    }
                }
                //row checker
                for(int ind_row(0); ind_row != nb_pred_row; ind_row++)
                {
                    int row0(g.get_pred_row(col1,row1,ind_row));

                    if(lb_mat.get_node_row(row0,col1,row1) > 0.)
                    {
                    	if(is_brk_cst_row_node(row0, col1, row1, sol) != 0)
                    	{
                            nz_sgv_r[nz_r][0] = row0;
                            nz_sgv_r[nz_r][1] = col1;
                            nz_sgv_r[nz_r][2] = row1;
                            nz_sgv_r[nz_r][3] = -1;
                            ++nz_r;
                        }
                    }
                }
            }
        }
    }

    //Step 2 : Find equations : - lambda (X - Sum Y) > 0  (lambda should be increased : maximization)
    //That is x_ik=0 and sum y_rsik=1
    int t1 = 0;
    int t2 = 0;
    int t3 = 0;
    int if1 = 0;
    int if2 = 0;
    int if3 = 0;
    for(int col1(0); col1 != nb_col; ++col1)
    {
        int row1(sol[col1]);
        if(row1 >= 0)
        {
            int nb_next_row(g.get_nb_next_row(col1,row1));
            int nb_next_col(g.get_nb_next_col(col1,row1));

            if(nb_next_col >= 1 && nb_next_row >= 1)
            {
                //Find activated edges which tail is not activated
                for(int ind_col2(0); ind_col2 != nb_next_col; ++ind_col2)
                {
                    int col2(g.get_next_col(col1,row1,ind_col2));
                    int arc_end(dp.get_arc_out(col1,row1,ind_col2));
                    int row2(sol[col2]);

                    if(arc_end != row2 && arc_end >= 0)
                    {
                        // col2.arc_end is not selected, x = 0,
                        // but edge col1,row1,col2,arc_end is selected, ie Y = 1
                        // Find equations containing this Y
                        int nb_pred_col = g.get_nb_pred_col(col2,arc_end);
                        int nb_pred_row = g.get_nb_pred_row(col2,arc_end);
                        int ind_col1 = g.get_ind_col_col(col2, col1);
                        int ind_row1 = g.get_ind_row_row(arc_end, row1);

                        if(ind_col1 == -1 || ind_row1 == -1)
                        {
                            cout << "fatal error in sub_gradient norm\n";
                            exit(0);
                        }

                        //lambda_col_node checker
                        if(ind_row1 > 0 && ind_row1 < nb_pred_row - 1)
                        {
                            //only 1 broken lambda_col_node constraint
                            nz_sgv_c[nz_c][0] = col1;
                            nz_sgv_c[nz_c][1] = col2;
                            nz_sgv_c[nz_c][2] = arc_end;
                            nz_sgv_c[nz_c][3] = 1;
                            ++nz_c;
                        }
                        else if(ind_row1 == 0)
                        {
                            //lambda_col_node from 0 to ind_col1
                            for(int ind_col0(0); ind_col0 != ind_col1 + 1; ++ind_col0)
                            {
                                int col0 = g.get_pred_col(col2,arc_end,ind_col0);
                                nz_sgv_c[nz_c][0] = col0;
                                nz_sgv_c[nz_c][1] = col2;
                                nz_sgv_c[nz_c][2] = arc_end;
                                nz_sgv_c[nz_c][3] = 1;
                                ++nz_c;
                            }
                        }
                        else
                        {
                            //lambda_col_node from nb_pred_col to ind_col1
                            for(int ind_col0(ind_col1); ind_col0 != nb_pred_col; ++ind_col0)
                            {
                                int col0 = g.get_pred_col(col2,arc_end,ind_col0);
                                nz_sgv_c[nz_c][0] = col0;
                                nz_sgv_c[nz_c][1] = col2;
                                nz_sgv_c[nz_c][2] = arc_end;
                                nz_sgv_c[nz_c][3] = 1;
                                ++nz_c;
                            }
                        }

                        //lambda_row_node checker
                        if(ind_col1 > 0 && ind_col1 < nb_pred_col - 1)
                        {
                        	++if1;
                            //only 1 broken lambda_row_node constraint
                            nz_sgv_r[nz_r][0] = row1;
                            nz_sgv_r[nz_r][1] = col2;
                            nz_sgv_r[nz_r][2] = arc_end;
                            nz_sgv_r[nz_r][3] = 1;
                            ++nz_r;
                            ++t1;

                        }
                        else if(ind_col1 == 0)
                        {
                        	++if2;
                            //lambda_row_node from 0 to ind_row1
                            for(int ind_row0(0); ind_row0 != ind_row1 + 1; ++ind_row0)
                            {
                                int row0 = g.get_pred_row(col2,arc_end,ind_row0);
                                nz_sgv_r[nz_r][0] = row0;
                                nz_sgv_r[nz_r][1] = col2;
                                nz_sgv_r[nz_r][2] = arc_end;
                                nz_sgv_r[nz_r][3] = 1;
                                ++nz_r;
                                ++t2;
                            }
                        }
                        else
                        {
                        	++if3;
                            //lambda_row_node from nb_pred_row to ind_row1
                            for(int ind_row0(ind_row1); ind_row0 != nb_pred_row; ++ind_row0)
                            {
                                int row0 = g.get_pred_row(col2,arc_end,ind_row0);
                                nz_sgv_r[nz_r][0] = row0;
                                nz_sgv_r[nz_r][1] = col2;
                                nz_sgv_r[nz_r][2] = arc_end;
                                nz_sgv_r[nz_r][3] = 1;
                                ++nz_r;
                                ++t3;
                            }
                        }
                    }
                }
            }
        }
    }

    //Step 3 : Find for activation constraints:
    //Equations : - lambda (1-x_ik+y_rsik-x_rs) > 0  (lambda should be increased : maximization)
    //This means: x_ik=1, y_rsik=0, x_rs=1
    for(int col0(0); col0 != nb_col; ++col0)
    {
        int row0(sol[col0]);
        if(row0 >= 0)
        {
            for(int col1(col0+1); col1 != nb_col; ++col1)
            {
            	int row1(sol[col1]);

            	int ind_col = g.get_ind_col_col(col0,col1);
            	if(ind_col<0)
            	{
            		continue;
            	}

            	if(row1 >= 0 && dp.get_arc_out(col0,row0,ind_col) != row1 && g.get_edge_coef(col0,row0,col1,row1) > 0.0) //>=
                {
                    nz_sgv_r_act[nz_r_act][0] = row0;
                    nz_sgv_r_act[nz_r_act][1] = col1;
                    nz_sgv_r_act[nz_r_act][2] = row1;
                    nz_sgv_r_act[nz_r_act][3] = 1;
                    ++nz_r_act;
                }
            }
        }
    }

    //Step 4: Find activation constraints:
    //Equations : - lambda (1-x_ik+y_rsik-x_rs) < 0  (lambda should be decreased : maximization)
    //This means: x_ik=0, y_rsik=0, x_rs=0 or x_ik=0, y_rsik=1, x_rs=1

    int ind_node1;
    int nbIter = g.get_nb_node();
    #pragma omp parallel for default(shared) schedule(dynamic) num_threads(NUM_THREADS)
    for(ind_node1=0; ind_node1 < nbIter; ++ind_node1)
    {
    	int col1;
    	int row1;
    	g.get_node_from_ind(ind_node1,col1,row1);

        int nb_pred_col(g.get_nb_pred_col(col1,row1));
        int nb_pred_row(g.get_nb_pred_row(col1,row1));

        //find which equations has a Sum Y-X = 0
        if(nb_pred_col >=1 && nb_pred_row >=1)
        {
            //row checker
            for(int ind_row(0); ind_row != nb_pred_row; ind_row++)
            {
                int row0(g.get_pred_row(col1,row1,ind_row));
                if(lb_mat.get_node_row_act(row0,col1,row1) > 0.0)
                {
                	if(is_brk_cst_row_node_act(row0, col1, row1, sol) != 0)
                	{
                		int nz_r_act_tmp;
						#pragma omp critical(incs)
                		{
                			nz_r_act_tmp = nz_r_act++;
                		}
			
                    	nz_sgv_r_act[nz_r_act_tmp][0] = row0;
                    	nz_sgv_r_act[nz_r_act_tmp][1] = col1;
                    	nz_sgv_r_act[nz_r_act_tmp][2] = row1;
                    	nz_sgv_r_act[nz_r_act_tmp][3] = -1;
                    	//++nz_r_act;
                    }
                }
            }
        }
    }

    lb_mat.set_nz_c(nz_c);
    lb_mat.set_nz_r(nz_r);
    lb_mat.set_nz_r_act(nz_r_act);

    
    return(nz_c + nz_r + nz_r_act);
}

/**************************************************
 * UPDATE_LAMBDA                                  *
 **************************************************/
void problem_apurva :: update_lambda(double step)
{
    int col0(0), col1(0), row0(0), row1(0), sens(0);

    int ** nz_sgv_c = lb_mat.get_nz_sgv_c();
    int ** nz_sgv_r = lb_mat.get_nz_sgv_r();
    int ** nz_sgv_r_act = lb_mat.get_nz_sgv_r_act();
    int nz_c = lb_mat.get_nz_c();
    int nz_r = lb_mat.get_nz_r();
    int nz_r_act = lb_mat.get_nz_r_act();

    
    //Updating lambda_node_col
    int i;
    #pragma omp parallel for schedule(dynamic), default(shared), num_threads(NUM_THREADS), private(col0, col1, row1, sens)
    for(i = 0; i< nz_c; ++i)
    {
        col0 = nz_sgv_c[i][0];
        col1 = nz_sgv_c[i][1];
        row1 = nz_sgv_c[i][2];
        sens = nz_sgv_c[i][3];

        if(sens == 1)
        {
            lb_mat.add_node_col(col0, col1, row1, step);
        }
        else
        {
            lb_mat.sub_node_col(col0, col1, row1, step);
        }
    }
    //Updating lambda_node_row
    #pragma omp parallel for schedule(dynamic), default(shared), num_threads(NUM_THREADS), private(row0, col1, row1, sens)
    for(i = 0; i< nz_r; ++i)
    {
        row0 = nz_sgv_r[i][0];
        col1 = nz_sgv_r[i][1];
        row1 = nz_sgv_r[i][2];
        sens = nz_sgv_r[i][3];

        if(sens == 1)
        {
            lb_mat.add_node_row(row0, col1, row1, step);
        }
        else
        {
            lb_mat.sub_node_row(row0, col1, row1, step);
        }
    }
    //Updating for activation constraints: lambda_node_row
    #pragma omp parallel for schedule(dynamic), default(shared), num_threads(NUM_THREADS), private(col0, col1, row1, sens)
    for(int i = 0; i< nz_r_act; ++i)
    {
        row0 = nz_sgv_r_act[i][0];
        col1 = nz_sgv_r_act[i][1];
        row1 = nz_sgv_r_act[i][2];
        sens = nz_sgv_r_act[i][3];

        if(sens == 1)
        {
            lb_mat.add_node_row_act(row0, col1, row1, step);
        }
        else
        {
            lb_mat.sub_node_row_act(row0, col1, row1, step);
        }
    }
}


/**************************************************
 * LR_SGD_SOLVE                                   *
 **************************************************/
void problem_apurva :: lr_sgd_solve(parameters & params)
{
    double step;

    double gamma(params.gamma);
    double theta(params.theta);
    double stepsize_min(params.stepsize_min);
    int nb_non_increasing_steps_max(params.nb_non_increasing_steps_max);
    int cnt_non_increas(0);
    int cnt_increas(0);
    int cnt_break(0);

    cout.precision(12); // utilité ??

    long tic_per_sec = sysconf(_SC_CLK_TCK);
    struct tms start, end;
    int tic1(times(&start));

    //Prohibit edge going out range
    lb_mat.prohibit(g, lo, up);

    int tic2 = times(&end);
    solve_time = ((double)tic2 - (double)tic1) / (double)tic_per_sec;

    while(status == 1)
    {
    	tic1 = times(&start);

        ++iter;
        dp.fill(g,lb_mat,lo,up);

        double current_lb(dp.solve_w_gapcosts(g,solution,lb_mat,lo,up));

        double current_ub(get_correct_value(solution));

        int sub_gr_norm(get_nz_sub_gr(solution));


		double ub_score = -2*current_lb/(self1+self2);
		double lb_score = 2*max(0.0,-current_ub)/(self1+self2);
    	

    	tic2 = times(&end);
    	solve_time += ((double)tic2 - (double)tic1) / (double)tic_per_sec;
    	tic1 = times(&start);


        if(lb != ub && -lb <= params.score_min)
        {
	    cout << -lb << " < " << params.score_min << endl;
            cout <<"Not similar: Stop now.\n";
	    status = NOT_SIMILAR;
        }
	else
	{
		if (current_lb > lb)    // if current_lb > lb then increase the number of improving iteration, and number of non-improving iteration = 0.
		{
		    lb = current_lb;
		    cnt_non_increas = 0;
		    cnt_break = 0;
		    cnt_increas++;
		}
		if (current_ub < ub)    // if current_ub < ub then increase the number of improving iteration, and number of non-improving iteration = 0.
		{
		    ub = current_ub;
		    cnt_non_increas = 0;
		    cnt_break = 0;
		    cnt_increas++;
		    for(int ii(0); ii < size ;ii++)
		    {
			 best_solution[ii] = solution[ii];
		    }
		}
		if (current_lb != lb && current_ub != ub)   // else, number of improving iteration = 0, and increase number of non-improving iteration
		{
		    cnt_non_increas++;
		    cnt_break++;
		    cnt_increas = 0;
		}
		if (cnt_non_increas > 5)   // after 20 non-improving iteration, reduce gamma
		{
		    gamma *= theta;
		    cnt_non_increas = 0;
		}
		if (cnt_increas > 5)       // after 20 improving iteration, increase gamma
		{
		    gamma /= theta;
		    cnt_increas = 0;
		}
		if (lb >= ub || ((int)lb >= ub && obj_is_int))
		{
			cout <<"Optimal. Stop now.\n";
		    status = OPTIMAL;
		}
		else if (ub-lb <= params.epsilon)
		{
			cout <<"Less than " << params.epsilon << " difference between upper and lower bound. Stop now.\n";
		    status = EPSILON;
		}
		else if(sub_gr_norm == 0)
		{
		    cout <<"Subgradient = 0. Stop now.\n";
		    status = APPROXIMATE;
		}
		else if (iter >= params.max_iteration)
		{
			cout <<"Maximum itrations reached. Stop now.\n";
		    status = MAX_ITER;
		}
		else if (lb >= params.limit_lb)
		{
			cout <<"Useless node. Stop now.\n";
		    status = APPROXIMATE;
		}
		else if(gamma <= stepsize_min)
		{
			cout <<"Stepsize less than " << stepsize_min << ". Stop now.\n";
		    status = APPROXIMATE;
		}
		else if(cnt_break > nb_non_increasing_steps_max)
		{
		    cout << "More than " << nb_non_increasing_steps_max <<" non increasing steps. Stop now.\n";
		    status = APPROXIMATE;
		}
		else if(solve_time >= params.my_time_limit)
		{
			cout <<"Time limit reached. Stop now.\n";
		    status = TIME_LIMIT;
		}
		else
		{
		    step = gamma * (current_ub - lb) / sub_gr_norm;
		    if(step <= 0.)
		    {
			cout << "Error in step value (negative step are not allowed)\n";
			status = APPROXIMATE;
		    }
		    update_lambda(step);
		}
	}
	cout << "UB-LB=" << ub - lb << endl;
    	tic2 = times(&end);
    	solve_time += ((double)tic2 - (double)tic1) / (double)tic_per_sec;

    }
    //Prohibited edges are reallowed
    lb_mat.allow(g, lo, up);

    if(obj_is_int)
    {
    	lb = (float)(int)lb;    //trick for apurva: rounding
    }

}


void problem_apurva :: recursive_brut_solve(int *brut_sol, int col, int cur_row, int nb_row, int &checked)
{
    /*
    *  End of recurence
    */
    if( cur_row >= nb_row || col >= size )
    {
        //evaluate solution :
        double value = get_correct_value(brut_sol);
        checked++;
        //
        if(value < ub)
        {
            cout << "brut_solve : new ub = " << value << ", sol checked = " << checked << "\n";
            //sauvegarder la solution et la borne

            ub = value;
        }
        if(checked % 100000 == 0)
        {
            cout << "brut_solve : " << checked << "solution checked\n";
        }
    }

    /*
    *  Else
    */
    else
    {
        int up_lim = up[col];
        int lo_lim = MAX(lo[col], cur_row);

        for(int row(lo_lim); row <= up_lim; ++row)
        {
            if(g.is_a_node(col, row))
            {
                brut_sol[col] = row;
                recursive_brut_solve(brut_sol, col+1, row+1, nb_row, checked);
            }
        }
        brut_sol[col] = -1;
        recursive_brut_solve(brut_sol, col+1, cur_row, nb_row, checked);
    }
}

/**************************************************
 * BRUTE_SOLVE                                    *
 **************************************************/
void problem_apurva :: brut_solve()
{
    int *brut_sol = new int[size];

    ub = 0.;
    lb = -INFINITY;
    int nb_row = g.get_nb_row();
    int checked = 0;

    double nb_pos(up[0] - lo[0] + 2);
    
    for(int i(1); i<size; i++)
    {
        nb_pos *= (up[i] - lo[i] + 2);
    }

    for(int i(0); i<size; i++)
    {
        cout << "lo = " << lo[i] << "up = " << up[i] << "\n";
    }
    
    cout << "nb combinaison a tester = " << nb_pos << "\n";
    
    for(int col(0); col < size; ++col)
    {
        brut_sol[col] = -1;
    }
    recursive_brut_solve(brut_sol, 0, -1, nb_row, checked);

    lb = ub;

    delete [] brut_sol;
}

problem * problem_apurva :: create_subproblem(int * lo1, int * up1)
{
    return( new problem_apurva((*this), lo1, up1));
}


void problem_apurva :: split(int * lo1, int * up1, int * lo2, int * up2)
{
    //int left;
    int best_col(-1), best_row(-1), min_area(0), max_area(0);
    int nb_col(g.get_nb_col());
    int nb_row(g.get_nb_row());

    int black[nb_col][nb_row];
    int white[nb_col][nb_row];
    
   
    /*****************************************
    * Step 1 - find the best splitting point *
    *****************************************/

    //Step 1.1 - Filling black values
    for(int row(nb_row -1); row >= 0; --row)
    {
        if(row > up[0])
        {
            black[0][row] = 0;
        }
        else
        {
            black[0][row] = up[0] - row;
        }
    }
    for(int col(1); col != nb_col; ++col)
    {
        for(int row(nb_row - 1); row >= 0; --row)
        {
            if(row > up[col])
            {
                black[col][row] = 0;
            }
            else
            {
                black[col][row] = (up[col] - row) + (black[col-1][row] + 1);
            }
        }
    }

    //Step 1.2 - filling white
    for(int row(0); row < nb_row; ++row)
    {
        if(row < lo[nb_col-1])
        {
            white[nb_col-1][row] = 0;
        }
        else
        {
            white[nb_col-1][row] = row - lo[nb_col-1] + 1;
        }
    }
    for(int col(nb_col-2); col >= 0; --col)
    {
        for(int row(0); row < nb_row; ++row)
        {
            if(row < lo[col])
            {
                white[col][row] = 0;
            }
            else
            {
                white[col][row] = (row - lo[col] + 1) + white[col+1][row];
            }
        }
    }

    //step 1.3 - finding best splitting point, maximising the minimum of black and white area
    for(int col(0); col != nb_col; ++col)
    {
        for(int row(lo[col]); row <= up[col]; ++row)
        {
            if(MIN(black[col][row], white[col][row]) > min_area)
            {
                best_col = col;
                best_row = row;
                min_area = MIN(black[col][row], white[col][row]);
                max_area = MAX(black[col][row], white[col][row]);
            }
            else if(MIN(black[col][row], white[col][row]) == min_area)
            {
                if(MAX(black[col][row], white[col][row]) > max_area)
                {
                    best_col = col;
                    best_row = row;
                    min_area = MIN(black[col][row], white[col][row]);
                    max_area = MAX(black[col][row], white[col][row]);
                }
            }
        }
    }

    if(best_col == -1 || best_row == -1 || min_area==0)
        //The problem cannot be splitted anymore
        return;

    /*****************************************************************
    * Step 2 - Computing the borders of white and black sub-problems *
    *****************************************************************/

    //Step 2.1 - Black sub-problem
    for(int col(0); col != nb_col; ++col)
    {
        //lower_limits don't change
        lo2[col] = lo[col];

        //but upper do
        if(col < best_col && up[col] > best_row -1)
            up2[col] = best_row -1;
        else if (col == best_col && up[col] > best_row)
            up2[col] = best_row;
        else
            up2[col] = up[col];
    }

    //Step 2.2 - White sub-problem
    for(int col(0); col != nb_col; ++col)
    {
        //upper_limits don't change
        up1[col] = up[col];

        //but lower do
        if(col >= best_col && lo[col] < best_row)
            lo1[col] = best_row;
        else
            lo1[col] = lo[col];
    }

    is_split = 1;
    return;
}
