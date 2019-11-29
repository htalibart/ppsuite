#include "main.h"
using namespace std;




graph_apurva :: graph_apurva(int nb_node1, int ** contacts1, int nb_node2, int ** contacts2, double (*s_vertex)(int,int), double (*s_edge)(int,int,int,int), int ** forbidden)
{
    edge_col = contacts2;
    edge_row = contacts1;

    nb_col=nb_node2;
    nb_row=nb_node1;

    score_vertex = s_vertex;
    score_edge = s_edge;

    nb_node = 0;
    node_to_ind = new int * [nb_node2];

    for(int col(0); col != nb_node2; ++col)
    {
        node_to_ind[col] = new int [nb_node1];
        for(int row(0); row != nb_node1; ++row)
        {
            if(!forbidden[col][row])
            {
                node_to_ind[col][row] = nb_node;
                ++nb_node;
            }
            else
                node_to_ind[col][row] = -1;
        }
    }

    ind_to_node = new int * [nb_node];
    for(int col(0); col != nb_node2; ++col)
    {
        for(int row(0); row != nb_node1; ++row)
        {
            int node_ind =  node_to_ind[col][row];
            if(node_ind >= 0)
            {
                ind_to_node[node_ind] = new int [2];
                ind_to_node[node_ind][0]=col;
                ind_to_node[node_ind][1]=row;
            }
        }
    }

    ind_col_col = new int * [nb_col];
    for(int col1(0); col1 != nb_col ; ++ col1)
    {
        ind_col_col[col1] = new int [nb_col];
        fill_n(ind_col_col[col1], nb_col, static_cast<int>( -1 ));
    }
    ind_row_row = new int * [nb_row];
    for(int row1(0); row1 != nb_row ; ++ row1)
    {
        ind_row_row[row1] = new int [nb_row];
        fill_n(ind_row_row[row1], nb_row, static_cast<int>( -1 ));
    }

    //Step 2 Next_col
    nb_next_col = new int [nb_col];
    next_col = new int * [nb_col];

    int temp_next_col[nb_col];

    for(int col1(0); col1 != nb_col; ++col1)
    {
        int nb_next(0);
        for(int col2(col1 + 1); col2 != nb_col; ++col2)
        {
            if(edge_col[col1][col2])
            {
                temp_next_col[nb_next] = col2;
                nb_next++;
            }
        }
        nb_next_col[col1] = nb_next;
        if(nb_next)
        {
            next_col[col1] = new int [nb_next];
            for(int n(0); n != nb_next; ++n)
            {
                next_col[col1][n] = temp_next_col[n];
                ind_col_col[col1][ next_col[col1][n] ] = n;
            }
        }
        else
            next_col[col1] = NULL;
    }

    //Step 3 Next_row
    nb_next_row = new int [nb_row];
    next_row = new int * [nb_row];

   int temp_next_row[nb_row];

    for(int row1(0); row1 != nb_row; ++row1)
    {
        int nb_next(0);
        for(int row2(row1 + 1); row2 != nb_row; ++row2)
        {
            if(edge_row[row1][row2])
            {
                temp_next_row[nb_next] = row2;
                nb_next++;
            }
        }
        nb_next_row[row1] = nb_next;
        if(nb_next)
        {
            next_row[row1] = new int [nb_next];
            for(int n(0); n != nb_next; ++n)
            {
                next_row[row1][n] = temp_next_row[n];
                ind_row_row[row1][ next_row[row1][n] ] = n;
            }
        }
        else
            next_row[row1] = NULL;
    }

    //Step 4 pred_col
    nb_pred_col = new int [nb_col];
    pred_col = new int * [nb_col];

    int temp_pred_col[nb_col];

    for(int col2(0); col2 != nb_col; ++col2)
    {
        int nb_pred(0);
        for(int col1(col2 - 1); col1 >= 0; --col1)
        {
            if(edge_col[col1][col2])
            {
                temp_pred_col[nb_pred] = col1;
                nb_pred++;
            }
        }
        nb_pred_col[col2] = nb_pred;
        if(nb_pred)
        {
            pred_col[col2] = new int [nb_pred];
            for(int n(0); n != nb_pred; ++n)
            {
                pred_col[col2][n] = temp_pred_col[n];
                ind_col_col[col2][ pred_col[col2][n] ] = n;
            }
        }
        else
            pred_col[col2] = NULL;
    }

    //Step 5 pred_row
    nb_pred_row = new int [nb_row];
    pred_row = new int * [nb_row];

    int temp_pred_row[nb_row];

    for(int row2(0); row2 != nb_row; ++row2)
    {
        int nb_pred(0);
        for(int row1(row2 - 1); row1 >= 0; --row1)
        {
            if(edge_row[row1][row2])
            {
                temp_pred_row[nb_pred] = row1;
                nb_pred++;
            }
        }
        nb_pred_row[row2] = nb_pred;
        if(nb_pred)
        {
            pred_row[row2] = new int [nb_pred];
            for(int n(0); n != nb_pred; ++n)
            {
                pred_row[row2][n] = temp_pred_row[n];
                ind_row_row[row2][ pred_row[row2][n] ] = n;
            }
        }
        else
            pred_row[row2] = NULL;
    }

}


