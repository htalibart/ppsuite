#ifndef LAMBDA_MAT_APURVA_1
#define LAMBDA_MAT_APURVA_1
using namespace std;


#define MAX(x,y) ((x<y)?y:x)
#define MAX_PRED_COL 10
#define MAX_PRED_ROW 10
#define MAX_BRK_CST 4000

/**
* Lambda_Mat_Apurva Class.
* Contain the lambda coeficients of the lagrangian relaxation of the A_purva model.
*/

class lambda_mat_apurva
{
    protected:
        /**
        * The Contact Map Alignment graph
        */
        graph_apurva & g;

        //Lambda_multiplier for each king of relaxed constraints

        /**
        * Lambda coeficients for relaxed column/node constraints
        */
        double ** lambda_node_col;
        /**
        * Lambda coeficients for relaxed row/node constraints
        */
        double ** lambda_node_row;
        /**
        * Lambda coeficients for relaxed row/node activation constraints
        */
        double ** lambda_node_row_act;

        //Lambda modified coeficient for a node and for en edge

        /**
        * Lambda modified coeficients for nodes
        */
        double * coef_lambda_node;
        /**
        * Lambda modified coeficients for edges
        */
        double ** coef_lambda_edge;

        //Non-Zero components of the subgradient vector

        /**
        * Non-Zero components of the subgradient vector for lambda column/node
        */
        int **nz_sgv_c;
        /**
        * Number of non-zero components of the subgradient vector for lambda column/node
        */
        int nz_c;
        /**
        * Non-Zero components of the subgradient vector for lambda row/node
        */
        int **nz_sgv_r;
        /**
        * Number of non-zero components of the subgradient vector for lambda row/node
        */
        int nz_r;
        /**
        * Non-Zero components of the subgradient vector for lambda row/node for activation constraints
        */
        int **nz_sgv_r_act;
        /**
        * Number of non-zero components of the subgradient vector for lambda row/node for activation constraints
        */
        int nz_r_act;

        /**
        * Table recordig which nodes were modified during last sgd iteration,
        * in order not to recompute the weight of all nodes...
        */
        int ** lambda_modified_nodes;

        /**
         * This sum is a term in LR(lambda) and needs to be added to the upper bound
         */
        double sum_lambda_act;


    public:

        /**
        * Constructors
        */
        lambda_mat_apurva(graph_apurva & g1):
            g(g1)
            {
                int nb_col(g1.get_nb_col());
                int nb_row(g1.get_nb_row());
                int nb_node(g1.get_nb_node());

                int size;
                int col1, row1, col2, row2;
                int nb_next_col, nb_next_row;

                lambda_node_col = new double* [nb_col];
                lambda_node_row = new double* [nb_row];
                lambda_node_row_act = new double* [nb_row];
                for(int col(0); col != nb_col; ++col)
                {
                    lambda_node_col[col] = new double [nb_node];
                    fill_n( lambda_node_col[col], nb_node, static_cast<double>( 0. ) );
                }
                for(int row(0); row != nb_row; ++row)
                {
                    lambda_node_row[row] = new double [nb_node];
                    fill_n( lambda_node_row[row], nb_node, static_cast<double>( 0. ) );
                }
                for(int row(0); row != nb_row; ++row)
                {
                    lambda_node_row_act[row] = new double [nb_node];
                    fill_n( lambda_node_row_act[row], nb_node, static_cast<double>( 0. ) );
                }

                coef_lambda_node = new double [nb_node];
                for(int cur_node(0); cur_node < nb_node; ++cur_node)
                {
                	int ncol;
                	int nrow;
                	g1.get_node_from_ind(cur_node,ncol,nrow);
                	coef_lambda_node[cur_node] = g1.get_node_coef(ncol,nrow);
                }
                coef_lambda_edge = new double* [nb_node];
                for(int node(0);node != nb_node; ++node)
                {
                    g1.get_node_from_ind(node, col1, row1);
                    nb_next_col = g1.get_nb_next_col(col1, row1);
                    nb_next_row = g1.get_nb_next_row(col1, row1);
                    size = nb_next_col * nb_next_row;
                    if(size > 0)
                    {
                        coef_lambda_edge[node] = new double [size];
                        for(int ncol(0); ncol < nb_next_col; ncol++)
                        {
                            col2 = g1.get_next_col(col1, row1, ncol);
                            for(int nrow(0); nrow < nb_next_row; nrow++)
                            {
                                row2 = g1.get_next_row(col1, row1, nrow);
                                coef_lambda_edge[node][ncol*nb_next_row + nrow] = g1.get_edge_coef(col1,row1,col2,row2);
                            }
                        }
                    }
                    else
                        coef_lambda_edge[node] = 0;
                }

                int max_broken_cstr = 100000;//MAX(nb_col,nb_row) * 7 *7;

                nz_sgv_c = new int* [ max_broken_cstr];
                for(int bc(0); bc != max_broken_cstr; ++bc)
                {
                    nz_sgv_c[bc] = new int [4];
                }
                nz_c = 0;

                nz_sgv_r = new int* [ max_broken_cstr ];
                for(int bc(0); bc != max_broken_cstr; ++bc)
                {
                    nz_sgv_r[bc] = new int [4];
                }
                nz_r = 0;

                nz_sgv_r_act = new int* [ max_broken_cstr ];
                for(int bc(0); bc != max_broken_cstr; ++bc)
                {
                    nz_sgv_r_act[bc] = new int [4];
                }
                nz_r_act = 0;

                lambda_modified_nodes = new int* [nb_col];
                for(int col(0); col < nb_col; ++col)
                {
                    lambda_modified_nodes[col] = new int [nb_row];
                    fill_n( lambda_modified_nodes[col], nb_row, static_cast<int>( 1 ) );
                }

                sum_lambda_act = 0.0;


            };

        /**
        * Destructor
        */
        virtual ~lambda_mat_apurva()
        {
            int nb_col(g.get_nb_col());
            int nb_row(g.get_nb_row());
            int nb_node(g.get_nb_node());
            int max_broken_cstr = 100000;//MAX(nb_col,nb_row) * 7 *7;

            for(int col(0); col != nb_col; ++col)
            {
                delete [] lambda_node_col[col];
            }
            delete [] lambda_node_col;

            for(int row(0); row != nb_row; ++row)
            {
                delete [] lambda_node_row[row];
            }
            delete [] lambda_node_row;

            for(int row(0); row != nb_row; ++row)
            {
                delete [] lambda_node_row_act[row];
            }
            delete [] lambda_node_row_act;

            delete [] coef_lambda_node;

            for(int node(0);node != nb_node; ++node)
            {
                delete [] coef_lambda_edge[node];
            }
            delete [] coef_lambda_edge;

            for(int bc(0); bc != max_broken_cstr; ++bc)
            {
                delete [] nz_sgv_c[bc];
            }
            delete [] nz_sgv_c;

            for(int bc(0); bc != max_broken_cstr; ++bc)
            {
                delete [] nz_sgv_r[bc];
            }
            delete [] nz_sgv_r;

            for(int bc(0); bc != max_broken_cstr; ++bc)
            {
                delete [] nz_sgv_r_act[bc];
            }
            delete [] nz_sgv_r_act;

            for(int col(0); col < nb_col; ++col)
            {
                delete [] lambda_modified_nodes[col];
            }
            delete [] lambda_modified_nodes;

        };

        /**
        * Accessors
        */

        /**
        * Return the lambda coeficients for the relaxed constraint between column col1 and node N_{col2,row2}.
        */
        inline double get_node_col(int col1, int col2, int row2)
        {
            //cout << "looking for node " << col2 << "." << row2 <<" with col " << col1 << "\n";
            int ind_node2 = g.get_ind_from_node(col2,row2);
            return(lambda_node_col[col1][ind_node2]);
        }
        /**
        * Return the lambda coeficients for the relaxed constraint between row row1 and node N_{col2,row2}.
        */
        inline double get_node_row(int row1, int col2, int row2)
        {
            //cout << "looking for node " << col2 << "." << row2 <<" with row " << row1 << "\n";
            int ind_node2 = g.get_ind_from_node(col2,row2);
            return(lambda_node_row[row1][ind_node2]);
        }
        /**
        * Return the lambda coeficients for the relaxed activation constraint between row row1 and node N_{col2,row2}.
        */
        inline double get_node_row_act(int row1, int col2, int row2)
        {
            //cout << "looking for node " << col2 << "." << row2 <<" with row " << row1 << "\n";
            int ind_node2 = g.get_ind_from_node(col2,row2);
            return(lambda_node_row_act[row1][ind_node2]);
        }

        /**
        * Modificators
        */

        virtual inline void clear_lambda(graph_apurva &g1)
        {
            int nb_col(g.get_nb_col());
            int nb_row(g.get_nb_row());
            int nb_node(g1.get_nb_node());
            int size;
            int col1, row1;

            for(int col(0); col != nb_col; ++col)
            {
                fill_n( lambda_node_col[col], nb_node, static_cast<double>( 0. ) );
            }
            for(int row(0); row != nb_row; ++row)
            {
                fill_n( lambda_node_row[row], nb_node, static_cast<double>( 0. ) );
            }

            for(int cur_node(0); cur_node < nb_node; ++cur_node)
            {
            	int ncol;
            	int nrow;
            	g1.get_node_from_ind(cur_node,ncol,nrow);
            	coef_lambda_node[cur_node] = g1.get_node_coef(ncol,nrow);
            }

            for(int node(0);node != nb_node; ++node)
            {
                g1.get_node_from_ind(node, col1, row1);
                size = g1.get_nb_next_col(col1, row1) * g1.get_nb_next_row(col1, row1);
                if(size > 0)
                {
                    fill_n( coef_lambda_edge[node], size, static_cast<double>( -1. ) );
                }
            }
        }
        
        
        /**
        * Prohibit edges going out of lower and upper limits
        */
        virtual inline void prohibit(graph_apurva &g, int *lo, int *up)
        {
            int nb_col(g.get_nb_col());
            int nb_row(g.get_nb_row());

                //prohibit edge out of lower and upper limit
            for(int col1(0); col1 != nb_col; ++col1)
            {
                for(int row1(0); row1 != nb_row; ++row1)
                {
                    //all node need to be recomputed
                    lambda_modified_nodes[col1][row1] = 1;

                    int ind_node1(g.get_ind_from_node(col1,row1));
                    if(ind_node1>=0)
                    {
                        int nb_next_col(g.get_nb_next_col(col1,row1));
                        int nb_next_row(g.get_nb_next_row(col1,row1));
                        if(nb_next_col >= 1 && nb_next_row >= 1)
                        {
                            for(int ind_col2(0); ind_col2 != nb_next_col; ++ind_col2)
                            {
                                int col2(g.get_next_col(col1,row1,ind_col2));
                                for(int ind_row2(0); ind_row2 != nb_next_row; ++ind_row2)
                                {
                                    int row2(g.get_next_row(col1,row1,ind_row2));
                                    if(row2 < lo[col2] || row2 > up[col2] || g.get_ind_from_node(col2,row2)<0)
                                    {
                                        coef_lambda_edge[ind_node1][ind_col2 * nb_next_row + ind_row2] += 1000.;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        /**
        * Allow edges forbidden by function prohibit
        */
        virtual inline void allow(graph_apurva &g, int *lo, int *up)
        {
            int nb_col(g.get_nb_col());
            int nb_row(g.get_nb_row());

                //prohibit edge out of lower and upper limit
            for(int col1(0); col1 != nb_col; ++col1)
            {
                for(int row1(0); row1 != nb_row; ++row1)
                {
                    int ind_node1(g.get_ind_from_node(col1,row1));
                    if(ind_node1>=0)
                    {
                        int nb_next_col(g.get_nb_next_col(col1,row1));
                        int nb_next_row(g.get_nb_next_row(col1,row1));
                        if(nb_next_col >= 1 && nb_next_row >= 1)
                        {
                            for(int ind_col2(0); ind_col2 != nb_next_col; ++ind_col2)
                            {
                                int col2(g.get_next_col(col1,row1,ind_col2));
                                for(int ind_row2(0); ind_row2 != nb_next_row; ++ind_row2)
                                {
                                    int row2(g.get_next_row(col1,row1,ind_row2));
                                    if(row2 < lo[col2] || row2 > up[col2] || g.get_ind_from_node(col2,row2)<0)
                                    {
                                        coef_lambda_edge[ind_node1][ind_col2 * nb_next_row + ind_row2] -= 1000.;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        /**
        * update lambda coeficients for constraint col1 / node N_{col2,row2}.
        * add step value to lambda coeficient.
        * This function also update the lambda modified coeficients for node and edges.
        */
        inline void add_node_col(int col1, int col2, int row2, double step)
        {

            int ind_node2 = g.get_ind_from_node(col2,row2);
            //update lambda value

       	    lambda_node_col[col1][ind_node2] += step;
	    
            //update node modifier
	    	#pragma omp critical(coef)
       	    {
            	coef_lambda_node[ind_node2] += step;
       	    }

            //update edge modifiers
            int nb_pred_row = g.get_nb_pred_row(col2,row2);
            int nb_pred_col = g.get_nb_pred_col(col2,row2);
            int ind_col1 = g.get_ind_col_col(col2,col1);

            //->up apurva part
            int row0 = g.get_pred_row(col2,row2,0);
            for(int ind_col0(nb_pred_col - 1); ind_col0 != ind_col1; --ind_col0)
            {
                int col0 = g.get_pred_col(col2, row2, ind_col0);
                int ind_node0 = g.get_ind_from_node(col0,row0);

                if(ind_node0 >= 0)
                {
                    int nb_next_row = g.get_nb_next_row(col0,row0);
                    int ind_col2 = g.get_ind_col_col(col0, col2);
                    int ind_row2 = g.get_ind_row_row(row0, row2);

                    if(ind_row2 == -1 || ind_col2 == -1)
                    {
                        cout << "Fatal error in lambda add_node_col\n";
                        exit(0);
                    }

                    coef_lambda_edge[ind_node0][ind_col2 * nb_next_row + ind_row2] += step;
                    
                    lambda_modified_nodes[col0][row0] = 1;
                    lambda_modified_nodes[col2][row2] = 1;

                }
            }

            //->col part
            for(int ind_row1(0); ind_row1 != nb_pred_row; ++ind_row1)
            {
                int row1 = g.get_pred_row(col2,row2,ind_row1);
                int ind_node1 = g.get_ind_from_node(col1,row1);

                if(ind_node1 >= 0)
                {
                    int nb_next_row = g.get_nb_next_row(col1,row1);
                    int ind_col2 = g.get_ind_col_col(col1, col2);
                    int ind_row2 = g.get_ind_row_row(row1, row2);

                    if(ind_row2 == -1 || ind_col2 == -1)
                    {
                        cout << "Fatal error in lambda add_node_col\n";
                        exit(0);
                    }

                    coef_lambda_edge[ind_node1][ind_col2 * nb_next_row + ind_row2] += step;
		    
                    lambda_modified_nodes[col1][row1] = 1;
                    lambda_modified_nodes[col2][row2] = 1;

                }
            }
            //->lo apurva part
            row0 = g.get_pred_row(col2,row2,nb_pred_row - 1);
            for(int ind_col0(0); ind_col0 != ind_col1; ++ind_col0)
            {
                int col0 = g.get_pred_col(col2, row2, ind_col0);
                int ind_node0 = g.get_ind_from_node(col0,row0);

                if(ind_node0 >= 0)
                {
                    int nb_next_row = g.get_nb_next_row(col0,row0);
                    int ind_col2 = g.get_ind_col_col(col0, col2);
                    int ind_row2 = g.get_ind_row_row(row0, row2);

                    if(ind_row2 == -1 || ind_col2 == -1)
                    {
                        cout << "Fatal error in lambda add_node_col\n";
                        exit(0);
                    }
                    	
                    coef_lambda_edge[ind_node0][ind_col2 * nb_next_row + ind_row2] += step;
		    
                    lambda_modified_nodes[col0][row0] = 1;
                    lambda_modified_nodes[col2][row2] = 1;
                }
            }

        }

        /**
        * update lambda coeficients for constraint row1 / node N_{col2,row2}.
        * add step value to lambda coeficient.
        * This function also update the lambda modified coeficients for node and edges.
        */
        inline void add_node_row(int row1, int col2, int row2, double step)
        {
            int ind_node2 = g.get_ind_from_node(col2,row2);
            //update lambda value
            lambda_node_row[row1][ind_node2] += step;
            //update node modifier
            #pragma omp critical(coef)
	    {
	    	coef_lambda_node[ind_node2] += step;
	    }
            //update edge modifiers
            int nb_pred_row = g.get_nb_pred_row(col2,row2);
            int nb_pred_col = g.get_nb_pred_col(col2,row2);
            int ind_row1 = g.get_ind_row_row(row2,row1);

            //->up apurva part
            int col0 = g.get_pred_col(col2,row2,nb_pred_col - 1);
            for(int ind_row0(0); ind_row0 != ind_row1; ++ind_row0)
            {
                int row0 = g.get_pred_row(col2, row2, ind_row0);
                int ind_node0 = g.get_ind_from_node(col0,row0);

                if(ind_node0 >= 0)
                {
                    int nb_next_row = g.get_nb_next_row(col0,row0);
                    int ind_col2 = g.get_ind_col_col(col0, col2);
                    int ind_row2 = g.get_ind_row_row(row0, row2);

                    if(ind_row2 == -1 || ind_col2 == -1)
                    {
                        cout << "Fatal error in lambda add_node_row\n";
                        exit(0);
                    }

                    coef_lambda_edge[ind_node0][ind_col2 * nb_next_row + ind_row2] += step;

                    lambda_modified_nodes[col0][row0] = 1;
                    lambda_modified_nodes[col2][row2] = 1;

                }
            }

            // -> row part
            for(int ind_col1(0); ind_col1 != nb_pred_col; ++ind_col1)
            {
                int col1 = g.get_pred_col(col2,row2,ind_col1);
                int ind_node1 = g.get_ind_from_node(col1,row1);
                if(ind_node1 >= 0)
                {
                    int nb_next_row = g.get_nb_next_row(col1,row1);
                    int ind_col2 = g.get_ind_col_col(col1, col2);
                    int ind_row2 = g.get_ind_row_row(row1, row2);

                    if(ind_row2 == -1 || ind_col2 == -1)
                    {
                        cout << "Fatal error in lambda add_node_row\n";
                        exit(0);
                    }

                    coef_lambda_edge[ind_node1][ind_col2 * nb_next_row + ind_row2] += step;

                    lambda_modified_nodes[col1][row1] = 1;
                    lambda_modified_nodes[col2][row2] = 1;

                }
            }

            //->lo apurva part
            col0 = g.get_pred_col(col2,row2,0);
            for(int ind_row0(nb_pred_row - 1); ind_row0 != ind_row1; --ind_row0)
            {
                int row0 = g.get_pred_row(col2, row2, ind_row0);
                int ind_node0 = g.get_ind_from_node(col0,row0);

                if(ind_node0 >= 0)
                {
                    int nb_next_row = g.get_nb_next_row(col0,row0);
                    int ind_col2 = g.get_ind_col_col(col0, col2);
                    int ind_row2 = g.get_ind_row_row(row0, row2);

                    if(ind_row2 == -1 || ind_col2 == -1)
                    {
                        cout << "Fatal error in lambda add_node_row\n";
                        exit(0);
                    }

                    coef_lambda_edge[ind_node0][ind_col2 * nb_next_row + ind_row2] += step;

                    lambda_modified_nodes[col0][row0] = 1;
                    lambda_modified_nodes[col2][row2] = 1;

                }
            }
        }

        /**
        * update lambda coeficients for activation constraint row1 / node N_{col2,row2}.
        * add step value to lambda coeficient.
        * This function also update the lambda modified coeficients for node and edges.
        */
        inline void add_node_row_act(int row1, int col2, int row2, double step)
        {

            int ind_node2 = g.get_ind_from_node(col2,row2);
            //update lambda value
            lambda_node_row_act[row1][ind_node2] += step;
            //update node modifier
            #pragma omp critical(coef)
	    {
	    	coef_lambda_node[ind_node2] -= step;
	    }
	    #pragma omp critical(sum)
	    {
            	sum_lambda_act += step;
	    }

            //update edge modifiers
            int nb_pred_row = g.get_nb_pred_row(col2,row2);
            int nb_pred_col = g.get_nb_pred_col(col2,row2);
            int ind_row1 = g.get_ind_row_row(row2,row1);

            //->up apurva part
            int col0 = g.get_pred_col(col2,row2,nb_pred_col - 1);
            for(int ind_row0(0); ind_row0 != ind_row1; ++ind_row0)
            {
                int row0 = g.get_pred_row(col2, row2, ind_row0);
                int ind_node0 = g.get_ind_from_node(col0,row0);

                if(ind_node0 >= 0)
                {
                    int nb_next_row = g.get_nb_next_row(col0,row0);
                    int ind_col2 = g.get_ind_col_col(col0, col2);
                    int ind_row2 = g.get_ind_row_row(row0, row2);

                    if(ind_row2 == -1 || ind_col2 == -1)
                    {
                        cout << "Fatal error in lambda add_node_row_act\n";
                        exit(0);
                    }

                    if(g.get_edge_coef(col0,row0,col2,row2) > 0.0)
                    {
                    	coef_lambda_edge[ind_node0][ind_col2 * nb_next_row + ind_row2] -= step;
                    	#pragma omp critical(coef)
                    	{
                    		coef_lambda_node[ind_node0] -= step;
                    	}
                    	lambda_modified_nodes[col0][row0] = 1;
                    	lambda_modified_nodes[col2][row2] = 1;
                    }
                }
            }

            // -> row part
            for(int ind_col1(0); ind_col1 != nb_pred_col; ++ind_col1)
            {
                int col1 = g.get_pred_col(col2,row2,ind_col1);
                int ind_node1 = g.get_ind_from_node(col1,row1);

                if(ind_node1 >= 0)
                {
                    int nb_next_row = g.get_nb_next_row(col1,row1);
                    int ind_col2 = g.get_ind_col_col(col1, col2);
                    int ind_row2 = g.get_ind_row_row(row1, row2);

                    if(ind_row2 == -1 || ind_col2 == -1)
                    {
                        cout << "Fatal error in lambda add_node_row_act\n";
                        exit(0);
                    }

                    if(g.get_edge_coef(col1,row1,col2,row2) >= 0.0)
                    {
                    	coef_lambda_edge[ind_node1][ind_col2 * nb_next_row + ind_row2] -= step;
						#pragma omp critical(coef)
                    	{
                    		coef_lambda_node[ind_node1] -= step;
                    	}
                    	lambda_modified_nodes[col1][row1] = 1;
                    	lambda_modified_nodes[col2][row2] = 1;
                    }
                }
            }


            //->lo apurva part
            col0 = g.get_pred_col(col2,row2,0);
            for(int ind_row0(nb_pred_row - 1); ind_row0 != ind_row1; --ind_row0)
            {
                int row0 = g.get_pred_row(col2, row2, ind_row0);
                int ind_node0 = g.get_ind_from_node(col0,row0);

                if(ind_node0 >= 0)
                {
                    int nb_next_row = g.get_nb_next_row(col0,row0);
                    int ind_col2 = g.get_ind_col_col(col0, col2);
                    int ind_row2 = g.get_ind_row_row(row0, row2);

                    if(ind_row2 == -1 || ind_col2 == -1)
                    {
                        cout << "Fatal error in lambda add_node_row_act\n";
                        exit(0);
                    }

                    if(g.get_edge_coef(col0,row0,col2,row2) > 0.0)
                    {
                    	coef_lambda_edge[ind_node0][ind_col2 * nb_next_row + ind_row2] -= step;
                    	#pragma omp critical(coef)
                    	{
                    		coef_lambda_node[ind_node0] -= step;
                    	}
                    	lambda_modified_nodes[col0][row0] = 1;
                    	lambda_modified_nodes[col2][row2] = 1;
                    }
                }
            }
        }


        /**
        * Same as add_node_col, but for substracting the step value.
        * It modified the step value to avoid negative lambda (lambda are in R+)
        */
        inline void sub_node_col(int col1, int col2, int row2, double step)
        {
            int ind_node2 = g.get_ind_from_node(col2,row2);
            double my_step(step);

            //Fix step value to avoid negative lagrangian multiplicator
            if(lambda_node_col[col1][ind_node2] - step <= 0.)
                my_step = lambda_node_col[col1][ind_node2];

            add_node_col(col1,col2,row2, -my_step);
        }

        /**
        * Same as add_node_row, but for substracting the step value.
        * It modified the step value to avoid negative lambda (lambda are in R+)
        */
        inline void sub_node_row(int row1, int col2, int row2, double step)
        {
            int ind_node2 = g.get_ind_from_node(col2,row2);
            double my_step(step);

            //Fix step value to avoid negative lagrangian multiplicator
            if(lambda_node_row[row1][ind_node2] - step <= 0.)
                my_step = lambda_node_row[row1][ind_node2];

            add_node_row(row1,col2,row2, -my_step);

        }

        /**
        * Same as add_node_row_act, but for substracting the step value.
        * It modified the step value to avoid negative lambda (lambda are in R+)
        */
        inline void sub_node_row_act(int row1, int col2, int row2, double step)
        {
            int ind_node2 = g.get_ind_from_node(col2,row2);
            double my_step(step);

            //Fix step value to avoid negative lagrangian multiplicator
            if(lambda_node_row_act[row1][ind_node2] - step <= 0.)
                my_step = lambda_node_row_act[row1][ind_node2];

            add_node_row_act(row1,col2,row2, -my_step);

        }

        /**
        * This should be forbidden !!!!
        */
        inline double** get_lnc()
        {
            return(lambda_node_col);
        };
        /**
        * This should be forbidden !!!!
        */
        inline double** get_lnr()
        {
            return(lambda_node_row);
        };
        /**
        * This should be forbidden !!!!
        */
        inline double** get_lnr_act()
        {
            return(lambda_node_row_act);
        };
        /**
        * This should be forbidden !!!!
        */
        inline double * get_coef_lambda_node()
        {
            return(coef_lambda_node);
        }
        /**
        * This should be forbidden !!!!
        */
        inline double ** get_coef_lambda_edge()
        {
            return(coef_lambda_edge);
        }
      
        int ** get_nz_sgv_c()
        {
            return(nz_sgv_c);
        }
        int ** get_nz_sgv_r()
        {
            return(nz_sgv_r);
        }
        int ** get_nz_sgv_r_act()
        {
            return(nz_sgv_r_act);
        }
        int get_nz_c()
        {
            return(nz_c);
        }
        int get_nz_r()
        {
            return(nz_r);
        }
        int get_nz_r_act()
        {
            return(nz_r_act);
        }
        void set_nz_c(int value)
        {
            nz_c = value;
        }
        void set_nz_r(int value)
        {
            nz_r = value;
        }
        void set_nz_r_act(int value)
        {
            nz_r_act = value;
        }
        int ** get_lambda_modified_nodes()
        {
            return(lambda_modified_nodes);
        }
        double get_sum_lambda_act()
        {
        	return sum_lambda_act;
        }

};
#endif
