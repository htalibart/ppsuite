#ifndef GRAPH_APURVA_1
#define GRAPH_APURVA_1

using namespace std;

#define MIN(x,y) ((x<y)?x:y)

class graph_apurva : public graph
{
    protected :

        int nb_node;
        int ** node_to_ind;
        int ** ind_to_node;

        int ** edge_col;
        int ** edge_row;

        // neighbourhood
        int * nb_next_col;
        int * nb_next_row;
        int * nb_pred_col;
        int * nb_pred_row;
        int ** next_col;
        int ** next_row;
        int ** pred_col;
        int ** pred_row;

        int ** ind_col_col;
        int ** ind_row_row;

        /*
         * The function pointer for the function assigning node weights
         * @author IW
         */
        double (*score_vertex)(int,int);

        /*
         * The function pointer for the function assigning edge weights
         * @author IW
         */
        double (*score_edge)(int,int,int,int);


    public :

        /* Constructors & Destructor */
        graph_apurva(int nb_node1, int ** contacts1, int nb_node2, int ** contacts2, double (*s_vertex)(int,int), double (*s_edge)(int,int,int,int), int ** forbidden);

        ~graph_apurva()
        {
            for(int col(0); col != nb_col; ++col)
            {
                delete [] node_to_ind[col];
            }
            delete [] node_to_ind;

            for(int node1(0); node1 != nb_node; ++node1)
            {
                delete [] ind_to_node[node1];
            }
            delete [] ind_to_node;

            for(int col1(0); col1 != nb_col ; ++ col1)
            {
                delete [] ind_col_col[col1];
            }
            delete [] ind_col_col;

            for(int row1(0); row1 != nb_row ; ++ row1)
            {
                delete [] ind_row_row[row1];
            }
            delete [] ind_row_row;

            delete [] nb_next_col;
            delete [] nb_next_row;
            delete [] nb_pred_col;
            delete [] nb_pred_row;

            for(int col1(0); col1 != nb_col; ++col1)
            {
                delete [] next_col[col1];
            }
            delete [] next_col;

            for(int row1(0); row1 != nb_row; ++row1)
            {
                delete [] next_row[row1];
            }
            delete [] next_row;

            for(int col2(0); col2 != nb_col; ++col2)
            {
                delete [] pred_col[col2];
            }
            delete [] pred_col;
            
            for(int row2(0); row2 != nb_row; ++row2)
            {
                delete [] pred_row[row2];
            }
            delete [] pred_row;
            
        };

        /* Accessors */
        inline int get_nb_node(){ return(nb_node); };
        inline int get_ind_from_node(int col, int row){ return(node_to_ind[col][row]); };
        inline void get_node_from_ind(int ind, int & col, int & row)
        {
            col = ind_to_node[ind][0];
            row = ind_to_node[ind][1];
        };

        /* Number of columns / rows that are successor / predecessor of a node */
        inline int get_nb_next_col(int col, int row){ return(nb_next_col[col]); };
        inline int get_nb_next_row(int col, int row){ return(nb_next_row[row]); };
        inline int get_nb_pred_col(int col, int row){ return(nb_pred_col[col]); };
        inline int get_nb_pred_row(int col, int row){ return(nb_pred_row[row]); };

        /* retrieve the column / row that is the nb^th successor / predecessor of a node */
        inline int get_next_col(int col, int row, int nb){ return(next_col[col][nb]); };
        inline int get_next_row(int col, int row, int nb){ return(next_row[row][nb]); };
        inline int get_pred_col(int col, int row, int nb){ return(pred_col[col][nb]); };
        inline int get_pred_row(int col, int row, int nb){ return(pred_row[row][nb]); };

        // return the coeficient of a node (0. for all node in the case of apurva)
        // this function doesn't check if col.row is valid node
        inline double get_node_coef(int col, int row)
        {
        	return score_vertex(col, row);
        };

        // return the coeficient of an edge (1. for all edges in the case of apurva)
        // this function doesn't check if <col1.row1, col2.row2> is valid edge
        inline double get_edge_coef(int col1, int row1, int col2, int row2)
        {

            int edge1 = edge_col[col1][col2];
            int edge2 = edge_row[row1][row2];
            double value = MIN(edge1, edge2);

            if(value != 0)
            {
            	return -1*score_edge(col1,row1,col2,row2);
            }
            else
            {
            	return 0;
            }
        };

        /* boolean function, true if col.row is a node, else false */
        inline int is_a_node(int col, int row){ return(node_to_ind[col][row] >= 0); };

        /* boolean function, true if <col1.row1, col2.row2> is an edge, else false */
        inline int is_an_edge(int col1, int row1, int col2, int row2)
        {
            return(edge_col[col1][col2] && edge_row[row1][row2] && node_to_ind[col1][row1] >= 0 && node_to_ind[col2][row2] >= 0);
        };

        /* for each couple of columns / rows, give the next/pred indices */
        inline int get_ind_col_col(int col1, int col2){ return(ind_col_col[col1][col2]); };
        inline int get_ind_row_row(int row1, int row2){ return(ind_row_row[row1][row2]); };

        /* Direct access to data... */
        /* Should be forbidden, for security, but needed for speeding up */
        inline int ** get_node_to_ind(){ return(node_to_ind); };
        inline int * get_nb_next_col(){ return(nb_next_col); };
        inline int * get_nb_pred_col(){ return(nb_pred_col); };
        inline int * get_nb_next_row(){ return(nb_next_row); };
        inline int * get_nb_pred_row(){ return(nb_pred_row); };
        inline int ** get_next_col(){ return(next_col); };
        inline int ** get_pred_col(){ return(pred_col); };
        inline int ** get_next_row(){ return(next_row); };
        inline int ** get_pred_row(){ return(pred_row); };
        inline int get_edge_col(int col1, int col2){ return(edge_col[col1][col2]);};
        inline int get_edge_row(int row1, int row2){ return(edge_row[row1][row2]);};

        /// return the function pointer to the scoring function for the vertices (IW)
        inline double (*get_vertex_scoringfunction())(int,int)
		{
        	return score_vertex;
		}

        /// return the function pointer to the scoring function for the edges (IW)
        inline double (*get_edge_scoringfunction())(int,int,int,int)
		{
        	return score_edge;
		}

};

#endif

