#ifndef GRAPH_1
#define GRAPH_1
using namespace std;

/**
*  VIRTUAL Class Graph
*
*/

class graph
{
    protected :

        /**
        * Number of columns.
        */
        int nb_col;

        /**
         * Number of rows.
        */
        int nb_row;

    public :

        /**
        * Constructors
        */

        /**
        * Build an empty graph (zero column and zero row).
        */
        graph(): nb_col(0), nb_row(0){};
        /**
        * Build a graph with nb_col columns and nb_row rows.
        */
        graph(int nb_col1, int nb_row1) : nb_col(nb_col1), nb_row(nb_row1) {};

        /**
        * Destructor
        */
        virtual ~graph(){};

        /**
        * Accessors
        */

        /**
        * Return column number
        */
        inline int get_nb_col(){return nb_col;}
        /**
        * Return row number
        */
        inline int get_nb_row(){return nb_row;}
        /**
        * Return coeficient S_{col,row} of node N_{col,row}
        */
        virtual inline float get_node_coef(int col, int row) = 0;
        /**
        * Return coeficient C_{col1,row1,col2,row2} of edge E_{col1,row1,col2,row2}
        */
        virtual inline float get_edge_coef(int col1, int row1, int col2, int row2) = 0;

};

#endif
