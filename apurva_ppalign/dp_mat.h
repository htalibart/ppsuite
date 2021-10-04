#ifndef DP_MAT1
#define DP_MAT1

using namespace std;

/**
* DP_MAT Class.
* Used for dynamic-programming.
*/
class dp_mat
{
    protected :

        int nb_col;
        int nb_row;

        /**
        * For each node N_{col_row} and each next column, indicates chosen next row
        * ie, gives all activated outgoing edge.
        */
        int *** dp_arc_out_col;
        /**
        * 2D table used to compute by dynamic programming the best set of outgoing edges.
        * Contains values of outgoing edges sets
        */
        double ** dp_arc;
        /**
        * 2D table used to compute by dynamic programming the best set of outgoing edges.
        * Contains moves made to obtain values in dp_arc.
        */
        int ** arc_from;


        /**
        * 2D table used to compute by dynamic programming the best set of nodes.
        * Contains values of path made in DP_score.
        */
        double ** dp_M;
        double ** dp_GA;
        double ** dp_GB;
        /**
        * 2D table used to compute by dynamic programming the best set of nodes.
        * Contains score of each node, based on there coeficient and best outgoing edge value.
        */
        double ** dp_score;
        /**
        * 2D table used to compute by dynamic programming the best set of nodes.
        * Contain moves made to obtain values in dp_score.
        */
        int ** dp_M_from;
        int ** dp_GA_from;
        int ** dp_GB_from;

        /**
         * The gap costs (>=0) (gap open and gap extension costs)
         */
        //double gap_open;
        //double gap_extend;
	

	float * insert_open_row;
	float * insert_open_col;
	float * insert_extend_row;
	float * insert_extend_col;

    public :

        /**
        * Constructors
        */
        dp_mat()
        {
            nb_col = 0;
            nb_row = 0;
            dp_arc_out_col = 0;
            dp_arc = 0;
            arc_from = 0;
            dp_M = 0;
            dp_score = 0;
            dp_M_from = 0;
        };
        /**
        * Destructor
        */
        virtual ~dp_mat()
        {
        };

        /**
        * For a given node N_{col1,row1} and a given next column (define by a next_col indice ind), return
        * the chosen row2, corresponding to the activated edge E_{col1,row1,col2,row2}.
        */
        virtual inline int get_arc_out(int col1, int row1, int ind)=0;
};

#endif

