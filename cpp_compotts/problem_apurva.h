#ifndef PROBLEM_APURVA_1
#define PROBLEM_APURVA_1
using namespace std;

/**
* Problem_apurva Class
* This is used to solve a Contact Map Overlap Maximisation problem,
* using the A_purva integer programming model and lagrangian relaxation approach.
*/
class problem_apurva : public problem
{
    private :

        /**
        * A Contact Map Alignment graph.
        */
        graph_apurva & g;
        /**
        * The Dynamic Programming object used to solve relaxed problem.
        */
        dp_mat_apurva & dp;
        /**
        * The Lambda coefficient associated to relaxed constraints.
        */
        lambda_mat_apurva & lb_mat;
        /**
         * The scores for comparing the proteins with themselves
         */
        double self1;
        double self2;
        double bound_ub_score;
        /**
         * A flag indicating whether the objective function has integer value,
         * in which case rounding can be used.
         */
        bool obj_is_int;

        /**
        * Return the value of the solution for the Primal problem.
        */
        double get_correct_value(int * sol);

        /**
        * Return the nb of non-zero component of the subgradient vector,
        * and the list of non-zero components is stored in nz_sgv
        */
        int get_nz_sub_gr(int * sol);

        /**
        * Update Lambda coeficients.
        */
        void update_lambda(double step);

        /**
        * return the value of the X_{col2,row2} - SUM Y_{col1,row1,col2,row2} col / node constraint.
        */
        inline int is_brk_cst_col_node(int col1, int col2, int row2, int * sol);
        /**
        * return the value of the X_{col2,row2} - SUM Y_{col1,row1,col2,row2} row / node constraint.
        */
        inline int is_brk_cst_row_node(int row1, int col2, int row2, int * sol);
        /**
        * return 1 if X_{col2,row2} - SUM Y_{col1,row1,col2,row2} = 0.
        */
        inline int is_brk_cst_row_node_act(int row1, int col2, int row2, int * sol);

        /**
        * return the value of the X_{col2,row2} - SUM Y_{col1,row1,col2,row2} col / node constraint.
        */
        inline int check_cst_col_node(int col1, int col2, int row2, int * sol);
        /**
        * return the value of the X_{col2,row2} - SUM Y_{col1,row1,col2,row2} row / node constraint.
        */
        inline int check_cst_row_node(int row1, int col2, int row2, int * sol);
        
        void recursive_brut_solve(int *sol, int col, int cur_row, int nb_row, int &checked);

    public :

        /**
        * Constructors : create a root problem.
        */
        problem_apurva(graph_apurva & g1, dp_mat_apurva & dp1, lambda_mat_apurva & lb1, double s1, double s2, double b_ub_score=-INFINITY, bool is_int=false) :
            problem(g1.get_nb_col(), g1.get_nb_row()),
            g(g1),
            dp(dp1),
            lb_mat(lb1),
        	self1(s1),
        	self2(s2),
        	bound_ub_score(b_ub_score),
        	obj_is_int(is_int)
            {cout << "creating problem apurva " << b_ub_score << endl;};
        /**
        * Constructors : create a sub problem of C, using lo1 and up1 as lower and upper limits.
        */
        problem_apurva(problem_apurva & c, int * lo1, int * up1) :
            problem(c,lo1,up1),
            g(c.get_graph()),
            dp(c.get_dp_mat_apurva()),
            lb_mat(c.get_lambda_mat()),
        	self1(c.get_self1()),
        	self2(c.get_self2()),
        	bound_ub_score(c.get_bound_ub_score()),
        	obj_is_int(c.get_obj_is_int())
            {};
        /**
        * Destructor
        */
        virtual ~problem_apurva(){};

        /**
        * Accessors
        */

        /**
        * Return the Lambda coeficients.
        */
        inline lambda_mat_apurva & get_lambda_mat(){return lb_mat;};
        /**
        * Return the Contact Map alignment graph.
        */
        inline graph_apurva & get_graph(){return g;};
        /**
        * Return the Dynamic programming object.
        */
        inline dp_mat_apurva & get_dp_mat_apurva(){return dp;};
        /**
         * Return the self similarity for protein 1 and 2
         */
        double get_self1(){return self1;};
        double get_self2(){return self2;};
        double get_bound_ub_score(){return bound_ub_score;};
        /**
         * Return the flag that denotes whether the objective function has integer value
         */
        double get_obj_is_int(){return obj_is_int;};

        /**
        * Methods
        */

        /**
        * Solve the problem using lagrangian relaxation, with Dual problem solved using Sub-Gradient Descent.
        */
        virtual void lr_sgd_solve(parameters &params);
        /**
        * Brut force solver.
        */
        virtual void brut_solve();
        /**
        * Being given a lower and an upper limits, create the corresponding sub-problem.
        */
        virtual problem * create_subproblem(int * lo1, int * up1);
        /**
        * This function split a problem into two subproblem. In fact it generate the two corresponding sets of upper and lower limits.
        */
        virtual void split(int * lo1, int * up1, int * lo2, int * up2);
};

#endif
