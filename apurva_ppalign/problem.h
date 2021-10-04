#ifndef PROBLEM_1
#define PROBLEM_1
using namespace std;

/**
*  VIRTUAL Class Problem
*/

class problem
{
    protected :

        /**
        * Size of the problem, in terms of columns.
        */
        int size;
        /**
        * Upper limits of the problem.
        */
        int * up;
        /**
        * Lower limits of the problem.
        */
        int * lo;
        /**
        * Current feasible solution
        */
        int * solution;
	/**
	 * Insertions in current feasible solution
	 */
	int * solution_insert_before;
	/**
        * Current best feasible solution, associated with ub
        */
	int * best_solution;
	/*
	 * Insertion in current best feasible solution
	 */
	int * best_solution_insert_before;
        /**
        * Current best upperbound of the problem, associated with best solution
        */
        double ub;
        /**
        * Current best lowerbound of the problem
        */
        double lb;
        /**
        * Lowerbound of father problem.
        */
        double lb_father;
        /**
        * Upperbound of father problem.
        */
        double ub_father;

        /**
        * Status of the problem
        */
        int status;
        /**
        * Boolean, =1 if the problem was divided without problem
        */
        int is_split;
        /**
        * Time taken by the resolution process
        */
        double solve_time;
        /**
        * Number of iteration of the resolution process
        */
        int iter;

    public :

        /**
        * Constructors
        */

        /**
        * Build an empty problem.
        */
        problem() : size(0),ub(INFINITY),lb(-INFINITY),status(NOT_INITIALIZED),is_split(0) {};

        /**
        * Build a root problem for a graph of nb_col columns and nb_row rows.
        */
        problem(int nb_col, int nb_row):
            size(nb_col),
            ub(INFINITY),
            lb(-INFINITY),
            lb_father(-INFINITY),
            ub_father(INFINITY),
            status(INITIALIZED),
            is_split(0),
            solve_time(0),
            iter(0)
            {
                up = new int [nb_col];
                fill_n(up, nb_col, static_cast <int> (nb_row - 1));
                lo = new int [nb_col];
                fill_n(lo, nb_col, static_cast <int> (0));
                solution = new int [nb_col];
		best_solution = new int [nb_col];
                fill_n(solution, nb_col, static_cast <int> (-1));
                fill_n(best_solution, nb_col, static_cast <int> (-1));
                solution_insert_before = new int [nb_col+1];
		best_solution_insert_before = new int [nb_col+1];
                fill_n(solution_insert_before, nb_col+1, static_cast <int> (0));
                fill_n(best_solution_insert_before, nb_col, static_cast <int> (0));
            };

        /**
        * Build a sub-problem of problem p, using lo1 for lowerlimits and up1for upperlimits.
        */
        problem(problem & p, int * lo1, int * up1):
            size(p.get_size()),
            ub(INFINITY),
            lb(p.get_lb()),
            lb_father(p.get_lb()),
            ub_father(p.get_ub()),
            status(INITIALIZED),
            is_split(0),
            solve_time(0),
            iter(0)
            {
                up = new int [size];
                lo = new int [size];
                solution = new int [size];
                solution_insert_before = new int [size+1];
		best_solution = new int [size];
		best_solution_insert_before = new int [size+1];
                for(int i(0); i != size; ++i)
                {
                    up[i] = up1[i];
                    lo[i] = lo1[i];
                }
                fill_n(solution, size, static_cast <int> (-1));
                fill_n(solution_insert_before, size+1, static_cast <int> (0));
		fill_n(best_solution, size, static_cast <int> (-1));
		fill_n(best_solution_insert_before, size+1, static_cast <int> (0));
            };

        /**
        * Destructor
        */
        virtual ~problem()
        {
            delete [] up;
            delete [] lo;
            delete [] solution;
            delete [] solution_insert_before;
	    delete [] best_solution;
	    delete [] best_solution_insert_before;
        };

        /**
        * Accessors
        */

        /**
        * Return the size of the problem.
        */
        inline int get_size(){ return size; };
        /**
        * Return the upper limits of the problem.
        */
        inline int * get_up(){ return up; };
        /**
        * Return the lower limits of the problem.
        */
        inline int * get_lo(){ return lo; };
        /**
        * Return the best feasible solution of the problem.
        */
        inline int * get_solution(){ return best_solution; };
	/**
	 * Returns insertion positions for the best feasible solution of the problem
	 */
	inline int * get_solution_insert_before() { return best_solution_insert_before; };
        /**
        * Return the current upperbound of the problem.
        */
        inline double get_ub(){ return ub; };
        /**
        * Return the current lowerbound of the problem.
        */
        inline double get_lb(){ return lb; };
        /**
        * Return the lowerbound of the father problem.
        */
        inline double get_lb_father(){ return lb_father; };
        /**
        * Return the upperbound of the father problem.
        */
        inline double get_ub_father(){ return ub_father; };
        /**
        * Return the status of the problem (see symbols.h).
        */
        inline int get_status(){ return status; };
        /**
        * Boolean function, return 1 if the problem was successfully splitted, 0 else.
        */
        inline int is_splitted(){ return is_split; };
        /**
        * Return the time needed to solve the problem.
        */
        inline double get_solve_time(){ return solve_time; };
        /**
        * Return the number of iteration needed to solve the problem.
        */
        inline int get_iter(){ return iter; };

        /**
        * Methods
        */

        /**
        * Set the upperlimits of the problem.
        */
        inline void set_up(int col, int val){ up[col]=val; };
        /**
        * Set the lowerlimits of the problem.
        */
        inline void set_lo(int col, int val){ lo[col]=val; };
        /**
        * Problem Solver, using Lagrangian Relaxation, with Dual problem solved using SubGradient Descent.
        */
        virtual void lr_sgd_solve(parameters &) = 0;
        /**
        * Brut-Force Problem Solver.
        */
        virtual void brut_solve() = 0;
        /**
        * This function split a problem into two subproblem. In fact it generate the two corresponding sets of upper and lower limits.
        */
        virtual void split(int * lo1, int * up1, int * lo2, int * up2) = 0;
        /**
        * Being given a lower and an upper limits, create the corresponding sub-problem.
        */
        virtual problem * create_subproblem(int * lo1, int * up1) = 0;
};


#endif
