#ifndef __B_AND_B__
#define __B_AND_B__

using namespace std;

/**
* Service Class Branch_and_Bound.
* A branch_and_bound object apply a Branch and Bound strategy to solve
* a given problem object.
*
* A branch_and_bound object can only be used to solve one problem !
*/

class branch_and_bound
{
    private :

        /**
        * List of sub-problem, ordered by decreasing lowerbound.
        */
        list <problem*> problem_list;

        /**
        * Size of a solution (Columns number in the graph).
        */
        int solution_size;

        /**
        * Status of the branch_and_bound object.
        */
        int status;
        /**
        * Time taken by the resolution process.
        */
        double solve_time;
        /**
        * Should be removed, still possibly used in the branch and bound but not needed...
        */
        int max_visited_depth;
        /**
        * Number of node visited during the resolution process.
        */
        int visited_nodes;
        /**
        * Status of the solution (see symbols.h).
        */
        int solution_status;

        /**
        * best upperbound obtained during the Branch and Bound process.
        */
        double ub;
        /**
        * best lowerbound obtained during the Branch and Bound process.
        */
        double lb;
        /**
        * best feasible solution obtained during the Branch and Bound process.
        */
        int * best_solution;
        int * best_solution_insert_before;

        /**
        * private functions
        */

        /**
        * Insert a sub-problem into the sub-problem list, respecting the lowerbound order.
        */
        void insert_prb_in_list(problem * my_problem);
        /**
        * Remove from the sub-problem list all sub-problem with a lowerbound bigger than the current best upperbound :
        * ie, fathom sub-problems by value dominance.
        */
        int sparse_list();
        /**
        * Update the best upperbound value, taking into account the results found in local_problem.
        * If the local_problem contain a better feasible solution, also update solution.
        */
        void update_ub(problem * local_problem);
        /**
        * Update the best lowerbound value, taking into account the results found in local_problem
        */
        void update_lb(problem * local_problem);
        //void clear_list();

    public :

        /**
        * Constructors
        */
        branch_and_bound()
        {
            solution_size = 0;
            status = NOT_INITIALIZED;
            solution_status = NOT_INITIALIZED;
            solve_time = 0.;
            max_visited_depth = 0;
            visited_nodes = 0;
            ub =  INFINITY;
            lb = -INFINITY;
            best_solution = 0;
		best_solution_insert_before = 0;
        };

        /**
        * Destructor
        */
        ~branch_and_bound()
        {
            delete [] best_solution;
            delete [] best_solution_insert_before;
            //should clear the list also ;-)
        };

        /**
        * Accessor
        */

        /**
        * Return the best upperbound found during the Branch and Bound process
        */
        double get_ub();
        /**
        * Return the best lowerbound found during the Branch and Bound process
        */
        double get_lb();
        /**
        * Return the best feasible solution found during the Branch and Bound process
        */
        void get_solution(int * my_solution, int * my_solution_insert_before);
        /**
        *
        */
        int get_status();
        
        /**
        *
        **/
        int get_solution_status();
        
        /**
        * Return the number of sub-problem in the sub-problem list
        */
        int get_nb_prb_in_list();
        /**
        * Return the time needed by the branch and bound resolution process.
        */
        double get_solve_time();
        /**
        * Return the number of visited nodes during the branch and bound process.
        */
        int get_nb_visited_nodes();

        /**
        * Solve the root problem, using a Branch and Bound strategy, accordingly to the given parameters
        */
        void solve(problem & root, parameters & my_param);
        /**
        * Print some stats about the run...
        */
        void print_stats();

        /**
         * Set the ub in case a feasible (heuritic) solution is known
         * (ub and not lb, because a_purva minimizes instead of maximizes)
         */
        void set_ub(double bound);
};

#endif
