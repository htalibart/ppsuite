#ifndef PARAMETERS_1
#define PARAMETERS_1
using namespace std;

/**
* Parameters Class.
* Contain all parameters needed to run a Branch and bound solver.
*/

class parameters
{
    public :
        /**
        * Maximum iteration number of sub-gradient descent for each problem / sub-problem.
        */
        int max_iteration;
        /**
        * Maximum visited node number during the Branch and Bound process.
        */
        int max_node;
        /**
        * Maximum computation time allowed for the Branch and Bound process.
        */
        double time_limit;
        /**
        * Maximum computation time allowed for the Branch and Bound process within the CURRENT NODE.
        */
        double my_time_limit;
        /**
        * Relative gap value allowed. For exact solving, epsilon = 0 (default).
        */
        double epsilon;

        /**
        * The limit of the lower bound; used for pruning; if this is less than the value
        * of the best feasible solution found so far, computation for this B&B node can be stopped
        */
        double limit_lb;

        /**
        * Gamma parameters, needed for the computation of the lambda update step
        */
        double gamma;
        /**
        * Theta parameters, used to dynamicaly change Gamma during Sub-Gradient Descent process.
        */
        double theta;
	/**
        * Stepsize_min parameter, used during Sub-Gradient Descent process.
        */
        double stepsize_min;
	/**
        * Max nb of non-increasint steps during Subg-Gradient Descent process.
        */
	int nb_non_increasing_steps_max;



        /**
        * Default Constructor, usualy for root problem
        */
        parameters():
                   max_iteration(4000),
                   max_node(1000000),
                   time_limit(1800.0),
                   my_time_limit(1800.0),
                   epsilon(0.0),
                   limit_lb(INFINITY),
                   gamma(1.),
                   theta(0.9), /*0.5*/
		   stepsize_min(0.000000005),
		   nb_non_increasing_steps_max(500)
                   {};

        /**
        * Constructor not used anymore, should be removed...
        */
        parameters(int iter, int node, double time, double eps, double ub=INFINITY) :
                max_iteration(iter),
                max_node(node),
                time_limit(time),
                my_time_limit(1800.0),
                epsilon(eps),
                limit_lb(ub),
                gamma(1.),
                theta(0.5)
                {};

        /**
        * Destructor
        */
        ~parameters(){};
};

#endif
