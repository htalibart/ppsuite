#include "main.h"
using namespace std;
#include <unistd.h>

/****************************************
 * Accessors                            *
 ***************************************/
double branch_and_bound :: get_ub(){ return(ub); }
double branch_and_bound :: get_lb(){ return(lb); }
void branch_and_bound :: get_solution(int * my_solution, int * my_solution_insert_before)
{
    for(int s(0); s != solution_size; ++s)
        my_solution[s] = best_solution[s];

    for (int s(0); s != solution_size+1; ++s)
	    my_solution_insert_before[s] = best_solution_insert_before[s];
}

int branch_and_bound :: get_status(){ return(status); }
int branch_and_bound :: get_solution_status(){ return(solution_status); }
double branch_and_bound :: get_solve_time(){ return(solve_time); }
int branch_and_bound :: get_nb_prb_in_list(){ return(problem_list.size()); }
int branch_and_bound :: get_nb_visited_nodes(){ return(visited_nodes); }
void branch_and_bound :: set_ub(double bound){ ub = bound; }

/****************************************
 * Private functions                    *
 ***************************************/

void branch_and_bound :: insert_prb_in_list(problem * my_problem)
{
    int inserted(0);
    double new_lb(my_problem->get_lb_father());

    //easy case, empty list
    if(problem_list.empty())
        problem_list.push_front(my_problem);

    //general case
    else
    {
        list<problem*>::iterator it_list = problem_list.begin();

        while(it_list != problem_list.end())
        {
            double current_lb((*it_list)->get_lb_father());

            if(new_lb <= current_lb)
            {
                //place found
                problem_list.insert(it_list,my_problem);
                it_list=problem_list.end();
                inserted = 1;
            }
            else
                ++it_list;
        }
        if(inserted == 0) problem_list.push_back(my_problem);
    }
}

int branch_and_bound :: sparse_list()
{
    /* To be done */
    return(0);
}

void branch_and_bound :: update_ub(problem * local_problem)
{
    double local_ub(local_problem->get_ub());
    int * local_solution = local_problem->get_solution();
    int * local_solution_insert_before = local_problem->get_solution_insert_before();

    //If the subproblem contain a better solution
    if(local_ub < ub)
    {
        //Save the new best solution
        ub = local_ub;
        for(int s(0); s != solution_size; ++s)
	{
            best_solution[s] = local_solution[s];
	}
        for(int s(0); s != solution_size+1; ++s)
	{
	    best_solution_insert_before[s] = local_solution_insert_before[s];
	}
    }
}

void branch_and_bound :: update_lb(problem * local_problem)
{
    double local_lb(local_problem->get_lb());
    double stack_lb;     //LB in the problem_list

    //Update Lowerbound
    //IW
    if(!problem_list.empty())
    {
    	problem * best_problem = problem_list.front();
    	stack_lb = best_problem->get_lb_father();

    	if(stack_lb > ub)
    	{
    		cout << "Discard remaining " << problem_list.size() << " nodes because they are useless." << endl;
    		visited_nodes += problem_list.size();
    		while(problem_list.size()>0)
    		{
    			problem* del_prob= problem_list.front();
    			problem_list.pop_front();
    			delete del_prob;
    		}
    	}
    }

    if(problem_list.empty())
    {
    	if(local_lb > ub)         //fix : when the ub is not found on the last node
    	{
    		stack_lb = ub;
    	}
    	else
    	{
    		stack_lb = local_lb;
    	}
    }

    if(MIN(stack_lb, local_lb) > lb)
    {
        lb = MIN(stack_lb, local_lb);
    }
}





/****************************************
 * Public functions                     *
 ***************************************/

void branch_and_bound :: print_stats()
{;}

void branch_and_bound :: solve(problem & root, parameters & my_param)
{
    long tic_per_sec = sysconf(_SC_CLK_TCK);
    struct tms start, end;
    int tic1(times(&start));

    solution_size = root.get_size();
    best_solution = new int [solution_size];
    fill_n(best_solution, solution_size, static_cast<int>(-1));
    best_solution_insert_before = new int [solution_size+1];
    fill_n(best_solution_insert_before, solution_size+1, static_cast<int>(0));


    if(status == 3)
    {
        cout << "Branch and Bound already solved\n";
        exit(0);
    }

    /**************************************
    * Step 1 : Root                       *
    *                                     *
    **************************************/

    my_param.my_time_limit = my_param.time_limit;

    root.lr_sgd_solve(my_param);

    //update BandB informations
    ub = min(ub,root.get_ub());
    lb = root.get_lb();

    cout << "update BandB informations: ub="<< ub <<", lb="<<lb<<endl;

    int * root_solution = root.get_solution();
    for(int l(0); l != solution_size; ++l)
        best_solution[l] = root_solution[l];
    int * root_solution_insert_before = root.get_solution_insert_before();
    for(int l(0); l != solution_size+1; ++l)
        best_solution_insert_before[l] = root_solution_insert_before[l];

    int tic2 = times(&end);
    solve_time = ((double)tic2 - (double)tic1) / (double)tic_per_sec;

    if(root.get_status() == OPTIMAL || root.get_status() == EPSILON)
    {
        solution_status = status;
        return;
    }
    if(solve_time > my_param.time_limit)
    {
        solution_status = status;
        return;
    }
    if(my_param.max_node == 0)
    {
        solution_status = status;
        return;
    }
    if (root.get_status() == NOT_SIMILAR)
    {
	    solution_status = NOT_SIMILAR;
	    return;
    }

    /**************************************
    * Step 2 : Branch and Bound           *
    *                                     *
    **************************************/

    int * up1 = new int [solution_size];
    int * lo1 = new int [solution_size];
    int * up2 = new int [solution_size];
    int * lo2 = new int [solution_size];


    root.split(lo1, up1, lo2, up2);
    if(!root.is_splitted())
    {
        cout << "unable to split Root problem\n !";
        exit(0);
    }

    problem * subp1(root.create_subproblem(lo1, up1));
    problem * subp2(root.create_subproblem(lo2, up2));

    subp1->split(lo1, up1, lo2, up2);
    if(!subp1->is_splitted())
    {
        cout << "unable to split subp1 problem\n !";
        exit(0);
    }
    problem * subpA(subp1->create_subproblem(lo1, up1));
    problem * subpB(subp1->create_subproblem(lo2, up2));

    subp2->split(lo1, up1, lo2, up2);
    if(!subp2->is_splitted())
    {
        cout << "unable to split subp2 problem\n !";
        exit(0);
    }
    problem * subpC(subp2->create_subproblem(lo1, up1));
    problem * subpD(subp2->create_subproblem(lo2, up2));


    insert_prb_in_list(subpA);
    insert_prb_in_list(subpB);
    insert_prb_in_list(subpC);
    insert_prb_in_list(subpD);

    delete subp1;
    delete subp2;

    while(solution_status == 0)
    {
        /**************************************
        * Step 2.1 : Solve best sub-problem   *
        *                                     *
        **************************************/


	cout << "solution_status=" << solution_status << endl;

        my_param.limit_lb = ub;     //to avoid spending too much time into useless subproblems
	cout << "limit_lb was set to ub=" << ub << endl;
        my_param.my_time_limit = my_param.time_limit - solve_time;

        problem * current_problem(problem_list.front());
        problem_list.pop_front();
        current_problem->lr_sgd_solve(my_param);

        //update upper_bound (and possibly get a new best solution) and lower_bound
        update_ub(current_problem);
        update_lb(current_problem);
	//cout << "lb was updated to " << lb << endl;
	//cout << "ub was updated to " << ub << endl;

        /**************************************
        * Step 2.2 : Divide into subproblem   *
        *            (if needed)              *
        **************************************/
        if( (ub-lb) > my_param.epsilon && (current_problem->get_lb() < ub) && (current_problem->get_status() != OPTIMAL)  && (current_problem->get_status() != NOT_SIMILAR))
        {
	    cout << "Dividing into subproblems" << endl;
            current_problem->split(lo1, up1, lo2, up2);
            if(!current_problem->is_splitted())
            {
                cout << "unable to split sub-problem\n !";
                cout << "Launch brut force solver !\n";
                int tic3 = times(&end);
                solve_time = ((double)tic3 - (double)tic1) / (double)tic_per_sec;

                cout << "ub = " << -lb << " lb=" << -ub << " gap = " << (lb - ub)/ lb << " time_mip = " << solve_time << " n1= " << 0 << " n2= " << 0 << " narcs1= " << 0 << ", narcs2= " << 0 << ", iter= 0" << "\n";
                current_problem->brut_solve();
                update_ub(current_problem);
                update_lb(current_problem);
            }
            else
            {
             subp1 = current_problem->create_subproblem(lo1, up1);
             subp2 = current_problem->create_subproblem(lo2, up2);
             
             subp1->split(lo1, up1, lo2, up2);
             if(!subp1->is_splitted())
             {
                 cout << "unable to split subp1 problem\n !";
                 exit(0);
             }
             subpA = subp1->create_subproblem(lo1, up1);
             subpB = subp1->create_subproblem(lo2, up2);

             subp2->split(lo1, up1, lo2, up2);
             if(!subp2->is_splitted())
             {
                 cout << "unable to split subp2 problem\n !";
                 exit(0);
             }
             subpC = subp2->create_subproblem(lo1, up1);
	     
             subpD = subp2->create_subproblem(lo2, up2);


             insert_prb_in_list(subpA);
             insert_prb_in_list(subpB);
             insert_prb_in_list(subpC);
             insert_prb_in_list(subpD);

             delete subp1;
             delete subp2;
            }
        }
        

        delete current_problem;

        int tic3 = times(&end);
        solve_time = ((double)tic3 - (double)tic1) / (double)tic_per_sec;

        visited_nodes++;

        /**************************************
        * Step 2.3 : Check end conditions     *
        *                                     *
        **************************************/
        if(lb == ub)
        {
            lb = ub;
            cout << "Optimal Solution found\n";
            solution_status = OPTIMAL;
        }
        else if((ub-lb) <= my_param.epsilon)
        {
            cout << "Solution found is in Epsilon range\n";
            solution_status = EPSILON;
        }
        else if(solve_time >= my_param.time_limit)
        {
            cout << "Time-limit : Solve process ended\n";
            solution_status = APPROXIMATE;
        }
        else if(visited_nodes >= my_param.max_node)
        {
            cout << "Node-limit : Solve process ended\n";
            solution_status = APPROXIMATE;
        }
	else if (-lb < my_param.score_min)
	{
		std::cout << "NOT SIMILAR : " << -lb << " < " << my_param.score_min << std::endl;
		solution_status = NOT_SIMILAR;
	}
	
    }

    delete[] up1;
    delete[] up2;
    delete[] lo1;
    delete[] lo2;

    while(problem_list.size()>0)
    {
    	problem* del_prob= problem_list.front();
    	problem_list.pop_front();
    	delete del_prob;
    }
}











