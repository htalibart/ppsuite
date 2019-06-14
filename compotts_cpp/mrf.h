#ifndef CONTACT_MAP_PAUL
#define CONTACT_MAP_PAUL


/*
 *  mrf.h
 *  
 *
 */



#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<math.h>
#include<string.h>
#include <msgpack.hpp>
#include <boost/any.hpp>
#include "util.h"

using namespace std;


class MRF
{
private:	
	int q = 21; // nombre d'acides aminés
	int nb_nodes; // nombre de noeuds du MRF
	int nb_edges;
	double **** w; // les wijab
	double ** v; // les via
	int ** edges_map; // carte des arcs (1 si le coupling est jugé significatif, 0 sinon)
//	string three_to_one_letter(string);

public:
	MRF(){};
	~MRF() // à changer
	{
		for (int i=0; i<nb_nodes ; ++i)
		{
			delete[] v[i];
			for (int j=0; j<nb_nodes; j++)
			{
				for (int a=0; a<q; a++)
				{
					delete[] w[i][j][a];
				}
				delete[] w[i][j];
			}
			delete[] w[i];
			delete[] edges_map[i];
		}
		delete[] v;
		delete[] w;
		delete[] edges_map;
	};

	int load_mrf_from_msgpack(char* filename, bool use_w, char* threshold_choice, double threshold_arg, int rescale, double rescaling_coeff_v, double rescaling_coeff_w);
	void init_edges_map(bool use_w, char* threshold_choice, double threshold_arg);
	void display_MRF(); // for debugging purposes
	
	//there should be accessors
	int get_nb_nodes(){return(nb_nodes);};
	int get_nb_edges(){return(nb_edges);};
	double **** get_w(){return(w);};
	double ** get_v(){return(v);};
	int ** get_edges_map(){return(edges_map);};
	int get_nb_aa(){return q;};
};







void MRF::init_edges_map(bool use_w, char* threshold_choice, double threshold_arg)
{
	edges_map = new int*[nb_nodes];
	nb_edges=0;
	double w_threshold;
	if(strcmp(threshold_choice, "manual_threshold")==0)
	{
		w_threshold = threshold_arg;
	}
	else if(strcmp(threshold_choice, "mean_fraction")==0)
	{
		double mean_frobenius_norm = 0.0;
		for (int i=0; i<nb_nodes; i++)
		{
			for (int j=0; j<nb_nodes; j++)
			{
				mean_frobenius_norm+=frobenius_norm(w[i][j], q);
			}
		}
		mean_frobenius_norm = mean_frobenius_norm/(nb_nodes*nb_nodes);
		w_threshold = threshold_arg*mean_frobenius_norm;
	}
	for (int i=0; i<nb_nodes; i++)
	{
		edges_map[i] = new int[nb_nodes];
		if (use_w)
		{
			for (int j=0; j<nb_nodes; j++)
			{
				double wij = frobenius_norm(w[i][j], q);
				if (wij>w_threshold)
				{
					edges_map[i][j]=1;
				}
				else
				{
					edges_map[i][j]=0;
				}
				nb_edges+=edges_map[i][j];
			}
		}
		else
		{
			for (int j=0; j<nb_nodes; j++)
			{
				edges_map[i][j]=0;
			}

		}
	}
}


void MRF::display_MRF()
{
	std::cout << "DISPLAYING MRF" << std::endl;
	for (int i=0; i<nb_nodes; i++)
	{
		for (int a=0; a<q; a++)
		{
			std::cout << "v[" << i << "][" << a << "]=" << v[i][a] << ", ";
		}
		std::cout << std::endl;
	}

	for (int i=0; i<nb_nodes; i++)
	{
		for (int j=0; j<nb_nodes; j++)
		{
			std::cout << edges_map[i][j];
		}
		std::cout << endl;
	}
	std::cout << endl << endl;	

	/*for (int i=0; i<nb_nodes; i++)
	{
		for (int j=0; j<nb_nodes; j++)
		{
			double sumw = 0;
			for (int a=0; a<q; a++)
			{
				for (int b=0; b<q; b++)
				{
					sumw+=w[i][j][a][b];
				}
			}
			if ( (sumw!=0) && (i!=j) )
			{
				std::cout << "(" << i << "," << j << ")" << std::endl;
			}
		} 
	}*/
	
	
}


int MRF::load_mrf_from_msgpack(char* filename, bool use_w, char* threshold_choice, double threshold_arg, int rescale, double rescaling_coeff_v, double rescaling_coeff_w)
{
	// récupération de l'objet msgpack
	std::ifstream ifs(filename, std::ifstream::in);
        std::stringstream buffer;
        buffer << ifs.rdbuf();
        msgpack::unpacked msg;
        msgpack::unpack(msg, buffer.str().data(), buffer.str().size());
	msgpack::object obj = msg.get();
	// transformation en map {"meta" : [truc], "ncol" : [truc], etc.}
	msgpack::object_map obj_map = obj.via.map;

	// objets temporaires pour stocker les trucs pendant le parcours de la map
	std::vector<double> x_single_obj;
	msgpack::object_map obj_pair;

	// parcours de cette map
	int SIZE_MAP = 4;
	for (int i=0; i<SIZE_MAP; i++)
	{
		std::string map_key;
		(obj_map.ptr[i].key).convert(map_key);

		if (map_key=="ncol")
		{
			(obj_map.ptr[i].val).convert(nb_nodes); // on récupère le nombre de colonnes dans nbcol
		}

		else if (map_key=="x_single")
		{
			(obj_map.ptr[i].val).convert(x_single_obj); // on récupère un vecteur contenant les via dans une forme pas très pratique dans x_single_obj
			
		}

		else if (map_key=="x_pair")
		{
			obj_pair = (obj_map.ptr[i].val).via.map; // on récupère une map contenant les wijab dans une forme pas encore pratique dans obj_pair
		}
	}



	// on met les via récupérés dans v
	v = new double*[nb_nodes];
	for (int k=0; k<nb_nodes; k++)
	{
		v[k] = new double[q]; 
		std::vector<double> x_single_k;
		for (int a=0; a<q-1; a++)
		{
			v[k][a] = x_single_obj[k*(q-1)+a];
		}
		v[k][q-1] = 0;
	}


	// on met les wijab récupérés dans w
	
	w = new double***[nb_nodes];
	for (int i=0; i<nb_nodes; i++)
	{
		w[i] = new double**[nb_nodes];
		for (int j=0; j<nb_nodes; j++)
		{
			w[i][j] = new double*[q];
			for (int a=0; a<q; a++)
			{
				w[i][j][a] = new double[q];
				for (int b=0; b<q; b++)
				{
					w[i][j][a][b] = 0;
				}
			}
		}
	}
	for (int k=0; k<obj_pair.size; k++)
	{
		int ind_i;
		int ind_j;
		std::vector<double> x;
		msgpack::object_map obj_x = (obj_pair.ptr[k].val).via.map;
		(obj_x.ptr[0].val).convert(ind_i);
		try
		{
			(obj_x.ptr[1].val).convert(ind_j);
		}
		catch(const std::bad_cast &e)
		{
			(obj_x.ptr[1].val).convert(x);
		}
		try
		{
			(obj_x.ptr[2].val).convert(x);
		}
		catch(const std::bad_cast &e)
		{
			(obj_x.ptr[2].val).convert(ind_j);
		}
		for (int a=0; a<q; a++)
		{
			for (int b=0; b<q; b++)
			{
				w[ind_i][ind_j][a][b] = x[a*q+b];
				w[ind_j][ind_i][a][b] = x[a*q+b];
			}
		}
	}

	if (rescale)
	{
		cout << "rescale" << endl;
		for (int i=0; i<nb_nodes; i++)
		{
			rescale_vi(v[i], q, rescaling_coeff_v);
			for (int j=0; j<nb_nodes; j++)
			{
				rescale_wij(w[i][j], q, rescaling_coeff_w);
			}
		}
	}

	init_edges_map(use_w, threshold_choice, threshold_arg);
	//display_MRF();
	
	return(0);
}


#endif

