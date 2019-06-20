/*
 * @file Domains.h
 * @brief Runs puu, a program for domain identification, reads the domains from
 * the puu result file and can adjust contact maps such that distances between
 * different domains are not considered (set to 0)
 * @author: Inken Wohlers (Inken.Wohlers@cwi.nl)
 */

#ifndef DOMAINS_H_
#define DOMAINS_H_

#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ostream>
#include <values.h>
#include <math.h>
#include <list>


using namespace std;

typedef vector<vector<pair<int,int> > > boundaries_vec;
typedef vector<vector<pair<int,int> > >::iterator boundaries_vec_it;
typedef vector<pair<int,int> >::iterator segment_it;

class Domains{

private:

	/// The pdb file
	string pdb_file;

	char chain;

	/*
	 * The location of the DaliLite Bin directory
	 * in which the dssp and the puu binary are located
	 */
	string dali_home;

	/// Length of the protein
	int n;

	/// The number of remaining distances (non-zero entries in cm)
	int num_dist;

	/// Contact map from protein
	int** cm;

	boundaries_vec boundaries;


private:

	/*
	 * Runs dssp;
	 * Uses dssp output file and pdb file in order to
	 * run puu (program for domain identification)
	 * Reads the domain boundaries from the puu result file and writes them to boundaries
	 */
	void RunPuu();

	/*
	 * Given an arbitrary contact map from initialization, set all entries from this
	 * contact map to 0 if the distances correspond to inter-domain distances (distances
	 * between two different domains)
	 */
	void SparsifyCm();

	/// Goes over cm and counts the number of non-zero entries
	void CompNumdist();

public:

	/// Default constructor
	Domains();

	/// Constructor
	Domains(string inputfile, char ch, string dalih, int len, int** contact_map);

	/// Destruktor
	~Domains();

	/// Copy constructor
	Domains(Domains const& dom);

	/// Assignment operator
	Domains& operator=(Domains const& dom);

	/// Returns the length of the protein
	int GetLen() const;

	/// Returns the path to theDaliLite Bin directory in which are binaries "dssp" and "puu"
	string GetDalihome() const;

	/// Returns the path and filename of the input file (pdb or dssp)
	string GetPdbfile() const;

	/*
	 * Returns the contact map of the protein;
	 * 0-entries are inter-domain distances
	 * 1-entries are intra-domain distances
	 */
	int** GetCm() const;

	/*
	 * Returns the domain boundaries;
	 * First dimension is the number of domains
	 * Second dimension is the number of segments for each corresponding domain
	 */
	boundaries_vec GetBoundaries() const;

	/// Returns the number of non-zero entries in cm, i.e. the number of remaining distances
	int GetNumdist() const;

};



#endif /* DOMAINS_H_ */

