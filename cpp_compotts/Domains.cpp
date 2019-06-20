#include "Domains.h"


Domains::Domains()
{
	pdb_file = "";
	chain = 'A';
	dali_home = "";
	n = -1;
	num_dist = -1;
	cm = NULL;
}

Domains::Domains(string pdbfile, char ch, string dali, int len, int** contact_map)
{
	pdb_file = pdbfile;
	chain = ch;
	dali_home = dali;
	n = len;
	num_dist = -1;
	cm = contact_map;
	RunPuu();
	SparsifyCm();
	CompNumdist();
}

Domains::~Domains()
{
}

Domains::Domains(Domains const& d)
{
	pdb_file = d.pdb_file;
	chain = d.chain;
	dali_home = d.dali_home;
	n = d.n;
	num_dist = d.num_dist;
	cm = d.cm;
	boundaries = d.boundaries;
}

Domains& Domains::operator=(Domains const& d)
{
	pdb_file = d.pdb_file;
	chain = d.chain;
	dali_home = d.dali_home;
	n = d.n;
	num_dist = d.num_dist;
	cm = d.cm;
	boundaries = d.boundaries;
	return *this;
}

void Domains::RunPuu()
{

	// Execute dssp
	/*
	string command = dali_home+"/dssp "+pdb_file+" >"+pdb_file+".dssp 2>/dev/null";
	int ret = system(command.c_str());
	if(ret != 0)
	{
		cout << "Some problem executing command " << command;
		cout << " during dssp execution!"<< endl;
		exit(1);
	}
	*/

	// Execute puu
	string command = "perl "+dali_home+"/forwarder.pl xxxx"+chain+" "+pdb_file+" "+pdb_file+".dssp END | "+dali_home+"/puu >/dev/null";
	int ret = system(command.c_str());
	if(ret != 0)
	{
		cout << "Some problem executing command " << command;
		cout << " during puu execution!"<< endl;
		exit(1);
	}

	// Extract the domain boundaries
	ifstream rsl_file("domains.puu");
	bool start_reading = false;
	string line;
	vector<pair<int,int> > doms;
	while(getline(rsl_file,line))
	{
		if(start_reading)
		{
			if(line[0] == '/') // completely finished
			{
				boundaries.push_back(doms);
				break;
			}
			else if(line.substr(0,3) != "   ") // new domain
			{
				if(doms.size() != 0)
				{
					boundaries.push_back(doms);
				}
				doms.clear();
				doms.push_back(make_pair(atof(line.substr(51,4).c_str()),atof(line.substr(58,4).c_str())));
			}
			else // old domain continues
			{
				doms.push_back(make_pair(atof(line.substr(51,4).c_str()),atof(line.substr(58,4).c_str())));
			}

		}
		if(line[0] == '#') // start reading at next line
		{
			start_reading = true;
		}
	}

	int max_len = 0;
	for(boundaries_vec_it it = boundaries.begin(); it != boundaries.end(); ++it)
	{
		for(segment_it it_seg = it->begin(); it_seg != it->end(); ++it_seg)
		{
			if(it_seg->second > max_len)
			{
				max_len = it_seg->second;
			}
		}
	}

	if( max_len != n)
	{
		cout << "puu length " << max_len << " differs from length n=" << n << endl;
        n = -1;
        return;
		//exit(1);
	}


}

void Domains::SparsifyCm()
{

	for(boundaries_vec_it it1 = boundaries.begin(); it1 != boundaries.end(); ++it1)
	{
		for(boundaries_vec_it it2 = it1+1; it2 != boundaries.end(); ++it2)
		{
			for(segment_it it1_seg = it1->begin(); it1_seg != it1->end(); ++it1_seg)
			{
				for(segment_it it2_seg = it2->begin(); it2_seg != it2->end(); ++it2_seg)
				{
					for(int i = it1_seg->first; i <= it1_seg->second; ++i)
					{
						for(int j = it2_seg->first; j <= it2_seg->second; ++j)
						{
							cm[i-1][j-1] = 0;
							cm[j-1][i-1] = 0;
						}
					}
				}
			}
		}
	}
}

void Domains::CompNumdist()
{
	num_dist = 0;
	for(int i=0; i<n; ++i)
	{
		for(int j=i+1; j<n; ++j)
		{
			if(cm[i][j] == 1)
			{
				++num_dist;
			}
		}
	}
}

int Domains::GetLen() const
{
	return n;
}

string Domains::GetPdbfile() const
{
	return pdb_file;
}

int** Domains::GetCm() const
{
	return cm;
}

string Domains::GetDalihome() const
{
	return dali_home;
}

boundaries_vec Domains::GetBoundaries() const
{
	return boundaries;
}

int Domains::GetNumdist() const
{
	return num_dist;
}
