#include "main_infer_insertion_penalties.h"

double compute_No(double* delta_ins, int pos, int L, int nseq)
{
	double No=0;
	for (int n=0; n<nseq; n++)
	{
		No+=(delta_ins[n*(L+1)+pos]>0);
	}
	return No;

}

double compute_Nt(double* delta_ins, int pos, int L, int nseq)
{
	double Nt=0;
	for (int n=0; n<nseq; n++)
	{
		Nt+=delta_ins[n*(L+1)+pos];
	}
	return Nt;

}



void write_insertion_penalties_in_file(double* insertion_penalties_open, double* insertion_penalties_extend, int L, char* output_fname)
{
	ofstream output_file;
  	output_file.open(output_fname);

	for (int pos=0; pos<L+1; pos++)
	{
		output_file << insertion_penalties_open[pos] << '\t' << insertion_penalties_extend[pos] << endl;
	}	
	output_file.close();
}	



extern "C" int call_infer_insertion_penalties_from_python(double* delta_ins, int L, int nseq, int maxit_infer_insertions, double tol_infer_insertions, double learning_coeff_insertions, double freq_insert_min, double pc_insertions_tau, int delta_n_max, double expected_No, double expected_Nt, char* output_fname)
{
	int status = 0;


	double* insertion_penalties_open = new double[L+1];
	double* insertion_penalties_extend = new double[L+1];

	/*for (int n=0; n<nseq; n++)
	{
		for (int pos=0; pos<L+1; pos++)
		{

			cout << delta_ins[n*(L+1)+pos] << " ";
		}
		cout << endl;
	}*/


	// gap open and extends for a position without gap (stored to save time)
	double le0=-1;
	double lo0=-1;


	double dLdlo;
	double dLdle;


	cout << "nseq=" << nseq << endl;

	for (int pos=0; pos<L+1; pos++)
	{
		double eps=1e8;
		double lo=1.0;
		double le=1.0;
		int it=1;

		double No = compute_No(delta_ins, pos, L, nseq);
		double Nt = compute_Nt(delta_ins, pos, L, nseq);


		bool is_pos_without_gap;

		//if (No==0)
		if (No<freq_insert_min*nseq)
		{
			No=freq_insert_min*nseq;
			Nt=freq_insert_min*nseq;
			is_pos_without_gap=true;
		}
		else
		{
			is_pos_without_gap=false;
		}


		// pseudo-counts
		No = (1-pc_insertions_tau)*No+pc_insertions_tau*expected_No;
		Nt = (1-pc_insertions_tau)*Nt+pc_insertions_tau*expected_Nt;


		double fo = No/nseq;
		double ft = Nt/nseq;



		if ((is_pos_without_gap) && (lo0!=-1))
		{
			lo = lo0;
			le = le0;
		}
		else
		{
			while ( (eps > tol_infer_insertions) && (it < maxit_infer_insertions) )
			{
				//dLdle = No-Nt+nseq*(exp(-lo-le)/((1-exp(-le)+exp(-lo))))-2*le;
				//dLdle = fo-ft+exp(-lo-le)/((1-exp(-le)+exp(-lo)))-2*le/nseq;
				dLdle = ( exp(-lo - le) *pow(1.0 - exp(-le),-2) )/(1.0 + exp(-lo) * pow(1.0-exp(-le),-1)) - ft + fo - (2.0/nseq)*le;
				le = le + learning_coeff_insertions * dLdle;
				//dLdlo = -No+nseq*(exp(-lo)/(1-exp(-le)+exp(-lo)))-2*lo;
				//dLdlo = -fo+exp(-lo)/(1-exp(-le)+exp(-lo))-2*lo/nseq;
				dLdlo = ( exp(-lo) * pow(1.0-exp(-le),-1) ) / (1.0 + exp(-lo)* pow(1.0 -exp(-le),-1 )) - fo - (2.0/nseq)*lo;
				lo = lo + learning_coeff_insertions * dLdlo;
				eps = MAX(abs(dLdle), abs(dLdlo));
				it+=1;

			}
			cout << "pos=" << pos << ", No=" << No << ", Nt=" << Nt << ", fo=" << fo <<", ft=" << ft << ", lo=" << lo << ", le=" << le << ",it=" << it << ",eps=" << eps << endl;

			if (is_pos_without_gap)
			{
				lo0 = lo;
				le0 = le;
			}

		}

		insertion_penalties_open[pos]=lo;
		insertion_penalties_extend[pos]=le;


	}


	write_insertion_penalties_in_file(insertion_penalties_open, insertion_penalties_extend, L, output_fname);



    	delete[] insertion_penalties_open;
    	delete[] insertion_penalties_extend;
	
	return status;
}
