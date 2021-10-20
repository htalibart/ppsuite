#include "cpp_rescale.h"



void expbw(float* w_flat, int len_w_flat, float beta_softmax_w)
{
	for (int k=0; k<len_w_flat; k++)
	{
		w_flat[k] = exp(beta_softmax_w * w_flat[k]);
	}
}


float get_sum_over_a_b(float* w_flat, int L, int q, int i, int j)
{
	float sum=0;
	for (int a=0; a<q; a++)
	{
		for (int b=0; b<q; b++)
		{
			sum+=w_flat[b+a*q+j*q*q+i*q*q*L];
		}
	}
	return sum;
}


void expbwij_to_fij(float* w_flat, int L, int q, int i, int j, float t, float sum_expbw)
{
	for (int a=0; a<q; a++)
	{
		for (int b=0; b<q; b++)
		{
			w_flat[b+a*q+j*q*q+i*q*q*L] = log((1-t)*w_flat[b+a*q+j*q*q+i*q*q*L]/sum_expbw + t/(q*q));
		}
	}
}


extern "C" int cpp_rescale_w(float* w_flat, float* w_rescaled_flat, int L, int q, float w_rescaling_tau, float beta_softmax_w)
{
	int status = 0;
	
	int len_w_flat = L*L*q*q;
	expbw(w_flat, len_w_flat, beta_softmax_w);

	for (int i=0; i<L-1; i++)
	{
		for (int j=i+1; j<L; j++)
		{
			float sum_expbw	= get_sum_over_a_b(w_flat, L, q, i, j);
			expbwij_to_fij(w_flat, L, q, i, j, w_rescaling_tau, sum_expbw);
			float sum_fij = get_sum_over_a_b(w_flat, L, q, i, j);
			for (int a=0; a<q; a++)
			{
				for (int b=0; b<q; b++)
				{
					w_rescaled_flat[b+a*q+j*q*q+i*q*q*L] = (1/beta_softmax_w)*(w_flat[b+a*q+j*q*q+i*q*q*L] - sum_fij/(q*q));
					w_rescaled_flat[b+a*q+i*q*q+j*q*q*L] = w_rescaled_flat[b+a*q+j*q*q+i*q*q*L];
				}
			}

		}
	}


	return status;
}
