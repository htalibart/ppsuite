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
					//w_rescaled_flat[a+b*q+i*q*q+j*q*q*L] = w_rescaled_flat[b+a*q+j*q*q+i*q*q*L];
					w_rescaled_flat[b+a*q+i*q*q+j*q*q*L] = w_rescaled_flat[b+a*q+j*q*q+i*q*q*L];
				}
			}

		}
	}


	return status;
}


float get_sum_over_a(float* v_flat, int q, int i)
{
	float sum=0;
	for (int a=0; a<q; a++)
	{
		sum+=v_flat[a+i*q];
	}
	return sum;
}

void expvi_to_fi(float* v_flat, int q, int i, float t, float sum_expvi)
{
	for (int a=0; a<q; a++)
	{
		v_flat[a+i*q] = log((1-t)*v_flat[a+i*q]/sum_expvi + t/q);
	}
}


extern "C" int cpp_rescale_v(float* v_flat, float* v_rescaled_flat, int L, int q, float v_rescaling_tau)
{
	int status = 0;
	
	int len_v_flat = L*q;
	expbw(v_flat, len_v_flat, 1);

	for (int i=0; i<L; i++)
	{
		float sum_expvi = get_sum_over_a(v_flat, q, i);
		expvi_to_fi(v_flat, q, i, v_rescaling_tau, sum_expvi);
		float sum_fi = get_sum_over_a(v_flat, q, i);
		for (int a=0; a<q; a++)
		{
			v_rescaled_flat[a+i*q] = v_flat[a+i*q]-sum_fi/q;
		}
	}

	return status;
}


