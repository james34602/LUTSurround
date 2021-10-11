#include "codelet.h"
void DFT512(float *A, const float *sinTab)
{
	int i, j, theta;
	float alpha, beta, beta1, beta2, y0, y1, y2, y3;
	for (i = 0; i < 512; i += 4)
	{
		alpha = A[i];
		beta = A[i + 1];
		beta1 = A[i + 2];
		beta2 = A[i + 3];
		y0 = alpha + beta;
		y1 = alpha - beta;
		y2 = beta1 + beta2;
		y3 = beta1 - beta2;
		A[i] = y0 + y2;
		A[i + 2] = y0 - y2;
		A[i + 1] = y1 + y3;
		A[i + 3] = y1 - y3;
	}
	for (i = 0; i < 512; i += 8)
	{
		alpha = A[i];
		beta = A[i + 4];
		A[i] = alpha + beta;
		A[i + 4] = alpha - beta;
		alpha = A[i + 2];
		beta = A[i + 6];
		A[i + 2] = alpha + beta;
		A[i + 6] = alpha - beta;
		alpha = A[i + 1];
		beta1 = 0.70710678118654752440084436210485f*(A[i + 5] + A[i + 7]);
		beta2 = 0.70710678118654752440084436210485f*(A[i + 5] - A[i + 7]);
		A[i + 1] = alpha + beta1;
		A[i + 5] = alpha - beta1;
		alpha = A[i + 3];
		A[i + 3] = alpha + beta2;
		A[i + 7] = alpha - beta2;
	}
	for (i = 0; i < 512; i += 16)
	{
		theta = 32;
		alpha = A[i];
		beta = A[i + 8];
		A[i] = alpha + beta;
		A[i + 8] = alpha - beta;
		alpha = A[i + 4];
		beta = A[i + 12];
		A[i + 4] = alpha + beta;
		A[i + 12] = alpha - beta;
		for (j = 1; j < 4; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 128];
			alpha = A[i + j];
			beta = A[i - j + 8];
			beta1 = A[i + j + 8] * y1 + A[i - j + 16] * y0;
			beta2 = A[i + j + 8] * y0 - A[i - j + 16] * y1;
			theta += 32;
			A[i + j] = alpha + beta1;
			A[i + j + 8] = alpha - beta1;
			A[i - j + 8] = beta + beta2;
			A[i - j + 16] = beta - beta2;
		}
	}
	for (i = 0; i < 512; i += 32)
	{
		theta = 16;
		alpha = A[i];
		beta = A[i + 16];
		A[i] = alpha + beta;
		A[i + 16] = alpha - beta;
		alpha = A[i + 8];
		beta = A[i + 24];
		A[i + 8] = alpha + beta;
		A[i + 24] = alpha - beta;
		for (j = 1; j < 8; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 128];
			alpha = A[i + j];
			beta = A[i - j + 16];
			beta1 = A[i + j + 16] * y1 + A[i - j + 32] * y0;
			beta2 = A[i + j + 16] * y0 - A[i - j + 32] * y1;
			theta += 16;
			A[i + j] = alpha + beta1;
			A[i + j + 16] = alpha - beta1;
			A[i - j + 16] = beta + beta2;
			A[i - j + 32] = beta - beta2;
		}
	}
	for (i = 0; i < 512; i += 64)
	{
		theta = 8;
		alpha = A[i];
		beta = A[i + 32];
		A[i] = alpha + beta;
		A[i + 32] = alpha - beta;
		alpha = A[i + 16];
		beta = A[i + 48];
		A[i + 16] = alpha + beta;
		A[i + 48] = alpha - beta;
		for (j = 1; j < 16; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 128];
			alpha = A[i + j];
			beta = A[i - j + 32];
			beta1 = A[i + j + 32] * y1 + A[i - j + 64] * y0;
			beta2 = A[i + j + 32] * y0 - A[i - j + 64] * y1;
			theta += 8;
			A[i + j] = alpha + beta1;
			A[i + j + 32] = alpha - beta1;
			A[i - j + 32] = beta + beta2;
			A[i - j + 64] = beta - beta2;
		}
	}
	for (i = 0; i < 512; i += 128)
	{
		theta = 4;
		alpha = A[i];
		beta = A[i + 64];
		A[i] = alpha + beta;
		A[i + 64] = alpha - beta;
		alpha = A[i + 32];
		beta = A[i + 96];
		A[i + 32] = alpha + beta;
		A[i + 96] = alpha - beta;
		for (j = 1; j < 32; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 128];
			alpha = A[i + j];
			beta = A[i - j + 64];
			beta1 = A[i + j + 64] * y1 + A[i - j + 128] * y0;
			beta2 = A[i + j + 64] * y0 - A[i - j + 128] * y1;
			theta += 4;
			A[i + j] = alpha + beta1;
			A[i + j + 64] = alpha - beta1;
			A[i - j + 64] = beta + beta2;
			A[i - j + 128] = beta - beta2;
		}
	}
	for (i = 0; i < 512; i += 256)
	{
		theta = 2;
		alpha = A[i];
		beta = A[i + 128];
		A[i] = alpha + beta;
		A[i + 128] = alpha - beta;
		alpha = A[i + 64];
		beta = A[i + 192];
		A[i + 64] = alpha + beta;
		A[i + 192] = alpha - beta;
		for (j = 1; j < 64; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 128];
			alpha = A[i + j];
			beta = A[i - j + 128];
			beta1 = A[i + j + 128] * y1 + A[i - j + 256] * y0;
			beta2 = A[i + j + 128] * y0 - A[i - j + 256] * y1;
			theta += 2;
			A[i + j] = alpha + beta1;
			A[i + j + 128] = alpha - beta1;
			A[i - j + 128] = beta + beta2;
			A[i - j + 256] = beta - beta2;
		}
	}
	theta = 1;
	alpha = A[0];
	beta = A[256];
	A[0] = alpha + beta;
	A[256] = alpha - beta;
	alpha = A[128];
	beta = A[384];
	A[128] = alpha + beta;
	A[384] = alpha - beta;
	for (j = 1; j < 128; j++)
	{
		y0 = sinTab[theta];
		y1 = sinTab[theta + 128];
		alpha = A[j];
		beta = A[-j + 256];
		beta1 = A[j + 256] * y1 + A[-j + 512] * y0;
		beta2 = A[j + 256] * y0 - A[-j + 512] * y1;
		theta += 1;
		A[j] = alpha + beta1;
		A[j + 256] = alpha - beta1;
		A[-j + 256] = beta + beta2;
		A[-j + 512] = beta - beta2;
	}
}
void DFT1024(float *A, const float *sinTab)
{
	int i, j, theta;
	float alpha, beta, beta1, beta2, y0, y1, y2, y3;
	for (i = 0; i < 1024; i += 4)
	{
		alpha = A[i];
		beta = A[i + 1];
		beta1 = A[i + 2];
		beta2 = A[i + 3];
		y0 = alpha + beta;
		y1 = alpha - beta;
		y2 = beta1 + beta2;
		y3 = beta1 - beta2;
		A[i] = y0 + y2;
		A[i + 2] = y0 - y2;
		A[i + 1] = y1 + y3;
		A[i + 3] = y1 - y3;
	}
	for (i = 0; i < 1024; i += 8)
	{
		alpha = A[i];
		beta = A[i + 4];
		A[i] = alpha + beta;
		A[i + 4] = alpha - beta;
		alpha = A[i + 2];
		beta = A[i + 6];
		A[i + 2] = alpha + beta;
		A[i + 6] = alpha - beta;
		alpha = A[i + 1];
		beta1 = 0.70710678118654752440084436210485f*(A[i + 5] + A[i + 7]);
		beta2 = 0.70710678118654752440084436210485f*(A[i + 5] - A[i + 7]);
		A[i + 1] = alpha + beta1;
		A[i + 5] = alpha - beta1;
		alpha = A[i + 3];
		A[i + 3] = alpha + beta2;
		A[i + 7] = alpha - beta2;
	}
	for (i = 0; i < 1024; i += 16)
	{
		theta = 64;
		alpha = A[i];
		beta = A[i + 8];
		A[i] = alpha + beta;
		A[i + 8] = alpha - beta;
		alpha = A[i + 4];
		beta = A[i + 12];
		A[i + 4] = alpha + beta;
		A[i + 12] = alpha - beta;
		for (j = 1; j < 4; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 256];
			alpha = A[i + j];
			beta = A[i - j + 8];
			beta1 = A[i + j + 8] * y1 + A[i - j + 16] * y0;
			beta2 = A[i + j + 8] * y0 - A[i - j + 16] * y1;
			theta += 64;
			A[i + j] = alpha + beta1;
			A[i + j + 8] = alpha - beta1;
			A[i - j + 8] = beta + beta2;
			A[i - j + 16] = beta - beta2;
		}
	}
	for (i = 0; i < 1024; i += 32)
	{
		theta = 32;
		alpha = A[i];
		beta = A[i + 16];
		A[i] = alpha + beta;
		A[i + 16] = alpha - beta;
		alpha = A[i + 8];
		beta = A[i + 24];
		A[i + 8] = alpha + beta;
		A[i + 24] = alpha - beta;
		for (j = 1; j < 8; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 256];
			alpha = A[i + j];
			beta = A[i - j + 16];
			beta1 = A[i + j + 16] * y1 + A[i - j + 32] * y0;
			beta2 = A[i + j + 16] * y0 - A[i - j + 32] * y1;
			theta += 32;
			A[i + j] = alpha + beta1;
			A[i + j + 16] = alpha - beta1;
			A[i - j + 16] = beta + beta2;
			A[i - j + 32] = beta - beta2;
		}
	}
	for (i = 0; i < 1024; i += 64)
	{
		theta = 16;
		alpha = A[i];
		beta = A[i + 32];
		A[i] = alpha + beta;
		A[i + 32] = alpha - beta;
		alpha = A[i + 16];
		beta = A[i + 48];
		A[i + 16] = alpha + beta;
		A[i + 48] = alpha - beta;
		for (j = 1; j < 16; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 256];
			alpha = A[i + j];
			beta = A[i - j + 32];
			beta1 = A[i + j + 32] * y1 + A[i - j + 64] * y0;
			beta2 = A[i + j + 32] * y0 - A[i - j + 64] * y1;
			theta += 16;
			A[i + j] = alpha + beta1;
			A[i + j + 32] = alpha - beta1;
			A[i - j + 32] = beta + beta2;
			A[i - j + 64] = beta - beta2;
		}
	}
	for (i = 0; i < 1024; i += 128)
	{
		theta = 8;
		alpha = A[i];
		beta = A[i + 64];
		A[i] = alpha + beta;
		A[i + 64] = alpha - beta;
		alpha = A[i + 32];
		beta = A[i + 96];
		A[i + 32] = alpha + beta;
		A[i + 96] = alpha - beta;
		for (j = 1; j < 32; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 256];
			alpha = A[i + j];
			beta = A[i - j + 64];
			beta1 = A[i + j + 64] * y1 + A[i - j + 128] * y0;
			beta2 = A[i + j + 64] * y0 - A[i - j + 128] * y1;
			theta += 8;
			A[i + j] = alpha + beta1;
			A[i + j + 64] = alpha - beta1;
			A[i - j + 64] = beta + beta2;
			A[i - j + 128] = beta - beta2;
		}
	}
	for (i = 0; i < 1024; i += 256)
	{
		theta = 4;
		alpha = A[i];
		beta = A[i + 128];
		A[i] = alpha + beta;
		A[i + 128] = alpha - beta;
		alpha = A[i + 64];
		beta = A[i + 192];
		A[i + 64] = alpha + beta;
		A[i + 192] = alpha - beta;
		for (j = 1; j < 64; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 256];
			alpha = A[i + j];
			beta = A[i - j + 128];
			beta1 = A[i + j + 128] * y1 + A[i - j + 256] * y0;
			beta2 = A[i + j + 128] * y0 - A[i - j + 256] * y1;
			theta += 4;
			A[i + j] = alpha + beta1;
			A[i + j + 128] = alpha - beta1;
			A[i - j + 128] = beta + beta2;
			A[i - j + 256] = beta - beta2;
		}
	}
	for (i = 0; i < 1024; i += 512)
	{
		theta = 2;
		alpha = A[i];
		beta = A[i + 256];
		A[i] = alpha + beta;
		A[i + 256] = alpha - beta;
		alpha = A[i + 128];
		beta = A[i + 384];
		A[i + 128] = alpha + beta;
		A[i + 384] = alpha - beta;
		for (j = 1; j < 128; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 256];
			alpha = A[i + j];
			beta = A[i - j + 256];
			beta1 = A[i + j + 256] * y1 + A[i - j + 512] * y0;
			beta2 = A[i + j + 256] * y0 - A[i - j + 512] * y1;
			theta += 2;
			A[i + j] = alpha + beta1;
			A[i + j + 256] = alpha - beta1;
			A[i - j + 256] = beta + beta2;
			A[i - j + 512] = beta - beta2;
		}
	}
	theta = 1;
	alpha = A[0];
	beta = A[512];
	A[0] = alpha + beta;
	A[512] = alpha - beta;
	alpha = A[256];
	beta = A[768];
	A[256] = alpha + beta;
	A[768] = alpha - beta;
	for (j = 1; j < 256; j++)
	{
		y0 = sinTab[theta];
		y1 = sinTab[theta + 256];
		alpha = A[j];
		beta = A[-j + 512];
		beta1 = A[j + 512] * y1 + A[-j + 1024] * y0;
		beta2 = A[j + 512] * y0 - A[-j + 1024] * y1;
		theta += 1;
		A[j] = alpha + beta1;
		A[j + 512] = alpha - beta1;
		A[-j + 512] = beta + beta2;
		A[-j + 1024] = beta - beta2;
	}
}
void DFT2048(float *A, const float *sinTab)
{
	int i, j, theta;
	float alpha, beta, beta1, beta2, y0, y1, y2, y3;
	for (i = 0; i < 2048; i += 4)
	{
		alpha = A[i];
		beta = A[i + 1];
		beta1 = A[i + 2];
		beta2 = A[i + 3];
		y0 = alpha + beta;
		y1 = alpha - beta;
		y2 = beta1 + beta2;
		y3 = beta1 - beta2;
		A[i] = y0 + y2;
		A[i + 2] = y0 - y2;
		A[i + 1] = y1 + y3;
		A[i + 3] = y1 - y3;
	}
	for (i = 0; i < 2048; i += 8)
	{
		alpha = A[i];
		beta = A[i + 4];
		A[i] = alpha + beta;
		A[i + 4] = alpha - beta;
		alpha = A[i + 2];
		beta = A[i + 6];
		A[i + 2] = alpha + beta;
		A[i + 6] = alpha - beta;
		alpha = A[i + 1];
		beta1 = 0.70710678118654752440084436210485f*(A[i + 5] + A[i + 7]);
		beta2 = 0.70710678118654752440084436210485f*(A[i + 5] - A[i + 7]);
		A[i + 1] = alpha + beta1;
		A[i + 5] = alpha - beta1;
		alpha = A[i + 3];
		A[i + 3] = alpha + beta2;
		A[i + 7] = alpha - beta2;
	}
	for (i = 0; i < 2048; i += 16)
	{
		theta = 128;
		alpha = A[i];
		beta = A[i + 8];
		A[i] = alpha + beta;
		A[i + 8] = alpha - beta;
		alpha = A[i + 4];
		beta = A[i + 12];
		A[i + 4] = alpha + beta;
		A[i + 12] = alpha - beta;
		for (j = 1; j < 4; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 512];
			alpha = A[i + j];
			beta = A[i - j + 8];
			beta1 = A[i + j + 8] * y1 + A[i - j + 16] * y0;
			beta2 = A[i + j + 8] * y0 - A[i - j + 16] * y1;
			theta += 128;
			A[i + j] = alpha + beta1;
			A[i + j + 8] = alpha - beta1;
			A[i - j + 8] = beta + beta2;
			A[i - j + 16] = beta - beta2;
		}
	}
	for (i = 0; i < 2048; i += 32)
	{
		theta = 64;
		alpha = A[i];
		beta = A[i + 16];
		A[i] = alpha + beta;
		A[i + 16] = alpha - beta;
		alpha = A[i + 8];
		beta = A[i + 24];
		A[i + 8] = alpha + beta;
		A[i + 24] = alpha - beta;
		for (j = 1; j < 8; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 512];
			alpha = A[i + j];
			beta = A[i - j + 16];
			beta1 = A[i + j + 16] * y1 + A[i - j + 32] * y0;
			beta2 = A[i + j + 16] * y0 - A[i - j + 32] * y1;
			theta += 64;
			A[i + j] = alpha + beta1;
			A[i + j + 16] = alpha - beta1;
			A[i - j + 16] = beta + beta2;
			A[i - j + 32] = beta - beta2;
		}
	}
	for (i = 0; i < 2048; i += 64)
	{
		theta = 32;
		alpha = A[i];
		beta = A[i + 32];
		A[i] = alpha + beta;
		A[i + 32] = alpha - beta;
		alpha = A[i + 16];
		beta = A[i + 48];
		A[i + 16] = alpha + beta;
		A[i + 48] = alpha - beta;
		for (j = 1; j < 16; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 512];
			alpha = A[i + j];
			beta = A[i - j + 32];
			beta1 = A[i + j + 32] * y1 + A[i - j + 64] * y0;
			beta2 = A[i + j + 32] * y0 - A[i - j + 64] * y1;
			theta += 32;
			A[i + j] = alpha + beta1;
			A[i + j + 32] = alpha - beta1;
			A[i - j + 32] = beta + beta2;
			A[i - j + 64] = beta - beta2;
		}
	}
	for (i = 0; i < 2048; i += 128)
	{
		theta = 16;
		alpha = A[i];
		beta = A[i + 64];
		A[i] = alpha + beta;
		A[i + 64] = alpha - beta;
		alpha = A[i + 32];
		beta = A[i + 96];
		A[i + 32] = alpha + beta;
		A[i + 96] = alpha - beta;
		for (j = 1; j < 32; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 512];
			alpha = A[i + j];
			beta = A[i - j + 64];
			beta1 = A[i + j + 64] * y1 + A[i - j + 128] * y0;
			beta2 = A[i + j + 64] * y0 - A[i - j + 128] * y1;
			theta += 16;
			A[i + j] = alpha + beta1;
			A[i + j + 64] = alpha - beta1;
			A[i - j + 64] = beta + beta2;
			A[i - j + 128] = beta - beta2;
		}
	}
	for (i = 0; i < 2048; i += 256)
	{
		theta = 8;
		alpha = A[i];
		beta = A[i + 128];
		A[i] = alpha + beta;
		A[i + 128] = alpha - beta;
		alpha = A[i + 64];
		beta = A[i + 192];
		A[i + 64] = alpha + beta;
		A[i + 192] = alpha - beta;
		for (j = 1; j < 64; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 512];
			alpha = A[i + j];
			beta = A[i - j + 128];
			beta1 = A[i + j + 128] * y1 + A[i - j + 256] * y0;
			beta2 = A[i + j + 128] * y0 - A[i - j + 256] * y1;
			theta += 8;
			A[i + j] = alpha + beta1;
			A[i + j + 128] = alpha - beta1;
			A[i - j + 128] = beta + beta2;
			A[i - j + 256] = beta - beta2;
		}
	}
	for (i = 0; i < 2048; i += 512)
	{
		theta = 4;
		alpha = A[i];
		beta = A[i + 256];
		A[i] = alpha + beta;
		A[i + 256] = alpha - beta;
		alpha = A[i + 128];
		beta = A[i + 384];
		A[i + 128] = alpha + beta;
		A[i + 384] = alpha - beta;
		for (j = 1; j < 128; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 512];
			alpha = A[i + j];
			beta = A[i - j + 256];
			beta1 = A[i + j + 256] * y1 + A[i - j + 512] * y0;
			beta2 = A[i + j + 256] * y0 - A[i - j + 512] * y1;
			theta += 4;
			A[i + j] = alpha + beta1;
			A[i + j + 256] = alpha - beta1;
			A[i - j + 256] = beta + beta2;
			A[i - j + 512] = beta - beta2;
		}
	}
	for (i = 0; i < 2048; i += 1024)
	{
		theta = 2;
		alpha = A[i];
		beta = A[i + 512];
		A[i] = alpha + beta;
		A[i + 512] = alpha - beta;
		alpha = A[i + 256];
		beta = A[i + 768];
		A[i + 256] = alpha + beta;
		A[i + 768] = alpha - beta;
		for (j = 1; j < 256; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 512];
			alpha = A[i + j];
			beta = A[i - j + 512];
			beta1 = A[i + j + 512] * y1 + A[i - j + 1024] * y0;
			beta2 = A[i + j + 512] * y0 - A[i - j + 1024] * y1;
			theta += 2;
			A[i + j] = alpha + beta1;
			A[i + j + 512] = alpha - beta1;
			A[i - j + 512] = beta + beta2;
			A[i - j + 1024] = beta - beta2;
		}
	}
	theta = 1;
	alpha = A[0];
	beta = A[1024];
	A[0] = alpha + beta;
	A[1024] = alpha - beta;
	alpha = A[512];
	beta = A[1536];
	A[512] = alpha + beta;
	A[1536] = alpha - beta;
	for (j = 1; j < 512; j++)
	{
		y0 = sinTab[theta];
		y1 = sinTab[theta + 512];
		alpha = A[j];
		beta = A[-j + 1024];
		beta1 = A[j + 1024] * y1 + A[-j + 2048] * y0;
		beta2 = A[j + 1024] * y0 - A[-j + 2048] * y1;
		theta += 1;
		A[j] = alpha + beta1;
		A[j + 1024] = alpha - beta1;
		A[-j + 1024] = beta + beta2;
		A[-j + 2048] = beta - beta2;
	}
}
void DFT4096(float *A, const float *sinTab)
{
	int i, j, theta;
	float alpha, beta, beta1, beta2, y0, y1, y2, y3;
	for (i = 0; i < 4096; i += 4)
	{
		alpha = A[i];
		beta = A[i + 1];
		beta1 = A[i + 2];
		beta2 = A[i + 3];
		y0 = alpha + beta;
		y1 = alpha - beta;
		y2 = beta1 + beta2;
		y3 = beta1 - beta2;
		A[i] = y0 + y2;
		A[i + 2] = y0 - y2;
		A[i + 1] = y1 + y3;
		A[i + 3] = y1 - y3;
	}
	for (i = 0; i < 4096; i += 8)
	{
		alpha = A[i];
		beta = A[i + 4];
		A[i] = alpha + beta;
		A[i + 4] = alpha - beta;
		alpha = A[i + 2];
		beta = A[i + 6];
		A[i + 2] = alpha + beta;
		A[i + 6] = alpha - beta;
		alpha = A[i + 1];
		beta1 = 0.70710678118654752440084436210485f*(A[i + 5] + A[i + 7]);
		beta2 = 0.70710678118654752440084436210485f*(A[i + 5] - A[i + 7]);
		A[i + 1] = alpha + beta1;
		A[i + 5] = alpha - beta1;
		alpha = A[i + 3];
		A[i + 3] = alpha + beta2;
		A[i + 7] = alpha - beta2;
	}
	for (i = 0; i < 4096; i += 16)
	{
		theta = 256;
		alpha = A[i];
		beta = A[i + 8];
		A[i] = alpha + beta;
		A[i + 8] = alpha - beta;
		alpha = A[i + 4];
		beta = A[i + 12];
		A[i + 4] = alpha + beta;
		A[i + 12] = alpha - beta;
		for (j = 1; j < 4; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 1024];
			alpha = A[i + j];
			beta = A[i - j + 8];
			beta1 = A[i + j + 8] * y1 + A[i - j + 16] * y0;
			beta2 = A[i + j + 8] * y0 - A[i - j + 16] * y1;
			theta += 256;
			A[i + j] = alpha + beta1;
			A[i + j + 8] = alpha - beta1;
			A[i - j + 8] = beta + beta2;
			A[i - j + 16] = beta - beta2;
		}
	}
	for (i = 0; i < 4096; i += 32)
	{
		theta = 128;
		alpha = A[i];
		beta = A[i + 16];
		A[i] = alpha + beta;
		A[i + 16] = alpha - beta;
		alpha = A[i + 8];
		beta = A[i + 24];
		A[i + 8] = alpha + beta;
		A[i + 24] = alpha - beta;
		for (j = 1; j < 8; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 1024];
			alpha = A[i + j];
			beta = A[i - j + 16];
			beta1 = A[i + j + 16] * y1 + A[i - j + 32] * y0;
			beta2 = A[i + j + 16] * y0 - A[i - j + 32] * y1;
			theta += 128;
			A[i + j] = alpha + beta1;
			A[i + j + 16] = alpha - beta1;
			A[i - j + 16] = beta + beta2;
			A[i - j + 32] = beta - beta2;
		}
	}
	for (i = 0; i < 4096; i += 64)
	{
		theta = 64;
		alpha = A[i];
		beta = A[i + 32];
		A[i] = alpha + beta;
		A[i + 32] = alpha - beta;
		alpha = A[i + 16];
		beta = A[i + 48];
		A[i + 16] = alpha + beta;
		A[i + 48] = alpha - beta;
		for (j = 1; j < 16; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 1024];
			alpha = A[i + j];
			beta = A[i - j + 32];
			beta1 = A[i + j + 32] * y1 + A[i - j + 64] * y0;
			beta2 = A[i + j + 32] * y0 - A[i - j + 64] * y1;
			theta += 64;
			A[i + j] = alpha + beta1;
			A[i + j + 32] = alpha - beta1;
			A[i - j + 32] = beta + beta2;
			A[i - j + 64] = beta - beta2;
		}
	}
	for (i = 0; i < 4096; i += 128)
	{
		theta = 32;
		alpha = A[i];
		beta = A[i + 64];
		A[i] = alpha + beta;
		A[i + 64] = alpha - beta;
		alpha = A[i + 32];
		beta = A[i + 96];
		A[i + 32] = alpha + beta;
		A[i + 96] = alpha - beta;
		for (j = 1; j < 32; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 1024];
			alpha = A[i + j];
			beta = A[i - j + 64];
			beta1 = A[i + j + 64] * y1 + A[i - j + 128] * y0;
			beta2 = A[i + j + 64] * y0 - A[i - j + 128] * y1;
			theta += 32;
			A[i + j] = alpha + beta1;
			A[i + j + 64] = alpha - beta1;
			A[i - j + 64] = beta + beta2;
			A[i - j + 128] = beta - beta2;
		}
	}
	for (i = 0; i < 4096; i += 256)
	{
		theta = 16;
		alpha = A[i];
		beta = A[i + 128];
		A[i] = alpha + beta;
		A[i + 128] = alpha - beta;
		alpha = A[i + 64];
		beta = A[i + 192];
		A[i + 64] = alpha + beta;
		A[i + 192] = alpha - beta;
		for (j = 1; j < 64; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 1024];
			alpha = A[i + j];
			beta = A[i - j + 128];
			beta1 = A[i + j + 128] * y1 + A[i - j + 256] * y0;
			beta2 = A[i + j + 128] * y0 - A[i - j + 256] * y1;
			theta += 16;
			A[i + j] = alpha + beta1;
			A[i + j + 128] = alpha - beta1;
			A[i - j + 128] = beta + beta2;
			A[i - j + 256] = beta - beta2;
		}
	}
	for (i = 0; i < 4096; i += 512)
	{
		theta = 8;
		alpha = A[i];
		beta = A[i + 256];
		A[i] = alpha + beta;
		A[i + 256] = alpha - beta;
		alpha = A[i + 128];
		beta = A[i + 384];
		A[i + 128] = alpha + beta;
		A[i + 384] = alpha - beta;
		for (j = 1; j < 128; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 1024];
			alpha = A[i + j];
			beta = A[i - j + 256];
			beta1 = A[i + j + 256] * y1 + A[i - j + 512] * y0;
			beta2 = A[i + j + 256] * y0 - A[i - j + 512] * y1;
			theta += 8;
			A[i + j] = alpha + beta1;
			A[i + j + 256] = alpha - beta1;
			A[i - j + 256] = beta + beta2;
			A[i - j + 512] = beta - beta2;
		}
	}
	for (i = 0; i < 4096; i += 1024)
	{
		theta = 4;
		alpha = A[i];
		beta = A[i + 512];
		A[i] = alpha + beta;
		A[i + 512] = alpha - beta;
		alpha = A[i + 256];
		beta = A[i + 768];
		A[i + 256] = alpha + beta;
		A[i + 768] = alpha - beta;
		for (j = 1; j < 256; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 1024];
			alpha = A[i + j];
			beta = A[i - j + 512];
			beta1 = A[i + j + 512] * y1 + A[i - j + 1024] * y0;
			beta2 = A[i + j + 512] * y0 - A[i - j + 1024] * y1;
			theta += 4;
			A[i + j] = alpha + beta1;
			A[i + j + 512] = alpha - beta1;
			A[i - j + 512] = beta + beta2;
			A[i - j + 1024] = beta - beta2;
		}
	}
	for (i = 0; i < 4096; i += 2048)
	{
		theta = 2;
		alpha = A[i];
		beta = A[i + 1024];
		A[i] = alpha + beta;
		A[i + 1024] = alpha - beta;
		alpha = A[i + 512];
		beta = A[i + 1536];
		A[i + 512] = alpha + beta;
		A[i + 1536] = alpha - beta;
		for (j = 1; j < 512; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 1024];
			alpha = A[i + j];
			beta = A[i - j + 1024];
			beta1 = A[i + j + 1024] * y1 + A[i - j + 2048] * y0;
			beta2 = A[i + j + 1024] * y0 - A[i - j + 2048] * y1;
			theta += 2;
			A[i + j] = alpha + beta1;
			A[i + j + 1024] = alpha - beta1;
			A[i - j + 1024] = beta + beta2;
			A[i - j + 2048] = beta - beta2;
		}
	}
	theta = 1;
	alpha = A[0];
	beta = A[2048];
	A[0] = alpha + beta;
	A[2048] = alpha - beta;
	alpha = A[1024];
	beta = A[3072];
	A[1024] = alpha + beta;
	A[3072] = alpha - beta;
	for (j = 1; j < 1024; j++)
	{
		y0 = sinTab[theta];
		y1 = sinTab[theta + 1024];
		alpha = A[j];
		beta = A[-j + 2048];
		beta1 = A[j + 2048] * y1 + A[-j + 4096] * y0;
		beta2 = A[j + 2048] * y0 - A[-j + 4096] * y1;
		theta += 1;
		A[j] = alpha + beta1;
		A[j + 2048] = alpha - beta1;
		A[-j + 2048] = beta + beta2;
		A[-j + 4096] = beta - beta2;
	}
}
void DFT8192(float *A, const float *sinTab)
{
	int i, j, theta;
	float alpha, beta, beta1, beta2, y0, y1, y2, y3;
	for (i = 0; i < 8192; i += 4)
	{
		alpha = A[i];
		beta = A[i + 1];
		beta1 = A[i + 2];
		beta2 = A[i + 3];
		y0 = alpha + beta;
		y1 = alpha - beta;
		y2 = beta1 + beta2;
		y3 = beta1 - beta2;
		A[i] = y0 + y2;
		A[i + 2] = y0 - y2;
		A[i + 1] = y1 + y3;
		A[i + 3] = y1 - y3;
	}
	for (i = 0; i < 8192; i += 8)
	{
		alpha = A[i];
		beta = A[i + 4];
		A[i] = alpha + beta;
		A[i + 4] = alpha - beta;
		alpha = A[i + 2];
		beta = A[i + 6];
		A[i + 2] = alpha + beta;
		A[i + 6] = alpha - beta;
		alpha = A[i + 1];
		beta1 = 0.70710678118654752440084436210485f*(A[i + 5] + A[i + 7]);
		beta2 = 0.70710678118654752440084436210485f*(A[i + 5] - A[i + 7]);
		A[i + 1] = alpha + beta1;
		A[i + 5] = alpha - beta1;
		alpha = A[i + 3];
		A[i + 3] = alpha + beta2;
		A[i + 7] = alpha - beta2;
	}
	for (i = 0; i < 8192; i += 16)
	{
		theta = 512;
		alpha = A[i];
		beta = A[i + 8];
		A[i] = alpha + beta;
		A[i + 8] = alpha - beta;
		alpha = A[i + 4];
		beta = A[i + 12];
		A[i + 4] = alpha + beta;
		A[i + 12] = alpha - beta;
		for (j = 1; j < 4; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 2048];
			alpha = A[i + j];
			beta = A[i - j + 8];
			beta1 = A[i + j + 8] * y1 + A[i - j + 16] * y0;
			beta2 = A[i + j + 8] * y0 - A[i - j + 16] * y1;
			theta += 512;
			A[i + j] = alpha + beta1;
			A[i + j + 8] = alpha - beta1;
			A[i - j + 8] = beta + beta2;
			A[i - j + 16] = beta - beta2;
		}
	}
	for (i = 0; i < 8192; i += 32)
	{
		theta = 256;
		alpha = A[i];
		beta = A[i + 16];
		A[i] = alpha + beta;
		A[i + 16] = alpha - beta;
		alpha = A[i + 8];
		beta = A[i + 24];
		A[i + 8] = alpha + beta;
		A[i + 24] = alpha - beta;
		for (j = 1; j < 8; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 2048];
			alpha = A[i + j];
			beta = A[i - j + 16];
			beta1 = A[i + j + 16] * y1 + A[i - j + 32] * y0;
			beta2 = A[i + j + 16] * y0 - A[i - j + 32] * y1;
			theta += 256;
			A[i + j] = alpha + beta1;
			A[i + j + 16] = alpha - beta1;
			A[i - j + 16] = beta + beta2;
			A[i - j + 32] = beta - beta2;
		}
	}
	for (i = 0; i < 8192; i += 64)
	{
		theta = 128;
		alpha = A[i];
		beta = A[i + 32];
		A[i] = alpha + beta;
		A[i + 32] = alpha - beta;
		alpha = A[i + 16];
		beta = A[i + 48];
		A[i + 16] = alpha + beta;
		A[i + 48] = alpha - beta;
		for (j = 1; j < 16; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 2048];
			alpha = A[i + j];
			beta = A[i - j + 32];
			beta1 = A[i + j + 32] * y1 + A[i - j + 64] * y0;
			beta2 = A[i + j + 32] * y0 - A[i - j + 64] * y1;
			theta += 128;
			A[i + j] = alpha + beta1;
			A[i + j + 32] = alpha - beta1;
			A[i - j + 32] = beta + beta2;
			A[i - j + 64] = beta - beta2;
		}
	}
	for (i = 0; i < 8192; i += 128)
	{
		theta = 64;
		alpha = A[i];
		beta = A[i + 64];
		A[i] = alpha + beta;
		A[i + 64] = alpha - beta;
		alpha = A[i + 32];
		beta = A[i + 96];
		A[i + 32] = alpha + beta;
		A[i + 96] = alpha - beta;
		for (j = 1; j < 32; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 2048];
			alpha = A[i + j];
			beta = A[i - j + 64];
			beta1 = A[i + j + 64] * y1 + A[i - j + 128] * y0;
			beta2 = A[i + j + 64] * y0 - A[i - j + 128] * y1;
			theta += 64;
			A[i + j] = alpha + beta1;
			A[i + j + 64] = alpha - beta1;
			A[i - j + 64] = beta + beta2;
			A[i - j + 128] = beta - beta2;
		}
	}
	for (i = 0; i < 8192; i += 256)
	{
		theta = 32;
		alpha = A[i];
		beta = A[i + 128];
		A[i] = alpha + beta;
		A[i + 128] = alpha - beta;
		alpha = A[i + 64];
		beta = A[i + 192];
		A[i + 64] = alpha + beta;
		A[i + 192] = alpha - beta;
		for (j = 1; j < 64; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 2048];
			alpha = A[i + j];
			beta = A[i - j + 128];
			beta1 = A[i + j + 128] * y1 + A[i - j + 256] * y0;
			beta2 = A[i + j + 128] * y0 - A[i - j + 256] * y1;
			theta += 32;
			A[i + j] = alpha + beta1;
			A[i + j + 128] = alpha - beta1;
			A[i - j + 128] = beta + beta2;
			A[i - j + 256] = beta - beta2;
		}
	}
	for (i = 0; i < 8192; i += 512)
	{
		theta = 16;
		alpha = A[i];
		beta = A[i + 256];
		A[i] = alpha + beta;
		A[i + 256] = alpha - beta;
		alpha = A[i + 128];
		beta = A[i + 384];
		A[i + 128] = alpha + beta;
		A[i + 384] = alpha - beta;
		for (j = 1; j < 128; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 2048];
			alpha = A[i + j];
			beta = A[i - j + 256];
			beta1 = A[i + j + 256] * y1 + A[i - j + 512] * y0;
			beta2 = A[i + j + 256] * y0 - A[i - j + 512] * y1;
			theta += 16;
			A[i + j] = alpha + beta1;
			A[i + j + 256] = alpha - beta1;
			A[i - j + 256] = beta + beta2;
			A[i - j + 512] = beta - beta2;
		}
	}
	for (i = 0; i < 8192; i += 1024)
	{
		theta = 8;
		alpha = A[i];
		beta = A[i + 512];
		A[i] = alpha + beta;
		A[i + 512] = alpha - beta;
		alpha = A[i + 256];
		beta = A[i + 768];
		A[i + 256] = alpha + beta;
		A[i + 768] = alpha - beta;
		for (j = 1; j < 256; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 2048];
			alpha = A[i + j];
			beta = A[i - j + 512];
			beta1 = A[i + j + 512] * y1 + A[i - j + 1024] * y0;
			beta2 = A[i + j + 512] * y0 - A[i - j + 1024] * y1;
			theta += 8;
			A[i + j] = alpha + beta1;
			A[i + j + 512] = alpha - beta1;
			A[i - j + 512] = beta + beta2;
			A[i - j + 1024] = beta - beta2;
		}
	}
	for (i = 0; i < 8192; i += 2048)
	{
		theta = 4;
		alpha = A[i];
		beta = A[i + 1024];
		A[i] = alpha + beta;
		A[i + 1024] = alpha - beta;
		alpha = A[i + 512];
		beta = A[i + 1536];
		A[i + 512] = alpha + beta;
		A[i + 1536] = alpha - beta;
		for (j = 1; j < 512; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 2048];
			alpha = A[i + j];
			beta = A[i - j + 1024];
			beta1 = A[i + j + 1024] * y1 + A[i - j + 2048] * y0;
			beta2 = A[i + j + 1024] * y0 - A[i - j + 2048] * y1;
			theta += 4;
			A[i + j] = alpha + beta1;
			A[i + j + 1024] = alpha - beta1;
			A[i - j + 1024] = beta + beta2;
			A[i - j + 2048] = beta - beta2;
		}
	}
	for (i = 0; i < 8192; i += 4096)
	{
		theta = 2;
		alpha = A[i];
		beta = A[i + 2048];
		A[i] = alpha + beta;
		A[i + 2048] = alpha - beta;
		alpha = A[i + 1024];
		beta = A[i + 3072];
		A[i + 1024] = alpha + beta;
		A[i + 3072] = alpha - beta;
		for (j = 1; j < 1024; j++)
		{
			y0 = sinTab[theta];
			y1 = sinTab[theta + 2048];
			alpha = A[i + j];
			beta = A[i - j + 2048];
			beta1 = A[i + j + 2048] * y1 + A[i - j + 4096] * y0;
			beta2 = A[i + j + 2048] * y0 - A[i - j + 4096] * y1;
			theta += 2;
			A[i + j] = alpha + beta1;
			A[i + j + 2048] = alpha - beta1;
			A[i - j + 2048] = beta + beta2;
			A[i - j + 4096] = beta - beta2;
		}
	}
	theta = 1;
	alpha = A[0];
	beta = A[4096];
	A[0] = alpha + beta;
	A[4096] = alpha - beta;
	alpha = A[2048];
	beta = A[6144];
	A[2048] = alpha + beta;
	A[6144] = alpha - beta;
	for (j = 1; j < 2048; j++)
	{
		y0 = sinTab[theta];
		y1 = sinTab[theta + 2048];
		alpha = A[j];
		beta = A[-j + 4096];
		beta1 = A[j + 4096] * y1 + A[-j + 8192] * y0;
		beta2 = A[j + 4096] * y0 - A[-j + 8192] * y1;
		theta += 1;
		A[j] = alpha + beta1;
		A[j + 4096] = alpha - beta1;
		A[-j + 4096] = beta + beta2;
		A[-j + 8192] = beta - beta2;
	}
}