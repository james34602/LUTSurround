#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "LUTSurround.h"
#include "codelet.h"
#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_2PI
#define M_2PI (M_PI * 2.0)
#endif
static inline float clampFLT(float x)
{
	return max(-1.0f, min(1.0f, x));
}
static inline int map_to_grid(float *x, int gridSize)
{
	float gp = ((*x + 1.0f)*0.5f)*((float)gridSize - 1.0f);
	float i = min((float)gridSize - 2.0f, floorf(gp));
	*x = gp - i;
	return (int)i;
}
static inline int map_to_grid2(float *x)
{
	float gp = (*x + 1.0f) * 127.5f;
	float i = floorf(gp);
	if (i > 254.0f)
		i = 254.0f;
	*x = gp - i;
	return (int)i;
}
static inline  float map(const float x, const float in_min, const float in_max, const float out_min, const float out_max)
{
	return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}
#include "tbl.h"
static inline void sincos_fast(float x, float g, float *pS, float *pC)
{
	const float cosOff4LUT[] = { 0x1.000000p+00,  0x1.6A09E6p-01,  0x0.000000p+00, -0x1.6A09E6p-01, -0x1.000000p+00, -0x1.6A09E6p-01,  0x0.000000p+00,  0x1.6A09E6p-01 };
	int     m, ms, mc;
	float   xI, xR, xR2;
	float   c, s, cy, sy;
	// Cody & Waite's range reduction Algorithm, [-pi/4, pi/4]
	xI = floorf(x * 0x1.45F306p+00 + 0.5);              // This is 4/pi.
	xR = (x - xI * 0x1.920000p-01) - xI * 0x1.FB5444p-13; // This is pi/4 in two parts per C&W.
	m = (int)xI;
	xR2 = xR * xR;
	// Find cosine & sine index for angle offsets indices
	mc = (m) & 0x7;     // two's complement permits upper modulus for negative numbers =P
	ms = (m + 6) & 0x7;     // phase correction for sine
	// Find cosine & sine
	cy = cosOff4LUT[mc];     // Load angle offset neighborhood cosine value 
	sy = cosOff4LUT[ms];     // Load angle offset neighborhood sine value
	//c = 0xf.ff79fp-4 + xR2 * (-0x7.e58e9p-4);               // TOL = 1.2786e-4
	c = 0xf.ffffdp-4 + xR2 * (-0x7.ffebep-4 + xR2 * 0xa.956a9p-8);  // TOL = 1.7882e-7
	//s = xR * (0xf.ffbf7p-4 + xR2 * (-0x2.a41d0cp-4));   // TOL = 4.835251e-6
	s = xR * (0xf.fffffp-4 + xR2 * (-0x2.aaa65cp-4 + xR2 * 0x2.1ea25p-8));  // TOL = 1.1841e-8
	*pC = clampFLT(g * (c * cy - s * sy));
	*pS = clampFLT(g * (c * sy + s * cy));
}
static inline void cartesianPrecomputer(float arg1, float arg2, float refangle, float separation, float *x, float *y)
{
	int u = map_to_grid2(&arg1);
	int v = map_to_grid2(&arg2);
	float panningPlaneMagnitude = ((1.0f - arg1) * (1.0f - arg2)) * angleMapper[0 * 256 * 256 + u * 256 + v] + (arg1 * (1.0f - arg2)) * angleMapper[0 * 256 * 256 + (u + 1) * 256 + v] + ((1.0f - arg1) * arg2) * angleMapper[0 * 256 * 256 + u * 256 + v + 1] + (arg1 * arg2) * angleMapper[0 * 256 * 256 + (u + 1) * 256 + v + 1];
	float panningPlanePhase = ((1.0f - arg1) * (1.0f - arg2)) * angleMapper[1 * 256 * 256 + u * 256 + v] + (arg1 * (1.0f - arg2)) * angleMapper[1 * 256 * 256 + (u + 1) * 256 + v] + ((1.0f - arg1) * arg2) * angleMapper[1 * 256 * 256 + u * 256 + v + 1] + (arg1 * arg2) * angleMapper[1 * 256 * 256 + (u + 1) * 256 + v + 1];
	float xAxis = map(panningPlanePhase, -M_PI, M_PI, -1.0f, 1.0f);
	arg2 = map(refangle, 0.0f, 2.0f * M_PI, -1.0f, 1.0f);
	u = map_to_grid2(&xAxis);
	v = map_to_grid2(&arg2);
	panningPlanePhase = ((1.0f - xAxis) * (1.0f - arg2)) * angleMapper[2 * 256 * 256 + u * 256 + v] + (xAxis * (1.0f - arg2)) * angleMapper[2 * 256 * 256 + (u + 1) * 256 + v] + ((1.0f - xAxis) * arg2) * angleMapper[2 * 256 * 256 + u * 256 + v + 1] + (xAxis * arg2) * angleMapper[2 * 256 * 256 + (u + 1) * 256 + v + 1];
	xAxis = map(panningPlaneMagnitude, 0.0, sqrt(2.0), -1.0f, 1.0f);
	arg2 = map(separation, -0.3f, 0.3f, -1.0f, 1.0f);
	u = map_to_grid2(&xAxis);
	v = map_to_grid2(&arg2);
	panningPlaneMagnitude = ((1.0f - xAxis) * (1.0f - arg2)) * angleMapper[3 * 256 * 256 + u * 256 + v] + (xAxis * (1.0f - arg2)) * angleMapper[3 * 256 * 256 + (u + 1) * 256 + v] + ((1.0f - xAxis) * arg2) * angleMapper[3 * 256 * 256 + u * 256 + v + 1] + (xAxis * arg2) * angleMapper[3 * 256 * 256 + (u + 1) * 256 + v + 1];
	sincos_fast(panningPlanePhase, panningPlaneMagnitude, x, y);
}
static inline void cartesianMap(float *magPhaseDiff2Cartesian, float xAxis, float arg, float *x, float *y)
{
	float yAxis = arg * (2.0f / M_PI) - 1.0f;
	int q = map_to_grid(&xAxis, TBL_SIZE);
	int p = map_to_grid(&yAxis, TBL_SIZE);
	*x = ((1.0f - xAxis) * (1.0f - yAxis)) * magPhaseDiff2Cartesian[0 * TBL_SIZE * TBL_SIZE + q * TBL_SIZE + p] + (xAxis * (1.0f - yAxis)) * magPhaseDiff2Cartesian[0 * TBL_SIZE * TBL_SIZE + (q + 1) * TBL_SIZE + p] + ((1.0f - xAxis) * yAxis) * magPhaseDiff2Cartesian[0 * TBL_SIZE * TBL_SIZE + q * TBL_SIZE + p + 1] + (xAxis * yAxis) * magPhaseDiff2Cartesian[0 * TBL_SIZE * TBL_SIZE + (q + 1) * TBL_SIZE + p + 1];
	*y = ((1.0f - xAxis) * (1.0f - yAxis)) * magPhaseDiff2Cartesian[1 * TBL_SIZE * TBL_SIZE + q * TBL_SIZE + p] + (xAxis * (1.0f - yAxis)) * magPhaseDiff2Cartesian[1 * TBL_SIZE * TBL_SIZE + (q + 1) * TBL_SIZE + p] + ((1.0f - xAxis) * yAxis) * magPhaseDiff2Cartesian[1 * TBL_SIZE * TBL_SIZE + q * TBL_SIZE + p + 1] + (xAxis * yAxis) * magPhaseDiff2Cartesian[1 * TBL_SIZE * TBL_SIZE + (q + 1) * TBL_SIZE + p + 1];
}
#define mainLoop \
symIdx = msr->fftLen - i; \
bitRevFwd = msr->mBitRev[i]; \
bitRevSym = msr->mBitRev[symIdx]; \
lR = msr->mTempLBuffer[i] + msr->mTempLBuffer[symIdx]; \
float lI = msr->mTempLBuffer[i] - msr->mTempLBuffer[symIdx]; \
rR = msr->mTempRBuffer[i] + msr->mTempRBuffer[symIdx]; \
float rI = msr->mTempRBuffer[i] - msr->mTempRBuffer[symIdx]; \
magnitudeLeft = hypotf(lR, lI); \
magnitudeRight = hypotf(rR, rI); \
float minMag = min(magnitudeLeft, magnitudeRight); \
float ratMinLeft = minMag / (magnitudeLeft + FLT_EPSILON); \
if (fabsf(ratMinLeft - msr->smoothedFunc[8][i]) < msr->smoothing) \
	msr->smoothedFunc[8][i] = ratMinLeft; \
else if (ratMinLeft > msr->smoothedFunc[8][i]) \
msr->smoothedFunc[8][i] += msr->smoothing; \
else \
msr->smoothedFunc[8][i] -= msr->smoothing; \
float ratMinRight = minMag / (magnitudeRight + FLT_EPSILON); \
if (fabsf(ratMinRight - msr->smoothedFunc[7][i]) < msr->smoothing) \
	msr->smoothedFunc[7][i] = ratMinRight; \
else if (ratMinRight > msr->smoothedFunc[7][i]) \
msr->smoothedFunc[7][i] += msr->smoothing; \
else \
msr->smoothedFunc[7][i] -= msr->smoothing; \
magnitudeSum = magnitudeLeft + magnitudeRight; \
magDiff = clampFLT((magnitudeSum < DBL_EPSILON) ? 0.0f : (magnitudeRight - magnitudeLeft) / magnitudeSum); \
phaseL = atan2f(lI, lR); \
phaseR = atan2f(rI, rR); \
phaseDiff = fabsf(phaseR - phaseL); \
if (phaseDiff > (float)M_PI) \
phaseDiff = (float)M_2PI - phaseDiff; \
cartesianMap(msr->magPhaseDiff2Cartesian, magDiff, phaseDiff, &x, &y); \
leftCosineTerm = cosf(phaseL); \
float leftSineTerm = sinf(phaseL); \
rightCosineTerm = cosf(phaseR); \
float rightSineTerm = sinf(phaseR); \
float centrePhase = atan2f(lI + rI, lR + rR); \
centreCosineTerm = cosf(centrePhase); \
float centreSineTerm = sinf(centrePhase); \
u = map_to_grid(&x, 21); \
v = map_to_grid(&y, 21); \
magSqrt = hypotf(magnitudeLeft, magnitudeRight); \
opLFwd = (lR + lI) - (rR + rI) * msr->smoothedFunc[7][i]; \
float opLSym = (lR - lI) - (rR - rI) * msr->smoothedFunc[7][i]; \
opRFwd = (rR + rI) - (lR + lI) * msr->smoothedFunc[8][i]; \
float opRSym = (rR - rI) - (lR - lI) * msr->smoothedFunc[8][i];

#define lerpCompCplx(idx) \
a = &chmaskPtr[idx * 21 * 21]; \
gf = ((1.0f - x)*(1.0f - y)*a[u * 21 + v] + x * (1.0f - y)*a[(u + 1) * 21 + v] + (1.0f - x)*y*a[u * 21 + v + 1] + x * y*a[(u + 1) * 21 + v + 1]); \
difGf = gf - msr->smoothedFunc[idx][i]; \
if (difGf > msr->smoothing) \
gf = msr->smoothedFunc[idx][i] + msr->smoothing; \
if (difGf < -msr->smoothing) \
gf = msr->smoothedFunc[idx][i] - msr->smoothing; \
msr->smoothedFunc[idx][i] = gf; \
magnitude = magSqrt * msr->smoothedFunc[idx][i];
static inline void LLPAMSProcessNPR(LUTSurround *msr, const float * const chmaskPtr, const unsigned int numOut)
{
	unsigned int i, j;
	for (i = 0; i < msr->frameLen; ++i)
	{
		unsigned int k = i + msr->mInputPos;
		if (k >= msr->frameLen)
			k = k - msr->frameLen;
		const float w = msr->analysisWnd[i];
		msr->mTempLBuffer[msr->mBitRev[i]] = msr->mInput[0][k] * w;
		msr->mTempRBuffer[msr->mBitRev[i]] = msr->mInput[1][k] * w;
	}
	for (i = msr->frameLen; i < msr->fftLen; ++i)
	{
		msr->mTempLBuffer[msr->mBitRev[i]] = 0;
		msr->mTempRBuffer[msr->mBitRev[i]] = 0;
	}
	msr->fft(msr->mTempLBuffer, msr->mSineTab);
	msr->fft(msr->mTempRBuffer, msr->mSineTab);
	const float *a;
	float gf, difGf, magnitude, real, imag, leftCosineTerm, rightCosineTerm, centreCosineTerm;
	unsigned int symIdx, bitRevFwd, bitRevSym;
	float lR = msr->mTempLBuffer[0] * 2.0f;
	float rR = msr->mTempRBuffer[0] * 2.0f;
	float magnitudeLeft = fabsf(lR);
	float magnitudeRight = fabsf(rR);
	float magnitudeSum = magnitudeLeft + magnitudeRight;
	float magDiff = clampFLT((magnitudeSum < FLT_EPSILON) ? 0.0f : ((magnitudeRight - magnitudeLeft) / magnitudeSum));
	float phaseL = lR < 0.0f ? M_PI : 0.0f;
	float phaseR = rR < 0.0f ? M_PI : 0.0f;
	float phaseDiff = fabsf(phaseR - phaseL);
	if (phaseDiff > M_PI)
		phaseDiff = (float)M_2PI - phaseDiff;
	float x, y;
	cartesianMap(msr->magPhaseDiff2Cartesian, magDiff, phaseDiff, &x, &y);
	leftCosineTerm = cosf(phaseL);
	rightCosineTerm = cosf(phaseR);
	centreCosineTerm = (lR + rR) < 0.0f ? -1.0f : 1.0f;
	int u = map_to_grid(&x, 21), v = map_to_grid(&y, 21);
	float magSqrt = hypotf(magnitudeLeft, magnitudeRight);
	float opLFwd = lR - rR * msr->smoothedFunc[7][0];
	float opRFwd = rR - lR * msr->smoothedFunc[8][0];
	if (numOut == 4) // LCR
	{
		lerpCompCplx(0);
		msr->timeDomainOut[0][0] = (magnitude * leftCosineTerm);
		lerpCompCplx(1);
		msr->timeDomainOut[1][0] = (magnitude * rightCosineTerm);
		lerpCompCplx(2);
		msr->timeDomainOut[2][0] = (magnitude * centreCosineTerm);
		for (i = 1; i < msr->procUpTo; i++)
		{
			mainLoop;
			lerpCompCplx(0);
			real = magnitude * leftCosineTerm;
			imag = magnitude * leftSineTerm;
			msr->timeDomainOut[0][bitRevFwd] = (real + imag);
			msr->timeDomainOut[0][bitRevSym] = (real - imag);
			lerpCompCplx(1);
			real = magnitude * rightCosineTerm;
			imag = magnitude * rightSineTerm;
			msr->timeDomainOut[1][bitRevFwd] = (real + imag);
			msr->timeDomainOut[1][bitRevSym] = (real - imag);
			lerpCompCplx(2);
			real = magnitude * centreCosineTerm;
			imag = magnitude * centreSineTerm;
			msr->timeDomainOut[2][bitRevFwd] = (real + imag);
			msr->timeDomainOut[2][bitRevSym] = (real - imag);
		}
	}
	else if (numOut == 5) // Quadraphonic
	{
		float pmix = msr->mix * 0.5f;
		float minusPmix = 1.0f - pmix;
		lerpCompCplx(0);
		msr->timeDomainOut[0][0] = (magnitude * leftCosineTerm);
		lerpCompCplx(1);
		msr->timeDomainOut[1][0] = (magnitude * rightCosineTerm);
		lerpCompCplx(2);
		msr->timeDomainOut[2][0] = (magnitude * leftCosineTerm) * minusPmix + opLFwd * pmix;
		lerpCompCplx(3);
		msr->timeDomainOut[3][0] = (magnitude * rightCosineTerm) * minusPmix + opRFwd * pmix;
		for (i = 1; i < msr->procUpTo; i++)
		{
			mainLoop;
			lerpCompCplx(0);
			real = magnitude * leftCosineTerm;
			imag = magnitude * leftSineTerm;
			msr->timeDomainOut[0][bitRevFwd] = (real + imag);
			msr->timeDomainOut[0][bitRevSym] = (real - imag);
			lerpCompCplx(1);
			real = magnitude * rightCosineTerm;
			imag = magnitude * rightSineTerm;
			msr->timeDomainOut[1][bitRevFwd] = (real + imag);
			msr->timeDomainOut[1][bitRevSym] = (real - imag);
			lerpCompCplx(2);
			real = magnitude * leftCosineTerm;
			imag = magnitude * leftSineTerm;
			msr->timeDomainOut[2][bitRevFwd] = (real + imag) * minusPmix + opLFwd * pmix;
			msr->timeDomainOut[2][bitRevSym] = (real - imag) * minusPmix + opLFwd * pmix;
			lerpCompCplx(3);
			real = magnitude * rightCosineTerm;
			imag = magnitude * rightSineTerm;
			msr->timeDomainOut[3][bitRevFwd] = (real + imag) * minusPmix + opRFwd * pmix;
			msr->timeDomainOut[3][bitRevSym] = (real - imag) * minusPmix + opRFwd * pmix;
		}
	}
	else if (numOut == 6) // 5.1
	{
		float pmix = msr->mix;
		float minusPmix = 1.0f - pmix;
		lerpCompCplx(0);
		msr->timeDomainOut[0][0] = (magnitude * leftCosineTerm) * 1.2f;
		lerpCompCplx(1);
		msr->timeDomainOut[1][0] = (magnitude * rightCosineTerm) * 1.2f;
		lerpCompCplx(2);
		msr->timeDomainOut[2][0] = (magnitude * centreCosineTerm) * 1.2f;
		lerpCompCplx(3);
		msr->timeDomainOut[3][0] = (magnitude * leftCosineTerm) * minusPmix + opLFwd * pmix;
		lerpCompCplx(4);
		msr->timeDomainOut[4][0] = (magnitude * rightCosineTerm) * minusPmix + opRFwd * pmix;
		for (i = 1; i < msr->procUpTo; i++)
		{
			mainLoop;
			lerpCompCplx(0);
			real = magnitude * leftCosineTerm;
			imag = magnitude * leftSineTerm;
			msr->timeDomainOut[0][bitRevFwd] = (real + imag) * 1.2f;
			msr->timeDomainOut[0][bitRevSym] = (real - imag) * 1.2f;
			lerpCompCplx(1);
			real = magnitude * rightCosineTerm;
			imag = magnitude * rightSineTerm;
			msr->timeDomainOut[1][bitRevFwd] = (real + imag) * 1.2f;
			msr->timeDomainOut[1][bitRevSym] = (real - imag) * 1.2f;
			lerpCompCplx(2);
			real = magnitude * centreCosineTerm;
			imag = magnitude * centreSineTerm;
			msr->timeDomainOut[2][bitRevFwd] = (real + imag) * 1.2f;
			msr->timeDomainOut[2][bitRevSym] = (real - imag) * 1.2f;
			lerpCompCplx(3);
			real = magnitude * leftCosineTerm;
			imag = magnitude * leftSineTerm;
			msr->timeDomainOut[3][bitRevFwd] = (real + imag) * minusPmix + opLFwd * pmix;
			msr->timeDomainOut[3][bitRevSym] = (real - imag) * minusPmix + opLSym * pmix;
			lerpCompCplx(4);
			real = magnitude * rightCosineTerm;
			imag = magnitude * rightSineTerm;
			msr->timeDomainOut[4][bitRevFwd] = (real + imag) * minusPmix + opRFwd * pmix;
			msr->timeDomainOut[4][bitRevSym] = (real - imag) * minusPmix + opRSym * pmix;
		}
	}
	else if (numOut == 8) // 7.1
	{
		float pmix = msr->mix * 0.2f;
		float minusPmix = 1.0f - pmix;
		float bpmix = msr->mix * (1.0f - 0.2f);
		float minusBpmix = 1.0f - bpmix;
		lerpCompCplx(0);
		msr->timeDomainOut[0][0] = (magnitude * leftCosineTerm) * 1.2f;
		lerpCompCplx(1);
		msr->timeDomainOut[1][0] = (magnitude * rightCosineTerm) * 1.2f;
		lerpCompCplx(2);
		msr->timeDomainOut[2][0] = (magnitude * centreCosineTerm) * 1.2f;
		lerpCompCplx(3);
		msr->timeDomainOut[3][0] = (magnitude * leftCosineTerm) * minusPmix + opLFwd * pmix;
		lerpCompCplx(4);
		msr->timeDomainOut[4][0] = (magnitude * rightCosineTerm) * minusPmix + opRFwd * pmix;
		lerpCompCplx(5);
		msr->timeDomainOut[5][0] = (magnitude * leftCosineTerm) * minusBpmix + opLFwd * bpmix;
		lerpCompCplx(6);
		msr->timeDomainOut[6][0] = (magnitude * rightCosineTerm) * minusBpmix + opRFwd * bpmix;
		for (i = 1; i < msr->procUpTo; i++)
		{
			mainLoop;
			lerpCompCplx(0);
			real = magnitude * leftCosineTerm;
			imag = magnitude * leftSineTerm;
			msr->timeDomainOut[0][bitRevFwd] = (real + imag) * 1.2f;
			msr->timeDomainOut[0][bitRevSym] = (real - imag) * 1.2f;
			lerpCompCplx(1);
			real = magnitude * rightCosineTerm;
			imag = magnitude * rightSineTerm;
			msr->timeDomainOut[1][bitRevFwd] = (real + imag) * 1.2f;
			msr->timeDomainOut[1][bitRevSym] = (real - imag) * 1.2f;
			lerpCompCplx(2);
			real = magnitude * centreCosineTerm;
			imag = magnitude * centreSineTerm;
			msr->timeDomainOut[2][bitRevFwd] = (real + imag) * 1.2f;
			msr->timeDomainOut[2][bitRevSym] = (real - imag) * 1.2f;
			lerpCompCplx(3);
			real = magnitude * leftCosineTerm;
			imag = magnitude * leftSineTerm;
			msr->timeDomainOut[3][bitRevFwd] = (real + imag) * minusPmix + opLFwd * pmix;
			msr->timeDomainOut[3][bitRevSym] = (real - imag) * minusPmix + opLSym * pmix;
			lerpCompCplx(4);
			real = magnitude * rightCosineTerm;
			imag = magnitude * rightSineTerm;
			msr->timeDomainOut[4][bitRevFwd] = (real + imag) * minusPmix + opRFwd * pmix;
			msr->timeDomainOut[4][bitRevSym] = (real - imag) * minusPmix + opRSym * pmix;
			lerpCompCplx(5);
			real = magnitude * leftCosineTerm;
			imag = magnitude * leftSineTerm;
			msr->timeDomainOut[5][bitRevFwd] = (real + imag) * minusBpmix + opLFwd * bpmix;
			msr->timeDomainOut[5][bitRevSym] = (real - imag) * minusBpmix + opLSym * bpmix;
			lerpCompCplx(6);
			real = magnitude * rightCosineTerm;
			imag = magnitude * rightSineTerm;
			msr->timeDomainOut[6][bitRevFwd] = (real + imag) * minusBpmix + opRFwd * bpmix;
			msr->timeDomainOut[6][bitRevSym] = (real - imag) * minusBpmix + opRSym * bpmix;
		}
	}
	for (i = msr->procUpTo; i < msr->halfLen; i++)
	{
		symIdx = msr->fftLen - i;
		bitRevFwd = msr->mBitRev[i];
		bitRevSym = msr->mBitRev[symIdx];
		msr->timeDomainOut[0][bitRevSym] = msr->mTempLBuffer[i];
		msr->timeDomainOut[0][bitRevFwd] = msr->mTempLBuffer[symIdx];
		msr->timeDomainOut[1][bitRevSym] = msr->mTempRBuffer[i];
		msr->timeDomainOut[1][bitRevFwd] = msr->mTempRBuffer[symIdx];
		for (unsigned int c = 2; c < numOut - 1; c++)
			msr->timeDomainOut[c][bitRevSym] = msr->timeDomainOut[c][bitRevFwd] = 0;
	}
	for (j = 0; j < numOut - 1; j++)
		msr->fft(msr->timeDomainOut[j], msr->mSineTab);
	msr->mOutputBufferCount++;
	if (msr->mOutputBufferCount > MAX_OUTPUT_BUFFERS)
		return;
	float *outBuffer = msr->mOutputBuffer[msr->mOutputBufferCount - 1];
	for (i = 0; i < msr->ovpLen; ++i)
	{
		// Subwoofer index
		unsigned int k = i + msr->mInputPos + msr->smpShift;
		if (k >= msr->frameLen)
			k = k - msr->frameLen;
		outBuffer[0] = msr->mOverlapStage2dash[0][i] + (msr->timeDomainOut[0][i + msr->smpShift] * msr->synthesisWnd[i]) + msr->mInputSub[0][k] * msr->minuscrossBass;
		msr->mOverlapStage2dash[0][i] = (msr->timeDomainOut[0][msr->smpShift + msr->ovpLen + i] * msr->synthesisWnd[i + msr->ovpLen]);
		outBuffer[1] = msr->mOverlapStage2dash[1][i] + (msr->timeDomainOut[1][i + msr->smpShift] * msr->synthesisWnd[i]) + msr->mInputSub[1][k] * msr->minuscrossBass;
		msr->mOverlapStage2dash[1][i] = (msr->timeDomainOut[1][msr->smpShift + msr->ovpLen + i] * msr->synthesisWnd[i + msr->ovpLen]);
		for (j = 2; j < numOut - 1; j++)
		{
			outBuffer[j] = msr->mOverlapStage2dash[j][i] + (msr->timeDomainOut[j][i + msr->smpShift] * msr->synthesisWnd[i]);
			msr->mOverlapStage2dash[j][i] = (msr->timeDomainOut[j][msr->smpShift + msr->ovpLen + i] * msr->synthesisWnd[i + msr->ovpLen]);
		}
		outBuffer[numOut - 1] = (msr->mInputSub[0][k] + msr->mInputSub[1][k]) * msr->crossBass;
		outBuffer += numOut;
	}
	msr->mInputSamplesNeeded = msr->ovpLen;
}
void LUTSurroundRefreshParameter(LUTSurround *msr)
{
	if (msr->spread < 0.3f)
		msr->spread = 0.3f;
	if (msr->spread > 2.0f * M_PI - 0.3f)
		msr->spread = 2.0f * M_PI - 0.3f;
	if (msr->separation < -0.3f)
		msr->separation = -0.3f;
	if (msr->separation > 0.3f)
		msr->separation = 0.3f;
	if (msr->spatialEnhancement < 0.0f)
		msr->spatialEnhancement = 0.0f;
	if (msr->spatialEnhancement > 2.0f)
		msr->spatialEnhancement = 2.0f;
	float x, y;
	for (int i = 0; i < TBL_SIZE; i++)
	{
		for (int j = 0; j < TBL_SIZE; j++)
		{
			cartesianPrecomputer(msr->XYaxis[i], msr->XYaxis[j], msr->spread, msr->separation, &x, &y);
			// Spatial enhancement
			x = clampFLT(x * (msr->spatialEnhancement + msr->spatialEnhancement * y + 1.0f - y) * 0.5f);
			msr->magPhaseDiff2Cartesian[0 * TBL_SIZE * TBL_SIZE + i * TBL_SIZE + j] = x;
			msr->magPhaseDiff2Cartesian[1 * TBL_SIZE * TBL_SIZE + i * TBL_SIZE + j] = y;
		}
	}
}
void LUTSurroundSetSmoothing(LUTSurround *msr, float sm)
{
	if (sm < 10.0f)
		sm = 10.0f;
	if (sm > 60.0f)
		sm = 60.0f;
	msr->smoothing = sm / msr->fs * (msr->frameLen / ANALYSIS_OVERLAP);
	msr->mix = 1.0f - tanhf(0.3f * powf(sm, 0.5f));
}
void getAsymmetricWindow(float *analysisWnd, float *synthesisWnd, int k, int m, int p, double freq_temporal)
{
	int i;
	if ((k / m) < 4)
		freq_temporal = 1.0f;
	if (freq_temporal > 9.0f)
		freq_temporal = 9.0f;
	memset(synthesisWnd, 0, k * sizeof(float));
	int n = ((k - m) << 1) + 2;
	for (i = 0; i < k - m; ++i)
		analysisWnd[i] = (float)pow(0.5 * (1.0 - cos(2.0 * M_PI * (i + 1.0) / (double)n)), freq_temporal);
	n = (m << 1) + 2;
	if (freq_temporal > 1.02)
		freq_temporal = 1.02;
	for (i = k - m; i < k; ++i)
		analysisWnd[i] = (float)pow(sqrt(0.5 * (1.0 - cos(2.0 * M_PI * ((m + i - (k - m)) + 1.0) / (double)n))), freq_temporal);
	n = m << 1;
	for (i = k - (m << 1); i < k; ++i)
		synthesisWnd[i] = (float)(0.5 * (1.0 - cos(2.0 * M_PI * (double)(i - (k - (m << 1))) / (double)n))) / analysisWnd[i];
	// Pre-shift window function
	for (i = 0; i < k - p; i++)
		synthesisWnd[i] = synthesisWnd[i + p];
}
size_t fast_upper_bound3(const float *vec, size_t n, float *value)
{
	size_t size = n;
	size_t low = 0;
	while (size > 0)
	{
		size_t half = size / 2;
		size_t other_half = size - half;
		size_t probe = low + half;
		size_t other_low = low + other_half;
		size = half;
		low = (*value >= vec[probe]) ? other_low : low;
	}
	return low;
}
static float lerp1DNoExtrapo(float val, const float *x, const float *y, int n)
{
	if (val <= x[0])
		return y[0];
	if (val >= x[n - 1])
		return y[n - 1];
	size_t j = fast_upper_bound3(x, n, &val);
	return ((val - x[j - 1]) / (x[j] - x[j - 1])) * (y[j] - y[j - 1]) + y[j - 1]; // Interpolation
}
static unsigned int closestInteger(unsigned int a, unsigned int b)
{
	unsigned int c1 = a - (a % b);
	unsigned int c2 = (a + b) - (a % b);
	if (a - c1 > c2 - a)
		return c2;
	else
		return c1;
}
static unsigned int LLIntegerLog2(unsigned int v)
{
	unsigned int i = 0;
	while (v > 1)
	{
		++i;
		v >>= 1;
	}
	return i;
}
static unsigned LLRevBits(unsigned int x, unsigned int bits)
{
	unsigned int y = 0;
	while (bits--)
	{
		y = (y + y) + (x & 1);
		x >>= 1;
	}
	return y;
}
static void LLbitReversalTbl(unsigned *dst, int n)
{
	unsigned int bits = LLIntegerLog2(n);
	for (int i = 0; i < n; ++i)
		dst[i] = LLRevBits(i, bits);
}
void InitLR2(SecondOrderCrossover *lr2, double fs, double fc)
{
	double fpi = M_PI * fc;
	double wc = 2 * fpi;
	double wc2 = wc * wc;
	double wc22 = 2 * wc2;
	double k = wc / tan(fpi / fs);
	double k2 = k * k;
	double k22 = 2 * k2;
	double wck2 = 2 * wc*k;
	double tmpk = (k2 + wc2 + wck2);
	lr2->b1 = (-k22 + wc22) / tmpk;
	lr2->b2 = (-wck2 + k2 + wc2) / tmpk;
	lr2->a1_lp = (wc22) / tmpk;
	lr2->a2_lp = (wc2) / tmpk;
	lr2->a1_hp = (-k22) / tmpk;
	lr2->a2_hp = (k2) / tmpk;
}
void ClearStateLR2(SecondOrderCrossover *lr2)
{
	lr2->lp1_xm0 = lr2->lp1_xm1 = lr2->lp2_xm0 = lr2->lp2_xm1 = lr2->hp1_xm0 = lr2->hp1_xm1 = lr2->hp2_xm0 = lr2->hp2_xm1 = 0.0;
}
void processLR2(SecondOrderCrossover *lr2, float *in1, float *in2, float *lp1, float *lp2, float *hp1, float *hp2)
{
	double lp1_out = lr2->a2_lp * *in1 + lr2->lp1_xm0;
	lr2->lp1_xm0 = lr2->a1_lp * *in1 - lr2->b1 * lp1_out + lr2->lp1_xm1;
	lr2->lp1_xm1 = lr2->a2_lp * *in1 - lr2->b2 * lp1_out;
	double lp2_out = lr2->a2_lp * *in2 + lr2->lp2_xm0;
	lr2->lp2_xm0 = lr2->a1_lp * *in2 - lr2->b1 * lp2_out + lr2->lp2_xm1;
	lr2->lp2_xm1 = lr2->a2_lp * *in2 - lr2->b2 * lp2_out;
	double hp1_out = lr2->a2_hp * *in1 + lr2->hp1_xm0;
	lr2->hp1_xm0 = lr2->a1_hp * *in1 - lr2->b1 * hp1_out + lr2->hp1_xm1;
	lr2->hp1_xm1 = lr2->a2_hp * *in1 - lr2->b2 * hp1_out;
	double hp2_out = lr2->a2_hp * *in2 + lr2->hp2_xm0;
	lr2->hp2_xm0 = lr2->a1_hp * *in2 - lr2->b1 * hp2_out + lr2->hp2_xm1;
	lr2->hp2_xm1 = lr2->a2_hp * *in2 - lr2->b2 * hp2_out;
	*lp1 = lp1_out;
	*lp2 = lp2_out;
	*hp1 = hp1_out;
	*hp2 = hp2_out;
}
void LUTSurroundSetSubwoofer(LUTSurround *msr, float fc)
{
	InitLR2(&msr->subBandSplit, msr->fs, fc < 40.0f ? 40.0 : (double)fc);
}
void LUTSurroundSetBassCross(LUTSurround *msr, float crossRat)
{
	if (crossRat < 0.0f)
		crossRat = 0.0f;
	if (crossRat > 1.0f)
		crossRat = 1.0f;
	msr->crossBass = crossRat;
	msr->minuscrossBass = 1.0f - msr->crossBass;
}
// Common output buffer macro switching
#define expandOutputBufferOp \
outSampleCount += copyCount; \
msr->mOutputReadSampleOffset += copyCount; \
if (msr->mOutputReadSampleOffset >= msr->ovpLen) \
{ \
	msr->mOutputBufferCount--; \
	msr->mOutputReadSampleOffset = 0; \
	if (msr->mOutputBufferCount > 0) { \
		unsigned int i; \
		float *moveToEnd = msr->mOutputBuffer[0]; \
		for (i = 1; i < MAX_OUTPUT_BUFFERS; i++) \
			msr->mOutputBuffer[i - 1] = msr->mOutputBuffer[i]; \
		msr->mOutputBuffer[MAX_OUTPUT_BUFFERS - 1] = 0; \
		for (i = 0; i < MAX_OUTPUT_BUFFERS; i++) \
		{ \
			if (!msr->mOutputBuffer[i]) \
			{ \
				msr->mOutputBuffer[i] = moveToEnd; \
				break; \
			} \
		} \
	} \
}
// Common input buffer macro for inputs
#define expandInputBufferOp \
unsigned int outSampleCount, maxOutSampleCount, copyCount; \
outSampleCount = 0; \
maxOutSampleCount = inSampleCount; \
float subTmp1, subTmp2; \
while (inSampleCount > 0) \
{ \
	copyCount = min(msr->mInputSamplesNeeded, inSampleCount); \
	float *sampDL = &msr->mInput[0][msr->mInputPos]; \
	float *sampDR = &msr->mInput[1][msr->mInputPos]; \
	float *sampDSubLeft = &msr->mInputSub[0][msr->mInputPos]; \
	float *sampDSubRight = &msr->mInputSub[1][msr->mInputPos]; \
	const float *bufBound = inLeft + copyCount; \
	while (inLeft < bufBound) \
	{ \
		processLR2(&msr->subBandSplit, inLeft, inRight, &subTmp1, &subTmp2, sampDL, sampDR); \
		*sampDSubLeft = -subTmp1; \
		*sampDSubRight = -subTmp2; \
		inLeft += 1; \
		inRight += 1; \
		sampDL += 1; \
		sampDR += 1; \
		sampDSubLeft += 1; \
		sampDSubRight += 1; \
	} \
	inSampleCount -= copyCount; \
	msr->mInputPos = msr->mInputPos + copyCount; \
	if (msr->mInputPos >= msr->frameLen) \
		msr->mInputPos = msr->mInputPos - msr->frameLen; \
	msr->mInputSamplesNeeded -= copyCount; \
	if (msr->mInputSamplesNeeded == 0) \
		LLPAMSProcessNPR(msr, chmaskPtr, numOut); \
}
// Macro for getting buffer(Multiplexed) pointer and deinterleaved output pointer
#define getBuf \
float *sampD = msr->mOutputBuffer[0]; \
copyCount = min(msr->ovpLen - msr->mOutputReadSampleOffset, maxOutSampleCount - outSampleCount); \
float *out = sampD + (msr->mOutputReadSampleOffset * numOut); \
float *max = components[0] + copyCount;
// Macro for incrementing output pointer
#define incrementPtr(idx) \
*components[idx] = (float)*out * 1.12202f; \
out++; \
components[idx]++;
void LUTSurroundProcessSamplesLCR(LUTSurround *msr, float *inLeft, float *inRight, unsigned int inSampleCount, float *components[8])
{
	const float * const chmaskPtr = LCR_mask;
	const unsigned int numOut = 4;
	expandInputBufferOp
	while ((msr->mOutputBufferCount > 0) && (outSampleCount < maxOutSampleCount))
	{
		getBuf;
		while (components[0] < max)
		{
			incrementPtr(0);
			incrementPtr(1);
			incrementPtr(2);
			incrementPtr(3);
		}
		expandOutputBufferOp
	}
}
void LUTSurroundProcessSamplesQuadra(LUTSurround *msr, float *inLeft, float *inRight, unsigned int inSampleCount, float *components[8])
{
	const float * const chmaskPtr = Quadraphonic_mask;
	const unsigned int numOut = 5;
	expandInputBufferOp
	while ((msr->mOutputBufferCount > 0) && (outSampleCount < maxOutSampleCount))
	{
		getBuf;
		while (components[0] < max)
		{
			incrementPtr(0);
			incrementPtr(1);
			incrementPtr(4);
			incrementPtr(5);
			incrementPtr(3);
		}
		expandOutputBufferOp
	}
}
void LUTSurroundProcessSamples5_1(LUTSurround *msr, float *inLeft, float *inRight, unsigned int inSampleCount, float *components[8])
{
	const float * const chmaskPtr = FivePt1_mask;
	const unsigned int numOut = 6;
	expandInputBufferOp
	while ((msr->mOutputBufferCount > 0) && (outSampleCount < maxOutSampleCount))
	{
		getBuf;
		while (components[0] < max)
		{
			incrementPtr(0);
			incrementPtr(1);
			incrementPtr(2);
			incrementPtr(4);
			incrementPtr(5);
			incrementPtr(3);
		}
		expandOutputBufferOp
	}
}
void LUTSurroundProcessSamples7_1(LUTSurround *msr, float *inLeft, float *inRight, unsigned int inSampleCount, float *components[8])
{
	const float * const chmaskPtr = SevenPt1_mask;
	const unsigned int numOut = 8;
	expandInputBufferOp
	while ((msr->mOutputBufferCount > 0) && (outSampleCount < maxOutSampleCount))
	{
		getBuf;
		while (components[0] < max)
		{
			incrementPtr(0);
			incrementPtr(1);
			incrementPtr(2);
			incrementPtr(6);
			incrementPtr(7);
			incrementPtr(4);
			incrementPtr(5);
			incrementPtr(3);
		}
		expandOutputBufferOp
	}
}
void LUTSurroundSetOutputChannel(LUTSurround *msr, unsigned int mode)
{
	if (mode == 0)
		msr->process = LUTSurroundProcessSamplesLCR;
	else if (mode == 1)
		msr->process = LUTSurroundProcessSamplesQuadra;
	else if (mode == 2)
		msr->process = LUTSurroundProcessSamples5_1;
	else if (mode == 3)
		msr->process = LUTSurroundProcessSamples7_1;
}
void LUTSurroundInit(LUTSurround *msr, float fs)
{
	memset(msr, 0, sizeof(LUTSurround));
	unsigned int i;
	msr->fs = fs;
	const float oX[6] = { 11025, 22050, 44100, 88200, 176400, 352800 };
	const float oY[6] = { 384, 768, 1536, 3072, 6144, 8192 };
	unsigned int proposedFrameLen = (unsigned int)lerp1DNoExtrapo(fs, oX, oY, 6);
	if (proposedFrameLen % ANALYSIS_OVERLAP != 0)
		proposedFrameLen = closestInteger(proposedFrameLen, ANALYSIS_OVERLAP);
	unsigned int nextPwr2 = (unsigned int)pow(2.0, ceil(log((double)proposedFrameLen) / log(2.0)));
	msr->frameLen = proposedFrameLen;
	msr->fftLen = nextPwr2;
	if (msr->fftLen == 512)
		msr->fft = DFT512;
	else if (msr->fftLen == 1024)
		msr->fft = DFT1024;
	else if (msr->fftLen == 2048)
		msr->fft = DFT2048;
	else if (msr->fftLen == 4096)
		msr->fft = DFT4096;
	else if (msr->fftLen == 8192)
		msr->fft = DFT8192;
	msr->minus_frameLen = msr->frameLen;
	msr->ovpLen = msr->frameLen / ANALYSIS_OVERLAP;
	msr->halfLen = (msr->fftLen >> 1) + 1;
	msr->smpShift = (msr->frameLen - (msr->ovpLen << 1));
	const float desiredProcessFreq = 24000.0f;
	unsigned int idx = (unsigned int)(desiredProcessFreq / (fs / (float)msr->fftLen)) + 1UL;
	if (idx > msr->halfLen)
		msr->procUpTo = msr->halfLen;
	else
		msr->procUpTo = idx;
	for (i = 0; i < MAX_OUTPUT_BUFFERS; i++)
		msr->mOutputBuffer[i] = msr->buffer[i];
	msr->mInputSamplesNeeded = 0;
	msr->mInputPos = 0;
	msr->mOutputBufferCount = 0;
	msr->mOutputReadSampleOffset = 0;
	LLbitReversalTbl(msr->mBitRev, msr->fftLen);
	const double twopi_over_n = 6.283185307179586476925286766559 / msr->fftLen;
	for (i = 0; i < msr->fftLen >> 1; ++i)
		msr->mSineTab[i] = (float)sin(twopi_over_n * i);
	getAsymmetricWindow(msr->analysisWnd, msr->synthesisWnd, msr->frameLen, msr->ovpLen, msr->smpShift, 1.0f);
	for (i = 0; i < msr->frameLen; i++)
		msr->analysisWnd[i] *= (1.0f / msr->fftLen) * 0.5f;
	msr->spread = M_PI / 2.0f; // 0 <-> 2pi
	msr->separation = 0.0f; // -0.3 <-> +0.3
	msr->spatialEnhancement = 1.0f;
#define linspace(vec, n, a, b) { double d = ((b) - (a)) / (double)((n) - 1); for (unsigned int i = 0; i < n; i++) vec[i] = (a) + i * d; }
	linspace(msr->XYaxis, TBL_SIZE, -1.0, 1.0);
	LUTSurroundRefreshParameter(msr);
	LUTSurroundSetSmoothing(msr, 15.0f);
	for (i = 0; i < 8; i++)
		for (unsigned int j = 0; j < HALFWNDLEN; j++)
			msr->smoothedFunc[i][j] = 0.5f;
	LUTSurroundSetSubwoofer(msr, 70.0f);
	LUTSurroundSetOutputChannel(msr, 0);
	LUTSurroundSetBassCross(msr, 0.8f);
	msr->setSubwoofer = LUTSurroundSetSubwoofer;
	msr->setBassCrossLR = LUTSurroundSetBassCross;
	msr->setSmoothing = LUTSurroundSetSmoothing;
	msr->setOutputMode = LUTSurroundSetOutputChannel;
	msr->refresh = LUTSurroundRefreshParameter;
}