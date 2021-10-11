#define FFTSIZE (8192)
#define ANALYSIS_OVERLAP 4
#define OVPSIZE (FFTSIZE / ANALYSIS_OVERLAP)
#define SAMPLESHIFT (FFTSIZE - (OVPSIZE << 1))
#define MINUSFFTSIZE (FFTSIZE - 1)
#define HALFWNDLEN ((FFTSIZE >> 1) + 1)
#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif
#define MAX_OUTPUT_BUFFERS 2
#define TBL_SIZE (64)
typedef struct
{
	double b1, b2, a1_lp, a2_lp, a1_hp, a2_hp;
	double lp1_xm0, lp1_xm1, lp2_xm0, lp2_xm1, hp1_xm0, hp1_xm1, hp2_xm0, hp2_xm1;
} SecondOrderCrossover;
typedef struct lutUpmix
{
	unsigned int fftLen, minus_frameLen, frameLen, ovpLen, halfLen, smpShift, procUpTo;
	void(*fft)(float*, const float*);
	unsigned int mBitRev[FFTSIZE];
	float 	mSineTab[FFTSIZE >> 1];
	float analysisWnd[FFTSIZE];
	float synthesisWnd[FFTSIZE];
	unsigned int  mOutputReadSampleOffset;
	unsigned int  mOutputBufferCount;
	unsigned int mInputSamplesNeeded;
	float 	*mOutputBuffer[MAX_OUTPUT_BUFFERS];
	float buffer[MAX_OUTPUT_BUFFERS][OVPSIZE * 8];
	unsigned int mInputPos;
	float 	mInput[2][FFTSIZE];
	float	mInputSub[2][FFTSIZE];
	float 	mOverlapStage2dash[7][OVPSIZE];
	float 	mTempLBuffer[FFTSIZE];
	float 	mTempRBuffer[FFTSIZE];
	float smoothedFunc[9][HALFWNDLEN];
	float timeDomainOut[7][FFTSIZE];
	float magPhaseDiff2Cartesian[2 * TBL_SIZE * TBL_SIZE];
	float spread, separation, spatialEnhancement;
	float XYaxis[TBL_SIZE];
	SecondOrderCrossover subBandSplit;
	float smoothing, mix, crossBass, minuscrossBass, fs;
	void(*setSubwoofer)(struct lutUpmix *, float);
	void(*setBassCrossLR)(struct lutUpmix *, float);
	void(*setSmoothing)(struct lutUpmix *, float);
	void(*setOutputMode)(struct lutUpmix *, unsigned int);
	void(*refresh)(struct lutUpmix *);
	void(*process)(struct lutUpmix *, float *, float *, unsigned int, float *[8]);
} LUTSurround;
void LUTSurroundInit(LUTSurround *msr, float fs);