#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
//#include <vld.h>
#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h"
const char *get_filename_ext(const char *filename)
{
	const char *dot = strrchr(filename, '.');
	if (!dot || dot == filename) return "";
	return dot + 1;
}
float* loadAudioFile(char *filename, unsigned int *fs, unsigned int *channels, drwav_uint64 *totalPCMFrameCount)
{
	const char *ext = get_filename_ext(filename);
	float *pSampleData = 0;
	if (!strncmp(ext, "wav", 5))
		pSampleData = drwav_open_file_and_read_pcm_frames_f32(filename, channels, fs, totalPCMFrameCount, 0);
	if (pSampleData == NULL)
	{
		printf("Error opening and reading WAV file");
		return 0;
	}
	// Sanity check
	if (*channels < 1)
	{
		printf("Invalid audio channels count");
		free(pSampleData);
		return 0;
	}
	if ((*totalPCMFrameCount <= 0) || (*totalPCMFrameCount <= 0))
	{
		printf("Invalid audio sample rate / frame count");
		free(pSampleData);
		return 0;
	}
	return pSampleData;
}
#include "LUTSurround.h"
//  Windows
#ifdef _WIN32
#include <Windows.h>
////////////////////////////////////////////////////////////////////
// Performance timer
double get_wall_time()
{
	LARGE_INTEGER time, freq;
	if (!QueryPerformanceFrequency(&freq))
		return 0;
	if (!QueryPerformanceCounter(&time))
		return 0;
	return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time()
{
	FILETIME a, b, c, d;
	if (GetProcessTimes(GetCurrentProcess(), &a, &b, &c, &d) != 0)
		return (double)(d.dwLowDateTime | ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
	else
		return 0;
}
#else
#include <time.h>
#include <sys/time.h>
double get_wall_time()
{
	struct timeval time;
	if (gettimeofday(&time, NULL))
		return 0;
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time()
{
	return (double)clock() / CLOCKS_PER_SEC;
}
#endif
////////////////////////////////////////////////////////////////////
char *inputString(FILE* fp, size_t size) {
	//The size is extended by the input with the value of the provisional
	char *str;
	int ch;
	size_t len = 0;
	str = (char*)realloc(NULL, sizeof(char)*size);//size is start size
	if (!str)return str;
	while (EOF != (ch = fgetc(fp)) && ch != '\n') {
		str[len++] = ch;
		if (len == size) {
			str = (char*)realloc(str, sizeof(char)*(size += 16));
			if (!str)return str;
		}
	}
	str[len++] = '\0';
	return (char*)realloc(str, sizeof(char)*len);
}
char *basename(char const *path)
{
#ifdef _MSC_VER
	char *s = (char*)strrchr(path, '\\');
#else
	char *s = strrchr(path, '/');
#endif
	if (!s)
		return strdup(path);
	else
		return strdup(s + 1);
}
void channel_join(float **chan_buffers, unsigned int num_channels, float *buffer, unsigned int num_frames)
{
	unsigned int i, samples = num_frames * num_channels;
	for (i = 0; i < samples; i++)
		buffer[i] = chan_buffers[i % num_channels][i / num_channels];
}
void channel_split(float *buffer, unsigned int num_frames, float **chan_buffers, unsigned int num_channels)
{
	unsigned int i, samples = num_frames * num_channels;
	for (i = 0; i < samples; i++)
		chan_buffers[i % num_channels][i / num_channels] = buffer[i];
}
int main(int argc, char *argv[])
{
	char *path = "test09.wav";
	char *filename = basename(path);
	unsigned int i, j;
	unsigned int channels;
	unsigned int sampleRate;
	drwav_uint64 totalPCMFrameCount;
	float* pFrameBuffer = loadAudioFile(path, &sampleRate, &channels, &totalPCMFrameCount);
	if (pFrameBuffer == NULL)
	{
		printf("Error opening and reading WAV file");
		return -1;
	}
	// Sanity check
	if (channels < 1)
	{
		printf("Invalid audio channels count");
		free(pFrameBuffer);
		return -1;
	}
	if ((totalPCMFrameCount <= 0) || (totalPCMFrameCount <= 0))
	{
		printf("Invalid audio sample rate / frame count");
		free(pFrameBuffer);
		return -1;
	}
	printf("[Status] Precomputing...\n");
	double wall0 = get_wall_time();
	double cpu0 = get_cpu_time();
	unsigned char outputMode = 2; // 0 -> LCR, 1 -> Quadraphonic, 2 -> 5.1, 3 -> 7.1
	unsigned int numCh[4] = { 4, 6, 6, 8 };
	LUTSurround msr;
	LUTSurroundInit(&msr, sampleRate);
	//msr.setSmoothing(&msr, 15.0f); // Default = 15, Recommended range = 12 - 20
	msr.setSubwoofer(&msr, 70.0f); // Default = 15, Recommended range = 30 - 120
	msr.setBassCrossLR(&msr, 0.8f); // Default = 0.8, Recommended range = 0.0 - 1.0
	msr.setOutputMode(&msr, outputMode);
	// You need to call of refresh() member after adjusting below parameters
	msr.spread = 1.57; // Default = 1.57, Recommended range 0.3 <-> 5.98
	msr.separation = 0.0; // Default = 0.0, Recommended range = -0.3 <-> 0.3
	msr.spatialEnhancement = 1.0; // Default = 1.0, Recommended range = 0.0 <-> 1.0
	msr.refresh(&msr);
	// 
	double wall1 = get_wall_time();
	double cpu1 = get_cpu_time();
	printf("[Info] Precomputation of filters take:\n[Info] %lf sec Wall time, %lf sec CPU time\n", wall1 - wall0, cpu1 - cpu0);
	printf("[Info] Samples: %I64d, Channels: %d\n[Status] Allocating buffer\n", totalPCMFrameCount, channels);
	unsigned int frameCountProvided = 128;
	unsigned int readcount = (unsigned int)ceil((float)totalPCMFrameCount / (float)frameCountProvided);
	unsigned int finalSize = frameCountProvided * readcount;
	const unsigned int targetChannels = numCh[outputMode];
	float **splittedBuffer = (float**)malloc(channels * sizeof(float*));
	float **outBuffer = (float**)malloc(targetChannels * sizeof(float*));
	for (i = 0; i < channels; i++)
	{
		splittedBuffer[i] = (float*)malloc(finalSize * sizeof(float));
		memset(splittedBuffer[i], 0, finalSize * sizeof(float));
	}
	for (i = 0; i < targetChannels; i++)
	{
		outBuffer[i] = (float*)malloc(finalSize * sizeof(float));
		memset(outBuffer[i], 0, finalSize * sizeof(float));
	}
	channel_split(pFrameBuffer, totalPCMFrameCount, splittedBuffer, channels);
	free(pFrameBuffer);
	printf("[Info] Processing...\n");
	float *ptr[8];
	wall1 = get_wall_time();
	cpu1 = get_cpu_time();
	for (i = 0; i < readcount; i++)
	{
		unsigned int pointerOffset = frameCountProvided * i;
		float *ptrInLeft = splittedBuffer[0] + pointerOffset;
		float *ptrInRight = splittedBuffer[1] + pointerOffset;
		for (j = 0; j < targetChannels; j++)
			ptr[j] = outBuffer[j] + pointerOffset;
		float *inputs[2] = { ptrInLeft, ptrInRight };
		unsigned int offset = 0;
		float* outputs[8];
		while (offset < frameCountProvided)
		{
			const int processing = min(frameCountProvided - offset, msr.ovpLen);
			for (j = 0; j < targetChannels; j++)
				outputs[j] = ptr[j] + offset;
			msr.process(&msr, inputs[0] + offset, inputs[1] + offset, processing, outputs);
			offset += processing;
		}
	}
	printf("[Info] Audio samples processing take:\n[Info] %lf sec Wall time, %lf sec CPU time\n", get_wall_time() - wall1, get_cpu_time() - cpu1);
	pFrameBuffer = (float*)malloc(targetChannels * finalSize * sizeof(float));
	channel_join(outBuffer, targetChannels, pFrameBuffer, finalSize);
	for (i = 0; i < channels; i++)
		free(splittedBuffer[i]);
	free(splittedBuffer);
	for (i = 0; i < targetChannels; i++)
		free(outBuffer[i]);
	free(outBuffer);
	unsigned int totalFrames = finalSize * targetChannels;
	size_t bufsz = snprintf(NULL, 0, "%s_Processed.wav", filename);
	char *filenameNew = (char*)malloc(bufsz + 1);
	snprintf(filenameNew, bufsz + 1, "%s_Processed.wav", filename);
	free(filename);
	drwav pWav;
	drwav_data_format format;
	format.container = drwav_container_riff;
	format.format = DR_WAVE_FORMAT_IEEE_FLOAT;
	format.channels = targetChannels;
	format.sampleRate = sampleRate;
	format.bitsPerSample = 32;
	unsigned int fail = drwav_init_file_write(&pWav, filenameNew, &format, 0);
	drwav_uint64 framesWritten = drwav_write_pcm_frames(&pWav, totalPCMFrameCount, pFrameBuffer);
	drwav_uninit(&pWav);
	free(filenameNew);
	free(pFrameBuffer);
	//system("pause");
	return 0;
}