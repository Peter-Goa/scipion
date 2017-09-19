/*
 * cuda_gpu_reconstruct_fourier.h
 *
 *  Created on: Aug 11, 2017
 *      Author: david
 */

#ifndef XMIPP_LIBRARIES_RECONSTRUCTION_CUDA_CUDA_GPU_RECONSTRUCT_FOURIER_H_
#define XMIPP_LIBRARIES_RECONSTRUCTION_CUDA_CUDA_GPU_RECONSTRUCT_FOURIER_H_

#include <data/projection_data.h>
#include <data/fourier_reconstruction_traverse_space.h>
#include <reconstruction_adapt_cuda/xmipp_gpu_utils.h>
#include <reconstruction_cuda/cuda_xmipp_utils.h>


//static ProjectionData* projData;



struct ProjectionDataGPU
{
	bool skip;
	float* img;
	float* CTF;
	float* modulator;
	int xSize;
	int ySize;
	int imgIndex;
	float weight;
	float localAInv[3][3];
	float localA[3][3];
public:
	ProjectionDataGPU() {
			setDefault();
		}
	ProjectionDataGPU(const ProjectionData& data) {
		skip = data.skip;
		if (skip) {
			setDefault();
			return;
		}
//		copy(*data.img, img);
//		if (NULL != data.CTF) {
			img = NULL;
//			copy(*data.CTF, CTF);
//		} else {
			CTF = NULL;
//		}
//		if (NULL != data.modulator) {
//			copy(*data.modulator, modulator);
//		} else {
			modulator = NULL;
//		}
//		xSize = data.img->getXSize();
//		ySize = data.img->getYSize();
		imgIndex = data.imgIndex;
		weight = data.weight;
		data.localAInv.convertTo(localAInv);
		data.localA.convertTo(localA);
	}
	/** Remove stored data and set to skip */
	void clean();
	void setDefault() {
		skip = true;
		img = NULL;
		CTF = NULL;
		modulator = NULL;
		xSize = 0;
		ySize = 0;
		imgIndex = -1;
		weight = 0;
	}
private:
	template<typename T, typename U>
	void copy(const Array2D<T>& from, U& to);
};


/** Method to allocate 3D array (continuous) of given size^3 */
float* allocateGPU(float*& where, int size, int typeSize);
void releaseGPU(float*& where);
void copyTempSpaces(std::complex<float>*** tempVol, float*** tempWeights,
		float* tempVolGPU, float* tempWeightsGPU,
		int size);



ProjectionDataGPU* prepareBuffer(GpuMultidimArrayAtGpu<float>& ffts,
		int sizeX, int sizeY, int paddedImgSize, float maxResolutionSqr, int, float*&);



void processBufferGPU(
		float* tempVolumeGPU,
				float* tempWeightsGPU,
				ProjectionData* data,int noOfPixel,int bufferSize,
				TraverseSpace* traverseSpaces, int noOfTransforms,
				int maxVolIndexX, int maxVolIndexYZ,
				bool useFast, float blobRadius,
				float iDeltaSqrt,
				float* blobTableSqrt, int blobTableSize,
				ProjectionDataGPU* images);
void getTempSpaces(int size, std::complex<float>***& volume, float***& tempWeights);
void copyBuffer(ProjectionData* data, int size);
void allocateTempSpaces(int size);

void a(std::complex<float>* a);
void b();
void c(std::complex<float>* a);

void moje();

static int num = 20;

#endif /* XMIPP_LIBRARIES_RECONSTRUCTION_CUDA_CUDA_GPU_RECONSTRUCT_FOURIER_H_ */