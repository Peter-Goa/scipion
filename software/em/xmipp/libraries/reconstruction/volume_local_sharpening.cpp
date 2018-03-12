/***************************************************************************
 *
 * Authors:    Jose Luis Vilas, 					  jlvilas@cnb.csic.es
 * 			   Carlos Oscar S. Sorzano            coss@cnb.csic.es (2016)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "volume_local_sharpening.h"
//#define DEBUG
//#define DEBUG_MASK
//#define DEBUG_FILTER
void ProgLocSharpening::readParams()
{
	fnVol = getParam("--vol");
	fnRes = getParam("--resolution_map");
	sampling = getDoubleParam("--sampling");
	lambda = getDoubleParam("-l");
	Niter = getIntParam("-n");
	fnOut = getParam("-o");
}


void ProgLocSharpening::defineParams()
{
	addUsageLine("This function performs local sharpening");
	addParamsLine("  --vol <vol_file=\"\">   : Input volume");
	addParamsLine("  --resolution_map <vol_file=\"\">: Resolution map");
	addParamsLine("  -o <output=\"Sharpening.vol\">: sharpening volume");
	addParamsLine("  --sampling <s=1>: sampling");
	addParamsLine("  -l <s=1>: regularization param");
	addParamsLine("  -n <s=5>: iteration");
}

void ProgLocSharpening::produceSideInfo()
{
	std::cout << "Starting..." << std::endl;
	Image<double> V;
    V.read(fnVol);
    V().setXmippOrigin();

    Image<double> resolutionVolume;
    resolutionVolume.read(fnRes);
    MultidimArray<double> resVol = resolutionVolume();
    resVol.setXmippOrigin();

	FourierTransformer transformer;
	MultidimArray<double> &inputVol = V();

	Vorig = inputVol;

	transformer.FourierTransform(inputVol, fftV);

	iu.initZeros(fftV);
	double uz, uy, ux, uz2, u2, uz2y2;
	long n=0;
	for(size_t k=0; k<ZSIZE(fftV); ++k)
	{
		FFT_IDX2DIGFREQ(k,ZSIZE(inputVol),uz);
		uz2=uz*uz;

		for(size_t i=0; i<YSIZE(fftV); ++i)
		{
			FFT_IDX2DIGFREQ(i,YSIZE(inputVol),uy);
			uz2y2=uz2+uy*uy;

			for(size_t j=0; j<XSIZE(fftV); ++j)
			{
				FFT_IDX2DIGFREQ(j,XSIZE(inputVol),ux);
				u2=uz2y2+ux*ux;
				if ((k != 0) || (i != 0) || (j != 0))
					DIRECT_MULTIDIM_ELEM(iu,n) =sqrt(u2);
				else
					DIRECT_MULTIDIM_ELEM(iu,n) = 1e-38;
				++n;
			}
		}
	}
	maxMinResolution(resVol, maxRes, minRes);
	std::cout << "maxRes = " << maxRes << "  minRes = " << minRes << std::endl;

	resVol2IdxVol(resVol, maxRes, idxVol, idxList);
	Image<int> filteredvolume;
	filteredvolume() = idxVol;
	filteredvolume.write(formatString("idxvolume.vol"));


	double u;
	int size_fourier = ZSIZE(fftV);
	freq_fourier.initZeros(size_fourier);

	int size = ZSIZE(inputVol);

	VEC_ELEM(freq_fourier,0) = 1e-38;

	for(size_t k=1; k<size_fourier; ++k)
	{
		FFT_IDX2DIGFREQ(k,size, u);
		VEC_ELEM(freq_fourier,k) = u;
//		std::cout << u << std::endl;
	}

}

void ProgLocSharpening::maxMinResolution(MultidimArray<double> &resVol,
										double &maxRes, double &minRes)
{
	// Count number of voxels with resolution
	size_t n=0;
	double lastMinRes=1e38, lastMaxRes=1e-38, value;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resVol)
	{
		value = DIRECT_MULTIDIM_ELEM(resVol, n);
		if (value>lastMaxRes)
			lastMaxRes = value;
		if (value<lastMinRes)
			lastMinRes = value;
	}

	maxRes = lastMaxRes;
	minRes = lastMinRes;
}


void ProgLocSharpening::lowPassFilterFunction(const MultidimArray< std::complex<double> > &myfftV,
		double w, double wL, MultidimArray<double> &filteredVol, int count)
{
	fftVfilter.initZeros(myfftV);
	size_t xdim, ydim, zdim, ndim;
	//amplitude.resizeNoCopy(filteredVol);
	Vorig.getDimensions(xdim, ydim, zdim, ndim);
	MultidimArray<double> testVol;
	testVol.initZeros(myfftV);
	filteredVol.resizeNoCopy(Vorig);
	//resizeNoCopy(sharpenedMap);

	// Filter the input volume and add it to amplitude
	long n=0;
	double ideltal=PI/(wL-w);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(myfftV)
	{
		double un=DIRECT_MULTIDIM_ELEM(iu,n);
		if (un>=w && un<wL)
		{
			//double H=0.5*(1+cos((un-w1)*ideltal));
			DIRECT_MULTIDIM_ELEM(fftVfilter, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
			DIRECT_MULTIDIM_ELEM(fftVfilter, n) *= 0.5*(1+cos((un-w)*ideltal));//H;
//			DIRECT_MULTIDIM_ELEM(testVol, n) = 0.5;
		} else if (un<w)
		{
			DIRECT_MULTIDIM_ELEM(fftVfilter, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
//			DIRECT_MULTIDIM_ELEM(testVol, n) = 1;
		}
	}

	transformer_inv.inverseFourierTransform(fftVfilter, filteredVol);
	#ifdef DEBUG
	Image<double> filteredvolume;
	filteredvolume() = filteredVol;
	filteredvolume.write(formatString("Volumen_filtrado_%i.vol", count));
	filteredvolume() = testVol;
	filteredvolume.write(formatString("filtro_%i.vol", count));
	#endif
}

void ProgLocSharpening::resVol2IdxVol(const MultidimArray<double> &resVol, double &maxres, MultidimArray<int> &idxVol,
		std::vector<int> &idxList)
{
	idxVol.initZeros(resVol);
	size_t xdim, ydim, zdim,Ndim;
	idxVol.getDimensions(xdim, ydim, zdim,Ndim);

	size_t N;
	N = xdim*ydim*zdim;

	MultidimArray<int> allIdx(N);

	long n=0;
	double resolution, freq;
	int idx;
	size_t sizeVol=ZSIZE(resVol);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resVol)
	{
		resolution=DIRECT_MULTIDIM_ELEM(resVol, n);
		//std::cout << "resolution = " << resolution << std::endl;
		if (resolution<2*sampling)
		{
//			std::cout << "resolution = " << resolution << std::endl;
			DIRECT_MULTIDIM_ELEM(idxVol, n) = -1;
			DIGFREQ2FFT_IDX(sampling/maxRes, sizeVol, idx);
			DIRECT_MULTIDIM_ELEM(allIdx, n) = idx;

		}
		else{
//			resolution = maxres;
		#ifdef DEBUG
		if (resolution>0)
			std::cout << "resolution = " << resolution << std::endl;
		#endif
		freq=sampling/(resolution);

		DIGFREQ2FFT_IDX(freq, sizeVol, idx);
		DIRECT_MULTIDIM_ELEM(idxVol, n) = idx;

		DIRECT_MULTIDIM_ELEM(allIdx, n) = idx;
		}
	}

	std::sort(&A1D_ELEM(allIdx,0),&A1D_ELEM(allIdx, N));

	int last_idx = -1;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(allIdx)
	{
		idx = DIRECT_MULTIDIM_ELEM(allIdx, n);
		if (idx>last_idx){
			last_idx = idx;
			idxList.push_back(idx);
//			std::cout << "idx = " << idx << std::endl;
		}
	}
}


/*
void ProgLocSharpening::resolution2eval(int &count_res, double step,
								double &resolution, double &last_resolution,
								double &freq, double &freqL,
								int &last_fourier_idx,
								bool &continueIter,	bool &breakIter,
								bool &doNextIteration)
{
	resolution = maxRes - count_res*step;
	freq = sampling/resolution;
	++count_res;

	double Nyquist = 2*sampling;
	double aux_frequency;
	int fourier_idx;

	DIGFREQ2FFT_IDX(freq, ZSIZE(VRiesz), fourier_idx);

//	std::cout << "Resolution = " << resolution << "   iter = " << count_res-1 << std::endl;
//	std::cout << "freq = " << freq << "   Fourier index = " << fourier_idx << std::endl;

	FFT_IDX2DIGFREQ(fourier_idx, ZSIZE(VRiesz), aux_frequency);

	freq = aux_frequency;

	if (fourier_idx == last_fourier_idx)
	{
//		std::cout << "entro en el if"  << std::endl;
		continueIter = true;
		return;
	}

	last_fourier_idx = fourier_idx;
	resolution = sampling/aux_frequency;


	if (count_res == 0)
		last_resolution = resolution;

	if ( ( resolution<Nyquist ))// || (resolution > last_resolution) )
	{
		//std::cout << "Nyquist limit reached" << std::endl;
		breakIter = true;
		return;
	}
	if ( ( freq>0.49))// || (resolution > last_resolution) )
	{
		std::cout << "Nyquist limit reached" << std::endl;
		doNextIteration = false;
		return;
	}

	freqL = sampling/(resolution + step);

	int fourier_idx_2;

	DIGFREQ2FFT_IDX(freqL, ZSIZE(VRiesz), fourier_idx_2);
	double caca, caca2;

	if (fourier_idx_2 == fourier_idx)
	{
		if (fourier_idx > 0){
			//std::cout << " index low =  " << (fourier_idx - 1) << std::endl;
			FFT_IDX2DIGFREQ(fourier_idx - 1, ZSIZE(VRiesz), freqL);
		}
		else{
			freqL = sampling/(resolution + step);
		}
	}
}
*/

void ProgLocSharpening::run()
{
	produceSideInfo();

	double lowestFreq = sampling/maxRes;
	double highestFreq = sampling/minRes;
	int lowestIdx;
	MultidimArray<double> filteredVol, Vrest, lastVrest, Vlocalfilt, lastVrest_aux;
	Vlocalfilt.initZeros(Vorig);

	std::cout << "reoslution list = "<< idxList[0] << std::endl;

	DIGFREQ2FFT_IDX(lowestFreq, ZSIZE(fftV), lowestIdx);

	std::cout << "lowestIdx = " << lowestIdx << "  lowestFreq= "  << lowestFreq << std::endl;

	double freq;
	int idx = lowestIdx-1;
	bool Nextiter = true;


	size_t size_list;
	size_list = idxList.size();

	lastVrest = Vorig;

    #ifdef DEBUG_FILTER
	lowPassFilterFunction(fftV, 0.1214, 0.1414, filteredVol, idx);
	Image<double> qqqq;
	qqqq() = filteredVol;
	qqqq.write(formatString("localfilter.vol"));
    #endif


    for (size_t i = 0; i<Niter; ++i)
	{
    	std::cout << "----------------Iteration " << i << "----------------" << std::endl;
    		Vlocalfilt=Vorig;
			for (size_t k = 0; k<size_list; ++k)
			{
				idx = idxList[k];
				FFT_IDX2DIGFREQ(idx, ZSIZE(fftV), freq);

				if (freq<0)
					Nextiter = false;
				else
				{
					if (freq>highestFreq)
						break;

					double wL = freq+0.02;
					std::cout << "i=" << idx << " freq = " << sampling/freq << std::endl;
					lowPassFilterFunction(fftV, freq, wL, filteredVol, idx);

					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(filteredVol)
					{
						if (DIRECT_MULTIDIM_ELEM(idxVol, n) == idx)
							DIRECT_MULTIDIM_ELEM(Vlocalfilt, n) = DIRECT_MULTIDIM_ELEM(filteredVol, n);
					}
				}
				if (Nextiter == false)
					break;
			}


			///////
			/*
			Vlocalfilt = Vorig-Vlocalfilt;
			FFT_IDX2DIGFREQ(idx, ZSIZE(fftV), freq);
			double wL = freq+0.02;

			FourierTransformer transformer;
			MultidimArray<double> &inputVol = Vlocalfilt;
			transformer.FourierTransform(inputVol, fftV);

			lowPassFilterFunction(fftV, idxList[0], wL, filteredVol, idx);

			sharpenedMap=lastVrest+lambda*(filteredVol);
			*/
			////////////
			sharpenedMap=lastVrest+lambda*(Vorig-Vlocalfilt);
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sharpenedMap)
			{
				if (DIRECT_MULTIDIM_ELEM(sharpenedMap, n) < 0)
					DIRECT_MULTIDIM_ELEM(sharpenedMap, n) = 0;
			}

			lastVrest = sharpenedMap;

			Image<double> filteredvolume0;
			filteredvolume0() = Vlocalfilt;
			filteredvolume0.write(formatString("localfilter_%i.vol",i));
			filteredvolume0() = sharpenedMap;
			filteredvolume0.write(formatString("sharpenedMap_%i.vol", i));

			lastVrest_aux = lastVrest;
			FourierTransformer transformer;
			transformer.FourierTransform(lastVrest_aux, fftV);

	}
		Image<double> filteredvolume;
		filteredvolume() = sharpenedMap;
		filteredvolume.write(formatString("sharpenedMap.vol"));


}
