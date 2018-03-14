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




//    Image<double> save;
//    save()=resVol;
//    save.write("PPP.vol");


	FourierTransformer transformer;
	MultidimArray<double> &inputVol = V();

//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(inputVol)
//	{
//		if (DIRECT_MULTIDIM_ELEM(inputVol, n) < 0)
//			DIRECT_MULTIDIM_ELEM(inputVol, n) = 0;
//	}

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

    Image<double> resolutionVolume;
    resolutionVolume.read(fnRes);

    resVol = resolutionVolume();

	maxMinResolution(resVol, maxRes, minRes);
	std::cout << "maxRes = " << maxRes << "  minRes = " << minRes << std::endl;

//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resVol)
//	{
//		if (DIRECT_MULTIDIM_ELEM(resVol, n) < 0)
//			DIRECT_MULTIDIM_ELEM(resVol, n) = 0;
//		if (DIRECT_MULTIDIM_ELEM(resVol, n) < 2*sampling)
//			DIRECT_MULTIDIM_ELEM(resVol, n) = 2*sampling;
//	}

	resVol.setXmippOrigin();

	FourierFilter Filter;
	Filter.FilterShape=REALGAUSSIAN;
	Filter.FilterBand=LOWPASS;
	Filter.w1=1;
	Filter.apply(resVol);




	VsoftMask.initZeros(resVol);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resVol)
		if (DIRECT_MULTIDIM_ELEM(resVol,n)>=2*sampling)
			DIRECT_MULTIDIM_ELEM(VsoftMask,n)=1;

    Filter.w1=1;
    VsoftMask.setXmippOrigin();
    Filter.apply(VsoftMask);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(VsoftMask)
		if (DIRECT_MULTIDIM_ELEM(VsoftMask,n)<0)
			DIRECT_MULTIDIM_ELEM(VsoftMask,n)=0;

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


//void ProgLocSharpening::lowPassFilterFunction(const MultidimArray< std::complex<double> > &myfftV,
//		double w, double wL, MultidimArray<double> &filteredVol, int count)
//{
//	fftVfilter.initZeros(myfftV);
//	size_t xdim, ydim, zdim, ndim;
//	//amplitude.resizeNoCopy(filteredVol);
//	Vorig.getDimensions(xdim, ydim, zdim, ndim);
//	MultidimArray<double> testVol;
//	testVol.initZeros(myfftV);
//	filteredVol.resizeNoCopy(Vorig);
//	//resizeNoCopy(sharpenedMap);
//
//	// Filter the input volume and add it to amplitude
//	long n=0;
//	double ideltal=PI/(wL-w);
//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(myfftV)
//	{
//		double un=DIRECT_MULTIDIM_ELEM(iu,n);
//		if (un>=w && un<wL)
//		{
//			//double H=0.5*(1+cos((un-w1)*ideltal));
//			DIRECT_MULTIDIM_ELEM(fftVfilter, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
//			DIRECT_MULTIDIM_ELEM(fftVfilter, n) *= 0.5*(1+cos((un-w)*ideltal));//H;
////			DIRECT_MULTIDIM_ELEM(testVol, n) = 0.5;
//		} else if (un<w)
//		{
//			DIRECT_MULTIDIM_ELEM(fftVfilter, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
////			DIRECT_MULTIDIM_ELEM(testVol, n) = 1;
//		}
//	}
//
//	transformer_inv.inverseFourierTransform(fftVfilter, filteredVol);
//	#ifdef DEBUG
//	Image<double> filteredvolume;
//	filteredvolume() = filteredVol;
//	filteredvolume.write(formatString("Volumen_filtrado_%i.vol", count));
//	filteredvolume() = testVol;
//	filteredvolume.write(formatString("filtro_%i.vol", count));
//	#endif
//}


void ProgLocSharpening::bandPassFilterFunction(const MultidimArray< std::complex<double> > &myfftV,
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

	double delta = wL-w;
	double w_inf = w-delta;
	// Filter the input volume and add it to amplitude
	long n=0;
	double ideltal=PI/(delta);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(myfftV)
	{
		double un=DIRECT_MULTIDIM_ELEM(iu,n);
		if (un>=w && un<=wL)
		{
			//double H=0.5*(1+cos((un-w1)*ideltal));
			DIRECT_MULTIDIM_ELEM(fftVfilter, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
			DIRECT_MULTIDIM_ELEM(fftVfilter, n) *= 0.5*(1+cos((un-w)*ideltal));//H;
//			DIRECT_MULTIDIM_ELEM(testVol, n) = 0.5;
		} else if (un<=w && un>=w_inf)
		{
			//double H=0.5*(1+cos((un-w1)*ideltal));
			DIRECT_MULTIDIM_ELEM(fftVfilter, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
			DIRECT_MULTIDIM_ELEM(fftVfilter, n) *= 0.5*(1+cos((un-w)*ideltal));//H;
//			DIRECT_MULTIDIM_ELEM(testVol, n) = 0.5;
		}
	}

	transformer_inv.inverseFourierTransform(fftVfilter, filteredVol);
	#ifdef DEBUG_FILTER
	Image<double> filteredvolume;
	filteredvolume() = filteredVol;
	filteredvolume.write(formatString("Volumen_filtrado_%i.vol", count));
//	filteredvolume() = testVol;
//	filteredvolume.write(formatString("filtro_%i.vol", count));
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

void ProgLocSharpening::estimateS(MultidimArray <double> vol)
{
	MultidimArray<double> &mS=vol;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mS)
	if (DIRECT_MULTIDIM_ELEM(idxVol,n)<0)
		DIRECT_MULTIDIM_ELEM(mS,n)=0;

	// Apply positivity
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mS)
	if (DIRECT_MULTIDIM_ELEM(mS,n)<0)
		DIRECT_MULTIDIM_ELEM(mS,n)=0;

	// Filter S
	FourierTransformer transformer;
	MultidimArray < std::complex <double> > fVol;
	transformer.FourierTransform(mS,fVol,false);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fVol)
		if (DIRECT_MULTIDIM_ELEM(iu,n)>0.25)
			DIRECT_MULTIDIM_ELEM(fVol,n)=0.0;
	transformer.inverseFourierTransform();

	// Calculate HS CDF in real space
	long N = 0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(idxVol)
	{
	if (DIRECT_MULTIDIM_ELEM(idxVol,n)>0)
		++N;
	}
	MultidimArray<double> aux;
	aux.resizeNoCopy(N);
	size_t idx=0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mS)
	if (DIRECT_MULTIDIM_ELEM(idxVol,n)>0)
		DIRECT_MULTIDIM_ELEM(aux,idx++)=DIRECT_MULTIDIM_ELEM(mS,n)*DIRECT_MULTIDIM_ELEM(mS,n);

	cdfS.calculateCDF(aux);
}


void ProgLocSharpening::significanceRealSpace(const MultidimArray<double> &Vi, MultidimArray<double> Vol2 ,MultidimArray<double> &Vir, int count)
{
	// Calculate N=Vi-H*S and its energy CDF

	const MultidimArray<double> &mS=Vol2;
	long N = 0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(idxVol)
		{
		if (DIRECT_MULTIDIM_ELEM(idxVol,n)>0)
			++N;
		}
	MultidimArray<double> mN(N);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mN)
	{
		if (DIRECT_MULTIDIM_ELEM(idxVol,n)>0)
		{
		double diff=DIRECT_MULTIDIM_ELEM(Vi,n)-DIRECT_MULTIDIM_ELEM(mS,n);
		DIRECT_MULTIDIM_ELEM(mN,n)=diff*diff;
		}
	}
	CDF cdfN;
	cdfN.calculateCDF(mN);

	MultidimArray <double> mask_smoothing;
	mask_smoothing.initZeros(idxVol);
	// Mask the input volume with a mask value that is proportional to the probability of being larger than noise
	Vir.resizeNoCopy(Vi);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mS)
	{
		double e=DIRECT_MULTIDIM_ELEM(Vi,n)*DIRECT_MULTIDIM_ELEM(Vi,n);
		double pN=cdfN.getProbability(e);
		if (pN<1)
		{
			double pS=cdfS.getProbability(e);
			double pp=pS*pN;
			DIRECT_MULTIDIM_ELEM(Vir,n)=pp*DIRECT_MULTIDIM_ELEM(Vi,n);
			DIRECT_MULTIDIM_ELEM(mask_smoothing,n) = pp;
		}
		else
		{
			DIRECT_MULTIDIM_ELEM(Vir,n)=DIRECT_MULTIDIM_ELEM(Vi,n);
			DIRECT_MULTIDIM_ELEM(mask_smoothing,n) = 1.0;
		}
	}
	Image<double> filteredvolume;
	filteredvolume() = mask_smoothing;
	filteredvolume.write(formatString("mask_%i.vol", count));
}


void ProgLocSharpening::run()
{
	produceSideInfo();

	double lowestFreq = sampling/maxRes;
	double highestFreq = sampling/minRes;
	int lowestIdx;
	MultidimArray<double> filteredVol, Vrest, lastVrest, Vlocalfilt, lastVrest_aux;
	Vlocalfilt.initZeros(Vorig);
	lastVrest=Vorig;

	std::cout << "reoslution list = "<< idxList[0] << std::endl;

	DIGFREQ2FFT_IDX(lowestFreq, ZSIZE(fftV), lowestIdx);

	std::cout << "lowestIdx = " << lowestIdx << "  lowestFreq= "  << lowestFreq << std::endl;

	double freq;
	int idx = lowestIdx-1;

	double step = 0.2;
	double lastResolution=1e38;

	size_t size_list;
	size_list = idxList.size();

	MultidimArray<double> Vdiff, lastfreqVol, lastweight, weight;
	Vdiff.initZeros(Vorig);
	lastfreqVol.initZeros(Vorig);
	weight.initZeros(Vorig);
	lastweight.initZeros(Vorig);

	minRes = 2*sampling;

	int lastidx = -1;

    for (size_t i = 0; i<Niter; ++i)
	{
    	lastfreqVol.initZeros(Vorig);
    	std::cout << "----------------Iteration " << i << "----------------" << std::endl;
//    		Vlocalfilt=Vorig;
    		Vlocalfilt.initZeros(Vorig);
			for (double res = minRes; res<maxRes; res+=step)
			{
//				idx = idxList[k];
//				FFT_IDX2DIGFREQ(idx, ZSIZE(fftV), freq);

				double freq = sampling/res;

				DIGFREQ2FFT_IDX(freq, ZSIZE(fftV), idx);

				if (idx == lastidx)
					continue;

				if ((fabs(lastResolution-res))<step)
				{
					continue;
				}

				double wL = freq+0.02;
				std::cout << " res = " << res << std::endl;
				bandPassFilterFunction(fftV, freq, wL, filteredVol, idx);

				std::cout << "freq  = " << freq << std::endl;

				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(filteredVol)
				{
					double res = DIRECT_MULTIDIM_ELEM(resVol, n)+1e-38;
					double freq_map = sampling/res;


					DIRECT_MULTIDIM_ELEM(weight, n) = exp(-0.5*(freq-freq_map)*(freq-freq_map));
					DIRECT_MULTIDIM_ELEM(filteredVol, n) *= DIRECT_MULTIDIM_ELEM(weight, n);

//					if (DIRECT_MULTIDIM_ELEM(resVol, n)<2*sampling)
//					{
//						std::cout << "res_resVol  = " << DIRECT_MULTIDIM_ELEM(resVol, n) << std::endl;
//						std::cout << "res = " << res << std::endl;
//						std::cout << "freq_map  = " << freq_map << std::endl;
//						std::cout << "weight  = " << DIRECT_MULTIDIM_ELEM(weight, n) << std::endl;
//						//exit(0);
//					}
				}


				lastfreqVol += filteredVol;
				lastweight += weight;
				lastResolution = res;
				lastidx = idx;
			}
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(filteredVol)
			{
				DIRECT_MULTIDIM_ELEM(lastfreqVol, n) = DIRECT_MULTIDIM_ELEM(lastfreqVol, n)/DIRECT_MULTIDIM_ELEM(lastweight, n);
			}

			Image<double> save;
			save() = lastfreqVol;
			save.write(fnOut);

			exit(0);
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


//			Vlocalfilt=(1-VsoftMask)*Vorig+VsoftMask*Vlocalfilt;

			Vdiff=Vorig-Vlocalfilt;

		    FourierFilter Filter;
		    Filter.FilterShape=REALGAUSSIAN;
		    Filter.FilterBand=LOWPASS;
		    Filter.w1=1;
		    Vdiff.setXmippOrigin();
		    Filter.apply(Vdiff);

			sharpenedMap=(1-VsoftMask)*Vorig+VsoftMask*(lastVrest+lambda*(Vdiff));



			Image<double> kk;
			kk() = Vdiff;
			kk.write(formatString("diff_%i.vol",i));



//			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sharpenedMap)
//			{
//				if ((DIRECT_MULTIDIM_ELEM(sharpenedMap, n) < 0) || (DIRECT_MULTIDIM_ELEM(idxVol, n)<0))
//					DIRECT_MULTIDIM_ELEM(sharpenedMap, n) = 0;
//			}
//
//			estimateS(sharpenedMap);
//			significanceRealSpace(sharpenedMap, Vorig, sharpenedMap, i);


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
