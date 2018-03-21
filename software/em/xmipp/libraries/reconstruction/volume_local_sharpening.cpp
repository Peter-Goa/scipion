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

    Image<double> resolutionVolume;
    resolutionVolume.read(fnRes);

    resVol = resolutionVolume();

	maxMinResolution(resVol, maxRes, minRes);
	std::cout << "maxRes = " << maxRes << "  minRes = " << minRes << std::endl;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resVol)
	{
		if (DIRECT_MULTIDIM_ELEM(resVol, n) < 2*sampling)
			DIRECT_MULTIDIM_ELEM(resVol, n) = 2*sampling;
	}

	resVol.setXmippOrigin();

	FourierFilter Filter;
	Filter.FilterShape=REALGAUSSIAN;
	Filter.FilterBand=LOWPASS;
	Filter.w1=1;
	Filter.apply(resVol);
//
//	Image<double> filteredvolume;
//	filteredvolume() = Vorig;
//	filteredvolume.write(fnOut);
//
//exit(0);
//	VsoftMask.initZeros(resVol);
//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resVol)
//		if (DIRECT_MULTIDIM_ELEM(resVol,n)>=2*sampling)
//			DIRECT_MULTIDIM_ELEM(VsoftMask,n)=1;

//    Filter.w1=1;
//    VsoftMask.setXmippOrigin();
//    Filter.apply(VsoftMask);
//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(VsoftMask)
//		if (DIRECT_MULTIDIM_ELEM(VsoftMask,n)<0)
//			DIRECT_MULTIDIM_ELEM(VsoftMask,n)=0;
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

void ProgLocSharpening::bandPassFilterFunction(const MultidimArray< std::complex<double> > &myfftV,
		double w, double step, MultidimArray<double> &filteredVol, int count)
{
	fftVfilter.initZeros(myfftV);
	size_t xdim, ydim, zdim, ndim;
	//amplitude.resizeNoCopy(filteredVol);
	Vorig.getDimensions(xdim, ydim, zdim, ndim);
	MultidimArray<double> testVol;
	testVol.initZeros(myfftV);
	filteredVol.resizeNoCopy(Vorig);
	//resizeNoCopy(sharpenedMap);

	//double delta = wL-w;
	//double w_inf = w-delta;
	double w_inf = w;
	double w_sup = w+step;
	// Filter the input volume and add it to amplitude
	long n=0;
//	double ideltal=PI/(delta);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(myfftV)
	{
		double un=DIRECT_MULTIDIM_ELEM(iu,n);

		if (un<w_sup && un>=w_inf)
		{
			DIRECT_MULTIDIM_ELEM(fftVfilter, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
            DIRECT_MULTIDIM_ELEM(testVol, n) = 1.0;
		}
//		if (un>=w && un<=wL)
//		{
//			//double H=0.5*(1+cos((un-w1)*ideltal));
//			DIRECT_MULTIDIM_ELEM(fftVfilter, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
//			DIRECT_MULTIDIM_ELEM(fftVfilter, n) *= 0.5*(1+cos((un-w)*ideltal));//H;
//			//DIRECT_MULTIDIM_ELEM(testVol, n) = 1.0;
//		} else if (un<=w && un>=w_inf)
//		{
//			//double H=0.5*(1+cos((un-w1)*ideltal));
//			DIRECT_MULTIDIM_ELEM(fftVfilter, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
//			DIRECT_MULTIDIM_ELEM(fftVfilter, n) *= 0.5*(1+cos((un-w)*ideltal));//H;
//			//DIRECT_MULTIDIM_ELEM(testVol, n) = 1.0;
//		}
	}

	transformer_inv.inverseFourierTransform(fftVfilter, filteredVol);
	#ifdef DEBUG_FILTER
	Image<double> filteredvolume;
	filteredvolume() = filteredVol;
	filteredvolume.write(formatString("Volumen_filtrado_%i.vol", count));
	filteredvolume() = testVol;
	filteredvolume.write(formatString("filtro_%i.vol", count));
	#endif
}


void ProgLocSharpening::localfiltering(MultidimArray< std::complex<double> > &myfftV,
										MultidimArray<double> &localfilteredVol,
										double &minFreq, double &maxFreq, double &step)
{
	MultidimArray<double> lastfreqVol, filteredVol, lastweight, weight;
	localfilteredVol.initZeros(Vorig);
	weight.initZeros(Vorig);
	lastweight.initZeros(Vorig);

	double lastResolution=1e38;
	int idx=1, lastidx = -1;

	//First band
	double freq0=0.0;
	double step0 = minFreq;
	bandPassFilterFunction(myfftV, freq0, step0, filteredVol, idx);
	localfilteredVol = filteredVol;

	for (double freq = minFreq; freq<=maxFreq; freq+=step)
	{
		double res = sampling/freq;
		double freqstep=sampling/(res-0.2);
               step=freqstep-freq;
//		DIGFREQ2FFT_IDX(freq, ZSIZE(myfftV), idx);
//
//		if (idx == lastidx)
//			continue;
//		if ((fabs(lastResolution-res))<step)
//			continue;

		//double wL = sampling/(res - step);
		//double wL = freq+0.001;
		std::cout << " freqmin = " << freq <<  " freqmax = " << freq+step << "  resol= " << res << std::endl;

		bandPassFilterFunction(myfftV, freq, step, filteredVol, idx);

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(filteredVol)
		{
			double resOk = DIRECT_MULTIDIM_ELEM(resVol, n)+1e-38;
			//double freq_map = sampling/res;

			if (res<resOk)
				DIRECT_MULTIDIM_ELEM(filteredVol, n)=0;

			if (DIRECT_MULTIDIM_ELEM(filteredVol, n)<0)
				DIRECT_MULTIDIM_ELEM(filteredVol, n)=0;


//			DIRECT_MULTIDIM_ELEM(weight, n) = (exp(-0.5*(freq-freq_map)*(freq-freq_map)));
//			DIRECT_MULTIDIM_ELEM(filteredVol, n) *= DIRECT_MULTIDIM_ELEM(weight, n);
		}

		localfilteredVol += filteredVol;
//		lastweight += weight;
//		lastResolution = res;
//		lastidx = idx;
		++idx;
	}

//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(filteredVol)
//	{
//		DIRECT_MULTIDIM_ELEM(localfilteredVol, n) = DIRECT_MULTIDIM_ELEM(localfilteredVol, n)/DIRECT_MULTIDIM_ELEM(lastweight, n);
//	}
}

void ProgLocSharpening::run()
{
	produceSideInfo();

	MultidimArray<double> auxVol;


	double freq;
	int idx;

	double step = 0.01;
	double lastResolution=1e38;

	MultidimArray<double> operatedfiltered, Vk, filteredVol;

	minRes = 2*sampling;
    maxFreq=0.5;
    minFreq=sampling/maxRes;

	std::cout << "Resolutions between " << minRes << " and " << maxRes << std::endl;

	int lastidx = -1;

	double aux_maxres = 20.0;
//	localfiltering(fftV, filteredVol, minRes, maxRes, step);
//
//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(filteredVol)
//	{
//		if (DIRECT_MULTIDIM_ELEM(resVol, n) <= 2*sampling)
//			DIRECT_MULTIDIM_ELEM(filteredVol, n) = 0;
////		if (DIRECT_MULTIDIM_ELEM(filteredVol, n) <0)
////			DIRECT_MULTIDIM_ELEM(filteredVol, n) = 0;
//	}
//
//
//
//
//
//	Image<double> save;
//	save() = filteredVol;
//	save.write("filtrado.vol");
//
//
//	exit(0);

	FourierTransformer transformer;

	//Vorig = filteredVol;

	filteredVol = Vorig;
    for (size_t i = 0; i<Niter; ++i)
	{
    	std::cout << "----------------Iteration " << i << "----------------" << std::endl;
    	auxVol = filteredVol;
    	transformer.FourierTransform(auxVol, fftV);

    	localfiltering(fftV, operatedfiltered, minFreq, maxFreq, step);

//    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(operatedfiltered)
//		{
//    		if (DIRECT_MULTIDIM_ELEM(operatedfiltered, n) <0)
//    			DIRECT_MULTIDIM_ELEM(operatedfiltered, n) = 0;
//		}

		Image<double> save;
		save() = operatedfiltered;
		save.write(formatString("sharpenedMapA_%i.vol", i));

		filteredVol = Vorig-operatedfiltered;

		save() = filteredVol;
		save.write(formatString("sharpenedMapB_%i.vol", i));

		////Second operator
    	transformer.FourierTransform(filteredVol, fftV);
    	localfiltering(fftV, filteredVol, minFreq, maxFreq, step);

		save() = filteredVol;
		save.write(formatString("sharpenedMapC_%i.vol", i));

		if (i == 0)
			Vk = Vorig;
		else
			Vk = sharpenedMap;

		sharpenedMap=Vk+lambda*(filteredVol);

//		Image<double> save;
//		save() = sharpenedMap;
//		save.write(fnOut);



		Image<double> filteredvolume0;
		filteredvolume0() = sharpenedMap;
		filteredvolume0.write(formatString("sharpenedMapD_%i.vol", i));

		filteredVol = sharpenedMap;
	}
	Image<double> filteredvolume;
	filteredvolume() = sharpenedMap;
	filteredvolume.write(fnOut);

		//*/
		////////////

//
////			Vlocalfilt=(1-VsoftMask)*Vorig+VsoftMask*Vlocalfilt;
//
//			Vdiff=Vorig-Vlocalfilt;
//
//		    FourierFilter Filter;
//		    Filter.FilterShape=REALGAUSSIAN;
//		    Filter.FilterBand=LOWPASS;
//		    Filter.w1=1;
//		    Vdiff.setXmippOrigin();
//		    Filter.apply(Vdiff);
//
//			sharpenedMap=(1-VsoftMask)*Vorig+VsoftMask*(lastVrest+lambda*(Vdiff));
//
//
//
//			Image<double> kk;
//			kk() = Vdiff;
//			kk.write(formatString("diff_%i.vol",i));
//
//
//
////			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sharpenedMap)
////			{
////				if ((DIRECT_MULTIDIM_ELEM(sharpenedMap, n) < 0) || (DIRECT_MULTIDIM_ELEM(idxVol, n)<0))
////					DIRECT_MULTIDIM_ELEM(sharpenedMap, n) = 0;
////			}
////
////			estimateS(sharpenedMap);
////			significanceRealSpace(sharpenedMap, Vorig, sharpenedMap, i);
//
//

//
//			Image<double> filteredvolume0;
//			filteredvolume0() = Vlocalfilt;
//			filteredvolume0.write(formatString("localfilter_%i.vol",i));
//			filteredvolume0() = sharpenedMap;
//			filteredvolume0.write(formatString("sharpenedMap_%i.vol", i));
//
//			lastVrest_aux = lastVrest;



}
