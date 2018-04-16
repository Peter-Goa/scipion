/***************************************************************************
 *
 * Authors:    Erney Ramirez-Aportela,                                    eramirez@cnb.csic.es
 *             Jose Luis Vilas,                                           jlvilas@cnb.csic.es
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
                //if (value<lastMinRes)
                if (value<lastMinRes && value>0)
                        lastMinRes = value;
        }

        maxRes = lastMaxRes;
        minRes = lastMinRes;
}

void ProgLocSharpening::bandPassFilterFunction(const MultidimArray< std::complex<double> > &myfftV,
                double w, double wL, MultidimArray<double> &filteredVol, int count)
{
        fftVfilter.initZeros(myfftV);
        size_t xdim, ydim, zdim, ndim;
        Vorig.getDimensions(xdim, ydim, zdim, ndim);
        filteredVol.resizeNoCopy(Vorig);

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
                        DIRECT_MULTIDIM_ELEM(fftVfilter, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
                        DIRECT_MULTIDIM_ELEM(fftVfilter, n) *= 0.5*(1+cos((un-w)*ideltal));//H;
                } else if (un<=w && un>=w_inf)
                {
                        DIRECT_MULTIDIM_ELEM(fftVfilter, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
                        DIRECT_MULTIDIM_ELEM(fftVfilter, n) *= 0.5*(1+cos((un-w)*ideltal));//H;
                }
        }

        transformer_inv.inverseFourierTransform(fftVfilter, filteredVol);

        #ifdef DEBUG_FILTER
        Image<double> filteredvolume;
        filteredvolume() = filteredVol;
        filteredvolume.write(formatString("Volumen_filtrado_%i.vol", count));
        #endif
}

void ProgLocSharpening::localfiltering(MultidimArray< std::complex<double> > &myfftV,
                                       MultidimArray<double> &localfilteredVol,
                                       double &minRes, double &maxRes, double &step)
{
        MultidimArray<double> lastfreqVol, filteredVol, lastweight, weight;
        localfilteredVol.initZeros(Vorig);
        weight.initZeros(Vorig);
        lastweight.initZeros(Vorig);

        double freq, lastResolution=1e38;
        int idx, lastidx = -1;

        for (double res = minRes; res<maxRes; res+=step)
        {
                freq = sampling/res;

                DIGFREQ2FFT_IDX(freq, ZSIZE(myfftV), idx);

                if (idx == lastidx)
                {
//                    std::cout << "idx = " << idx << std::endl;
                        continue;
                }
//                else
//                	 std::cout << "Esta idx si entra = " << idx << std::endl;

//                if ((fabs(lastResolution-res))<step)
//                {
//                	    std::cout << "res = " << res << "  sampling/step = " << sampling/step << std::endl;
//                	    std::cout << "lastResolution = " << lastResolution << "  res = " << res << "  step = " << step << std::endl;
//                        continue;
//                }

                double wL = sampling/(res - step);
                //double wL = freq+0.02;
                //std::cout << " resolution = " << res << "  freq = " << freq << std::endl;

                bandPassFilterFunction(myfftV, freq, wL, filteredVol, idx);

                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(filteredVol)
                {
                        double res = DIRECT_MULTIDIM_ELEM(resVol, n)+1e-38;
                        double freq_map = sampling/res;

                        DIRECT_MULTIDIM_ELEM(weight, n) = (exp(-0.5*(freq-freq_map)*(freq-freq_map)));
                        DIRECT_MULTIDIM_ELEM(filteredVol, n) *= DIRECT_MULTIDIM_ELEM(weight, n);
                }

                localfilteredVol += filteredVol;
                lastweight += weight;
                lastResolution = res;
                lastidx = idx;
        }

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(filteredVol)
        {
         DIRECT_MULTIDIM_ELEM(localfilteredVol, n) = DIRECT_MULTIDIM_ELEM(localfilteredVol, n)/DIRECT_MULTIDIM_ELEM(lastweight, n);
        }
}

void ProgLocSharpening::run()
{
        produceSideInfo();

        MultidimArray<double> auxVol;
        double lastnorm=1;
        double freq;
        int idx;

        double step = 0.2;
        double lastResolution=1e38;

        MultidimArray<double> operatedfiltered, Vk, filteredVol;

        minRes = 2*sampling;
        //maxRes = 20;//Esto solo para este caso

        std::cout << "Resolutions between " << minRes << " and " << maxRes << std::endl;

        int lastidx = -1;

        FourierTransformer transformer;

        filteredVol = Vorig;
      	sharpenedMap.resizeNoCopy(Vorig);

    for (size_t i = 0; i<Niter; ++i)
        {
        std::cout << "----------------Iteration " << i << "----------------" << std::endl;
        auxVol = filteredVol;
        transformer.FourierTransform(auxVol, fftV);

        localfiltering(fftV, operatedfiltered, minRes, maxRes, step);

                Image<double> save;
//                save() = operatedfiltered;
//                save.write(formatString("sharpenedMapA_%i.vol", i));

                filteredVol = Vorig;
        		filteredVol -= operatedfiltered;

        		//calculate norm^2 del Vorig
        		if (i==0)
        		{
            		double norm=0;
            		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Vorig)
            		{
            			 norm +=(DIRECT_MULTIDIM_ELEM(Vorig,n)*DIRECT_MULTIDIM_ELEM(Vorig,n));
            		}
            		norm=sqrt(norm);
                    std::cout << "norma del original  " << norm << std::endl;
        		}



        		//calculate norm^2
        		double norm=0;
        		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(operatedfiltered)
        		{
        			 norm +=(DIRECT_MULTIDIM_ELEM(operatedfiltered,n)*DIRECT_MULTIDIM_ELEM(operatedfiltered,n));
        		}
        		norm=sqrt(norm);

                double porc=lastnorm*100/norm;
                std::cout << "norma " << norm << " porciento " << porc << std::endl;
                lastnorm=norm;



//                save() = filteredVol;
//                save.write(formatString("sharpenedMapB_%i.vol", i));

                ////Second operator
        transformer.FourierTransform(filteredVol, fftV);
        localfiltering(fftV, filteredVol, minRes, maxRes, step);

//                save() = filteredVol;
//                save.write(formatString("sharpenedMapC_%i.vol", i));

                if (i == 0)
                        Vk = Vorig;
                else
                        Vk = sharpenedMap;

                //sharpenedMap=Vk+lambda*(filteredVol);
        		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sharpenedMap)
                {
        			DIRECT_MULTIDIM_ELEM(sharpenedMap,n)=DIRECT_MULTIDIM_ELEM(Vk,n)+
        			                     lambda*DIRECT_MULTIDIM_ELEM(filteredVol,n);
                }

                Image<double> filteredvolume0;
                filteredvolume0() = sharpenedMap;
                filteredvolume0.write(formatString("sharpenedMapD_%i.vol", i));

                filteredVol = sharpenedMap;
        }
        Image<double> filteredvolume;
        filteredvolume() = sharpenedMap;
        filteredvolume.write(fnOut);


}





























//
//
///***************************************************************************
// *
// * Authors:    Erney Ramirez , 				     	  eramirez@cnb.csic.es
// * 			   Jose Luis Vilas                        jlvilas@cnb.csic.es
// *
// * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
// *
// * This program is free software; you can redistribute it and/or modify
// * it under the terms of the GNU General Public License as published by
// * the Free Software Foundation; either version 2 of the License, or
// * (at your option) any later version.
// *
// * This program is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// * GNU General Public License for more details.
// *
// * You should have received a copy of the GNU General Public License
// * along with this program; if not, write to the Free Software
// * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
// * 02111-1307  USA
// *
// *  All comments concerning this program package may be sent to the
// *  e-mail address 'xmipp@cnb.csic.es'
// ***************************************************************************/
//
//#include "volume_local_sharpening.h"
////#define DEBUG
////#define DEBUG_MASK
////#define DEBUG_FILTER
//void ProgLocSharpening::readParams()
//{
//	fnVol = getParam("--vol");
//	fnRes = getParam("--resolution_map");
//	sampling = getDoubleParam("--sampling");
//	lambda = getDoubleParam("-l");
//	Niter = getIntParam("-n");
//	fnOut = getParam("-o");
//}
//
//void ProgLocSharpening::defineParams()
//{
//	addUsageLine("This function performs local sharpening");
//	addParamsLine("  --vol <vol_file=\"\">   : Input volume");
//	addParamsLine("  --resolution_map <vol_file=\"\">: Resolution map");
//	addParamsLine("  -o <output=\"Sharpening.vol\">: sharpening volume");
//	addParamsLine("  --sampling <s=1>: sampling");
//	addParamsLine("  -l <s=1>: regularization param");
//	addParamsLine("  -n <s=5>: iteration");
//}
//
//void ProgLocSharpening::produceSideInfo()
//{
//	std::cout << "Starting..." << std::endl;
//	Image<double> V;
//    V.read(fnVol);
//    V().setXmippOrigin();
//
//	FourierTransformer transformer;
//	MultidimArray<double> &inputVol = V();
//
//	Vorig = inputVol;
//
//	transformer.FourierTransform(inputVol, fftV);
//
//	iu.initZeros(fftV);
//	double uz, uy, ux, uz2, u2, uz2y2;
//	long n=0;
//	for(size_t k=0; k<ZSIZE(fftV); ++k)
//	{
//		FFT_IDX2DIGFREQ(k,ZSIZE(inputVol),uz);
//		uz2=uz*uz;
//
//		for(size_t i=0; i<YSIZE(fftV); ++i)
//		{
//			FFT_IDX2DIGFREQ(i,YSIZE(inputVol),uy);
//			uz2y2=uz2+uy*uy;
//
//			for(size_t j=0; j<XSIZE(fftV); ++j)
//			{
//				FFT_IDX2DIGFREQ(j,XSIZE(inputVol),ux);
//				u2=uz2y2+ux*ux;
//				if ((k != 0) || (i != 0) || (j != 0))
//					DIRECT_MULTIDIM_ELEM(iu,n) =sqrt(u2);
//				else
//					DIRECT_MULTIDIM_ELEM(iu,n) = 1e-38;
//				++n;
//			}
//		}
//	}
//
//    Image<double> resolutionVolume;
//    resolutionVolume.read(fnRes);
//    resVol = resolutionVolume();
//
//	resVol.setXmippOrigin();
//
//	maxMinResolution(resVol, maxRes, minRes);
//	std::cout << "minRes = " << minRes << "  maxRes = " << maxRes << std::endl;
//
//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resVol)
//	{
//		if (DIRECT_MULTIDIM_ELEM(resVol, n) < 2*sampling)
//			DIRECT_MULTIDIM_ELEM(resVol, n) = 2*sampling;
//	}
//
////	resVol.setXmippOrigin();
//
//	FourierFilter Filter;
//	Filter.FilterShape=REALGAUSSIAN;
//	Filter.FilterBand=LOWPASS;
//	Filter.w1=1;
////	Filter.apply(resVol);
//
//}
//
//void ProgLocSharpening::maxMinResolution(MultidimArray<double> &resVol,
//										double &maxRes, double &minRes)
//{
//	// Count number of voxels with resolution
//	size_t n=0;
//	double lastMinRes=1e38, lastMaxRes=1e-38, value;
//
//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(resVol)
//	{
//		value = DIRECT_MULTIDIM_ELEM(resVol, n);
////		if (value!=0)
////			{std::cout << " valor de n = " << n <<  "  resol= " << value << std::endl;}
////		    exit(0);}
//		if (value>lastMaxRes)
//			lastMaxRes = value;
//		if (value<lastMinRes)
//			lastMinRes = value;
//	}
//
//	maxRes = lastMaxRes;
//	minRes = lastMinRes;
//}
//
//void ProgLocSharpening::bandPassFilterFunction(const MultidimArray< std::complex<double> > &myfftV,
//		double w, double step, MultidimArray<double> &filteredVol, int count)
//{
//	fftVfilter.initZeros(myfftV);
//	size_t xdim, ydim, zdim, ndim;
//	Vorig.getDimensions(xdim, ydim, zdim, ndim);
//	filteredVol.resizeNoCopy(Vorig);
//
//	double w_inf = w;
//	double w_sup = w+step;
//	// Filter the input volume
//	long n=0;
//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(myfftV)
//	{
//		double un=DIRECT_MULTIDIM_ELEM(iu,n);
//
//		if (un<w_sup && un>=w_inf)
//		{
//			DIRECT_MULTIDIM_ELEM(fftVfilter, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
//		}
//	}
//
//	transformer_inv.inverseFourierTransform(fftVfilter, filteredVol);
//
//	#ifdef DEBUG_FILTER
//	Image<double> filteredvolume;
//	filteredvolume() = filteredVol;
//	filteredvolume.write(formatString("Volumen_filtrado_%i.vol", count));
//	#endif
//}
//
////void ProgLocSharpening::sameEnergy(MultidimArray< std::complex<double> > &myfftV,
////										MultidimArray<double> &localfilteredVol,
////										double &minFreq, double &maxFreq, double &step)
////{
////	std::cout << "caca" << std::endl;
////	amplitudeMonogenicSignalBP(fftV, freq, freqL, amplitudeMS, iter);
////}
//
////void ProgLocSharpening::amplitudeMonogenicSignalBP(MultidimArray< std::complex<double> > &myfftV,
////		double w, double step, MultidimArray<double> &amplitude, int count)
////{
////	double u;
////	Matrix1D<double> freq_fourier;
////	int size_fourier = ZSIZE(myfftV);
////	freq_fourier.initZeros(size_fourier);
////
////	int size = ZSIZE(Vorig);
////
////	VEC_ELEM(freq_fourier,0) = 1e-38;
////	for(size_t k=0; k<size_fourier; ++k)
////	{
////		FFT_IDX2DIGFREQ(k,size, u);
////		VEC_ELEM(freq_fourier,k) = u;
////	}
////
////	MultidimArray< std::complex<double> > fftVRiesz, fftVRiesz_aux;
////	MultidimArray<double> VRiesz;
////	fftVRiesz.initZeros(myfftV);
////	amplitude.resizeNoCopy(Vorig);
////	fftVRiesz_aux.initZeros(myfftV);
////	std::complex<double> J(0,1);
////
////	// Filter the input volume and add it to amplitude
////	long n=0;
////	double w_inf = w;
////	double w_sup = w+step;
////
////	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(myfftV)
////	{
////		double un=DIRECT_MULTIDIM_ELEM(iu,n);
////		if (un<w_sup && un>=w_inf)
////		{
////			DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
////			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
////			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
////			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= 1/un;
////		}
////	}
////	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
////
////	#ifdef DEBUG
////	Image<double> filteredvolume;
////	filteredvolume = VRiesz;
////	filteredvolume.write(formatString("Volumen_filtrado_%i.vol", count));
////	#endif
////
////	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
////		DIRECT_MULTIDIM_ELEM(amplitude,n)=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
////
////	// Calculate first component of Riesz vector
////	fftVRiesz.initZeros(myfftV);
////	double uz, uy, ux;
////	n=0;
////
////	for(size_t k=0; k<ZSIZE(myfftV); ++k)
////	{
////		for(size_t i=0; i<YSIZE(myfftV); ++i)
////		{
////			for(size_t j=0; j<XSIZE(myfftV); ++j)
////			{
////				ux = VEC_ELEM(freq_fourier,j);
////				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = ux*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
////				++n;
////			}
////		}
////	}
////	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
////	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
////		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
////
////	// Calculate second and third component of Riesz vector
////	fftVRiesz.initZeros(myfftV);
////	n=0;
////
////	for(size_t k=0; k<ZSIZE(myfftV); ++k)
////	{
////		uz = VEC_ELEM(freq_fourier,k);
////		for(size_t i=0; i<YSIZE(myfftV); ++i)
////		{
////			uy = VEC_ELEM(freq_fourier,i);
////			for(size_t j=0; j<XSIZE(myfftV); ++j)
////			{
////				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = uy*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
////				DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = uz*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
////				++n;
////			}
////		}
////	}
////	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
////	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
////		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
////
////
////	transformer_inv.inverseFourierTransform(fftVRiesz_aux, VRiesz);
////	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
////	{
////		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
////		DIRECT_MULTIDIM_ELEM(amplitude,n)=sqrt(DIRECT_MULTIDIM_ELEM(amplitude,n));
////	}
////
////	FileName iternumber;
////	Image<double> saveImg;
////	saveImg = amplitude;
////	iternumber = formatString("_Amplitude_%i.vol", count);
////	saveImg.write(fnDebug+iternumber);
////	saveImg.clear();
////
////
//////
////	// Low pass filter the monogenic amplitude
////	FourierFilter lowPassFilter;
////	lowPassFilter.w1 = w;
////	amplitude.setXmippOrigin();
////	lowPassFilter.applyMaskSpace(amplitude);
////
////	#ifdef DEBUG
////	saveImg2 = amplitude;
////	FileName fnSaveImg2;
////	if (fnDebug.c_str() != "")
////	{
////		iternumber = formatString("_Filtered_Amplitude_%i.vol", count);
////		saveImg2.write(fnDebug+iternumber);
////	}
////	saveImg2.clear();
////	#endif // DEBUG
////}
//
//void ProgLocSharpening::localfiltering(MultidimArray< std::complex<double> > &myfftV,
//										MultidimArray<double> &localfilteredVol,
//										double &minFreq, double &maxFreq, double &step)
//{
//	MultidimArray<double>  filteredVol;
//	localfilteredVol.initZeros(Vorig);
//
//	double lastResolution=1e38;
//	int idx=1, lastidx = -1;
//
//	//First band
//	double freq0=0.0;
//	double step0 = minFreq;
//	bandPassFilterFunction(myfftV, freq0, step0, filteredVol, idx);
//
////	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(filteredVol)//Esto es nuevo
////	{
////		if (DIRECT_MULTIDIM_ELEM(filteredVol, n)<0)
////			DIRECT_MULTIDIM_ELEM(filteredVol, n)=0;
////	}
//	localfilteredVol = filteredVol;
//
//	for (double freq = minFreq; freq<=maxFreq; freq+=step)
//	{
//		double res = sampling/freq;
//		double freqstep=sampling/(res-0.2);
//               step=freqstep-freq;
//
//		std::cout << " freqmin = " << freq <<  " freqmax = " << freq+step << "  resol= " << res << std::endl;
//
//		bandPassFilterFunction(myfftV, freq, step, filteredVol, idx);
//
//		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(filteredVol)
//		{
//			double resOk = DIRECT_MULTIDIM_ELEM(resVol, n);//+1e-38;
//
//			if (res<resOk)
//				DIRECT_MULTIDIM_ELEM(filteredVol, n)=0;
//
////	        if (DIRECT_MULTIDIM_ELEM(filteredVol, n)<0)
////				DIRECT_MULTIDIM_ELEM(filteredVol, n)=0;
//
//		}
////		Image<double> save;
////		save() = filteredVol;
////		save.write(formatString("filtro_%f.vol", freq));
//
//		localfilteredVol += filteredVol;
//
//		++idx;
//	}
//}
//
//
//
//void ProgLocSharpening::run()
//{
//	produceSideInfo();
//
//	MultidimArray<double> auxVol;
//
//	double step;
//
//	MultidimArray<double> operatedfiltered, Vk, filteredVol;
//
//	minRes = 2*sampling;
//    maxFreq=0.5;
//    minFreq=sampling/maxRes;
//
//	std::cout << "Resolutions between " << minRes << " and " << maxRes << std::endl;
//
//	FourierTransformer transformer;
//
//	filteredVol = Vorig;
//	sharpenedMap.resizeNoCopy(Vorig);
//
//
//    for (size_t i = 0; i<Niter; ++i)
//	{
//    	std::cout << "----------------Iteration " << i << "----------------" << std::endl;
//    	auxVol = filteredVol;
//    	transformer.FourierTransform(auxVol, fftV);
//
//
//    	   localfiltering(fftV, operatedfiltered, minFreq, maxFreq, step);
//
////    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(operatedfiltered)
////		{
////    		if (DIRECT_MULTIDIM_ELEM(operatedfiltered, n) <0)
////    			DIRECT_MULTIDIM_ELEM(operatedfiltered, n) = 0;
////		}
//
//
//		Image<double> save;
//		save() = operatedfiltered;
//		save.write(formatString("sharpenedMapA_%i.vol", i));
//
//		filteredVol = Vorig;
//		filteredVol -= operatedfiltered;
//
//      //Image<double> save;
//		save() = filteredVol;
//		save.write(formatString("sharpenedMapB_%i.vol", i));
//
//		////Second operator
//    	transformer.FourierTransform(filteredVol, fftV);
//    	localfiltering(fftV, filteredVol, minFreq, maxFreq, step);
//
////    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(filteredVol)
////		{
////    		if (DIRECT_MULTIDIM_ELEM(filteredVol, n) <0)
////    			DIRECT_MULTIDIM_ELEM(filteredVol, n) = 0;
////		}
//
//
//   	//Image<double> save;
//		save() = filteredVol;
//		save.write(formatString("sharpenedMapC_%i.vol", i));
//
//		if (i == 0)
//			Vk = Vorig;
//		else
//			Vk = sharpenedMap;
//
//		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sharpenedMap)
//			DIRECT_MULTIDIM_ELEM(sharpenedMap,n)=DIRECT_MULTIDIM_ELEM(Vk,n)+
//			   lambda*DIRECT_MULTIDIM_ELEM(filteredVol,n);
//
////    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sharpenedMap)
////		{
////    		if (DIRECT_MULTIDIM_ELEM(sharpenedMap, n) <0)
////    			DIRECT_MULTIDIM_ELEM(sharpenedMap, n) = 0;
////		}
//
//		Image<double> filteredvolume0;
//		filteredvolume0() = sharpenedMap;
//		filteredvolume0.write(formatString("sharpenedMapD_%i.vol", i));
//
//		filteredVol = sharpenedMap;
//	}
//	Image<double> filteredvolume;
//	filteredvolume() = sharpenedMap;
//	filteredvolume.write(fnOut);
//
//}
