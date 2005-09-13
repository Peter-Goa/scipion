/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/

#ifndef _FOURIER_FILTER_HH
   #define _FOURIER_FILTER_HH

#include "../CTF.hh"
#include <XmippData/xmippMatrices3D.hh>
#include <XmippData/xmippFFT.hh>

/**@name Fourier Masks */
//@{
/** Filter class for Fourier space.

   Example of use for highpass filtering
   \begin{verbatim}
      ImageXmipp I("image.xmp");
      FourierMask Filter;
      Filter.FilterShape=RAISED_COSINE;
      Filter.FilterBand=HIGHPASS;
      Filter.w1=w_cutoff;
      Filter.raised_w=slope;
      I().set_Xmipp_origin();
      Filter.generate_mask(I());
      Filter.apply_mask_Space(I());
      I.write("filtered_image.xmp");
   \end{verbatim}

   Example of use reading a mask from file
   \begin{verbatim}
      ImageXmipp I("image.xmp");
      I().set_Xmipp_origin();
      FourierMask Filter;
      Filter.read_mask("mask.fft");
      Filter.apply_mask_Space(I());
      I.write("filtered_image.xmp");
   \end{verbatim}
*/
class FourierMask {
public:
   #define RAISED_COSINE 1
   #define FROM_FILE     6
   /** Shape of the decay in the filter.
      Valid types are RAISED_COSINE, FROM_FILE. */
   int FilterShape;
   
   #define LOWPASS       1
   #define HIGHPASS      2
   #define BANDPASS      3
   #define STOPBAND      4
   #define CTF           5
   #define WEDGE         7
   /** Pass band. LOWPASS, HIGHPASS, BANDPASS, STOPBAND, CTF, WEDGE, FROM_FILE */
   int FilterBand;
   
   /** Cut frequency for Low and High pass filters, first freq for bandpass.
       Normalized to 1/2*/
   double w1;
   
   /** Second frequency for bandpass and stopband. Normalized to 1/2 */
   double w2;
   
   /** Pixels around the central frequency for the raised cosine */
   double raised_w;

   /** File from which the mask is read */
   FileName fn_mask;

   /** CTF parameters. */
   XmippCTF ctf;
   
   /** Mask1D */
   matrix1D< complex<double> > mask1D;
   
   /** Mask2D */
   matrix2D< complex<double> > mask2D;
   
   /** Mask3D */
   matrix3D< complex<double> > mask3D;
public:
   /** Empty constructor */
   FourierMask() {clear();}

   /** Destructor */
   ~FourierMask() {clear();}

   /** Assignment */
   FourierMask & operator = (const FourierMask &F);

   /** Clear */
   void clear();

   /** Read parameters from command line.
       If a CTF description file is provided it is read. */
   void read(int argc, char **argv);

   /** Show. */
   void show();

   /** Usage. */
   void usage();

   /** Generate mask for a resized image.
   It is supposed that the image is already resized and with its logical
   origin set. If the filter is a CTF it must be already read and prepared
   in the ctf variable.
   
   The dimension is 1, 2 or 3 depending it is a signal,
   an image or a volume.
   
   This function cannot be used to generate CTF masks.*/
   template <class T>
   void generate_mask(T &v);

   /** Read mask from file. */
   void read_mask(const FileName &fn);

   /** Generate CTF mask for images */
   template <class T>
   void generate_CTF_mask(T &v);

   /** Flip the phase of an already generated 2D mask.
       Those frequencies that have negative amplitude are flipped. */
   void correct_phase();

   /** Save mask as a text file.
       Indicate which is the dimension od the mask to save */
   void write_mask(const FileName &fn, int dim);

   /** Save amplitude as a text file, Indicate
       which is the dimension of the mask to save. */
   void write_amplitude(const FileName &fn, int dim,
      bool do_not_center=FALSE);

   /** Apply mask (argument is in Fourier space).
       It should have been already generated. The given image is modified.
       An exception is thrown if the mask do not fit the size and shape of
       the */
   void apply_mask_Fourier(matrix1D< complex<double> > &v);

   /** Apply mask in 2D. */
   void apply_mask_Fourier(matrix2D< complex<double> > &v);

   /** Apply mask in 3D. */
   void apply_mask_Fourier(matrix3D< complex<double> > &v);

   /** Apply mask (argument is in real space)..
       It doesn't need to have a mask already generated. If the mask is equal
       to NULL, the apropiate mask is generated if not, nothing is done with
       the mask except filtering. The given image is modified. */
   void apply_mask_Space(matrix1D<double> &v);

   /** Apply mask in 2D. */
   void apply_mask_Space(matrix2D<double> &v);

   /** Apply mask in 3D. */
   void apply_mask_Space(matrix3D<double> &v);

   /** Resize fourier mask to a desired scale. */
   void resize_mask(int Ydim, int Xdim);

   /** Resize fourier mask to a desired scale. */
   void resize_mask(int Zdim, int Ydim, int Xdim);

   /** Mask power. Return the power of the Fourier Image contained in mask 
   within the given frequencies. */
   double mask2D_power(double wmin=0, double wmax=1);
};
//@}
#endif
