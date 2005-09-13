/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2002)
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

#include "../Prog_make_spectra.hh"
#include <XmippData/xmippArgs.hh>

// Empty constructor -------------------------------------------------------
Prog_make_spectra_prm::Prog_make_spectra_prm(): Prog_parameters() {
   each_image_produces_an_output=false;
   rot_spt.x0=rot_spt.y0=-1;
}

// Read from command line --------------------------------------------------
void Prog_make_spectra_prm::read(int argc, char **argv) {
   Prog_parameters::read(argc,argv);
   fn_out=get_param(argc,argv,"-o");
   rot_spt.read(argc,argv);
   produce_side_info();
}

// Produce side info -------------------------------------------------------
void Prog_make_spectra_prm::produce_side_info() {
}

// Show --------------------------------------------------------------------
void Prog_make_spectra_prm::show() {
   Prog_parameters::show();
   show_specific();
}

void Prog_make_spectra_prm::show_specific() {
   cout << "Output file: " << fn_out << endl;
   cout << rot_spt << endl;
}

// Usage -------------------------------------------------------------------
void Prog_make_spectra_prm::usage() {
   Prog_parameters::usage();
   usage_specific();
}

void Prog_make_spectra_prm::usage_specific() {
   cerr << "   -o <fn_out>                 : Output file with the spectra\n"
   ;
   rot_spt.usage();
}

// Process an image --------------------------------------------------------
void Prog_make_spectra_prm::process_img(ImageXmipp &img) {
   rot_spt.compute_rotational_spectrum(img(),rot_spt.rl, rot_spt.rh,
      rot_spt.dr, rot_spt.rh-rot_spt.rl);
   Harmonics.push_back(rot_spt.rot_spectrum);
   Img_name.push_back(img.name());
}

// Finish processing -------------------------------------------------------
void Prog_make_spectra_prm::finish_processing() {
   ofstream fh_out;
   fh_out.open(fn_out.c_str());
   
   if (!fh_out)
      REPORT_ERROR(1,(string)"Prog_make_spectra_prm::finish_processing: "
         "Cannot open"+fn_out+" for output");
   if (Harmonics.size()!=0) {
      fh_out << XSIZE(Harmonics[0]) << " " << Harmonics.size() << endl;
      int imax=Harmonics.size();
      for (int i=0; i<imax; i++) {
         double norm=Harmonics[i].sum()/100.0;
	 for (int j=0; j<XSIZE(Harmonics[i]); j++)
            fh_out << FtoA(Harmonics[i](j)/norm,6,4) << " ";
	 fh_out << Img_name[i] << endl;
      }
   }
   fh_out.close();
}
