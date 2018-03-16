# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from pyworkflow import VERSION_1_1
from pyworkflow.protocol.params import (PointerParam, StringParam, 
                                        BooleanParam, FloatParam, IntParam, LEVEL_ADVANCED)
from pyworkflow.em.protocol.protocol_3d import ProtAnalysis3D
from pyworkflow.object import Float
from pyworkflow.em import ImageHandler
from pyworkflow.utils import getExt
from pyworkflow.em.data import Volume
import numpy as np
import pyworkflow.em.metadata as md

CHIMERA_RESOLUTION_VOL = 'MG_Chimera_resolution.vol'



class XmippProtLocSharp(ProtAnalysis3D):
    """    
    Given a map the protocol assigns local resolutions to each voxel of the map.
    """
    _label = 'local Sharpening'
    _lastUpdateVersion = VERSION_1_1
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
        self.min_res_init = Float() 
        self.max_res_init = Float()
       
    
    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input Map", important=True,
                      help='Select a volume for sharpening.')

        form.addParam('resolutionVolume', PointerParam, pointerClass='Volume',
                      label="Resolution Map", important=True,
                      help='Select a local resolution map.')

        form.addParam('iterations', IntParam, default=5, 
                      expertLevel=LEVEL_ADVANCED,
                      label="Iterations",
                      help='Number of iterations.')
        
        form.addParam('const', FloatParam, default=5, 
                      expertLevel=LEVEL_ADVANCED,
                      label="lambda",
                      help='Regularization Param.')
  
    # --------------------------- INSERT steps functions --------------------------------------------


    def _insertAllSteps(self):
            # Convert input into xmipp Metadata format
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('sharpenStep')
        self._insertFunctionStep('createOutputStep')


    def convertInputStep(self):
        """ Read the input volume.
        """

        self.volFn = self.inputVolume.get().getFileName()
        self.resFn = self.resolutionVolume.get().getFileName()      
        extVol = getExt(self.volFn)
        extRes = getExt(self.resFn)        
        if (extVol == '.mrc') or (extVol == '.map'):
            self.volFn = self.volFn + ':mrc'
        if (extRes == '.mrc') or (extRes == '.map'):
            self.resFn = self.resFn + ':mrc'
            
    def  sharpenStep(self):   
           
        params = ' --vol %s' % self.volFn
        params += ' --resolution_map %s' % self.resFn
        params += ' --sampling %f' % self.inputVolume.get().getSamplingRate()
        params += ' -l %f' % self.const
        params += ' -n %i' % self.iterations
        params += ' -o %s' % self._getExtraPath('filtered.vol')
        
        self.runJob('xmipp_volume_local_sharpening', params)


    def createOutputStep(self):
        volume=Volume()
        volume.setFileName(self._getExtraPath('filtered.vol'))
        volume.setSamplingRate(self.inputVolume.get().getSamplingRate())
#         readSetOfVolumes(volume_path, self.volumesSet)
        self._defineOutputs(sharpened_map=volume)
        self._defineSourceRelation(self.inputVolume, volume)
            


    # --------------------------- INFO functions ------------------------------

    def _methods(self):
        messages = []
        if hasattr(self, 'sharpened_map'):
            messages.append(
                'Information about the method/article in ' + MONORES_METHOD_URL)
        return messages
    
    def _summary(self):
        summary = []
        summary.append("Highest resolution jejejeje")
        return summary

    def _citations(self):
        return ['Vilas2017']

