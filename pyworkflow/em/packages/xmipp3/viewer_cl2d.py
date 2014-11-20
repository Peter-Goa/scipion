# **************************************************************************
# *
# * Authors:     L. del Cano (ldelcano@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
from pyworkflow.em.packages.xmipp3.convert import readSetOfClasses
"""
This module implement the wrappers aroung Xmipp CL2D protocol
visualization program.
"""
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from protocol_cl2d import XmippProtCL2D
from pyworkflow.gui.text import *
from pyworkflow.gui.dialog import showError, showWarning
from pyworkflow.protocol.params import HiddenBooleanParam
import glob

CLASSES = 0
CLASS_CORES = 1
CLASS_STABLE_CORES = 2   
CLASS_CHOICES = ['Classes', 'Class Cores', 'Class Stable Cores']

LEVEL_LAST = 0
LEVEL_SEL = 1
LEVEL_CHOICES = ['last', 'selection']
        
        
class XmippCL2DViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _label = 'viewer cl2d'
    _targets = [XmippProtCL2D]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def __init__(self, **args):
        ProtocolViewer.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('classesToShow', EnumParam, choices=CLASS_CHOICES,
                      label="What to show", default=CLASS_CORES,
                      display=EnumParam.DISPLAY_LIST)
        form.addParam('doShowLastLevel', EnumParam, default=LEVEL_LAST, 
                      choices=LEVEL_CHOICES,
                      display=EnumParam.DISPLAY_LIST,
                      label="Level to visualize")     
        form.addParam('showSeveralLevels', StringParam, default='',
              label='Levels selection', condition='doShowLastLevel==%d' % LEVEL_SEL,
              help='Specify a  list of levels like: 0,1,3 or 0-3 ')    
        form.addParam('doShowClassHierarchy', HiddenBooleanParam, default=False,
                      label="Visualize class hierarchy.")      
    
    def _getVisualizeDict(self):
        return {'doShowClassHierarchy': self._viewClassHierarchy,
                'doShowLastLevel': self._viewLevelFiles}        

    def _getSubset(self):
        fnSubset = ""
        classesToShow = self.classesToShow.get()
        if classesToShow == CLASSES:
            fnSubset = ""
        elif classesToShow == CLASS_CORES:
            fnSubset = "_core"
        elif classesToShow == CLASS_STABLE_CORES:
            fnSubset = "_stable_core" 
        return fnSubset
                        
    def _viewClassHierarchy(self, e=None):
        fnSubset = self._getSubset()
        fnHierarchy = self.protocol._getExtraPath("classes%s_hierarchy.txt" % fnSubset)
        if os.path.exists(fnHierarchy):
            return [self.textView([fnHierarchy])] 
            
    def _getInputParticles(self):
        return self.protocol.inputParticles.get()
    
    def _getClassesSqlite(self, blockName, fn):
        """ Read the classes from Xmipp metadata and write as sqlite file. """
        fnSqlite = fn + '.sqlite'
        
        if not os.path.exists(fnSqlite):
            classesSet = SetOfClasses2D(filename=fnSqlite)
            classesSet.setImages(self._getInputParticles())
            readSetOfClasses(classesSet, fn, blockName)
            classesSet.write()
        return fnSqlite
    
    def _viewClasses(self, blockName, fn):
        fnSqlite = self._getClassesSqlite(blockName, fn)
        return ClassesView(self._project.getName(), self.protocol.strId(), fnSqlite, 
                           other=self._getInputParticles().strId())
    
    def _viewLevelFiles(self, e=None):
        fnSubset = self._getSubset()
        views = []
        errors = []
        #obj = getattr(self.protocol, "outputClasses" + fnSubset)
        levelFiles = self.protocol._getLevelMdFiles(fnSubset)
        
        if levelFiles:
            levelFiles.sort()
            
            if self.doShowLastLevel == LEVEL_LAST:
                views.append(self._viewClasses("classes", levelFiles[-1]))
            else:
                if self.showSeveralLevels.empty():
                    errors.append('Please select the levels that you want to visualize.')
                else:
                    listOfLevels = []
                    try:
                        listOfLevels = self._getListFromRangeString(self.showSeveralLevels.get())
                    except Exception, ex:
                        errors.append('Invalid levels range.')
                        
                    files = []
                    for level in listOfLevels:
                        fn = self.protocol._getExtraPath("level_%02d/level_classes%s.xmd"%(level,fnSubset))
                        if os.path.exists(fn):
                            files.append(("classes_sorted", fn))
                        else:
                            errors.append('Level %s does not exist.' % level)
                    
                    for block, fn in files:                        
                        views.append(self._viewClasses(block, fn))
        
        self.errorList(errors, views)
            
        return views        

