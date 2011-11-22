#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based particle picking
#
# Author: Carlos Oscar Sorzano, August, 2011
#
#------------------------------------------------------------------------------------------------

from config_protocols import protDict
from protlib_base import *
from protlib_utils import runJob
from protlib_filesystem import createLink
import xmipp
from glob import glob
from os.path import exists, join

# Create a GUI automatically from a selfile of micrographs
class ProtParticlePicking(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.particle_pick.name, scriptname, project)
        self.Import = "from protocol_particle_pick import *"
        self.micrographSelfile = os.path.join(getWorkingDirFromRunName(self.ImportRun), "micrographs.sel")

    def defineSteps(self):
        self.Db.insertStep('createLink',execution_mode=SqliteDb.EXEC_MAINLOOP,
                           source=self.micrographSelfile,dest=self.workingDirPath("micrographs.sel"))
        self.Db.insertStep('launchParticlePickingGUI',execution_mode=SqliteDb.EXEC_ALWAYS,
                           MicrographSelfile=self.micrographSelfile, WorkingDir=self.WorkingDir,
                           AutomaticPicking=self.AutomaticPicking, NumberOfThreads=self.NumberOfThreads,
                           Fast=self.Fast,InCore=self.InCore)       

    def summary(self):
        summary = []
        
        mD = xmipp.MetaData(self.micrographSelfile)
        isPairList = mD.containsLabel(xmipp.MDL_ASSOCIATED_IMAGE1) and not xmipp.FileName(self.micrographSelfile).isStar1()

        if isPairList: 
            suffix = "tilt pairs"
        else: 
            suffix = "micrographs"
        
        
        summary=["Input: <%s> with <%u> %s" % (self.micrographSelfile, mD.size(), suffix)]        
        total_manual = 0
        N_manual = 0
        total_auto = 0
        N_auto = 0
        Nblock = {}
        for posfile in glob(self.workingDirPath('*.pos')):
            blockList = xmipp.getBlocksInMetaDataFile(posfile)
            if 'families' in blockList:
                blockList.remove('families')
            particles = 0
            for block in blockList:
                mD = xmipp.MetaData(block+"@"+posfile);
                mD.removeDisabled();
                Nparticles = mD.size()
                particles += Nparticles
                if block in Nblock.keys():
                    Nblock[block] += Nparticles
                else:
                    Nblock[block] = Nparticles
            if particles > 0:
                if 'auto' in posfile:
                    total_auto += particles
                    N_auto += 1
                else:
                    total_manual += particles
                    N_manual += 1
        summary.append("Number of particles manually picked: <%d> (from <%d> micrographs)" % (total_manual, N_manual))
        if N_auto > 0:
            summary.append("Number of particles automatically picked: <%d> (from <%d> micrographs)" % (total_auto, N_auto))
        fnFamilies = self.workingDirPath("families.xmd")
        if exists(fnFamilies):
            mD = xmipp.MetaData(fnFamilies)
            Nfamilies = mD.size()
            if Nfamilies>1:
                summary.append("Number of families: <%u>" % Nfamilies)
        for block in Nblock.keys():
            summary.append("Family <%s>: <%u> particles" % (block, Nblock[block]))
        return summary
    
    def validate(self):
        errors = []
        # Check that there is a valid list of micrographs
        if not exists(self.micrographSelfile)>0:
            errors.append("Cannot find "+self.micrographSelfile)
            return errors
        # Check that all micrographs exist
        errMsg = ""
        mD = xmipp.MetaData(self.micrographSelfile)
        isPairList = mD.containsLabel(xmipp.MDL_ASSOCIATED_IMAGE1) and not xmipp.FileName(self.micrographSelfile).isStar1()
        for id in mD:
            micrograph = mD.getValue(xmipp.MDL_IMAGE,id)
            if not exists(micrograph):
                errMsg += "  " + micrograph
            if isPairList:
                micrograph = mD.getValue(xmipp.MDL_ASSOCIATED_IMAGE1,id)
                if not exists(micrograph):
                    errMsg += "  " + micrograph
        
        if len(errMsg):
            errors.append("Cannot find the following micrographs: " + errMsg)
                
        # Check that automatic particle picking is not for tilted
        if isPairList and self.AutomaticPicking:
            errors.append("Automatic particle picking cannot be done on tilt pairs")
        
        return errors
    
    def visualize(self):
        launchParticlePickingGUI(None, self.micrographSelfile, self.WorkingDir)

# Execute protocol in the working directory
def launchParticlePickingGUI(log,MicrographSelfile,WorkingDir,
                             AutomaticPicking=False, NumberOfThreads=1, Fast=True, InCore=False):
    mode = "manual"
    if AutomaticPicking:
        mode = "supervised %(NumberOfThreads)d %(Fast)s %(InCore)s"
    params = "-i %(MicrographSelfile)s -o %(WorkingDir)s --mode " + mode 
    runJob(log,"xmipp_micrograph_particle_picking", params % locals(), RunInBackground=True)

#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtParticlePicking)
