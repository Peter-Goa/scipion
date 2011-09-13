#!/usr/bin/env python
'''
#/***************************************************************************
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
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
# *  e-mail address 'xmipp@cnb.csic.es'
# ***************************************************************************
 '''
 
import os
import sys
import shutil
import ConfigParser
from config_protocols import projectDefaults, sections, protDict
from protlib_sql import SqliteDb, XmippProjectDb, XmippProtocolDb
from protlib_utils import XmippLog, loadModule, reportError, getScriptPrefix, failStr
from protlib_filesystem import deleteFiles

class XmippProject():
    def __init__(self, projectDir=None):
        if not projectDir:
            self.projectDir = os.getcwd()
        else:
            self.projectDir = projectDir
        self.cfgName = projectDefaults['Cfg']
        self.dbName =  projectDefaults['Db']
        self.logsDir = projectDefaults['LogsDir']
        self.runsDir = projectDefaults['RunsDir']
        
    
    def exists(self):
        ''' a project exists if the data base can be opened and directories Logs and Runs
        exists'''
        from protlib_sql import existsDB
        status = True
        if not os.path.exists(self.logsDir):
            status = False
        elif not os.path.exists(self.runsDir):
            status = False
        elif not os.path.exists(self.cfgName):
            status = False
        elif not existsDB(self.dbName):
            status = False
        return status
    
    def create(self):
        os.chdir(self.projectDir)
        print "Creating new project on directory: '%s'" % self.projectDir
        #==== CREATE CONFIG file
        self.config = ConfigParser.RawConfigParser()            
        self.config.add_section('project')
        self.config.set('project', 'projectdir', self.projectDir)
        self.config.set('project', 'systemflavour', '')
        self.SystemFlavour = self.config.get('project', 'systemflavour')

        self.writeConfig()
        #===== CREATE LOG AND RUN directories
        if not os.path.exists(self.logsDir):
            os.mkdir(self.logsDir)
        if not os.path.exists(self.runsDir):
            os.mkdir(self.runsDir)
        #===== CREATE DATABASE
        self.projectDb  = XmippProjectDb(self.dbName)
        #===== POPULATE SOME TABLES
        groupName = ""
        for section, groupList in sections:
            for group in groupList:
                groupName = group[0]
                prots = group[1:]
                self.projectDb.insertGroup(groupName)
                for p in prots:
                    self.projectDb.insertProtocol(groupName, p)
        #Hard coded insertion of xmipp_program protocol
        #this is an special case of protocols
        self.projectDb.insertProtocol(protDict.xmipp_program.title, protDict.xmipp_program.name)
        # commit changes
        self.projectDb.connection.commit()
        
    def load(self):
        self.config = ConfigParser.RawConfigParser()
        self.config.read(self.cfgName)
        # Load database
        self.projectDb = XmippProjectDb(self.dbName)
        self.SystemFlavour = self.config.get('project', 'systemflavour')
        

    def writeConfig(self):
        with open(self.cfgName, 'wb') as configfile:
            self.config.write(configfile) 
            
    def clean(self):
        ''' delete files related with a project'''
        if os.path.exists(self.logsDir):
            shutil.rmtree(self.logsDir)
        if os.path.exists(self.runsDir):
            shutil.rmtree(self.runsDir)
        if os.path.exists(self.cfgName):
            os.remove(self.cfgName)
        if os.path.exists(self.dbName):
            os.remove(self.dbName)
            
        for p in protDict.values():
            if os.path.exists(p.dir):
                shutil.rmtree(p.dir)
        self.create()
            
    
    def cleanRun(self, run):
        script = run['script']
        toDelete = [script.replace(self.runsDir, self.logsDir).replace(".py", ext) for ext in ['.log', '.err', '.out']]
        deleteFiles(None, toDelete, False)
        workingDir = self.getWorkingDir(run['protocol_name'], run['run_name'])
        from distutils.dir_util import remove_tree
        if os.path.exists(workingDir):
            remove_tree(workingDir, True)

    def deleteRun(self, run):
        script = run['script']
        deleteFiles(None, [script, script+'c'], False)
        self.cleanRun(run)
        self.projectDb.deleteRun(run)
        
    def deleteTmpFiles(self):
        for section, groupList in sections:
            for group in groupList:
                groupName = group[0]
                runs = self.projectDb.selectRuns(groupName)
                for run in runs:
                    tmpDir = os.path.join(self.getWorkingDir(run['protocol_name'], run['run_name']), 'tmp')
                    if os.path.exists(tmpDir):
                        shutil.rmtree(tmpDir)

    def getWorkingDir(self,protocol_name,run_name):
        return os.path.join(protDict[protocol_name].dir,run_name)
            
    def getRunScriptFileName(self, protocol_name, runName):
        return os.path.join(self.runsDir, '%s_%s.py' % (protocol_name, runName))
    
    def createRunFromScript(self, protocol_name, script, prefix=None):
        if prefix: #Remove protocol_name from prefix if present
            prefix = prefix.replace(protocol_name+"_", '')      
        #suggest a new run_name  
        runName = self.projectDb.suggestRunName(protocol_name, prefix)
        dstAbsPath = self.getRunScriptFileName(protocol_name, runName)
        run = {
               'protocol_name':protocol_name, 
               'run_name': runName, 
               'script': dstAbsPath, 
               'comment': "",
               'source': script
               }
        return run
    
    def newProtocol(self, protocol_name):
        from protlib_filesystem import getXmippPath
        srcProtAbsPath = os.path.join(getXmippPath('protocols'),
                                       'protocol_%s_header.py' % protocol_name)
        return self.createRunFromScript(protocol_name, srcProtAbsPath)
    
    def copyProtocol(self, protocol_name, script):
        prefix, suffix  = getScriptPrefix(script)
        return self.createRunFromScript(protocol_name, script, prefix)
        
    def loadProtocol(self, protocol_name, runName):
        run = self.projectDb.selectRunByName(self, protocol_name, runName)
        if not run is None:
            run['source'] = run['script']
        return run
    
    def newOrLoadProtocol(self, protocol_name, runName):
        run = self.loadProtocol(self.project, protocol_name, runName)
        if run is None:
            run = self.newProtocol(protocol_name)
        return run
            
class XmippProtocol(object):
    '''This class will serve as base for all Xmipp Protocols'''
    def __init__(self, protocolName, scriptname, project):
        '''Basic constructor of the Protocol
        protocolName -- the name of the protocol, should be unique
        scriptname,     file containing the protocol instance
        runName      -- the name of the run,  should be unique for one protocol
        project      -- project instance
        '''
        #Import all variables in the protocol header
        self.Header = loadModule(scriptname)
        for k, v in self.Header.__dict__.iteritems():
            self.__dict__[k] = v
            
        self.DoParallel = 'NumberOfMpi' in dir() and self.NumberOfMpi > 1
        self.Name = protocolName
        self.scriptName = scriptname
        #A protocol must be able to find its own project
        self.project = project
        self.Import = '' # this can be used by database for import modules
        self.WorkingDir = project.getWorkingDir(protocolName, self.RunName)
        self.TmpDir = os.path.join(self.WorkingDir, 'tmp')
        self.projectDir = project.projectDir  
        #Setup the Log for the Protocol
        self.LogDir = project.logsDir
        runName = self.RunName
        self.uniquePrefix = "%(protocolName)s_%(runName)s" % locals()
        self.LogPrefix = os.path.join(self.LogDir, self.uniquePrefix)       
        self.Err = self.LogPrefix+".err"
        self.Out = self.LogPrefix+".out"
        self.LogFile = self.LogPrefix + ".log"
        self.SystemFlavour = project.SystemFlavour
        
    def getProjectId(self):
        pass
        #self.project = project.getId(launchDict['Projection Matching'],runName,)
        
    def workingDirPath(self, path):
        os.path.join(self.WorkingDir, path)

    def getRunState(self):
        return self.project.projectDb.getRunStateByName(self.Name, self.RunName)
    
    def validateBase(self):
        '''Validate if the protocols is ready to be run
        in this function will be implemented all common
        validations to all protocols and the function
        validate will be called for specific validation
        of each protocol'''
        errors=[]
        #check if there is a valid project, otherwise abort
        if not self.project.exists():
            errors.append("Not valid project available")
        # Check if there is runname
        if self.RunName == "":
            errors.append("No run name given")
            
        #Check that number of threads and mpi are int and greater than 0
        if 'NumberOfThreads' in dir(self):
            if self.NumberOfThreads<1:
                errors.append("Number of threads has to be >=1")
        if 'NumberOfMpi' in dir(self):
            if self.NumberOfMpi<1:
                errors.append("Number of MPI processes has to be >=1")
        
        #specific protocols validations
        errors += self.validate()
        
        return errors
    
    def validate(self):
        '''Validate if the protocols is ready to be run
        it may be redefined in derived protocol classes but do not forget
        to call the main class with
        super(ProtProjMatch, self).validate()
        '''
        return []
    
    def summary(self):
        '''Produces a summary with the most relevant information of the protocol run'''
        return []
    
    def visualize(self):
        '''Visualizes the results of this run'''
        pass
    
    def warningsBase(self):
        '''Output some warnings that can be errors and require user confirmation to proceed'''
        warningList=[]
        if self.Behavior=="Restart":
            warningList.append("Restarting a protocol will delete previous results")
        warningList += self.warnings()
        return warningList

    def warnings(self):
        '''Output some warnings that can be errors and require user confirmation to proceed'''
        return []
    
    def postRun(self):
        '''This function will be called after the run function is executed'''
        pass   
    
    def defineSteps(self):
        '''In this function the actions to be performed by the protocol will be added to db.
        each particular protocol need to add its specific actions. Thing to be performed before 
        "run" should be added here'''
        pass
    
    def runSetup(self, isMainLoop=True):

        #Redirecting standard output and error to files
        self.fOut = open(self.Out, 'a')
        self.fErr = open(self.Err, 'a')
        self.stderr = sys.stderr # backup stderr
        sys.stdout = self.fOut
        sys.stderr = self.fErr
        self.Log = XmippLog(self.LogFile)
        self.Db  = XmippProtocolDb(self, isMainLoop)

    def run(self):
        '''Run of the protocols
        if the other functions have been correctly implemented, this not need to be
        touched in derived class, since the run of protocols should be the same'''
        
        #Change to project dir
        os.chdir(self.projectDir)

        errors = self.validateBase()
        if len(errors) > 0:
            raise Exception('\n'.join(errors))
        
        #insert basic operations for all scripts
        if self.Behavior=="Restart":
            run = {
               'protocol_name':self.Name, 
               'run_name': self.RunName, 
               'script': self.scriptName, 
               'source': self.scriptName
               }
            self.project.cleanRun(run)
        #Remove temporaly files
        if os.path.exists(self.TmpDir):
            shutil.rmtree(self.TmpDir)
        #Initialization of log and db
        retcode = 0
        try:
            self.runSetup()
            self.Db.insertStep('createDir', path = self.WorkingDir)
            self.Db.insertStep('createDir', path = self.TmpDir)
            self.defineSteps()
            self.Db.runSteps()
            self.postRun()
        except Exception, e:
            retcode = 1;
            print >> sys.stderr, failStr("ERROR: %s" %  e)
            #THIS IS DURING DEVELOPMENT ONLY
            print >> self.stderr, failStr("ERROR: %s" %  e)
            import traceback
            traceback.print_exc(file=self.stderr)
            
        finally:
            self.fOut.close()
            self.fErr.close()  
                            
        return retcode
    
def command_line_options():
    '''process protocol command line'''
    import optparse
    parser = optparse.OptionParser()
    parser.add_option('-g', '--gui',
                              dest="gui",
                              default=False,
                              action="store_true",
                              help="use graphic interface to launch protocol "
                              )
    parser.add_option('-c', '--no_check',
                              dest="no_check",
                              default=False,
                              action="store_true",
                              help="do NOT check run checks before execute protocols"
                              )
    parser.add_option('-n', '--no_confirm',
                          dest="no_confirm",
                          default=False,
                          action="store_true",
                          help="do NOT ask confirmation for warnings"
                          )
        
    options = parser.parse_args()[0]
    return options

def getProtocolFromModule(script, project):
    mod = loadModule(script)
    from inspect import isclass
    for v in mod.__dict__.values():
        if isclass(v) and issubclass(v, XmippProtocol) and v != XmippProtocol:
            return v(script, project)
    reportError("Can load protocol from " + script)

def getExtendedRunName(run):
    ''' Return the extended run name, ie: protocol_runName '''
    return "%s_%s" % (run['protocol_name'], run['run_name'])

def splitExtendedRunName(extendedRunName):
    ''' Take an extended runname of the form protocol_runName and returns (protocol,runName).
    If no protocol can be matched, then None is returned '''
    from config_protocols import protDict
    runName=""
    protocolName=""
    bestLength=0
    for k in protDict.keys():
        if extendedRunName.startswith(k):
            len_k=len(k)
            if len_k>bestLength:
                runName=extendedRunName.replace(k+"_","")
                protocolName=k
                bestLength=len_k
    if bestLength==0:
        return None
    else:
        return (protocolName,runName)

def getWorkingDirFromRunName(extendedRunName):
    # The extended run name has the name of the protocol in front
    tuple=splitExtendedRunName(extendedRunName)
    if tuple is None:
        return None
    protocolName=tuple[0]
    runName=tuple[1]
    return os.path.join(protDict[protocolName].dir,runName)

def getScriptFromRunName(extendedRunName):
    return os.path.join(projectDefaults['RunsDir'],extendedRunName+".py")

def protocolMain(ProtocolClass, script=None):
    gui = False
    no_check = False
    no_confirm = False
    doRun = True
    
    if script is None:
        script  = sys.argv[0]
        options = command_line_options()
        gui = options.gui
        no_check = options.no_check
        no_confirm = options.no_confirm
    
    mod = loadModule(script)
    #init project
    project = XmippProject()
    #load project: read config file and open conection database
    project.load()
    #register run with runName, script, comment=''):
    p = ProtocolClass(script, project)
    run_id = project.projectDb.getRunId(p.Name, mod.RunName)
    
    #1) just call the gui    
    if gui:
        run = {
           'run_id': run_id,
           'protocol_name':p.Name, 
           'run_name': mod.RunName, 
           'script': script, 
           'comment': "",
           'source': script
           }
        from protlib_gui import ProtocolGUI 
        gui = ProtocolGUI()
        gui.createGUI(project, run)
        gui.fillGUI()
        gui.launchGUI()
    else:#2) Run from command line
        
        if no_check: 
            if not run_id:
                reportError("Protocol run '%s' has not been registered in project database" % mod.RunName)
        else:
            from protlib_utils import showWarnings
            if not showWarnings(p.warningsBase(), no_confirm):
                exit(0)
            
            if not run_id:
                _run={
                       'comment' : ''
                      ,'protocol_name': p.Name
                      ,'run_name' : mod.RunName
                      ,'script'  : script 
                      }
                project.projectDb.insertRun(_run)
                run_id = _run['run_id']
            else:
                _run = {'run_id': run_id}
                
            if 'SubmitToQueue' in dir(mod) and mod.SubmitToQueue:
                from protlib_utils import submitProtocol
                NumberOfThreads = 1
                if 'NumberOfThreads' in dir(mod):
                    NumberOfThreads=mod.NumberOfThreads
                pbsPid = submitProtocol(script,
                               jobId = p.uniquePrefix,
                               queueName = mod.QueueName,
                               nodes = mod.NumberOfMpi,
                               threads = NumberOfThreads,
                               hours = mod.QueueHours,
                               command = 'python %s --no_check' % script
                               )
                project.projectDb.updateRunState(SqliteDb.RUN_LAUNCHED, run_id)
                doRun = False
                _run['pid'] = os.getpid()
                _run['pid_type'] = pbsPid              
            else:
                _run['pid'] = os.getpid()
                _run['pid_type'] = SqliteDb.PID_POSIX
            
            print _run
            # Update run's process info in DB
            project.projectDb.updateRunPid(_run)
        
        if doRun:
            return p.run()
