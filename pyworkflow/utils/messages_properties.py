# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *              Adrian Quintana (aquintana@cnb.csic.es)
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
"""
This module defines the text used in the application.
"""

class Message():
    # Example Usage: 
    # MyMessage = Message()
    # print MyMessage.label

    # Header List
    VIEW_PROJECTS = 'Projects'
    VIEW_PROTOCOLS = 'Protocols'
    VIEW_DATA = 'Data'
    VIEW_HOSTS = 'Hosts'
    
    # Projects Template
    LABEL_PROJECTS='Projects'
    LABEL_CREATE_PROJECT='Create Project'
    TITLE_CREATE_PROJECT='Enter the project name'
    TITLE_CREATE_PROJECT_NAME='Project Name: '
    MESSAGE_CREATE_PROJECT=', Are you sure to *DELETE* the project and all its *DATA*?'
    LABEL_DELETE_PROJECT='Delete Project'
    TITLE_DELETE_PROJECT='Confirm project deletion'
    LABEL_MODIFIED='Modified: '
    
    # Project Content Template
    LABEL_PROJECT='Project '
    #-- Protocol Treeview --
    LABEL_WORKFLOW='Workflow View: '
    LABEL_PROTTREE_NONE='No athletes.'
    #-- Toolbar --
    LABEL_NEW='New'
    LABEL_NEW_ACTION='New     '
    LABEL_EDIT='Edit'
    LABEL_EDIT_ACTION='Edit     '
    LABEL_COPY='Copy'
    LABEL_COPY_ACTION='Copy   '
    LABEL_DELETE='Delete'
    LABEL_DELETE_ACTION='Delete    '
    LABEL_BROWSE='Browse'
    LABEL_BROWSE_ACTION='Browse '
    LABEL_STOP='Stop'
    LABEL_STOP_ACTION='Stop execution'
    LABEL_ANALYZE='Analyze Results'
    LABEL_TREE='Tree'
    LABEL_LIST='List'
    LABEL_REFRESH='Refresh'
    LABEL_DEFAULT = 'Default'
    LABEL_CONTINUE = 'Continue'
    LABEL_CONTINUE_ACTION = 'Approve continue'
    #-- Tabs --
    LABEL_DATA='Data'
    LABEL_SUMMARY='Summary'
    LABEL_INPUT='Input'
    LABEL_OUTPUT='Output'
    LABEL_METHODS='Methods'
    LABEL_LOGS_OUTPUT='Output Log'
    LABEL_LOGS_ERROR='Error Log'
    NO_INFO_SUMMARY='No summary information.'
    NO_INFO_METHODS='No methods information.'
    NO_INFO_LOGS='No logs information'
    NO_SAVE_SETTINGS='Error try to save settings.'
    
    # Protocol Form Template
    TITLE_NAME_RUN=' Protocol Run: '
    TITLE_RUN='Run'
    TITLE_RUN_NAME='Run name'
    VAR_RUN_NAME='runName'
    TITLE_LABEL='Label'
    LABEL_COMMENT='Describe your run here...'
    TITLE_COMMENT='Comment'
    TITLE_RUN_MODE='Run mode'
    VAR_RUN_MODE='runMode'
    LABEL_RUN_MODE_RESUME='resume'
    LABEL_RUN_MODE_RESTART='restart'
    TITLE_EXPERT='Expert Level'
    VAR_EXPERT='expertLevel'
    LABEL_EXPERT_NORMAL='Normal'
    LABEL_EXPERT_ADVANCE='Advanced'
    LABEL_EXPERT_EXPERT='Expert'
    TITLE_EXEC='Execution'
    TITLE_EXEC_HOST='Execution host'
    VAR_EXEC_HOST='hostName'
    TITLE_THREADS='Threads'
    VAR_THREADS='numberOfThreads'
    TITLE_MPI='MPI'
    VAR_MPI='numberOfMpi'
    TITLE_QUEUE='Launch to queue?'
    VAR_QUEUE='_useQueue'
    TITLE_BROWSE_DATA='Protocol data'
    LABEL_QUEUE_YES='Yes'
    LABEL_QUEUE_NO='No'
    LABEL_PARAM_YES='Yes'
    LABEL_PARAM_NO='No'
    LABEL_BUTTON_CLOSE='Close'
    LABEL_BUTTON_SAVE='Save'
    LABEL_BUTTON_EXEC='Execute'
    LABEL_BUTTON_VIS='Visualize'
    LABEL_BUTTON_WIZ='Wizard'
    LABEL_BUTTON_HELP='Help'
#    LABEL_BUTTON_VIS_ACTION='  Visualize  '
    LABEL_BUTTON_RETURN='Save'
    
    LABEL_PATTERN='Pattern'
    TEXT_PATTERN="The pattern (with wildcard expansion) of the files to import\n"+"For example: \n*data/particles/***.spi* \n or \n*~/Micrographs/mic/***.mrc*"
    ERROR_PATTERN_EMPTY='The *pattern* cannot be empty.'
    ERROR_PATTERN_FILES='There are no files matching the *pattern*'
    LABEL_CHECKSTACK='Check stack files?'
    LABEL_VOLTAGE='Microscope voltage (in kV)'
    LABEL_SPH_ABERRATION='Spherical aberration (in mm)'
    LABEL_AMPLITUDE='Amplitude Contrast'
    TEXT_AMPLITUDE='It should be a positive number, typically between 0.05 and 0.3.'
    
    LABEL_SAMP_MODE='Sampling rate mode'
    LABEL_SAMP_MODE_1='From image'
    LABEL_SAMP_MODE_2='From scanner'
    LABEL_SAMP_RATE='Sampling rate (A/px)'
    LABEL_MAGNI_RATE='Magnification rate'
    LABEL_SCANNED='Scanned pixel size'
    
    LABEL_IMPORT_VOL='Import volumes'
    ERROR_IMPORT_VOL='importVolumes:There is not filePaths matching pattern'
    
    LABEL_CTF_ESTI='CTF Estimation'
    LABEL_INPUT='Input'
    LABEL_INPUT_MIC='Input Micrographs'
    LABEL_INPUT_PART='Input Particles'
    LABEL_ALIG_PART='Write aligned particles?'
    TEXT_ALIG_PART='If set to *Yes*, the aligment will be applied to \n'+'input particles and a new aligned set will be created.'
    LABEL_LOW_RES='Lowest resolution'
    TEXT_LOW_RES='Give a value in digital frequency (i.e. between 0.0 and 0.5). '+'This cut-off prevents the typically peak at the center of the PSD '+'to interfere with CTF estimation. The default value is 0.05, but for '+'micrographs with a very fine sampling this may be lowered towards 0.0'
    LABEL_HIGH_RES='Highest resolution'
    TEXT_HIGH_RES='Give a value in digital frequency (i.e. between 0.0 and 0.5). '+'This cut-off prevents high-resolution terms where only noise exists '+'to interfere with CTF estimation. The default value is 0.35, but it should '+'be increased for micrographs with signals extending beyond this value. '+'However, if your micrographs extend further than 0.35, you should consider '+'sampling them at a finer rate.'
    LABEL_MIN_FOCUS='Minimum defocus to search (in microns)'
    TEXT_MIN_FOCUS='<Minimum defocus value (in microns) to include in defocus search. '+'Underfocus is represented by a positive number.'
    LABEL_MAX_FOCUS='Maximum defocus to search (in microns)'
    TEXT_MAX_FOCUS='Maximum defocus value (in microns) to include in defocus search. '+'Underfocus is represented by a positive number.'
    LABEL_WINDOW_SIZE='Window size'
    TEXT_WINDOW_SIZE='The PSD is estimated from small patches of this size. Bigger patches '+'allow identifying more details. However, since there are fewer windows, '+'estimations are noisier'
    
    TEXT_NO_INPUT_MIC='No *Input Micrographs* selected.'
    TEXT_NO_OUTPUT_CO='Output coordinates not ready yet.'
    ERROR_NO_EST_CTF='_estimateCTF should be implemented'
    
    
    TITLE_LAUNCHED='Success'
    LABEL_LAUNCHED='The protocol was launched successfuly.'
    LABEL_FOUND_ERROR='Errors found'
    TITLE_SAVED_FORM='Success'
    LABEL_SAVED_FORM='The protocol was saved successfuly.'
    TITLE_DELETE_FORM='Confirm DELETE'
    LABEL_DELETE_FORM='*ALL DATA* related to this _protocol run_ will be *DELETED*.\nDo you really want to continue?'    
    TITLE_STOP_FORM='Confirm STOP'
    LABEL_STOP_FORM='Do you really want to *STOP* this run?'
    
    NO_VIEWER_FOUND='There is not viewer for protocol:'
    
    # Hosts Template
    LABEL_HOSTS='Hosts'
    LABEL_HOSTS_ACTION=' Hosts '
    TITLE_HOST_DELETE_FORM='Confirm DELETE'
    LABEL_HOST_DELETE_FORM='*ALL CONFIGURATION* related to this _host_ will be *DELETED*.\nDo you really want to continue?'    
    BUTTON_TEST='Test configuration'
    TITLE_HOST_SAVE_FORM='Action SAVE'
    LABEL_HOST_SAVE_FORM_SUCESS='Host sucessfully saved. \n'
    LABEL_HOST_SAVE_FORM_FAIL='Error saving host. \n'
    LABEL_HOST_SAVE_FORM='*ALL CONFIGURATION* related to this _host_ will be *DELETED*.\nDo you really want to continue?'    
    LABEL_HOST_CONECTION_SUCESS='Connection SUCCEEDED to host'
    LABEL_HOST_CONECTION_FAIL='Connection FAILED to host'
    LABEL_HOST_CONECTION_TEST='Connection TEST'
    
    #SHOWJ_WEB
    SHOWJ_TITLE='Showj'
    
    LABEL_RESLICE='Reslice'
    RESLICE_Z_NEGATIVE='Z Negative (Front)'
    RESLICE_Y_NEGATIVE='Y Negative (Top)'
    RESLICE_Y_POSITIVE='Y Positive (Bottom)'
    RESLICE_X_NEGATIVE='X Negative (Left)'
    RESLICE_X_POSITIVE='X Positive (Right)'
    
    LABEL_COLS='Cols'
    LABEL_ROWS='Rows'
    
    LABEL_MIRRORY='Invert Y Axis'
    LABEL_APPLY_TRANSFORM='Apply Transform Matrix'
    LABEL_ONLY_SHIFTS='Only Shifts'
    LABEL_WRAP='Wrap'
    
    LABEL_BLOCK_SELECTION='Select Block'
    LABEL_LABEL_SELECTION='Select Label'
    LABEL_VOLUME_SELECTION='Select Volume'
    
    LABEL_ENABLE='Enable'
    LABEL_DISABLE='Disable'
    LABEL_SELECT_ALL='Select All'
    LABEL_SELECT_FROM='Select From Here'
    LABEL_SELECT_TO='Select To Here'
    
    LABEL_DISPLAY_TABLE_CONFIG='Display Table Configuration'
    
    LABEL_LABEL='Label'
    LABEL_VISIBLE='Visible'
    LABEL_RENDER='Render'
    
    LABEL_BUTTON_OK='Ok'
    LABEL_BUTTON_CANCEL='Cancel'
    
    LABEL_THRESHOLD='Threshold:'
    
    ERROR_DIMENSIONS='Incorrect table width or height: '
    ERROR_WEBGL='Your web browser does not support or is not configured for WebGL. See [[http://get.webgl.org/][WebGL Support]] for more information.'
    
class Icon():
    # Project Content Template
    RUNS_TREE = 'fa-sitemap.png'
    RUNS_LIST = 'fa-bars.png'
    ACTION_NEW = 'new_object.gif',
    ACTION_EDIT = 'fa-pencil.png'
    ACTION_COPY = 'fa-files-o.png'
    ACTION_DELETE = 'fa-trash-o.png'
    ACTION_REFRESH = 'fa-refresh.png'
    ACTION_STEPS = 'fa-folder-open.png'
    ACTION_TREE = None
    ACTION_LIST = 'fa-bars.png'
    ACTION_STOP = 'fa-stop.png'
    ACTION_CONTINUE = 'fa-play-circle-o.png'
    ACTION_RESULTS = 'fa-eye.png'
    ACTION_CLOSE = 'fa-times.png'
    ACTION_SAVE = 'fa-save.png'
    ACTION_VISUZALIZE = 'fa-eye.png'
    ACTION_WIZ = 'fa-magic.png'
    ACTION_HELP = 'fa-question-circle.png'
    ACTION_REFERENCES = 'fa-external-link.png'
    
    ACTION_SEARCH = 'fa-search.png'
    
    #Host template
    BUTTON_CLOSE = 'dialog_close.png'
    BUTTON_SAVE = 'filesave.png'
