#!/usr/bin/env python3

# win32com examples:
# https://mhammond.github.io/pywin32/html/com/win32com/HTML/docindex.html
import win32com.client


def check_dependencies():
    """Check environment dependances.  Might need these on install, not here"""
    import sys
    assert sys.platform == 'win32'
    
    # http://timgolden.me.uk/pywin32-docs/html/com/win32com/HTML/QuickStartClientCom.html
    # Use makepy.py -i to poke around in what might be useful
    try:
        # 'ASCOM Master Interfaces for .NET and COM' constants.
        # Use example: win32com.client.constants.shutterOpen
        win32com.client.gencache.EnsureModule('{76618F90-032F-4424-A680-802467A55742}', 0, 1, 0)
    except:
        log.info('ASCOM does not seem to be installed.  MaxIm/telescope control will not work.')
        return False
    try:
        # MaxIm constants.  The long string is the GUID of MaxIm found by makepy.py
        win32com.client.gencache.EnsureModule('{B4955EC7-F7F2-11D2-AA9C-444553540000}', 0, 1, 0)
    except:
        log.info('MaxIm not found.  MaxIm/telescope control will not work.')
        return False


DEFAULT_EXPTIME = 1
DEFAULT_FILT = 0
DEFAULT_CENT_TOL = 5   # Pixels
DEFAULT_GUIDER_EXPTIME = 1 # s 1 for night, 0.2 for day

# In principle, it is possible to use only MaxIm guiding stuff for
# this, alleviating the need for us to connect directly to the
# telescope.  In practice, with a GEM and for setting DEC conveniently
# in the guider, MaxImControl really needs to be connected to the scope.
# Since unexpected behavior may result if there is not a hard failure
# when a GEM is not connected, indicate here whether or not you want
# that hard failure
TELESCOPE_MUST_BE_CONNECTABLE = True

# These are necessary for GEMs because MaxIm does not reveal the
# contents of the Camera Control -> Guide Tab -> Settings dialog box
# -> Advanced Tab -> Guider Motor Control radio buttons to the
# scripting interface.  As explained in guider_motor_reverse_setup,
# when using ACP or otherwise not having MaxIm connected to the
# telescope, we need to manage motor reversal ourselves
GUIDER_MOTOR_CONTROL_REVERSE_X = True
GUIDER_MOTOR_CONTROL_REVERSE_Y = False
# Misalignment in deg
GUIDER_CAL_ASTROMETRY_MAX_MISALIGNMENT = 10

# --> hack
HORIZON_LIMIT = 8.5

# --> hack
# --> I may improve this location or the technique of message passing
hostname = socket.gethostname()
if hostname == "greyhound" or hostname == "gigabyte":
    RAW_DATA_ROOT = r'\\snipe\data\io\IoIO\raw'
    DEFAULT_TELESCOPE = 'ScopeSim.Telescope'
elif socket.gethostname() == "IoIO1U1":
    RAW_DATA_ROOT = r'C:\Users\PLANETARY SCIENCE\Desktop\IoIO\data'
    # --> Eventually, it would be nice to have this in a chooser
    DEFAULT_TELESCOPE = 'AstroPhysicsV2.Telescope'
# For weather synchronization with ACP
ACPUtil = 'ACP.Util'
DEFAULT_GUIDE_BOX_COMMAND_FILE = os.path.join(RAW_DATA_ROOT, 'GuideBoxCommand.txt')
DEFAULT_GUIDE_BOX_LOG_FILE = os.path.join(RAW_DATA_ROOT, 'GuideBoxLog.txt')

RUN_LEVEL_MAIN_ASTROMETRY = os.path.join(
    raw_data_root, '2020-09_Astrometry/Main_Astrometry_East_of_Pier.fit')
    #raw_data_root, '2020-03_Astrometry/Main_Astrometry_East_of_Pier.fit')
    #raw_data_root, '2019-04_Astrometry/Main_Astrometry_East_of_Pier.fit')
    #raw_data_root, '2019-04_Astrometry/Main_Astrometry_West_of_Pier.fit')
    #raw_data_root, '2019-02_Astrometry/PinPointSolutionEastofPier.fit')
    #raw_data_root, '2019-02_Astrometry/PinPointSolutionWestofPier.fit')
    #raw_data_root, '2018-04_Astrometry/PinPointSolutionEastofPier.fit')

RUN_LEVEL_GUIDER_ASTROMETRY = os.path.join(
    raw_data_root, '2020-09_Astrometry/Guider_Astrometry_East_of_Pier.fit')
    #raw_data_root, '2020-03_Astrometry/Guider_Astrometry_East_of_Pier.fit')
    #raw_data_root, '2019-04_Astrometry/Guider_Astrometry_West_of_Pier.fit')
    #raw_data_root, '2019-02_Astrometry/GuiderPinPointSolutionWestofPier.fit')
    #raw_data_root, '2019-02_Astrometry/GuiderPinPointSolutionEastofPier.fit')    
    #raw_data_root, '2018-04_Astrometry/GuiderPinPointSolutionWestofPier.fit')
    #raw_data_root, '2018-01_Astrometry//GuiderPinPointSolutionEastofPier.fit')

class Astrometries:
    # --> Development left off here.  I might want to take some hints
    # --> from progess I made in IoIO.photometry and friends
    def __init__(self,
                 
