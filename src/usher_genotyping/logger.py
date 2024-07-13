import logging
import os
import sys

def setup_logger():
    """
    set up the logger
    log file is either taken from ENV var USHER_GENOTYPING_LOGFILE or will be 
    set to file 'usher_genotyping.log' in current dir
    
    Once setup, the logger is global, so this function does not return anything
    """
    if os.getenv('USHER_GENOTYPING_LOGFILE'):
        logfile = os.getenv('USHER_GENOTYPING_LOGFILE')
    else:
        logfile = 'usher_genotyping.log'

    logging.basicConfig(
        filename= logfile,
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%d/%m/%Y %H:%M:%S')
    
def log( level:str, msg:str ):
    """
    If we have a logger attached, write to the logger, otherwise print the msg to stdout.
    
    Parameters:
        level: str, required
          one of these allowed log levels: 
            - critical
            - error
            - warning
            - info
            - debug
        msg: str, required
          the log message. For 'critical' level, the message can be an empty string in which case a default
          message is provided.
    """
    if not msg:
        if level =='critical':
            msg = '!!! Code was terminated due to an EXCEPTION! If a traceback is available, it will be shown here'
        else:
            msg = '(no log message provided)'
        
    try:
        if sys.exc_info() != (None, None, None):
            getattr( logging, level)( msg, exc_info = sys.exc_info() )
        else:
            getattr( logging, level)( msg )
    except:
        print( msg )