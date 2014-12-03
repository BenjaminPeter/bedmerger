import logging

from bedmerger.parameters import Parameters
from bedmerger import merge


def setup_logging(logfile, loglevel):
    """ initializes basic logging"""
    
    if logfile is not None:
        logging.basicConfig(filename=logfile, level=loglevel)
        logging.debug("set logfile to %s, loglevel to %s",
                logfile, loglevel)
    else:
        logging.basicConfig(level=loglevel)
        logging.debug("no logfile set, loglevel to %s",
                loglevel)
    
    return None


if __name__ == "__main__":
    params = Parameters.from_command_line()
    params.setup()
    params.sanity_checks()
    logger = setup_logging(params.logfile, params.loglevel)
    merge.merge(params)
