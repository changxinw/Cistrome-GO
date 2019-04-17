import os
import sys
from corelib import *

def optValidate(optparser):
    """
    :param optparser:
    :return: Validated option object
    """
    options = optparser
    if not os.path.exists(options.bed):
        Info("ERROR: Bed file not exist.")
        sys.exit(1)
    if options.peakn == "all":
        pass
    if options.peakn.endswith("%"):
        options.peakn = float(options.peakn[:-1]) / 100
    options.prefix = "%s/%s" % (options.output, options.name)
    if not os.path.exists(options.output):
        os.makedirs(options.output)
    try:
        if options.expr:
            pass
        else:
            options.expr = False
    except:
        options.expr = False
    return options
