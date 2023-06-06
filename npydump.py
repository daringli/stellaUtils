#!/usr/bin/env python

import sys
import numpy as np
if len(sys.argv[1]) > 1:
    for f in sys.argv[1:]:
        print(np.load(f))
