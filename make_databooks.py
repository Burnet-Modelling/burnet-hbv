import atomica as at
import sciris as sc
import pandas as pd
import numpy as np
import hbv_vimc as hbv
from hbv_vimc.utils import databook_gen 
    

if __name__ == '__main__':
   #__spec__=None
   sc.parallelize(databook_gen, hbv.missing.keys(), die=False)



