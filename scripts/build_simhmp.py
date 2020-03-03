
#!/usr/bin/env python


from subprocess import Popen, PIPE
import pandas as pd
import sys
import argparse
import numpy as np

def get_args():
    parser=argparse.ArgumentParser(description="""Program description""")
    parser.add_argument("infile",type=str,help="""The input vcf file""")
    parser.add_argument("outfile",type=str,help="""The output vcf file""")
    parser.add_argument("markerfile",type=str,help=""""File with list of marker postions""")
    args = parser.parse_args()
    return args





if __name__ == "__main__":
    build_hmp()
