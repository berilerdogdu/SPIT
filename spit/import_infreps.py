import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import os
import subprocess

def r_package_installed(package_name):
    try:
        importr(package_name)
        return True
    except:
        return False

def call_import_infreps(args):
    try:
        import rpy2
    except ImportError:
        subprocess.run(['pip', 'install', 'rpy2'])
    # List of R libraries to install and load
    library_names = ["tidyverse", "progress"]
    
    if not r_package_installed("tximport"):
        if not r_package_installed("BiocManager"):
            robjects.r('install.packages("BiocManager")')
        robjects.r('BiocManager::install("tximport")')
    robjects.r('library("tximport")')
    
    for library_name in library_names:
        if not r_package_installed(library_name):
            robjects.r(f'install.packages("{library_name}")')
        robjects.r(f'library("{library_name}", character.only = TRUE)')
    directory_path = os.path.join(args.O, "SPIT_analysis")
    infReps_directory_path = os.path.join(directory_path, "infReps")
    if os.path.exists(infReps_directory_path) == False:
        os.mkdir(infReps_directory_path)
    robjects.r.assign("directory_path", directory_path)
    current_script_dir = os.path.dirname(os.path.abspath(__file__))
    tximport_path = os.path.join(current_script_dir, 'r_scripts/tximport_quant_files.R')
    robjects.r.source(tximport_path)
    import_infreps = robjects.r['import_infreps']
    import_infreps(args.l, args.m, args.quant_type, args.quant_path, args.i)
