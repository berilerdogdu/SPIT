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

def call_hclust(args):
    try:
        import rpy2
    except ImportError:
        subprocess.run(['pip', 'install', 'rpy2'])
    # List of R libraries to install and load
    library_names = ["tidyverse", "ggplot2", "tidyr", "grid", "RColorBrewer", "gplots", "viridis"]

    for library_name in library_names:
        if not r_package_installed(library_name):
            robjects.r(f'install.packages("{library_name}")')
        robjects.r(f'library("{library_name}", character.only = TRUE)')
    
    directory_path = os.path.join(args.O, "SPIT_analysis")
    robjects.r.assign("directory_path", directory_path)
    current_script_dir = os.path.dirname(os.path.abspath(__file__))
    hclust_path = os.path.join(current_script_dir, 'r_scripts/hclust.R')
    robjects.r.source(hclust_path)
    perform_hclust = robjects.r['perform_hclust']
    perform_hclust(args.l, args.include_shared_dtu, args.color_palette, args.color_covariate)
