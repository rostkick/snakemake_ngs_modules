import sys


from modules.scripts.functions import *
from modules.scripts.params_builder import NGSData

args = sys.argv
config_path = os.path.abspath(args[args.index("--configfile") + 1])


ngs = NGSData(config['grm_dir'], config['tmr_dir'])
print(ngs.long_df)

