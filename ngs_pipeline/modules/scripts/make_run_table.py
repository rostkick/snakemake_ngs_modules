#%%
import pandas as pd
from glob import glob
import os
import argparse
import re
# %%
path_to_grm = '/mnt/rskitchenko/projects/medulo_part3/data/grm'
path_to_tmr = '/mnt/rskitchenko/projects/medulo_part3/data/tmr'
# %%
g_grm = glob(path_to_grm + '/*.fastq.gz')
g_tmr = glob(path_to_tmr + '/*.fastq.gz')
# %%
l_grm = [i for i in g_grm]
l_tmr = [i for i in g_tmr]
# %%
df = pd.DataFrame({'fastq': l_grm + l_tmr})
df
# %%
df['base_fastq'] = df['fastq'].apply(lambda x: os.path.basename(x))
# %%
pattern = r'^[^_]+_[^_]+'
pattern = '(?P<sample>' + pattern + ')'
print(pattern)
df['base_fastq'].str.extractall(pattern).droplevel(1)
# %%
