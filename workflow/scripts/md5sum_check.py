import hashlib
import os
import pandas
import pathlib
import sys

#%%############################################################################
#
# arguments
#
###############################################################################

try:
    query_path = sys.argv[1]
    md5sum_txt_path = sys.argv[2]
    if len(sys.argv) > 4:
        print("""Two many arguments!""")
        sys.exit(1)
except IndexError:
    print("""Argument missing!""")
    sys.exit(1)

#%%
md5sum_df = pandas.read_csv(md5sum_txt_path, delim_whitespace=True, names=['md5', 'filename'])
filename = os.path.basename(query_path)

#%%
md5_server = md5sum_df.loc[md5sum_df['filename'] == filename, 'md5'].values[0]
md5_local = hashlib.md5(pathlib.Path(query_path).read_bytes()).hexdigest()

if md5_local != md5_server:
    raise ValueError('Different MD5 between server and local file: ' + query_path)
