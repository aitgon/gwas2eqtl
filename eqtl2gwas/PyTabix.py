import subprocess
# from subprocess import Popen, PIPE

def bgzip(filename):
    """Call bgzip to compress a file."""
    subprocess.run(['bgzip', '-f', filename])

def tabix_index(filename, preset=None, chrom=None, start=None, end=None, skip=0, comment="#"):
    command_lst = ['tabix', '-S', str(skip), '-c', comment]
    if not preset is None:
        command_lst.append('-p')
        command_lst.append(preset)
    if not chrom is None:
        command_lst.append('-s')
        command_lst.append(str(chrom))
    if not start is None:
        command_lst.append('-b')
        command_lst.append(str(start))
    if not end is None:
        command_lst.append('-e')
        command_lst.append(str(end))
    command_lst.append(filename)
    """Call tabix to create an index for a bgzip-compressed file."""
    subprocess.run(command_lst)

def tabix_query(filename, chrom, start, end):
    """Call tabix and generate an array of strings for each line it returns."""
    query = '{}:{}-{}'.format(chrom, start, end)
    process = Popen(['tabix', '-f', filename, query], stdout=PIPE)
    for line in process.stdout:
        yield line.strip().split()
