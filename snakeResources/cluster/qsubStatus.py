#!/ifs/scratch/c2b2/ac_lab/jk3755/atac/conda/envs/snakemake/bin/python3
import sys
import re
import subprocess
from subprocess import Popen, PIPE
import xml
import xml.etree.cElementTree as ET

# setup
jobid = sys.argv[1]
regex = re.compile(r' +')

# checking for running job with qstat
try:
    p = Popen(['qstat', '-xml', '-j', jobid], stdout=PIPE)
    output, err = p.communicate()
    root = ET.ElementTree(ET.fromstring(output.decode())).getroot()
    print("running")
except KeyboardInterrupt:
    print("failed")
except (subprocess.CalledProcessError, xml.etree.ElementTree.ParseError) as e:
    # if not running, use qacct to check for exit status
    try:
        p = Popen(['qacct', '-j', jobid], stdout=PIPE)
        output, err = p.communicate()
        for x in output.decode().split('\n'):
            y = regex.split(x.strip())
            if y[0] == 'exit_status':
                if y[1] == '0':
                    print("success")
                else:
                    print("failed")
                break
    except (subprocess.CalledProcessError, IndexError, KeyboardInterrupt) as e:
        print("failed")