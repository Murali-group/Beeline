"""
This file contains the functions I use for submitting jobs to the baobab cluster
More info about the baobab cluster is available here: https://github.com/Murali-group/utils/wiki/Baobab-info
"""

import os
import re
import socket
import subprocess
#from tqdm import tqdm


def copyToBaobabNodes(file_to_copy):
    """ 
    Copy a given file to each of the 6 baobab nodes. 
    Keeps the same filepath, but replaces /data to /localdisk which is unique to each of the nodes.
    Useful to speed up reading input files if many jobs are using the same file (such as an interactome)
    """
    print("Copying %s to the 6 baobab nodes" % (file_to_copy))
    print("\tTip: If you are being asked for a password, setup ssh keys to allow passwordless ssh and scp")
    ip_template = "192.168.200.%s"
    # use the localdisk storage on the baobab nodes to speed up read/write times
    copy_to_dir = re.sub("^/data", "/localdisk", os.path.dirname(os.path.abspath(file_to_copy))) 
    # loop through the 6 nodes
    #for i in tqdm(range(1,7)):
    for i in range(1,7):
        command = "ssh %s 'mkdir -p %s'; scp %s %s:%s" % (ip_template%i, copy_to_dir, os.path.abspath(file_to_copy), ip_template%i, copy_to_dir)
        runCommandOnBaobab(command)


def runCommandOnBaobab(command):
    """ 
    Run a given command on baobab. 
    Useful for submitting a job to the baobab cluster or copying a file from baobab to each of the nodes 
    """
    print("Running: %s" % (command))
    if 'baobab' in socket.gethostname():
        subprocess.check_call(command, shell=True)
    else:
        command = "ssh -t baobab.cbb.lan \"%s\"" % (command)
        subprocess.check_call(command, shell=True)


def submitQsubFile(qsub_file):
    """ Submit a qsub file to the baobab cluster using the qsub command
    """
    command = "qsub " + qsub_file
    runCommandOnBaobab(command)


def writeQsubFile(jobs, qsub_file, submit=False, log_file=None, err_log_file=None, name=None, nodes=1, ppn=24, walltime='10:00:00'):
    """ Function to write a qsub bash script which can be submitted to the baobab cluster.
    *jobs*: a list or set of commands/jobs to be run in this job. 
            This could include a 'cd' to move to your project directory, or 'export PYTHONPATH' to setup environmental variables for example
    *qsub_file*: path/to/file to write. Really just a bash script with special headers recognized by PBS. 
    *submit*: option to submit the written qsub file to the baobab cluster
    *log_file*: file which will contain the stdout output of the submitted qsub file. If None, -out.log will be appended to qsub_file
    *err_log_file*: file which will contain the stderr output of the submitted qsub file. If None, -err.log will be appended to qsub_file
    *name*: name to give the job 
    *nodes*: total number of nodes you need
    *ppn*: processors per node that you will need. Max is 24
    *walltime*: amount of time your job will be allowed before being forcefully removed. 'HH:MM:SS'
    """
    if log_file is None:
        std_out_log = "%s-out.log" % qsub_file
        std_err_log = "%s-err.log" % qsub_file
    else:
        std_out_log = log_file
        std_err_log = log_file
    # start the qsub file 
    with open(qsub_file, 'w') as out:
        out.write('#PBS -l nodes=%d:ppn=%d,walltime=%s\n' % (nodes, ppn, walltime))
        # set the job name
        out.write('#PBS -N %s\n' % (name))
        out.write('#PBS -o %s\n' % (std_out_log))
        out.write('#PBS -e %s\n' % (std_err_log))
        out.write('#####################################################################################\n')
        out.write('echo "Job Started at: `date`"\n')
        # write each job, as well as an echo (print) statement of the job to be run to the qsub file
        out.write('\n'.join(['echo """%s"""\n%s' % (cmd, cmd) for cmd in jobs]) + '\n')
        out.write('echo "Job Ended at: `date`"\n')
        # TODO some kind of email or notification if any of the jobs failed

    if submit:
        submitQsubFile(qsub_file)
