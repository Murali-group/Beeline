import os
import subprocess
import pandas as pd
from pathlib import Path
import numpy as np
import socket
from src.utils import baobab_utils


# these must be in the right order!
params_order = [
    'lambda', 'dT', 'num_lags', 'kernel_width',
    'prob_zero_removal', 'prob_remove_samples',
    'family', 'num_replicates',
]
default_params = {
    'lambda': '0.01',
    'dT': '10',
    'num_lags': '5',
    'kernel_width': '4',
    'prob_zero_removal': '0',
    'prob_remove_samples': '0.2',
    'family': 'gaussian',
    'num_replicates': '2',
}


def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for SCINGE.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("SCINGE").exists():
        print("Input folder for SCINGE does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("SCINGE").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("SCINGE/ExpressionData.csv").exists():
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
        newExpressionData = ExpressionData.T.copy()
        PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
        newExpressionData['PseudoTime'] = PTData['PseudoTime']
        newExpressionData.to_csv(RunnerObj.inputDir.joinpath("SCINGE/ExpressionData.csv"),
                             sep = ',', header  = True, index = False)

    setup_params(RunnerObj)


def setup_params(RunnerObj):
    # if the parameters aren't specified, then use default parameters
    # TODO allow passing in multiple sets of hyperparameters
    params = RunnerObj.params
    for param, val in default_params.items():
        if param not in params:
            params[param] = val
    # check if specific combinations of dT and num_lags are specified
    if 'dT_num_lags' in params:
        params['dT'], params['num_lags'] = params['dT_num_lags'].split(',')
    params_str = build_params_str(params)
    RunnerObj.params_str = params_str
    RunnerObj.params = params

    # the final file is written here:
    outDir = str(RunnerObj.inputDir).replace("inputs/","outputs/")+"/SCINGE/"
    RunnerObj.final_ranked_edges = "%s/%s-rankedEdges.csv" % (outDir, RunnerObj.params_str)


def build_params_str(params):
    params_str = "l%s-dT%s-nl%s-kw%s-pz%s-pr%s-nr%s" % (
        params['lambda'], params['dT'], params['num_lags'], 
        params['kernel_width'], params['prob_zero_removal'], 
        params['prob_remove_samples'], params['num_replicates'],
        )
    return params_str


def run(RunnerObj):
    '''
    Function to run SCINGE algorithm
    '''
    # input path on docker
    inputPath = "data/" + str(RunnerObj.inputDir).split("RNMethods/")[1] + \
                    "/SCINGE/ExpressionData.csv"

    params = RunnerObj.params
    params_str_to_run = ' '.join(str(params[p]) for p in params_order) 

    # make output dirs if they do not exist:
    outDir = str(RunnerObj.inputDir).replace("inputs/","outputs/")+"/SCINGE/"
    os.makedirs(outDir, exist_ok = True)

    # if this has already been run, then skip it
    if params.get('forced') is False and \
            os.path.isfile(RunnerObj.final_ranked_edges):
        print("%s already exists. Set forced=True to overwrite" % (RunnerObj.final_ranked_edges))
        return 'already_exists'

    # also add the parameters to the output dir
    outDir += RunnerObj.params_str

    if params.get('docker') is True:
        os.makedirs(outDir, exist_ok = True)
        outPath = "data/" + str(outDir) 
        cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/runSCINGE/data/ scinge:base /bin/sh -c \"time -v -o', "data/" + str(outDir) + 'time.txt', './runSCINGE ',
                            inputPath, outPath, params_str_to_run, '\"'])
    else:
        outPath = os.path.abspath(outDir) 
        if 'baobab' in socket.gethostname():
            # write the localdisk to speed up file IO
            outPath = outPath.replace("/data/","/localdisk/")
        os.makedirs(outPath, exist_ok = True)
        scinge_path = "Algorithms/SCINGE/runSCINGE/runSCINGE"
        # input path
        input_path = "%s/SCINGE/ExpressionData.csv " % (str(RunnerObj.inputDir)) 
        scinge_command = "time %s %s" % (
            ' '.join(os.path.abspath(p) for p in [scinge_path, input_path, outPath]), params_str)
        ml_lib_command = create_ml_lib_command()
        # if this is on a cluser (e.g., baobab), write a qsub file and submit the job
        if 'qsub' in params and params['qsub'] is True:
            qsub_file = "%s/cmd.qsub" % (outDir)
            name = "scinge-%s" % (RunnerObj.params_str)
            jobs = [ml_lib_command, scinge_command]
            # TODO make the nodes, ppn and walltime parameters(?)
            baobab_utils.writeQsubFile(
                jobs, qsub_file, name=name, nodes=1, ppn=1, walltime='10:00:00')
            print(scinge_command)
            cmdToRun = "qsub %s" % (qsub_file)
        else:
            # move to the output dir because each run makes some temporary files that will be copied if we stay in the base dir
            curr_dir = os.getcwd()
            os.chdir(outPath)
            # otherwise just run like normal
            cmdToRun = "%s\n%s" % (ml_lib_command, scinge_command)

    print(cmdToRun)
    # also print the parameters
    print("\tParameters: %s" % (', '.join("%s: %s" % (p, str(params[p])) for p in params_order)))

    #os.system(cmdToRun)
    subprocess.check_call(cmdToRun, shell=True)

    # move back to the base dir
    os.chdir(curr_dir)
    # delete the copy of the expression data
    if os.path.isfile("%s/ExpressionData.mat" % (outPath)):
        os.remove("%s/ExpressionData.mat" % (outPath))


def create_ml_lib_command():
    # these are the bash commands I had to run for the ML binary to run
    #MCROOT="/data/jeff-law/tools/matlab/matlab-v96/v96";
    #LD_LIBRARY_PATH="$MCROOT/runtime/glnxa64:$MCROOT/bin/glnxa64:$MCROOT/sys/os/glnxa64/:$MCROOT/sys/opengl/lib/glnxa64:"; export LD_LIBRARY_PATH;
    #LD_PRELOAD="$MCROOT/sys/os/glnxa64/libstdc++.so.6"; export LD_PRELOAD
    command = '\n'.join([
        'MCROOT="/data/jeff-law/tools/matlab/matlab-v96/v96";',
        'LD_LIBRARY_PATH="$MCROOT/runtime/glnxa64:$MCROOT/bin/glnxa64:$MCROOT/sys/os/glnxa64/:$MCROOT/sys/opengl/lib/glnxa64:"; export LD_LIBRARY_PATH;',
        'LD_PRELOAD="$MCROOT/sys/os/glnxa64/libstdc++.so.6"; export LD_PRELOAD;',
        ])
    return command


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SCINGE.
    '''
    # Quit if output directory does not exist
    outDir = str(RunnerObj.inputDir).replace("inputs/","outputs/")+"/SCINGE/"
    outPath = "%s/%s/" % (outDir, RunnerObj.params_str)
    if 'baobab' in socket.gethostname():
        # write the localdisk to speed up file IO
        outPath = outPath.replace("/data/","/localdisk/")
    out_file = outPath+'SCINGE_Ranked_Edge_List.txt'
    if not Path(out_file).exists():
        print(out_file + ' does not exist, skipping...')
        return

    # Read output
    OutDF = pd.read_csv(out_file, sep = '\t', header = 0)

    # write this to the normal directory
    print("\twriting processed rankedEdges to %s" % (RunnerObj.final_ranked_edges))
    with open(RunnerObj.final_ranked_edges,'w') as out:
        out.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

        for idx, row in OutDF.iterrows():
            out.write('\t'.join([row['Regulator'],row['Target'],str(row['SCINGE_Score'])])+'\n')

    if os.path.isdir(outPath) and RunnerObj.params.get('cleanup') is True:
        command = "rm -r %s" % (outPath)
        print("\tremoving temp files (%s)" % (command))
        os.system(command)

    print("-"*50)
    print("")

