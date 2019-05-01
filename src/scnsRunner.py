import os
from pathlib import Path
import pandas as pd

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for SCNS.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    
    if not RunnerObj.inputDir.joinpath("SCNS").exists():
        print("Input folder for SCNS does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("SCNS").mkdir(exist_ok = False)
        
    if not RunnerObj.inputDir.joinpath("SCNS/ExpressionData.csv").exists():
        # input file
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)

        # Convert input expression to boolean
        # If  the gene's expression value is >= it's avg. expression across cells
        # it receieves a "True", else "False"
        BinExpression = ExpressionData.T >= ExpressionData.mean(axis = 'columns')
        BinExpression.drop_duplicates(inplace= True)
        # Write unique cells x genes output to a file
        BinExpression.to_csv(RunnerObj.inputDir.joinpath("SCNS/ExpressionData.csv")) 

        # Read PseudoTime file to figure out
        # initial and final states
        PseudoTimeDF = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                                              header = 0, index_col = 0)
        # Get the Time corresponding to cells in BinExpression dataframe
        # Identify cells in initial and final states from them
        # cells in initial states are the ones with the earliest time pts (<10th percentile)
        # cells in target states are the ones with the latest time pts (>90th percentile)
                                   
        StateDF = PseudoTimeDF.loc[BinExpression.index, 'Time']
        initialStates =  open(RunnerObj.inputDir.joinpath("SCNS/initial.txt"),'w')                           
        for ix in StateDF[StateDF <= StateDF.quantile(0.1)].index:
               initialStates.write(ix +'\n')
        initialStates.close()
                                   
        targetStates =  open(RunnerObj.inputDir.joinpath("SCNS/target.txt"),'w')                           
        for ix in StateDF[StateDF >= StateDF.quantile(0.9)].index:
               targetStates.write(ix +'\n')
        targetStates.close()
        
        
        parameters =  open(RunnerObj.inputDir.joinpath("SCNS/Parameters.csv"),'w')  
        parameters.write('Gene,MaxActivators,MaxRepressors,Threshold\n')
        for cols in BinExpression.columns:
               parameters.write(cols+',3,3,75'+'\n')
        parameters.close()
        
        # generate Edges file
        # two cells are connected by an edge
        # if they only differ in the boolean 
        # expression of one gene
        States =  open(RunnerObj.inputDir.joinpath("SCNS/Edges.csv"),'w') 
        States.write('Gene,StateA,StateB\n')
        colNames =  BinExpression.columns
        for idx_i, row_i in BinExpression.iterrows():
            for idx_j, row_j in BinExpression.iterrows():
                if list(row_i == row_j).count(False) == 1:
                    States.write(str(colNames[list(row_i == row_j).index(False)])+','
                                 +str(idx_i)+','+str(idx_j)+'\n')
        States.close()

def run(RunnerObj):
    '''
    Function to run SCNS algorithm
    '''

    inputPath = "data/" + str(RunnerObj.inputDir).split("RNMethods/")[1] + "/SCNS/"
                    
    
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SCNS/"
    os.makedirs(outDir, exist_ok = True)
    
    outPath = "data/" +  str(outDir)
    cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/SCNS-Toolkit/SynthesisEngine/data/', 
                         'scns:base /bin/sh -c \"time -v -o', "data/" + str(outDir) + 'time.txt',
                         'mono SynthesisEngine.exe', inputPath+'ExpressionData.csv',
                          inputPath+'Edges.csv',  inputPath+'Parameters.csv',
                          inputPath+'initial.txt',  inputPath+'target.txt',
                          outPath, '\"'])

    print(cmdToRun)
    os.system(cmdToRun)
                                   
def parseOutput(Runner):
    '''
    Function to parse output from SCNS
    '''
    print("HA!")
    