import os


def get_output_path(runner_obj, algo_name):
    return "outputs/"+"/".join(str(runner_obj.inputDir).split("inputs" + os.sep)[1].split(os.sep))+algo_name
