import sys


def separator():
    if sys.platform.startswith("linux"):
        return "/"
    elif sys.platform.startswith("win"):
        return "\\"
    else:  # Leave optional cases for os (darwin) or cygwin (cygwin)
        return "/"


def get_output_path(runner_obj, algo_name):
    sep = separator()
    return "outputs/"+"/".join(str(runner_obj.inputDir).split("inputs" + sep)[1].split(sep))+algo_name
