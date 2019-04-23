import os
import matplotlib.pyplot as plt


class TimerObj:
    def __init__(self, algo_name, dataset, time_path):
        self.dataset = dataset
        self.time_path = time_path+"/time.txt"
        self.algo_name = algo_name


    def parse_time(self):
        """
       gets the user time given by the time command, which is stored along with the output when the algorithm runs,
       in a file called time.txt
       :return:
       float containing the time this object took to run on the dataset
        """
        try:
            with open(self.time_path, "r+") as f:
                lines = f.readlines()
                line = lines[1]
                time_val = float(line.split()[-1])

        except FileNotFoundError:
            print("Time file not found, setting time value to 0")
            time_val = 0.0
        except ValueError:
            print("Algorithm running failed, setting time value to 0")
            time_val = 0.0

        return time_val
    def get_time(self):
        try:
            with open(self.time_path, "r+") as f:
                lines = f.readlines()
                line = lines[1]
                time_val = float(line.split()[-1])

        except FileNotFoundError:
            print("Time file not found, setting time value to 0")
            time_val = 0.0
        except ValueError:
            print("Algorithm running failed, setting time value to 0")
            time_val = 0.0

        return time_val

    def get_dataset_name(self):
        return self.dataset

    def _get_time_path(self):
        return self.time_path

    def _get_algo_name(self):
        return self.algo_name


    def make_csv(self, algo_times):
        pass

    def plot_times(self, csv_file):
        # labels = []
        # vals = []
        # for dir in dirs:
        #     if os.path.isdir(path+dir):
        #         timefile = path+dir+"/time.txt"
        #         t = get_time(timefile)
        #         labels.append(dir)
        #         vals.append(t)
        # x_labels = [4*i for i in range(len(labels))]
        # plt.bar(x_labels, vals, align='center', alpha=0.5)
        # plt.xticks(x_labels, labels)
        # plt.xlabel('Algorithm')
        # plt.ylabel('Time')
        # plt.title('Time Taken for different algorithms')
        # plt.show()
        pass

