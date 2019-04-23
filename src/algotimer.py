import os
import matplotlib.pyplot as plt


class TimerObj:
    def __init__(self, input_path, time_path):
        self.input_path = input_path
        self.time_path = time_path+"/time.txt"

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
        return self.input_path

        labels = []
        vals = []
        for dir in dirs:
            if os.path.isdir(path+dir):
                timefile = path+dir+"/time.txt"
                t = get_time(timefile)
                labels.append(dir)
                vals.append(t)
        x_labels = [4*i for i in range(len(labels))]
        plt.bar(x_labels, vals, align='center', alpha=0.5)
        plt.xticks(x_labels, labels)
        plt.xlabel('Algorithm')
        plt.ylabel('Time')
        plt.title('Time Taken for different algorithms')
        plt.show()

if __name__ == "__main__":
    main()
