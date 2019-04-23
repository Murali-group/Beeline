import os
import matplotlib.pyplot as plt

def get_time(timefile):
    with open(timefile, "r+") as f:
        lines = f.readlines()
        line = lines[1]
        return float(line.split()[-1])



def main():
    path = '../outputs/simulated/HSC/'

    dirs = os.listdir(path)

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
