from cProfile import label
from matplotlib import markers
import numpy as np
import matplotlib.pyplot as plt
import csv

from pytest import mark


def calc(x):
    return pow(2, x)-1


e_list = [0.01, 0.1, 1]
test_list = ['Income', 'Retail']


def work(name):
    plt.clf()
    x = np.array([calc(8), calc(9), calc(10), calc(11), calc(12), calc(13)])
    px = np.array([500, 2000, 4000, 8000])
    plt.xticks(px)
    plt.xlabel("Size of data(T)")
    plt.ylabel("mean square error")
    with open(name+'.csv')as f:
        f_csv = csv.reader(f)
        for row in f_csv:
            np_csv = []
            for data in row[1:]:
                np_csv.append(float(data))

            plt.plot(x, np.array(np_csv), marker='o', label=row[0])
    plt.legend()
    plt.savefig(name+'.png')


if __name__ == '__main__':
    for e in e_list:
        for test in test_list:
            work(str(int(e*100))+test)
    work("1lenRetail")
