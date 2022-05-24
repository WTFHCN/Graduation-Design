from cProfile import label
from matplotlib import markers
import numpy as np
import matplotlib.pyplot as plt
import csv

from pytest import mark


def calc(x):
    return pow(2, x)-1


if __name__ == '__main__':

    x = np.array([calc(8), calc(9), calc(10), calc(11), calc(12), calc(13)])
    px = np.array([700,  2000, 4000, 8000])
    plt.xticks(px)
    plt.xlabel("T")
    plt.ylabel("Squared Error")
    with open('AgeTest.csv')as f:
        f_csv = csv.reader(f)
        for row in f_csv:
            # print(row[1:])
            np_csv = []
            for data in row[1:]:
                np_csv.append(float(data))

            plt.plot(x, np.array(np_csv), marker='o', label=row[0])
    plt.legend()
    # plt.show()
    plt.savefig("2333.png")
