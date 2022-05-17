from cProfile import label
from matplotlib import markers
import numpy as np
import matplotlib.pyplot as plt
import csv

from pytest import mark


if __name__ == '__main__':

    x = np.array([8, 9, 10, 11, 12, 13])
    # plt.xticks(x)

    with open('AgeTest.csv')as f:
        f_csv = csv.reader(f)
        for row in f_csv:
            # print(row[1:])
            np_csv = []
            for data in row[1:]:
                np_csv.append(float(data))

            plt.plot(x, np.array(np_csv), marker='o', label=row[0])
    plt.legend()
    plt.show()
