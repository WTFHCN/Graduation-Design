from cProfile import label
from matplotlib import markers
import numpy as np
import matplotlib.pyplot as plt
import csv

from pytest import mark


if __name__ == '__main__':

    x = np.linspace(1, 10, 10)

    plt.xticks(x)

    with open('test.csv')as f:
        f_csv = csv.reader(f)
        for row in f_csv:
            # print(row[1:])
            np_csv = []
            for data in row[1:]:
                np_csv.append(float(data))
                if len(np_csv) == 10:
                    break
            plt.plot(x, np.array(np_csv), marker='o', label=row[0])
    plt.legend()
    plt.show()
