from cProfile import label
from matplotlib import markers
import numpy as np
import matplotlib.pyplot as plt
import csv

from pytest import mark
import matplotlib.colors as mcolors
colors = list(mcolors.TABLEAU_COLORS.keys())


def calc(x):
    return pow(2, x)-1


e_list = [0.01, 0.1, 1]
test_list = ['Income', 'Retail']


def draw_plot(name):
    plt.clf()
    x = np.array([calc(8), calc(9), calc(10), calc(11), calc(12), calc(13)])
    px = np.array([500, 2000, 4000, 8000])
    plt.xticks(px)
    plt.xlabel("Size of data(T)")
    plt.ylabel("Mean square error")
    # plt.ylabel("Sensitivity")
    with open(name+'.csv')as f:
        f_csv = csv.reader(f)
        for row in f_csv:
            np_csv = []
            for data in row[1:]:
                np_csv.append(float(data))

            plt.plot(x, np.array(np_csv), marker='o',  label=row[0])
    plt.legend()
    plt.savefig(name+'.png')


def draw_plot2(name):
    plt.clf()
    x = np.array([calc(8), calc(9), calc(10), calc(11), calc(12)])
    px = np.array([500, 1000, 2000, 4000])
    plt.xticks(px)
    plt.xlabel("Size of data(T)")
    plt.ylabel("Mean square error")
    # plt.ylabel("Sensitivity")
    with open(name+'.csv')as f:
        f_csv = csv.reader(f)
        for row in f_csv:
            np_csv = []
            for data in row[1:]:
                np_csv.append(float(data))
            print(np_csv)
            plt.plot(x, np.array(np_csv), marker='o',  label=row[0])
    plt.legend()
    plt.savefig(name+'.png')


def draw_bar(name):
    plt.clf()
    px = np.array([calc(9), calc(10), calc(11), calc(12)])
    x = np.arange(len(px))
    plt.xticks(x, labels=px)
    # plt.xticklabels(px)
    plt.xlabel("Size of data(T)")
    plt.ylabel("Time(s)")
    plt.ylim((0, 1000))
    width = 0.2
    with open(name+'.csv')as f:
        f_csv = csv.reader(f)
        pp = -1
        for row in f_csv:
            np_csv = []
            for data in row[1:]:
                np_csv.append(float(data))
            plt.bar(x+width*pp,
                    np.array(np_csv), label=row[0], width=width)
            pp += 1
    plt.legend()
    # plt.show()

    plt.savefig(name+'.png')


if __name__ == '__main__':
    # for e in e_list:
    #     for test in test_list:
    #         draw_plot(str(int(e*100))+test)
    # draw_plot("1lenRetail")
    # draw_bar('time10Income')
    draw_bar('mbmtime')
    draw_plot2('mbmerror')
