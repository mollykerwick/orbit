import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if __name__ == '__main__':
    # file = open('Earth_periodicity.txt', 'r')
    df = pd.read_csv('hw4_orbit_Earth.txt', sep=' ', names=['x','y'])
    df2 = pd.read_csv('hw4_orbit_Jupiter.txt', sep=' ', names=['x','y'])
    # plot.scatter('Earth_periodicity.txt')
    # df.plot(x=0,y=1,kind="scatter")
    # print(df.columns)
    ax = df.plot(x=df.columns[0],y=df.columns[1],kind='line')
    ax = df2.plot(x=df2.columns[0],y=df2.columns[1],kind='line',ax=ax)
    ax.plot(0,0,'o')
    ax.legend(['Earth', 'Jupiter'])
    plt.show()
