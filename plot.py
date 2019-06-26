import matplotlib.pyplot as plt
# import numpy as np
import pandas as pd

if __name__ == '__main__':
    df = pd.read_csv('hw4_orbit_Earth.txt', sep=' ', names=['x','y'])
    df2 = pd.read_csv('hw4_orbit_Jupiter.txt', sep=' ', names=['x','y'])
    
    ax = df.plot(x=df.columns[0],y=df.columns[1],kind='line')
    ax = df2.plot(x=df2.columns[0],y=df2.columns[1],kind='line',ax=ax)
    ax.plot(0,0,'o')
    ax.legend(['Earth', 'Jupiter'])
    plt.show()
