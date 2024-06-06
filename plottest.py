import matplotlib
import matplotlib.pyplot as plt
import numpy as np


matplotlib.use("TkAgg")
# added to display plot initialized from terminal using tkinter
# from: https://github.com/matplotlib/matplotlib/issues/13414


x = [1,2,3,4,5,6,7,8,9]
y = [2,4,6,7,8,3,4,5,6]

plt.plot(x, y)
plt.show()

