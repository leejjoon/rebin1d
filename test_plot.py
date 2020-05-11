import numpy as np
from intp_integrated import (AccumulatedInterpNN,
                             AccumulatedInterpDrizzle,
                             AccumulatedInterpLinear)

import matplotlib.pyplot as plt


def test_plot_accumulate_interp_nn():
    x = [0, 1, 2, 3, 4]
    y = [1, 3, 2, 4, 1.5]

    xx = np.arange(-1., 5., 0.05)

    intp = AccumulatedInterpNN(x, y)

    yy = intp.rebin(xx)

    plt.clf()
    plt.plot(x, y, "o")
    plt.plot(xx[:-1], yy/0.05)

    intp = AccumulatedInterpDrizzle(x, y, drizzle=0.6)

    yy = intp.rebin(xx)

    plt.plot(xx[:-1], yy/0.05)



def test_plot_accumulate_interp_linear():
    x = [0, 1, 2, 3, 4]
    y = [1, 3, 2, 4, 1.5]

    dx = 0.2
    xx = np.arange(0.01, 4., dx)

    intp = AccumulatedInterpLinear(x, y)

    yy = intp.rebin(xx) / dx
    xc = 0.5 * (xx[:-1] + xx[1:])

    plt.clf()
    plt.plot(x, y, "o")
    plt.plot(xc, yy)
