import numpy as np
from scipy.interpolate import interp1d


def _prepend_column(a, c):
    c = np.asarray(c)
    r = np.hstack([c[..., np.newaxis], a])
    return r


def _append_column(a, c):
    c = np.asarray(c)
    r = np.hstack([a, c[..., np.newaxis]])
    return r


def _pre_ap_pend_column(a, c_left, c_right):
    c_left = np.asarray(c_left)
    c_right = np.asarray(c_right)
    r = np.hstack([c_left[..., np.newaxis], a, c_right[..., np.newaxis]])
    return r


def test_append():
    a = np.arange(6).reshape((2, 3))
    r = _append_column(a, [0, 0])

    assert np.all(r[..., -1] == [0, 0])


def get_bin_boundaries(wv):
    """
    Given the bin centers, try to estimate bin boundaries.
    """
    wv = np.asarray(wv)

    # dww_ = (wv[..., 1:] - wv[..., :-1])
    # dww_left = dww_[0]
    # dww_right = dww_[-1]
    dww_left = wv[..., 1] - wv[..., 0]
    dww_right = wv[..., -1] - wv[..., -2]
    cww_ = .5*(wv[..., 1:] + wv[..., :-1])

    boundaries = _pre_ap_pend_column(cww_,
                                     cww_[..., 0] - dww_left,
                                     cww_[..., -1] + dww_right)
    dw = boundaries[1:] - boundaries[:-1]
    return boundaries, dw


class AccumulatedInterpBase(object):
    def rebin(self, new_xbins):
        v = self.integrate_to(new_xbins)
        return v[..., 1:] - v[..., :-1]


class AccumulatedInterpNN(AccumulatedInterpBase):

    def __init__(self, x, y):
        """
        x, y : flux density y at position x
        """

        x = np.asarray(x)
        y = np.asarray(y)

        if len(x) == len(y) + 1:
            b = x
            dw = b[1:] - b[:-1]
        else:
            b, dw = get_bin_boundaries(x)

        acc = np.add.accumulate(y * dw)

        self._intp = interp1d(b, np.hstack([[0], acc]),
                              kind="linear",
                              bounds_error=False,
                              fill_value=(0, acc[-1]))

    def integrate_to(self, x0, debug=False):
        return self._intp(x0)


class AccumulatedInterpDrizzle(AccumulatedInterpBase):

    def __init__(self, x, y, drizzle=0.8):
        assert drizzle < 1.

        x = np.asarray(x)
        y = np.asarray(y)
        dw_ = x[1:] - x[:-1]

        dw = 0.5 * (np.hstack([dw_[0], dw_]) +
                    np.hstack([dw_, dw_[-1]]))
        dw2 = 0.5 * drizzle * dw

        n = len(x)
        xx = np.empty(n*2, dtype="d")
        xx[::2] = x - dw2
        xx[1::2] = x + dw2

        # yy = np.zeros(n*2, dtype="d")
        # yy[1::2] = y / drizzle

        yy = np.empty(n*2, dtype="d")

        acc = np.add.accumulate(y * dw)

        yy[0] = 0
        yy[1::2] = acc
        yy[2::2] = acc[:-1]

        self._intp = interp1d(xx, yy,
                              kind="linear",
                              bounds_error=False,
                              fill_value=(0, yy[-1]))

    def integrate_to(self, x0, debug=False):
        return self._intp(x0)


def test_accumulate_interp_nn():
    x = [0, 1, 2]
    y = [1, 3, 2]

    intp = AccumulatedInterpNN(x, y)

    assert intp.rebin([0.5, 1.5])[0] == 3.
    assert intp.rebin([0., 1.])[0] == 2.

    x = [-0.5, 0.5, 1.5, 2.5]
    y = [1, 3, 2]

    intp = AccumulatedInterpNN(x, y)

    assert intp.rebin([0.5, 1.5])[0] == 3.
    assert intp.rebin([0., 1.])[0] == 2.


class AccumulatedInterpLinear(AccumulatedInterpBase):
    def __init__(self, x, y):
        """
        x, y : flux density y at position x
        """

        x = np.asarray(x)
        y = np.asarray(y)

        self._x = x
        self._y = y
        self._c_dx = (x[1:] - x[:-1])
        self._c_dy = (y[1:] - y[:-1])
        # _c_dy_dx = _c_dy / _c_dx
        self._c_dxdy = self._c_dx * self._c_dy

        self._c_y = 0.5 * (y[1:] + y[:-1])
        self._acc_y = np.concatenate([[0],
                                      np.add.accumulate(self._c_y * self._c_dx)])
        # self._acc_y = np.concatenate([[0],
        #                               np.add.accumulate(y[1:] * self._c_dx)])

    def integrate_to(self, x0, debug=False):

        i0 = self._x.searchsorted(x0, side="left")
        m = (0 < i0) & (i0 < len(self._x))
        if debug:
            print(m, np.all(m))

        if not np.all(m):
            raise RuntimeError("index out of [0, {n}]: {indices}"
                               .format(n=len(self._x),
                                       indices=np.array(x0)[~m]))

        u = (x0 - self._x[i0 - 1]) / self._c_dx[i0 - 1]
        if debug:
            print(i0, u)

        sa = self._acc_y[i0 - 1]
        sb = self._c_dx[i0 - 1] * self._y[i0 - 1] * u
        sc = .5 * self._c_dxdy[i0 - 1] * u**2

        if debug:
            print(sa, sb, sc)

        return sa + sb + sc



if False:

    # x = np.array([1, 2, 3, 4, 6, 8])
    x = np.array([0, 1, 2, 3, 4, 5, 6])
    y = np.array([0, 1, 1, 1, 2, 2, 2])

    # assert np.allclose(IntP.t(x, y, 1.6), 0.6)
    # assert np.allclose(IntP.t(x, y, 2.5), 1.625)

    intp = AccumulatedInterpLinear(x, y)

    t34 = intp.integrate_to([3, 4])

    v = intp.integrate_to([1, 2, 3, 4])
    dv = v[1:] - v[:-1]

    clf()
    plot(x, y)
    plot([1.5, 2.5, 3.5], dv, "o")

    xk = np.arange(1.5, 4.5, 0.1)
    dx = 1.
    kk = [intp.rebin([x1-dx, x1+dx]) for x1 in xk]

if False:
    b = np.array([1, 2, 3, 4, 5])
    # bb = .5 * (b[:-1] + b[1:])
    h = intp.rebin(b)

    v = intp.integrate_to([1.5, 2.5, 3.5, 4.5])
    v3 = intp.integrate_to(3, debug=True)
    v35 = intp.integrate_to(3.5, debug=True)

    vv = (v35 - v3)/0.5
    plot(3.25, vv, "*")

    # b = np.array([1.5, 2.5, 3.5, 4.5, 5.5])
    # bb = .5 * (b[:-1] + b[1:])
    # h = intp.rebin(b)
    # db = b[1:] - b[:-1]

    # clf()
    # plot(x, y)
    # plot(bb, h, "o")
