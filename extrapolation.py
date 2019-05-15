import numpy as np
import matplotlib.pyplot as plt


class Extrapolate:

    def __init__(self, amu, spacing):
        self.amu_mean = []
        self.amu_err = []
        for i in range(len(amu)):
            self.amu_mean.append(amu[i][0])
            self.amu_err.append(amu[i][1])
        self.amu_mean = np.array(self.amu_mean)
        self.amu_err = np.array(self.amu_err)

        self.cov = self.get_cov()
        self.cov_inv = np.linalg.inv(self.cov)
        self.spacing2 = np.array(spacing) ** 2.
        return

    def get_cov(self):
        cov = np.zeros((len(self.amu_err), len(self.amu_err)))
        for i in range(len(self.amu_err)):
            cov[i][i] = self.amu_err[i] ** 2.
        return cov

    def linear_fit(self):
        self.fit_a_b = np.polyfit(self.spacing2, self.amu_mean, 1)
        hessian = self.get_linear_hessian()
        hessian_inv = np.linalg.inv(hessian)
        self.fit_a_b_err = np.array([2. * hessian_inv[0][0], 2. * hessian_inv[1][1]]) ** (1. / 2.)
        return

    def get_linear_hessian(self):
        hessian = np.zeros((2, 2))
        for a in range(2):
            for b in range(2):
                for i in range(len(self.amu_err)):
                    pfpa = self.spacing2[i] if a == 0 else 1.
                    pfpb = self.spacing2[i] if b == 0 else 1.
                    hessian[a][b] += 2. * pfpa * self.cov_inv[i][i] * pfpb
        return hessian

    def plt_amu(self, fmt, color, label):
        plt.errorbar(self.spacing2, self.amu_mean, self.amu_err,
                     fmt=fmt, capsize=2.5, elinewidth=1.2, label=label, color=color)
        return

    def plt_extrapolation_point(self, fmt, color, label):
        plt.errorbar(0, self.fit_a_b[1], self.fit_a_b_err[1],
                     fmt=fmt, capsize=2.5, elinewidth=1.2, label=label, color=color)
        return

    def plt_extrapolation_linear(self, fmt, label):
        x = np.arange(0.0, 0.016, 0.0001)
        f = lambda x: 1.0 * self.fit_a_b[0] * x + self.fit_a_b[1]
        plt.plot(x, f(x), fmt, label=label)
        return


if __name__ == '__main__':
    lat_spacing = [0.05684, 0.08787, 0.12121]
    amu = [(622.4, 26.9),
           (594.4, 10.3),
           (563.9, 8.4)]
    print('amu')
    print(amu)

    delta_amu_nlo = [15.6, 6.9, 2.1]
    amu_nlo = []
    for i in range(len(amu)):
        amu_nlo.append((amu[i][0] + delta_amu_nlo[i], amu[i][1]))
    print('amu_nlo')
    print(amu_nlo)

    delta_amu_taste = [9.5, 34.2, 51.6]
    amu_nlo_taste = []
    for i in range(len(amu_nlo)):
        amu_nlo_taste.append((amu_nlo[i][0] + delta_amu_taste[i], amu_nlo[i][1]))
    print('amu_nlo_taste')
    print(amu_nlo_taste)

    delta_amu_pion = [-1.1, -5.7, -2.2]
    amu_nlo_taste_pion = []
    for i in range(len(amu_nlo)):
        amu_nlo_taste_pion.append((amu_nlo_taste[i][0] + delta_amu_pion[i], amu_nlo[i][1]))
    print('amu_nlo_taste_pion')
    print(amu_nlo_taste_pion)

    exp_amu = Extrapolate(amu, lat_spacing)
    exp_amu.plt_amu('s', 'blue', '')
    exp_amu.linear_fit()
    # exp_amu.plt_extrapolation_point('x', 'black', '')
    print('amu:', exp_amu.fit_a_b[1], exp_amu.fit_a_b_err[1])

    exp_amu_nlo = Extrapolate(amu_nlo, lat_spacing)
    exp_amu_nlo.plt_amu('x', 'red', '')
    exp_amu_nlo.linear_fit()
    exp_amu_nlo.plt_extrapolation_point('x', 'red', '')
    exp_amu_nlo.plt_extrapolation_linear('r-', '')
    print('amu nlo:', exp_amu_nlo.fit_a_b[1], exp_amu_nlo.fit_a_b_err[1])

    exp_amu_nlo_taste = Extrapolate(amu_nlo_taste, lat_spacing)
    exp_amu_nlo_taste.plt_amu('o', 'black', '')
    exp_amu_nlo_taste.linear_fit()
    exp_amu_nlo_taste.plt_extrapolation_point('x', 'black', '')
    exp_amu_nlo_taste.plt_extrapolation_linear('b-', '')
    print('amu nlo taste:', exp_amu_nlo_taste.fit_a_b[1], exp_amu_nlo_taste.fit_a_b_err[1])

    exp_amu_nlo_taste_pion = Extrapolate(amu_nlo_taste_pion, lat_spacing)
    exp_amu_nlo_taste_pion.plt_amu('^', 'green', '')
    exp_amu_nlo_taste_pion.linear_fit()
    exp_amu_nlo_taste_pion.plt_extrapolation_point('x', 'green', '')
    exp_amu_nlo_taste_pion.plt_extrapolation_linear('g-', '')
    print('amu nlo taste pion:', exp_amu_nlo_taste_pion.fit_a_b[1], exp_amu_nlo_taste_pion.fit_a_b_err[1])

    plt.xlim(-0.001, 0.016)
    plt.show()
