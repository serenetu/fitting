import numpy as np


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
        print(self.cov)
        self.cov_inv = np.linalg.inv(self.cov)
        print(self.cov_inv)
        self.spacing2 = np.array(spacing) ** 2.
        print(self.spacing2)
        return

    def get_cov(self):
        cov = np.zeros((len(self.amu_err), len(self.amu_err)))
        for i in range(len(self.amu_err)):
            cov[i][i] = self.amu_err[i] ** 2.
        return cov

    def linear_fit(self):
        self.fit_a_b = np.polyfit(self.spacing2, self.amu_mean, 1)
        print(self.fit_a_b)
        hessian = self.get_linear_hessian()
        hessian_inv = np.linalg.inv(hessian)
        self.fit_a_b_err = np.array([hessian_inv[0][0], hessian_inv[1][1]]) ** (1. / 2.)
        print(self.fit_a_b_err)
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


if __name__ == '__main__':
    lat_spacing = [0.05684, 0.08787, 0.12121]
    amu = [(622.4, 26.9),
           (594.4, 10.3),
           (563.9, 8.4)]

    exp = Extrapolate(amu, lat_spacing)
    exp.linear_fit()
