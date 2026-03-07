import numpy as np
from scipy.optimize import newton
from globals import *


class One_Dim_Ising:

    @staticmethod
    def Kprime(K):
        return .5 * np.log(np.cosh(2 * K))


class Two_Dim_Ising:
    @staticmethod
    def Kprime(K):
        return (3 / 8) * np.log(np.cosh(4 * K))

    @staticmethod
    def f(K):
        return 2 * (np.cosh(2 * K) ** .5) * 2 * (np.cosh(4 * K) ** (1 / 8))

    @staticmethod
    def zeta(K):
        Kprime = Two_Dim_Ising.Kprime(K)
        f = Two_Dim_Ising.f(K)
        fprime = Two_Dim_Ising.f(Kprime)
        return .5 * np.log(f * fprime)


class Two_Dim_Ising_Dbnet:

    def __init__(self, Pp=1/4, lam=1/8, eta=4):
        self.Pp = Pp
        self.lam = lam
        self.eta = eta

    def Kprime(self, K):
        Pm = 1 - self.Pp

        Denom = self.Pp / np.cosh(K * (1 + self.lam)) ** 4 + \
                Pm / np.cosh(K * (-1 + self.lam)) ** 4
        Num = self.Pp * np.exp(K * self.eta * (1 + 4 * self.lam)) /\
              np.cosh(K * (1 + self.lam))** 4 + \
              Pm * np.exp(K * self.eta * (-1 + 4 * self.lam)) /\
              np.cosh(K * (-1 + self.lam)) ** 4

        return (np.log(Num) - np.log(Denom)) / (4 * self.eta * self.lam)

class Two_Dim_Ising_Dbnet1:

    def __init__(self, mean_spin):
        self.mean_spin = mean_spin

    def Kprime_av(self, K):
        cosh_0 = np.cosh(2 * K)
        cosh_m = np.cosh(K * (self.mean_spin - 2))
        cosh_p = np.cosh(K * (self.mean_spin + 2))
        return K + (2 / AA) * np.log(cosh_0 ** 2 / (cosh_m * cosh_p))

    def Kprime_plus(self, K):
        cosh_0 = np.cosh(2 * K)
        cosh_p = np.cosh(K * (self.mean_spin + 2))
        return K * (1 + (4 / AA) * self.mean_spin) + \
            (4 / AA) * np.log(cosh_0 / cosh_p)

    def Kprime_minus(self, K):
        cosh_0 = np.cosh(2 * K)
        cosh_m = np.cosh(K * (self.mean_spin - 2))
        return K * (1 - (4 / AA) * self.mean_spin) + \
            (4 / AA) * np.log(cosh_0 / cosh_m)


class Two_Dim_Ising_Dbnet2:
    def __init__(self, prob_plus):
        self.prob_plus = prob_plus

    def f(self, Kprime):
        x = np.exp(AA * Kprime * K_CURIE) \
            / (16 * (np.cosh(TWO * Kprime * K_CURIE) ** 4))
        return x

    def der_f(self, Kprime):
        coeff = self.f(Kprime) * K_CURIE
        return coeff * (AA - 4 * TWO * np.tanh(TWO * Kprime * K_CURIE))

    def inv_f(self, y, Kprime_init=.1):
        x0 = Kprime_init
        return newton(
            func=lambda x: self.f(x) - y,
            x0=x0,
            fprime=self.der_f
        )

    def x_3(self, K):
        return np.exp((AA + 4) * K * K_CURIE) \
            / (16 * (np.cosh((TWO + 1) * K * K_CURIE) ** 4))

    def x_1(self, K):
        return np.exp((AA - 4) * K * K_CURIE) \
            / (16 * (np.cosh((TWO - 1) * K * K_CURIE) ** 4))

    def Kprime_plus(self, K):
        prob_minus = 1 - self.prob_plus
        y = self.prob_plus * self.x_3(K) + prob_minus * self.x_1(K)
        return self.inv_f(y)

    def Kprime_minus(self, K):
        prob_minus = 1 - self.prob_plus
        y = self.prob_plus * self.x_1(K) + prob_minus * self.x_3(K)
        return self.inv_f(y)


if __name__ == "__main__":
    def main1():
        K = 1
        for mean_spin in [-1, 1, -.3, .3, 0]:
            Kprime = Two_Dim_Ising_Dbnet1(mean_spin).Kprime_av(K)
            print(mean_spin, Kprime)


    def main2():
        def fun(x):
            return x ** 3 + x

        def dfun(x):
            return 3 * x ** 2 + 1

        def invfun(y, x0=20):
            return newton(
                func=lambda x: fun(x) - y,
                x0=x0,
                fprime=dfun
            )

        x = invfun(1)

        print("x =", x)
        print("Check:", fun(x))


    def main3(x):
        mom = Two_Dim_Ising_Dbnet2(0)
        print("f(0)", mom.f(x))


    # main1()
    # main2()
    main3(0)
