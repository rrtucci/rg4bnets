import numpy as np


class One_Dim_Ising:

    @staticmethod
    def K_prime(K):
        return .5 * np.log(np.cosh(2 * K))


class Two_Dim_Ising:
    @staticmethod
    def K_prime(K):
        return (3 / 8) * np.log(np.cosh(4 * K))

    @staticmethod
    def f(K):
        return 2 * (np.cosh(2 * K) ** .5) * 2 * (np.cosh(4 * K) ** (1 / 8))

    @staticmethod
    def zeta(K):
        Kprime = Two_Dim_Ising.K_prime(K)
        f = Two_Dim_Ising.f(K)
        fprime = Two_Dim_Ising.f(Kprime)
        return .5 * np.log(f * fprime)
