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


class Two_Dim_Ising_Dbnet:

    def __init__(self, mean_spin):
        self.mean_spin = mean_spin

    def K_prime_av(self, K):
        cosh_0 = np.cosh(2 * K)
        cosh_m = np.cosh(K * (self.mean_spin - 2))
        cosh_p = np.cosh(K * (self.mean_spin + 2))
        return K + 2 * .25 * np.log(cosh_0 ** 2 / (cosh_m * cosh_p))

    def K_prime_plus(self, K):
        cosh_0 = np.cosh(2 * K)
        cosh_p = np.cosh(K * (self.mean_spin + 2))
        return K * (1 + .5 * self.mean_spin) + 2 * .5 * np.log(cosh_0 / cosh_p)

    def K_prime_minus(self, K):
        cosh_0 = np.cosh(2 * K)
        cosh_m = np.cosh(K * (self.mean_spin - 2))
        return K * (1 - .5 * self.mean_spin) + 2 * .5 * np.log(cosh_0 / cosh_m)


if __name__ == "__main__":
    def main1():
        K = 1
        for mean_spin in [-1, 1, -.3, .3, 0]:
            K_prime = Two_Dim_Ising_Dbnet(mean_spin).K_prime_av(K)
            print(mean_spin, K_prime)


    main1()
