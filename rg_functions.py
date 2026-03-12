import numpy as np
from scipy.optimize import newton
from globals import *


class One_Dim_Ising:

    @staticmethod
    def Kprime(K):
        """
        This method outputs the RG transformation K'(K) for the standard 1-dim
        Ising model

        Parameters
        ----------
        K: float
            initial coupling constant (K = \beta *J)

        Returns
        -------
        Kprime: float
            the value of K' = final coupling constant
        """
        return .5 * np.log(np.cosh(2 * K))

    @staticmethod
    def zeta(K, tol=1e-8):
        """
        Partition function for 1-dim Ising Model
        (recursion formula from misc/rg-lecture)

        zeta(K') = 2*zeta(K) - ln f(K)
        Z(K=0)= 2^N
        zeta(0) = ln(Z)/N = ln(2)

        Parameters
        ----------
        K: float

        Returns
        -------
        float

        """
        Ks = [K]

        # RG flow forward until K ~ 0
        while Ks[-1] > tol:
            Ks.append(One_Dim_Ising.Kprime(Ks[-1]))

        # boundary condition
        zeta = np.log(2)

        # propagate backwards
        for K in reversed(Ks[:-1]):
            zeta = 0.5 * (zeta + 0.5 * np.log(np.cosh(2 * K)) + np.log(2))

        return zeta


class Two_Dim_Ising:
    @staticmethod
    def Kprime(K):
        """
        This method outputs the RG transformation K'(K) for the standard 2-dim
        Ising model

        Parameters
        ----------
        K: float
            initial coupling constant (K = \beta *J)

        Returns
        -------
        Kprime: float
            the value of K' = final coupling constant

        """
        return (3 / 8) * np.log(np.cosh(4 * K))

    @staticmethod
    def f(K):
        """
        This method calculates f(K), a function of K used in calculation of
        zeta(K) below

        Parameters
        ----------
        K: float

        Returns
        -------
        float

        """
        return 2 * (np.cosh(2 * K) ** .5) * (np.cosh(4 * K) ** (1 / 8))

    @staticmethod
    def zeta(K):
        """
        Partition function for 2-dim Ising Model (formula from misc/rg-lecture)

        Parameters
        ----------
        K: float

        Returns
        -------
        float

        """
        Kprime = Two_Dim_Ising.Kprime(K)
        f = Two_Dim_Ising.f(K)
        fprime = Two_Dim_Ising.f(Kprime)
        return .5 * np.log(f * fprime)


class Two_Dim_Ising_Dbnet:
    """
    This class calculates quantities related to the 2-dim Ising-dbnet model
    (see https://github.com/rrtucci/ising-dbnet) In particular,
    it calculates its real space RG transformation proposed in the
    white paper included in this repo.


    Attributes
    ----------
    Pp: float
        Probability $P_+=P(S_5=+1)$. Pm =  probability $P_-=P(S_5=-1)$
    eta: float
        sum of spins of 4 nearest neighbors to S_5. eta=4 if all 4 nearest
        neighbors are +1
    lam: float
        positive or negative float, proportional to <S_1>=the average spin for
        the next-to-next nearest neighbors of S_5

    """

    def __init__(self, Pp=1 / 4, lam=1 / 8, eta=4):
        """
        Constructor

        Parameters
        ----------
        Pp: float
        lam: float
        eta: float
        """
        self.Pp = Pp
        self.lam = lam
        self.eta = eta

    def Kprime(self, K):
        """
        This method outputs the RG transformation K'(K) for the 2-dim
        Ising-dbnet model

        Parameters
        ----------
        K: float
            initial coupling constant (K = \beta *J)

        Returns
        -------
        Kprime: float
            the value of K' = final coupling constant

        """
        Pm = 1 - self.Pp

        Denom = self.Pp / np.cosh(K * (1 + self.lam)) ** 4 + \
                Pm / np.cosh(K * (-1 + self.lam)) ** 4
        Num = self.Pp * np.exp(K * self.eta * (1 + 4 * self.lam)) / \
              np.cosh(K * (1 + self.lam)) ** 4 + \
              Pm * np.exp(K * self.eta * (-1 + 4 * self.lam)) / \
              np.cosh(K * (-1 + self.lam)) ** 4

        return (np.log(Num) - np.log(Denom)) / (4 * self.eta * self.lam)

# if __name__ == "__main__":
