import numpy as np
from scipy.optimize import newton
from globals import *


class One_Dim_Ising:

    @staticmethod
    def Kprime(K):
        """

        Parameters
        ----------
        K: float

        Returns
        -------
        float

        """
        return .5 * np.log(np.cosh(2 * K))


class Two_Dim_Ising:
    @staticmethod
    def Kprime(K):
        """

        Parameters
        ----------
        K: float

        Returns
        -------
        float

        """
        return (3 / 8) * np.log(np.cosh(4 * K))

    @staticmethod
    def f(K):
        """

        Parameters
        ----------
        K: float

        Returns
        -------
        float

        """
        return 2 * (np.cosh(2 * K) ** .5) * 2 * (np.cosh(4 * K) ** (1 / 8))

    @staticmethod
    def zeta(K):
        """

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
    Attributes
    ----------
    Pp: float
    eta: float
    lam: float

    """

    def __init__(self, Pp=1/4, lam=1/8, eta=4):
        """

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

        Parameters
        ----------
        K: float

        Returns
        -------
        float

        """
        Pm = 1 - self.Pp

        Denom = self.Pp / np.cosh(K * (1 + self.lam)) ** 4 + \
                Pm / np.cosh(K * (-1 + self.lam)) ** 4
        Num = self.Pp * np.exp(K * self.eta * (1 + 4 * self.lam)) /\
              np.cosh(K * (1 + self.lam))** 4 + \
              Pm * np.exp(K * self.eta * (-1 + 4 * self.lam)) /\
              np.cosh(K * (-1 + self.lam)) ** 4

        return (np.log(Num) - np.log(Denom)) / (4 * self.eta * self.lam)


#if __name__ == "__main__":

