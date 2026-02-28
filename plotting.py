import matplotlib.pyplot as plt
from rg_functions import *

def plot_function(fun,
                  xlabel= "x",
                  ylabel="y",
                  caption="",
                  x_min=-5,
                  x_max=5,
                  num_points=100):
    """
    Plots the function fun(x) over the interval [x_min, x_max].

    Parameters:
        fun         : callable, function of one variable
        x_min       : float, left endpoint
        x_max       : float, right endpoint
        num_points  : int, number of sample points
    """
    x = np.linspace(x_min, x_max, num_points)
    y = fun(x)

    plt.figure()
    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(caption)
    plt.grid(True)
    # Plot x = y in red
    plt.plot(x, x, 'r-', label="x = y")
    plt.show()


if __name__ == "__main__":
    def main1():
        plot_function(One_Dim_Ising.K_prime,
            xlabel="K",
            ylabel= "K'(K)",
            caption= "plot of K'(K) for 1-dim ising",
            x_min=0,
            x_max=2,
            num_points=100)

        plot_function(Two_Dim_Ising.K_prime,
            xlabel="K",
            ylabel= "K'(K)",
            caption= "plot of K'(K) for 2-dim ising",
            x_min=0,
            x_max=2,
            num_points=100)
    def main2():
        plot_function(Two_Dim_Ising.zeta,
            xlabel="K",
            ylabel= "zeta(K)",
            caption= "plot of zeta(K) for 2-dim ising",
            x_min=0,
            x_max=2,
            num_points=100)

    main1()
    #main2()