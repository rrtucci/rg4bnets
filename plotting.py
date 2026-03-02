import matplotlib.pyplot as plt
from rg_functions import *


def plot_functions(funs,
                   fun_names=None,
                   xlabel="x",
                   ylabel="y",
                   caption="",
                   x_min=-5,
                   x_max=5,
                   num_points=100):
    """
    Plots several functions fun(x) over the interval [x_min, x_max].

    Parameters:
        funs         : list of functions of one variable
        x_min       : float, left endpoint
        x_max       : float, right endpoint
        num_points  : int, number of sample points
    """
    # plt.close('all')
    x = np.linspace(x_min, x_max, num_points)
    fig, ax = plt.subplots()  # single frame / single axes

    for i, fun in enumerate(funs):
        y = fun(x)
        if not fun_names:
            label = f"fun_{i + 1}(x)"
        else:
            label = fun_names[i]
        ax.plot(x, y, label=label)

    # Add x = y in red
    ax.plot(x, x, 'r-', linewidth=2, label="x = y")

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(caption)
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    def main0():
        funs = [
            lambda x: np.sin(x),
            lambda x: np.cos(x),
            lambda x: 0.5 * x
        ]

        plot_functions(funs)


    def main1():
        plot_functions([One_Dim_Ising.K_prime],
                       xlabel="K",
                       ylabel="K'(K)",
                       caption="plot of K'(K) for 1-dim ising",
                       x_min=0,
                       x_max=2,
                       num_points=100)

        plot_functions([Two_Dim_Ising.K_prime],
                       xlabel="K",
                       ylabel="K'(K)",
                       caption="plot of K'(K) for 2-dim ising",
                       x_min=0,
                       x_max=2,
                       num_points=100)


    def main2():
        plot_functions([Two_Dim_Ising.zeta],
                       xlabel="K",
                       ylabel="zeta(K)",
                       caption="plot of zeta(K) for 2-dim ising",
                       x_min=0,
                       x_max=2,
                       num_points=100)


    def main3(mean_spin):
        mom = Two_Dim_Ising_Dbnet(mean_spin)
        funs = [mom.K_prime_minus, mom.K_prime_av, mom.K_prime_plus]
        fun_names = ["K' minus", "K' av", "K' plus"]

        plot_functions(funs,
                       fun_names=fun_names,
                       xlabel="K",
                       ylabel="K'(K)",
                       caption="plot of K'(K) for 2-dim ising-dbnet",
                       x_min=0,
                       x_max=2,
                       num_points=100)


    # main1()
    # main2()
    main3(1)
