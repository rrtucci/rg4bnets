import matplotlib.pyplot as plt
from rg_functions import *


def plot_functions(funs,
                   fun_names=None,
                   xlabel="x",
                   ylabel="y",
                   caption="",
                   x_min=-5,
                   x_max=5,
                   num_points=100,
                   do_vert_line = True,
                   vert_x=K_CURIE,
                   vert_label=r"$K = K_c$"):
    """

    Parameters
    ----------
    funs: list[function]
    fun_names: lis[str]
    xlabel: str
    ylabel: str
    caption: str
    x_min: int
    x_max: int
    num_points: int
    do_vert_line: bool
    vert_x: float
    vert_label: str

    Returns
    -------

    """
    # plt.close('all')
    xx = np.linspace(x_min, x_max, num_points)
    fig, ax = plt.subplots()  # single frame / single axes

    for i, ff in enumerate(funs):
        yy = [ff(x) for x in xx]
        if not fun_names:
            label = f"fun_{i + 1}(x)"
        else:
            label = fun_names[i]
        ax.plot(xx, yy, label=label)

    # Add x = y in red
    ax.plot(xx, xx, 'r-', linewidth=2, label="x = y")

    if do_vert_line:
        ax.axvline(vert_x, color='black', linestyle='--', label=vert_label)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(caption)
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":

    def main1():
        plot_functions([One_Dim_Ising.Kprime],
                       xlabel="K",
                       ylabel="K'(K)",
                       caption="plot of K'(K) for 1-dim ising",
                       x_min=0,
                       x_max=2,
                       num_points=100)

        plot_functions([Two_Dim_Ising.Kprime],
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

    def main3():
        mom = Two_Dim_Ising_Dbnet()
        plot_functions([Two_Dim_Ising.Kprime, mom.Kprime],
                       fun_names=["2-dim Ising", "2-dim Ising-dbnet"],
                       xlabel="K",
                       ylabel="K'(K)",
                       caption="plots of K'(K)",
                       x_min=0,
                       x_max=1,
                       num_points=100)


    # main1()
    # main2()
    main3()
