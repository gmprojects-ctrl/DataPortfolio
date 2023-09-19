""" My Code """
# solution to parts (b), (c)
# YOUR CODE HERE
a = 1.0e-11


def eqn2(x, y, energy):
    """Eqn 2"""
    w_, dw_ = y
    V_0 = 50 * e_el
    V_ = V_0 * ((x * x) / (a * a))
    f_0 = dw_
    f_1 = (-2.0 * m_el / (hbar * hbar)) * (energy * w_ - V_ * w_)
    # f_1 = -2m/h^2 ( Ew-Vw )
    return [f_0, f_1]


def solve2(energy):
    """Solve 2"""
    x_range = [-10 * a, 10 * a]
    x_eval = np.linspace(x_range[0], x_range[1], 1000)
    init = [0, 1]
    answer = solve_ivp(eqn2, x_range, init, args=(energy,), t_eval=x_eval)
    return answer.y[0, -1]


def solve3(energy):
    """Solve 3"""
    x_range = [-10 * a, 10 * a]
    x_eval = np.linspace(x_range[0], x_range[1], 1000)
    init = [0, 1]
    answer = solve_ivp(eqn2, x_range, init, args=(energy,), t_eval=x_eval)
    return answer


e_1 = bisect(solve2, 100 * e_el, 200 * e_el, xtol=1e-22)
e_2 = bisect(solve2, 200 * e_el, 1000 * e_el, xtol=1e-22)
e_3 = bisect(solve2, e_2, 900 * e_el, xtol=1e-22)
difference = e_2 - e_1


def result():
    """Result Function"""
    return difference / e_el


E_1 = solve3(e_1)
E_2 = solve3(e_2)
E_3 = solve3(e_3)
E_1y = E_1.y[0] / np.sqrt(trapezoid(np.abs(E_1.y[0]) ** 2, x=E_1.t))
E_2y = E_2.y[0] / np.sqrt(trapezoid(np.abs(E_2.y[0]) ** 2, x=E_2.t))
E_3y = E_3.y[0] / np.sqrt(trapezoid(np.abs(E_3.y[0]) ** 2, x=E_3.t))
plt.figure(figsize=(16, 9))
plt.title("Normalised Wave Forms", fontsize=18)
plt.grid(True)
plt.xlim(-5 * a, 5 * a)
plt.ylim(-2e5, 2e5)
plt.xlabel("x (10e-11)", fontsize=18)
plt.ylabel(r"$\omega(x)$", fontsize=18)
plt.plot(E_1.t, E_1y, label=f"When E = {round(e_1/e_el,3)}eV")
plt.plot(E_2.t, E_2y, label=f"When E = {round(e_2/e_el,3)}eV")
plt.plot(E_3.t, E_3y, label=f"When E = {round(e_3/e_el,3)}eV")
plt.legend(fontsize=18)
plt.show()
