import numpy as np
import scipy.integrate as integrate


def sol(s, c, fun0, B, F, G, D, m, Time_end):
    global h, w, ds, S, T, X, Y, \
        B11, B12, B21, B22, \
        c1, c2, fun_x0, fun_y0, \
        G1, G2, \
        F1, F2, \
        D1, D2, \
        s0, s_end, \
        Y_left, P2_left, T_left, \
        P1_right, X_right, T_right, \
        kr, kl, h_l, h_r

    s0 = s[0]
    s_end = s[1]

    c1 = c[0]
    c2 = c[1]

    fun_x0 = fun0[0]
    fun_y0 = fun0[1]

    B11 = B[0][0]
    B12 = B[0][1]
    B21 = B[1][0]
    B22 = B[1][1]

    F1 = F[0]
    F2 = F[1]

    G1 = G[0]
    G2 = G[1]

    D1 = D[0]
    D2 = D[1]


    ### Настройка метода данных для сетки и шагов    ---------------------------------
    ds = (s_end - s0) / m
    h = ds / (c1 + c2)  # высота уровня
    w = c2 * h  # постоянный сдвиг
    level = round(Time_end // h)  # кол-во уровней

    h_l = ds / c2  # Шаг времени 1 граничного условия
    h_r = ds / c1  # Шаг времени 2 граничного условия
    ### Конец настройки метода ------------------------------------------------------

    # Выделение памяти --------------------------------------------------------------
    X = np.zeros((level + 1, m))
    Y = np.zeros((level + 1, m))
    S = np.zeros((level + 1, m))
    T = np.zeros((level + 1, m))

    Y_left = np.zeros(round(level + 1))
    P2_left = np.zeros(round(level + 1))
    T_left = np.zeros(round(level + 1))

    P1_right = np.zeros(round(level + 1))
    X_right = np.zeros(round(level + 1))
    T_right = np.zeros(round(level + 1))

    s = np.zeros(2)
    t = np.zeros(2)
    x = np.zeros(2)
    y = np.zeros(2)
    ### Конец выделения памяти ------------------------------------------------------

    # Обнуление массивов (отслеживание ошибок) --------------------------------------
    X[:, :] = None
    Y[:, :] = None
    S[:, :] = None
    T[:, :] = None

    Y_left[:] = None
    P2_left[:] = None
    T_left[:] = None

    P1_right[:] = None
    X_right[:] = None
    T_right[:] = None
    # Конец обнуления массивов ------------------------------------------------------

    # Начальные условия -------------------------------------------------------------
    for i in range(m):
        X[0, i] = fun_x0(ds * i)
        Y[0, i] = fun_y0(ds * i)
        S[0, i] = ds * i
        T[0, i] = 0.0

    Y_left[0] = fun_y0(s0)
    P2_left[0] = 0.0
    T_left[0] = 0.0

    P1_right[0] = 0.0
    X_right[0] = fun_x0(s_end)
    T_right[0] = 0.0

    kl = 0
    kr = 0
    # Конец начальных условий --------------------------------------------------------

    # М Е Т О Д ======================================================================
    for l in range(1, level + 1):
        for i in range(m):
            if i == m - 1:
                if was_shift_right(l):
                    s[:] = [S[l - 1, i], s_end]
                    t[:] = [T[l - 1, i], T_right[kr]]
                    x[:] = [X[l - 1, i], X_right[kr]]
                    y[:] = [Y[l - 1, i], P1_right[kr]*G1(T_right[kr])/c2]

                    X[l, i], Y[l, i] = right_sol(i, l, t, s, x, y)
                    kr += 1

                else:
                    s[:] = [S[l - 1, i - 1], S[l - 1, i]]
                    t[:] = [T[l - 1, i - 1], T[l - 1, i]]
                    x[:] = [X[l - 1, i - 1], X[l - 1, i]]
                    y[:] = [Y[l - 1, i - 1], Y[l - 1, i]]

                    X[l, i], Y[l, i] = center_sol(i, l, t, s, x, y)
                    print_rez(S[l, i], T[l, i], [X[l, i], Y[l, i]])
            elif i == 0:
                if not was_shift_right(l):  # ТОЧНО not, проверил
                    s[:] = [s0, S[l - 1, i]]
                    t[:] = [T_left[kl], T[l - 1, i]]
                    x[:] = [P2_left[kl]*G2(T_left[kl])/c1, X[l - 1, i]]
                    y[:] = [Y_left[kl], Y[l - 1, i]]

                    X[l, i], Y[l, i] = left_sol(i, l, t, s, x, y)
                    kl += 1
                else:
                    s[:] = [S[l - 1, i], S[l - 1, i + 1]]
                    t[:] = [T[l - 1, i], T[l - 1, i + 1]]
                    x[:] = [X[l - 1, i], X[l - 1, i + 1]]
                    y[:] = [Y[l - 1, i], Y[l - 1, i + 1]]

                    X[l, i], Y[l, i] = center_sol(i, l, t, s, x, y)
            else:
                if was_shift_right(l):  # ТОЧНО НЕТ not, проверил
                    s[:] = [S[l - 1, i], S[l - 1, i + 1]]
                    t[:] = [T[l - 1, i], T[l - 1, i + 1]]
                    x[:] = [X[l - 1, i], X[l - 1, i + 1]]
                    y[:] = [Y[l - 1, i], Y[l - 1, i + 1]]
                else:
                    s[:] = [S[l - 1, i - 1], S[l - 1, i]]
                    t[:] = [T[l - 1, i - 1], T[l - 1, i]]
                    x[:] = [X[l - 1, i - 1], X[l - 1, i]]
                    y[:] = [Y[l - 1, i - 1], Y[l - 1, i]]

                X[l, i], Y[l, i] = center_sol(i, l, t, s, x, y)
    print(X, Y)
    rez = {
        "s": S[level, :],
        "t1": T[level, 0],
        "x(t1)": X[level, :],
        "y(t1)": Y[level, :],
        "p2(t)": P2_left,
        "y(s0, t)": Y_left,
        "t(s0)": T_left,
        "x(s1, t)": X_right,
        "p1(t)": P1_right,
        "t(s1)": T_right,
    }
    return rez


### Переход из моих координат в координаты задачи ---------------------------------
def get_st_point_level_i(level, i):
    t_rez = h * level
    rs = w * level % ds
    s_rez = ds * i + rs
    return s_rez, t_rez


### Сдвиг (решения проблемы моих координат) ---------------------------------------
def was_shift_right(l):
    return get_st_point_level_i(l, 0)[0] > get_st_point_level_i(l - 1, 0)[0]


### Красивый вывод для решателей --------------------------------------------------
def print_input(s, t, x, y, l, i):
    pass


#     print(f"| level: {l:.0f}, i: {i:.0f} |")
#     print(f'vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv')
#     print(f'input\t left \t right')
#     print(f's:\t {s[0]:.4f}\t {s[1]:.4f}')
#     print(f't:\t {t[0]:.4f}\t {t[1]:.4f}')
#     print(f'x:\t {x[0]:.4f}\t {x[1]:.4f}')
#     print(f'y:\t {y[0]:.4f}\t {y[1]:.4f}\n')


def print_rez(si, ti, rez):
    pass


#     print(f'R E Z ---------------')
#     print(f's:\t {si:.4f}')
#     print(f't:\t {ti:.4f}')
#     print(f'x:\t {rez[0]:.4f}')
#     print(f'y:\t {rez[1]:.4f}')
#     print(f'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')

### Общий случай решения и частный ------------------------------------------------


def sol_xy_(dt_array, X, Y, M1, M2, f):
    dt_l, dt_r = dt_array
    x_l, x_r = X
    y_l, y_r = Y
    f_l, f_r = f
    b = np.zeros(2)
    b[0] = (M1[0, 0] * x_r + M1[0, 1] * y_r + f_r) * dt_r / 2 + x_r
    b[1] = (M1[1, 0] * x_l + M1[1, 1] * y_l + f_l) * dt_l / 2 + y_l

    A = np.zeros((2, 2))
    A[0, 0] = M2[0, 0]
    A[0, 1] = M2[0, 1]
    A[1, 0] = M2[1, 0]
    A[1, 1] = M2[1, 1]

    rez = np.linalg.solve(A, b)
    return rez


def center_sol(i, l, t, s, x, y):
    print_input(s, t, x, y, l, i)

    S[l, i], T[l, i] = get_st_point_level_i(l, i)
    si = S[l, i]
    ti = T[l, i]
    hl = ti - t[0]
    hr = ti - t[1]

    M1 = np.zeros((2, 2))
    M2 = np.zeros((2, 2))
    M1[0, 0] = -B11(s[1], t[1])
    M1[0, 1] = -B21(s[1], t[1])
    M1[1, 0] = -B12(s[0], t[0])
    M1[1, 1] = -B22(s[0], t[0])
    M2[0, 0] = 1 + B11(si, ti) * hr / 2
    M2[0, 1] = 0 + B21(si, ti) * hr / 2
    M2[1, 0] = 0 + B12(si, ti) * hl / 2
    M2[1, 1] = 1 + B22(si, ti) * hl / 2
    fr = F1(s[1], t[1]) + F1(si, ti)
    fl = F2(s[0], t[0]) + F2(si, ti)
    rez = sol_xy_([hl, hr], x, y, M1, M2, [fl, fr])

    print_rez(si, ti, rez)
    return rez


def right_sol(i, l, t, s, x, y):
    """
    s[:] = [S[l-1, i], s_end]
    t[:] = [T[l-1, i], T_right[kr]]
    x[:] = [X[l-1, i], X_right[kr]]
    y[:] = [Y[l-1, i], Y_right[kr]]
    """
    print_input(s, t, x, y, l, i)
    if l == 1:
        x[1] = fun_x0(s_end)
        y[1] = fun_y0(s_end)
        print_input(s, t, x, y, l, i)
        rez = center_sol(i, l, t, s, x, y)
        T_right[kr + 1] = t[1]
        X_right[kr + 1] = x[1]
        P1_right[kr + 1] = y[1]
    else:
        T_right[kr + 1] = t[1] + h_r
        cf = np.zeros(2)
        cf[1] = ((s[0] + c1 * t[0]) - (s[1] - c2 * T_right[kr + 1])) / (c1 + c2)
        if t[0] != t[1]:
            cf[0] = (cf[1] - t[1]) / (t[0] - t[1])
            if cf[0] > 1:
                cf[0] = 1.
        else:
            cf[0] = 1.

        cf[1] = 1. - cf[0]

        tl = cf[1] * t[1] + cf[0] * t[0]
        sl = cf[1] * s[1] + cf[0] * s[0]
        xl = cf[1] * x[1] + cf[0] * x[0]
        yl = cf[1] * y[1] + cf[0] * y[0]

        hl = T_right[kr + 1] - tl
        M1 = np.zeros((2, 2))
        M2 = np.zeros((2, 2))
        M1[0, 0] = -c1
        M1[0, 1] = G1(T_right[kr])
        M1[1, 0] = -B12(sl, tl)
        M1[1, 1] = -B22(sl, tl)

        M2[0, 0] = hl / 2 * B12(s_end, T_right[kr + 1])
        M2[0, 1] = (1 + hl / 2 * B22(s_end, T_right[kr + 1])) / c2 * G1(T_right[kr + 1])
        M2[1, 0] = c1 * h_r/2
        M2[1, 1] = 1 - h_r/2 * G1(T_right[kr + 1])
        fr = D1(T_right[kr + 1]) + D1(T_right[kr])
        fl = F2(s_end, T_right[kr + 1]) + F2(sl, tl)
        X_right[kr + 1], P1_right[kr + 1] = sol_xy_([hl, h_r],
                                                    [xl, x[1]], [yl, y[1]], M1, M2, [fl, fr])

        t[1] = T_right[kr + 1]
        s[1] = s_end
        x[1] = X_right[kr + 1]
        y[1] = P1_right[kr + 1]*G1(t[1])/c2
        print_input(s, t, x, y, l, i)
        rez = center_sol(i, l, t, s, x, y)
    print_rez(S[l, i], T[l, i], rez)

    return rez


def left_sol(i, l, t, s, x, y):
    print_input(s, t, x, y, l, i)
    """
    s[:] = [s0,         S[l-1, i]]
    t[:] = [T_left[kl], T[l-1, i]]
    x[:] = [X_left[kl], X[l-1, i]]
    y[:] = [Y_left[kl], Y[l-1, i]]
    """
    cf = np.zeros(2)
    T_left[kl + 1] = t[0] + h_l
    cf[1] = ((s[0] + c1 * T_left[kl + 1]) - (s[1] - c2 * t[1])) / (c1 + c2)
    if t[1] != t[0]:
        cf[0] = (cf[1] - t[0]) / (t[1] - t[0])
        if cf[0] > 1:
            cf[0] = 1.
    else:
        cf[0] = 1.

    cf[1] = 1. - cf[0]

    tr = cf[1] * t[0] + cf[0] * t[1]
    sr = cf[1] * s[0] + cf[0] * s[1]
    xr = cf[1] * x[0] + cf[0] * x[1]
    yr = cf[1] * y[0] + cf[0] * y[1]

    # hl = T_right[kr + 1] - tl
    # M1 = np.zeros((2, 2))
    # M2 = np.zeros((2, 2))
    # M1[0, 0] = -c1
    # M1[0, 1] = G1(T_right[kr])
    # M1[1, 0] = -B12(sl, tl)
    # M1[1, 1] = -B22(sl, tl)
    #
    # M2[0, 0] = c1 * h_r / 2
    # M2[0, 1] = 1 - h_r / 2 * G1(T_right[kr + 1])
    # M2[1, 0] = hl / 2 * B12(s_end, T_right[kr + 1])
    # M2[1, 1] = (1 + hl / 2 * B22(s_end, T_right[kr + 1])) / c2 * G1(T_right[kr + 1])
    # fr = D1(T_right[kr + 1]) + D1(T_right[kr])
    # fl = F2(s_end, T_right[kr + 1]) + F2(sl, tl)
    # X_right[kr + 1], P1_right[kr + 1] = sol_xy_([hl, h_r],
    #                                             [xl, x[1]], [yl, y[1]], M1, M2, [fl, fr])

    hr = T_left[kl + 1] - tr
    M1 = np.zeros((2, 2))
    M2 = np.zeros((2, 2))
    M1[0, 0] = -B11(sr, tr)
    M1[0, 1] = -B21(sr, tr)
    M1[1, 0] = G2(T_left[kl])
    M1[1, 1] = -c2
    M2[0, 0] = (1 + hr / 2 * B11(s0, T_left[kl + 1])) / c1 * G2(T_left[kl + 1])
    M2[0, 1] = hr / 2 * B21(s0, T_left[kl + 1])
    M2[1, 0] = 1 - h_l / 2 * G2(T_left[kl + 1])
    M2[1, 1] = c2 * h_l/2
    fr = F1(s0, T_left[kl + 1]) + F1(sr, tr)
    fl = D2(T_left[kl + 1]) + D2(T_left[kl])
    P2_left[kl + 1], Y_left[kl + 1] = sol_xy_([h_l, hr],
                                              [x[0], xr], [y[0], yr], M1, M2, [fl, fr])
    t[0] = T_left[kl + 1]
    s[0] = s0
    x[0] = P2_left[kl + 1]*G2(t[0])/c1
    y[0] = Y_left[kl + 1]

    rez = center_sol(i, l, t, s, x, y)
    print_rez(S[l, i], T[l, i], rez)
    return rez
