import numpy as np


def basic_musk(I,
               K=3,
               X=0.1,
               delta_t=1):
    """
    基础MUSK方程
    :param I: 入流流量列表，每个元素代表一个时间步的入流流量
    :param K: 蓄量系数
    :param X: 流量比重系数
    :param delta_t:
    :return: 出流流量列表，河段蓄水量列表
    """
    O = [I[0]]  # 初始化出流流量列表，第一个值设为初始入流流量
    S = K * (X * I[0] + (1 - X) * O[0])  # 初始化河段蓄水量
    S_process = [S]

    for i in range(1, len(I)):
        dI_dt = (I[i] - I[i - 1]) / delta_t  # 前向差分近似计算入流流量对时间的导数

        # 根据非线性槽蓄方程和水量平衡方程推导的公式来计算下一个时间步的出流流量
        O_next = ((I[i] - K * X * I[i] * dI_dt) / (1 + K * (1 - X)))
        O.append(O_next)

        # 更新河段蓄水量
        S = K * (X * I[i] + (1 - X) * O[i])
        S_process.append(S)

    return O, S_process


def nonlinear_musk(I,
                   K=3,
                   X=0.1,
                   n=2,
                   delta_t=1):
    """
    基于非线性槽蓄方程的MUSK

    :param I: 入流流量列表，每个元素代表一个时间步的入流流量
    :param K: 蓄量系数
    :param X: 流量比重系数
    :param n: 非线性指数
    :param delta_t: 时间步长
    :return: 出流流量列表，河段蓄水量列表
    """
    O = [I[0]]  # 初始化出流流量列表，第一个值设为初始入流流量
    S = K * (X * I[0] ** n + (1 - X) * O[0] ** n)  # 初始化河段蓄水量
    S_process = [S]

    for i in range(1, len(I)):
        dI_dt = (I[i] - I[i - 1]) / delta_t  # 前向差分近似计算入流流量对时间的导数

        # 根据非线性槽蓄方程和水量平衡方程推导的公式来计算下一个时间步的出流流量
        O_next = ((I[i] - K * n * X * I[i] ** (n - 1) * dI_dt) / (1 + K * n * (1 - X) * O[i - 1] ** (n - 1)))
        O.append(O_next)

        # 更新河段蓄水量
        S = K * (X * I[i] ** n + (1 - X) * O[i] ** n)
        S_process.append(S)

    return O, S_process


def variable_musk(I,
                  K_values,
                  X_values,
                  delta_t=1):
    """
    基于变系数模型的MUSK

    :param I: 入流流量列表，每个元素代表一个时间步的入流流量
    :param K_values: 与入流流量对应的蓄量系数列表，每个元素对应一个时间步的蓄量系数
    :param X_values: 与入流流量对应的流量比重系数列表，每个元素对应一个时间步的流量比重系数
    :param delta_t: 时间步长
    :return: 出流流量列表，河段蓄水量列表
    """
    O = [I[0]]  # 初始化出流流量列表，第一个值设为初始入流流量
    S = K_values[0] * (X_values[0] * I[0] + (1 - X_values[0]) * O[0])  # 初始化河段蓄水量
    S_process = [S]

    for i in range(1, len(I)):
        dI_dt = (I[i] - I[i - 1]) / delta_t  # 前向差分近似计算入流流量对时间的导数

        # 计算变系数的变化率
        dK_dt = (K_values[i] - K_values[i - 1]) / delta_t
        dX_dt = (X_values[i] - X_values[i - 1]) / delta_t

        # 根据变系数模型的槽蓄方程和水量平衡方程推导的公式来计算下一个时间步的出流流量
        O_next = ((I[i] - (dK_dt * (X_values[i] * I[i] + (1 - X_values[i]) * O[i - 1])) -
                   K_values[i] * (dX_dt * I[i] + X_values[i] * dI_dt - dX_dt * O[i - 1])) /
                  (1 + K_values[i] * (1 - X_values[i]) * (O[i - 1] - O[i - 2]) / delta_t))

        O.append(O_next)

        # 更新河段蓄水量
        S = K_values[i] * (X_values[i] * I[i] + (1 - X_values[i]) * O[i])
        S_process.append(S)

    return O, S_process


def variable_nonlinear_musk(I,
                            n_values,
                            K=3,
                            X=0.1,
                            delta_t=1):
    """
    基于变指数非线性MUSK

    :param I: 入流流量列表，每个元素代表一个时间步的入流流量
    :param n_values: 与入流流量对应的指数值列表，每个元素对应一个时间步的指数值
    :param K: 蓄量系数
    :param X: 流量比重系数
    :param delta_t: 时间步长
    :return: 出流流量列表，河段蓄水量列表
    """
    O = [I[0]]  # 初始化出流流量列表，第一个值设为初始入流流量
    S = K * (X * I[0] ** n_values[0] + (1 - X) * O[0] ** n_values[0])  # 初始化河段蓄水量
    S_process = [S]

    for i in range(1, len(I)):
        dI_dt = (I[i] - I[i - 1]) / delta_t  # 前向差分近似计算入流流量对时间的导数

        # 根据变指数非线性槽蓄方程和水量平衡方程推导的公式来计算下一个时间步的出流流量
        O_next = ((I[i] - K * n_values[i] * X * I[i] ** (n_values[i] - 1) * dI_dt) /
                  (1 + K * n_values[i] * (1 - X) * O[i - 1] ** (n_values[i] - 1)))

        O.append(O_next)

        # 更新河段蓄水量
        S = K * (X * I[i] ** n_values[i] + (1 - X) * O[i] ** n_values[i])
        S_process.append(S)

    return O, S_process


def nonlinear_with_continues_variable_musk(I,
                                           K=3,
                                           X=0.1,
                                           A=0.1,
                                           B=0.4,
                                           C=980,
                                           D=4.3,
                                           delta_t=1):
    """
    @article{
        罗宇轩2021变参数非线性马斯京根分段演算模型研究与应用,
        title={变参数非线性马斯京根分段演算模型研究与应用},
        author={罗宇轩 and 陈华 and 林康聆 and 王俊 and 王金星},
        journal={人民长江},
        year={2021},
    }

    :param I: 入流流量列表，每个元素代表一个时间步的入流流量
    :param K: 蓄量系数
    :param X: 流量比重系数
    :param A: 计算指数N所需参数
    :param B: 计算指数N所需参数
    :param C: 计算指数N所需参数
    :param D: 入流指数，取正数
    :param delta_t: 时间步长
    :return: 出流流量列表，河段蓄水量列表
    """
    O = [I[0]]  # 初始化出流流量列表，第一个值设为初始入流流量

    S = K * (X * I[0] + (1 - X) * O[0])  # 初始化河段蓄水量
    S_process = [S]
    for i in range(1, len(I)):
        # 无量纲入流变量 u_t = I_t / I_max
        u_t = I[i] / np.max(I)

        # 根据公式4计算 k(u_t)
        k_u_t = K / (u_t ** D)

        # 根据公式3计算 β(u_t)
        beta_u_t = A + B * np.log(1 + C * u_t)

        dI_dt = (I[i] - I[i - 1]) / delta_t

        # 根据变参数非线性槽蓄方程和水量平衡方程推导的公式计算下一个时间步的出流流量
        O_next = np.real(((I[i] - k_u_t * X * dI_dt) / (1 + k_u_t * (1 - X))) ** (1 / beta_u_t + 0j))
        O.append(O_next)

        # 更新河段蓄水量
        S = k_u_t * (X * I[i] + (1 - X) * O[i]) ** beta_u_t
        S_process.append(S)

    return O, S_process
