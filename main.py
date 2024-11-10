import numpy as np
from musk_algorithm import *

if __name__ == '__main__':
    # 基础MUSK参数
    # 入流序列
    input_stream = [100, 200, 300, 400, 500,
                    500, 400, 300, 200, 100]
    # 入流步长
    delta_t = 1
    # 蓄量系数
    K = 3
    # 流量比重系数
    X = 0.1

    # 基础MUSK
    basic_output_stream, basic_S_process = basic_musk(input_stream, K, X, delta_t)
    print("Basic Output: {}".format(basic_output_stream))
    print("Basic S Process: {}".format(basic_S_process))

    # 非线性MUSK
    n = 2
    nonlinear_output_stream, nonlinear_S_process = nonlinear_musk(input_stream, K, X, n, delta_t)
    print("Nonlinear Output: {}".format(nonlinear_output_stream))
    print("Nonlinear S Process: {}".format(nonlinear_S_process))

    # 变系数MUSK
    K_values = [3, 3, 3, 3, 3,
                3, 3, 3, 3, 3]
    X_values = [0.1, 0.1, 0.1, 0.1, 0.1,
                0.1, 0.1, 0.1, 0.1, 0.1]
    variable_output_stream, variable_S_process = variable_musk(input_stream, K_values, X_values, delta_t)
    print("Variable Output: {}".format(variable_output_stream))
    print("Variable S Process: {}".format(variable_S_process))

    # 变指数非线性MUSK
    n_values = [2, 2, 2, 2, 2,
                2, 2, 2, 2, 2]
    variable_nonlinear_output_stream, variable_nonlinear_S_process = variable_nonlinear_musk(input_stream, n_values, K,
                                                                                             X, delta_t)

    # 变参数非线性MUSK
    A = 0.1
    B = 0.4
    C = 980
    D = 4.3
    nonlinear_with_continues_variable_output_stream, nonlinear_with_continues_variable_S_process = nonlinear_with_continues_variable_musk(
        input_stream, K, X,
        A, B, C, D,
        delta_t)
    print("Nonlinear with Continues Variable Output: {}".format(nonlinear_with_continues_variable_output_stream))
    print("Nonlinear with Continues Variable S Process: {}".format(nonlinear_with_continues_variable_S_process))
