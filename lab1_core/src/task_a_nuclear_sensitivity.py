import math

# 1. 实现rate_3alpha(T)：计算3-α反应率q(T)
def rate_3alpha(T):
    """
    输入：T (float) - 开尔文温度
    输出：q(T) (float) - 3-α反应率的温度相关部分
    """
    if T <= 0:
        raise ValueError("温度T必须大于0")
    T8 = T / 1e8  # 题目定义T8 = T / 10^8
    # 核心公式：q = 5.09e11 * T8^(-3) * exp(-44.027 / T8)
    q = 5.09e11 * (T8 ** (-3)) * math.exp(-44.027 / T8)
    return q

# 2. 实现finite_diff_dq_dT(T0, h)：前向差分近似dq/dT
def finite_diff_dq_dT(T0, h):
    """
    输入：T0 (float) - 参考温度；h (float) - 相对步长（默认1e-8）
    输出：dq/dT在T0处的前向差分近似值
    """
    if T0 <= 0:
        raise ValueError("参考温度T0必须大于0")
    if h <= 0:
        raise ValueError("步长h必须大于0")
    delta_T = h * T0  # 绝对步长 = h * T0（题目要求，不是h本身）
    q0 = rate_3alpha(T0)
    q1 = rate_3alpha(T0 + delta_T)
    dq_dT = (q1 - q0) / delta_T
    return dq_dT

# 3. 实现sensitivity_nu(T0, h)：计算温度敏感性指数ν
def sensitivity_nu(T0, h):
    """
    输入：T0 (float) - 参考温度；h (float) - 相对步长（默认1e-8）
    输出：ν(T0) (float) - 温度敏感性指数
    """
    q0 = rate_3alpha(T0)
    dq_dT = finite_diff_dq_dT(T0, h)
    nu = (T0 / q0) * dq_dT  # 核心公式：ν = (T/q) * dq/dT
    return nu

# 4. 实现nu_table(T_values, h)：批量计算多个温度点的ν
def nu_table(T_values, h):
    """
    输入：T_values (list) - 温度列表（单位K）；h (float) - 相对步长
    输出：list - [(T0, nu0), (T1, nu1), ...]
    """
    result = []
    for T in T_values:
        nu = sensitivity_nu(T, h)
        result.append( (T, nu) )
    return result

# ------------------- 主程序：运行计算并输出结果 -------------------
if __name__ == "__main__":
    # 题目要求的必算温度点
    T_list = [1.0e8, 2.5e8, 5.0e8, 1.0e9, 2.5e9, 5.0e9]
    h_default = 1e-8  # 题目建议的默认h值

    # 计算所有温度点的ν
    nu_results = nu_table(T_list, h_default)

    # 按题目要求的格式输出结果
    print("=== 3-α反应率温度敏感性指数计算结果 ===")
    for T, nu in nu_results:
        # 格式：1.000e+08 K : nu = 41.03
        print(f"{T:.3e} K : nu = {nu:.2f}")