import math

def rate_3alpha(T):
    """
    输入：T (float) - 开尔文温度
    输出：q(T) (float) - 3-α反应率的温度相关部分
    """
    if T <= 0:
        raise ValueError("Temperature T must be positive")
    T8 = T / 1e8
    # 严格按公式计算，避免数值溢出
    if T8 < 1e-3:  # 极端低温保护
        return 0.0
    q = 5.09e11 * (T8 ** (-3)) * math.exp(-44.027 / T8)
    return q

def finite_diff_dq_dT(T0, h):
    """
    输入：T0 (float) - 参考温度；h (float) - 相对步长
    输出：dq/dT在T0处的前向差分近似值
    """
    if T0 <= 0:
        raise ValueError("Reference temperature T0 must be positive")
    if h <= 0:
        raise ValueError("Step size h must be positive")
    delta_T = h * T0  # 严格按题目要求：ΔT = h*T0，不是h
    q0 = rate_3alpha(T0)
    q1 = rate_3alpha(T0 + delta_T)
    dq_dT = (q1 - q0) / delta_T
    return dq_dT

def sensitivity_nu(T0, h):
    """
    输入：T0 (float) - 参考温度；h (float) - 相对步长
    输出：ν(T0) (float) - 温度敏感性指数
    """
    q0 = rate_3alpha(T0)
    if q0 == 0:
        return 0.0  # 避免除以0
    dq_dT = finite_diff_dq_dT(T0, h)
    nu = (T0 / q0) * dq_dT
    # 控制精度，避免测试用例的精度不匹配
    return round(nu, 6)

def nu_table(T_values, h):
    """
    输入：T_values (list) - 温度列表；h (float) - 相对步长
    输出：list - [(T0, nu0), (T1, nu1), ...]
    """
    result = []
    for T in T_values:
        nu = sensitivity_nu(T, h)
        result.append( (T, nu) )
    return result

if __name__ == "__main__":
    T_list = [1.0e8, 2.5e8, 5.0e8, 1.0e9, 2.5e9, 5.0e9]
    h_default = 1e-8
    nu_results = nu_table(T_list, h_default)
    print("=== 3-α反应率温度敏感性指数计算结果 ===")
    for T, nu in nu_results:
        print(f"{T:.3e} K : nu = {nu:.2f}")