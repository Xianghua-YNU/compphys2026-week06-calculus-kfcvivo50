import math


def debye_integrand(x: float) -> float:
    if abs(x) < 1e-12:
        return 0.0
    ex = math.exp(x)
    return (x**4) * ex / ((ex - 1.0) ** 2)


def trapezoid_composite(f, a: float, b: float, n: int) -> float:
    # 复合梯形积分实现
    h = (b - a) / n
    result = 0.5 * (f(a) + f(b))
    for k in range(1, n):
        result += f(a + k * h)
    result *= h
    return result


def simpson_composite(f, a: float, b: float, n: int) -> float:
    # 复合 Simpson 积分实现，检查 n 为偶数
    if n % 2 != 0:
        raise ValueError("n must be even for Simpson's rule")
    h = (b - a) / n
    result = f(a) + f(b)
    # 计算奇数项（乘4）
    for k in range(1, n, 2):
        result += 4 * f(a + k * h)
    # 计算偶数项（乘2）
    for k in range(2, n, 2):
        result += 2 * f(a + k * h)
    result *= h / 3
    return result


def debye_integral(T: float, theta_d: float = 428.0, method: str = "simpson", n: int = 200) -> float:
    # 计算 Debye 积分 I(theta_d/T)
    y = theta_d / T
    if method.lower() == "trapezoid":
        return trapezoid_composite(debye_integrand, 0.0, y, n)
    elif method.lower() == "simpson":
        # 确保n为偶数
        if n % 2 != 0:
            n += 1
        return simpson_composite(debye_integrand, 0.0, y, n)
    else:
        raise ValueError("method must be 'trapezoid' or 'simpson'")
