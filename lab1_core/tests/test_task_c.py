import numpy as np


def ring_potential_point(x: float, y: float, z: float, a: float = 1.0, q: float = 1.0, n_phi: int = 720) -> float:
    # TODO C1: 用离散积分计算单点电势（已修复奇点）
    phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    
    # 【关键修复】：在分母的根号里加入极小值 eps=1e-12，防止除零
    # 这是天体物理和计算物理中处理奇点的标准方法
    eps = 1e-12
    denominator = np.sqrt(
        (x - a * np.cos(phi)) ** 2 +
        (y - a * np.sin(phi)) ** 2 +
        z ** 2 +
        eps  # 软化因子在此
    )
    
    integral_sum = np.sum(1.0 / denominator)
    potential = q * integral_sum / n_phi
    return potential


def ring_potential_grid(y_grid, z_grid, x0: float = 0.0, a: float = 1.0, q: float = 1.0, n_phi: int = 720):
    # TODO C2: 在 yz 网格上计算电势矩阵（已修复奇点+完全匹配测试）
    Y, Z = np.meshgrid(y_grid, z_grid, indexing='xy')
    
    phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    
    # 扩展维度
    Y_expand = Y[..., np.newaxis]
    Z_expand = Z[..., np.newaxis]
    phi_expand = phi[np.newaxis, np.newaxis, :]
    
    # 【关键修复】：同样加入软化因子 eps=1e-12
    eps = 1e-12
    denominator = np.sqrt(
        (x0 - a * np.cos(phi_expand)) ** 2 +
        (Y_expand - a * np.sin(phi_expand)) ** 2 +
        Z_expand ** 2 +
        eps
    )
    
    integral_sum = np.sum(1.0 / denominator, axis=-1)
    potential_grid = q * integral_sum / n_phi
    
    return potential_grid


def axis_potential_analytic(z: float, a: float = 1.0, q: float = 1.0) -> float:
    return q / np.sqrt(a * a + z * z)