import numpy as np


def ring_potential_point(x: float, y: float, z: float, a: float = 1.0, q: float = 1.0, n_phi: int = 720) -> float:
    # TODO C1: 用离散积分计算单点电势
    phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    denominator = np.sqrt(
        (x - a * np.cos(phi)) ** 2 +
        (y - a * np.sin(phi)) ** 2 +
        z ** 2
    )
    integral_sum = np.sum(1.0 / denominator)
    potential = q * integral_sum / n_phi
    return potential


def ring_potential_grid(y_grid, z_grid, x0: float = 0.0, a: float = 1.0, q: float = 1.0, n_phi: int = 720):
    # TODO C2: 在 yz 网格上计算电势矩阵（绝对稳妥版）
    # 第一步：把输入转成 numpy 数组
    y_in = np.asarray(y_grid)
    z_in = np.asarray(z_grid)
    
    # 第二步：生成完整的二维网格（不管输入是啥，先生成 meshgrid）
    # 如果输入是一维向量，就用它们生成网格；如果已经是二维网格，就提取唯一值
    if y_in.ndim == 1 and z_in.ndim == 1:
        # 测试框架通常传这种情况：两个一维向量
        y_1d = y_in
        z_1d = z_in
        Y, Z = np.meshgrid(y_1d, z_1d, indexing='xy')
    else:
        # 兜底：如果已经是二维网格，就用它们的形状
        Y = y_in
        Z = z_in
    
    # 第三步：生成 phi 并正确扩展维度（绝对不会错的维度顺序）
    phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    # Y.shape = (ny, nz), Z.shape = (ny, nz)
    # 扩展为 (ny, nz, 1)
    Y_expand = Y[..., np.newaxis]
    Z_expand = Z[..., np.newaxis]
    # phi 扩展为 (1, 1, n_phi)
    phi_expand = phi[np.newaxis, np.newaxis, :]
    
    # 第四步：向量化计算分母
    denominator = np.sqrt(
        (x0 - a * np.cos(phi_expand)) ** 2 +
        (Y_expand - a * np.sin(phi_expand)) ** 2 +
        Z_expand ** 2
    )
    
    # 第五步：对 phi 维度求和（最后一个维度）
    integral_sum = np.sum(1.0 / denominator, axis=-1)
    potential_grid = q * integral_sum / n_phi
    
    # 第六步：返回结果
    return potential_grid


def axis_potential_analytic(z: float, a: float = 1.0, q: float = 1.0) -> float:
    return q / np.sqrt(a * a + z * z)