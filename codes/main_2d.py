"""
境界要素法 (BEM) による2次元定常ポテンシャル流れ解析

このコードは、境界要素法を用いて2次元の定常ポテンシャル流れを解析します。
境界条件として、ノイマン条件（流速指定）とディリクレ条件（ポテンシャル指定）をサポートしています。

Classes
-------
- BoundaryType: 境界条件の種類を定義する列挙型
- BoundaryCondition: 境界要素の境界条件を表すデータクラス
- BoundaryConditions: 境界全体の情報と操作を提供するクラス
- BEMPotentialFlow: 境界要素法による定常ポテンシャル流れ解析を行うクラス
"""
import argparse
from dataclasses import dataclass
from enum import Enum, auto
from time import perf_counter
from typing import Optional, Final

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path  # ポリゴン内外判定用

_INV_2PI: Final[float] = 1.0 / (2.0 * np.pi)

def _normalize_angles(angles: np.ndarray) -> np.ndarray:
    """角度を -π〜π の範囲に正規化する (ベクトル化版)

    Parameters
    ----------
    angles : np.ndarray
        入力角度の配列 (ラジアン)

    Returns
    -------
    normalized_angles : np.ndarray
        正規化された角度の配列 (ラジアン)
    """
    normalized = angles.copy()
    normalized = (normalized + np.pi) % (2.0 * np.pi) - np.pi
    return normalized



class BoundaryType(Enum):
    """境界条件の種類
    NEUMANN: q (流速) が既知、φ (ポテンシャル) が未知
    DIRICHLET: φ (ポテンシャル) が既知、q (流速) が未知
    """
    NEUMANN = auto()
    """流速指定"""
    DIRICHLET = auto()
    """ポテンシャル指定"""

@dataclass
class BoundaryCondition:
    """境界要素の情報: 線分の境界条件

    Parameters
    ----------
    boundary_type : BoundaryType
        境界条件の種類 (NEUMANN または DIRICHLET)
    value : float
        境界条件の値 (流速またはポテンシャル値)

    Notes
    -----
    - 基本的には、流入部ではDIRICHLET条件で正の値を、流出部ではDIRICHLET条件で負の値を指定すること。
    - 壁面などの不透過境界ではNEUMANN条件で0を指定すること。
    """
    boundary_type: BoundaryType
    """境界条件の種類"""
    value: float
    """境界条件の値 (流速またはポテンシャル値)"""

    def __init__(self, boundary_type: BoundaryType, value: float):
        self.boundary_type = boundary_type
        self.value = value

class BoundaryConditions:
    """境界条件オブジェクト

    Parameters
    ----------
    vertices : np.ndarray, optional
        2次元頂点リスト (N x 2)、
        反時計回りに定義することを推奨
    conditions : list of BoundaryCondition, optional
        各辺の境界条件リスト (BoundaryConditionのリスト)
        verticesの数と同じ長さであること
    max_length : float, optional
        境界要素の最大長さ（細分割のため）

    Notes
    -----
    - 閉じたポリゴン内部を流体領域とみなし、各辺を流入・流出部、または不透過境界として定義する.
    - max_lengthが指定された場合、各辺はmax_length以下になるように細分割される.
    その場合、conditionsは各細分割要素に対して元の条件が適用される (複製される).
    """
    def __init__(self, vertices: np.ndarray,
                 conditions: list[BoundaryCondition],
                 max_length: Optional[float] = None):
        self.vertices: np.ndarray
        """2次元頂点リスト (N x 2)"""
        self._original_vertices: Optional[np.ndarray] = None
        """元の頂点リスト、分割していない場合はNone"""
        self.conditions: list[BoundaryCondition]
        """境界条件リスト (BoundaryConditionのリスト)"""

        self._segmentation(vertices, conditions, max_length)
        self.validate()

        self._normals: np.ndarray
        """各要素の法線ベクトル"""
        self._lengths: np.ndarray
        """各要素の長さ"""
        self._init_geometry()

    def original_vertices(self) -> np.ndarray:
        """元の頂点リストを返す"""
        if self._original_vertices is not None:
            return self._original_vertices

        # 分割していない場合は現在の頂点を返す
        return self.vertices

    def count(self) -> int:
        """境界要素の数を返す"""
        return self.vertices.shape[0]

    def endpoints(self, index: int) -> tuple[np.ndarray, np.ndarray]:
        """指定した要素の端点を返す"""
        p1 = self.vertices[index]
        p2 = self.vertices[(index + 1) % self.count()] # 次の頂点（最後は最初に戻る）
        return p1, p2
    def all_endpoints(self) -> tuple[np.ndarray, np.ndarray]:
        """全要素の端点を返す

        Returns
        -------
        np.ndarray, np.ndarray
            n個の始点と終点のリスト
        """
        p1s = self.vertices
        p2s = np.roll(self.vertices, -1, axis=0)
        return p1s, p2s

    def midpoint(self, index: int) -> np.ndarray:
        """指定したインデックスのセグメントの中点を返す

        Parameters
        ----------
        index : int
            セグメントのインデックス

        Returns
        -------
        midpoint : np.ndarray
            セグメントの中点座標 (x, y)
        """
        p1, p2 = self.endpoints(index)
        return (p1 + p2) / 2.0

    def normal(self, index: int) -> np.ndarray:
        """指定したインデックスのセグメントの法線ベクトルを返す

        Parameters
        ----------
        index : int
            セグメントのインデックス

        Returns
        -------
        normal : np.ndarray
            セグメントの法線ベクトル (nx, ny) (単位ベクトル)
        """
        return self._normals[index]
    def normals(self) -> np.ndarray:
        """全セグメントの法線ベクトルを返す"""
        return self._normals

    def tangent(self, index: int) -> np.ndarray:
        """指定したインデックスのセグメントの接線ベクトルを返す

        Parameters
        ----------
        index : int
            セグメントのインデックス

        Returns
        -------
        tangent : np.ndarray
            セグメントの接線ベクトル (tx, ty) (単位ベクトル)
        """
        n = self._normals[index]
        # 接線ベクトルは法線ベクトルを90度回転させたもの
        t = np.array([-n[1], n[0]])
        return t
    def tangents(self) -> np.ndarray:
        """全セグメントの接線ベクトルを返す"""
        normals = self._normals
        tangents = np.zeros_like(normals)
        tangents[:, 0] = -normals[:, 1]
        tangents[:, 1] = normals[:, 0]
        return tangents

    def length(self, index: int) -> float:
        """指定したインデックスのセグメントの長さを返す

        Parameters
        ----------
        index : int
            セグメントのインデックス

        Returns
        -------
        length : float
            セグメントの長さ
        """
        return self._lengths[index]
    def lengths(self) -> np.ndarray:
        """全セグメントの長さを返す"""
        return self._lengths

    def aabb(self) -> tuple[float, float, float, float]:
        """境界の軸平行境界ボックス (AABB) を返す

        Returns
        -------
        min_x : float
            境界ボックスの最小x座標
        max_x : float
            境界ボックスの最大x座標
        min_y : float
            境界ボックスの最小y座標
        max_y : float
            境界ボックスの最大y座標
        """
        min_x = np.min(self.vertices[:, 0])
        max_x = np.max(self.vertices[:, 0])
        min_y = np.min(self.vertices[:, 1])
        max_y = np.max(self.vertices[:, 1])
        return min_x, max_x, min_y, max_y

    def validate(self):
        """verticesとconditionsの整合性チェック"""
        # もしverticesが空ならスキップ
        if self.vertices.size == 0:
            return

        # verticesがNx2であること
        if self.vertices.ndim != 2 or self.vertices.shape[1] != 2:
            raise ValueError("vertices must be a Nx2 array.")
        n = self.vertices.shape[0]

        # conditionsの数が辺の数と一致すること
        if len(self.conditions) != n:
            raise ValueError("Number of boundary conditions must match number of edges."
                             f" Got {len(self.conditions)} conditions for {n} edges.")

    def _segmentation(
            self, vertices: np.ndarray,
            conditions: list[BoundaryCondition],
            max_length: Optional[float] = None):
        """境界の細分割処理"""
        if max_length is None:
            self.vertices = vertices
            self.conditions = conditions
            return
        self._original_vertices = vertices.copy()
        new_vertices = []
        new_conditions = []
        n = vertices.shape[0]

        for i in range(n):
            p1 = vertices[i]
            p2 = vertices[(i + 1) % n]
            cond = conditions[i]

            # セグメントの長さ
            seg_vec = p2 - p1
            seg_length = np.linalg.norm(seg_vec)

            # 分割数
            num_div = max(1, int(np.ceil(seg_length / max_length)))

            for j in range(num_div):
                t = j / num_div
                new_point = (1 - t) * p1 + t * p2
                new_vertices.append(new_point)
                new_conditions.append(cond)

        self.vertices = np.array(new_vertices)
        self.conditions = new_conditions

    def _init_geometry(self):
        """各要素の中点、法線、長さを計算"""
        n = self.count()
        self._normals = np.zeros((n, 2))
        self._lengths = np.zeros(n)

        for i in range(n):
            p1 = self.vertices[i]
            p2 = self.vertices[(i + 1) % n] # 次の頂点（最後は最初に戻る）

            # ベクトルと長さ
            dx = p2[0] - p1[0]
            dy = p2[1] - p1[1]
            length = np.sqrt(dx**2 + dy**2)
            self._lengths[i] = length

            # 法線ベクトル（反時計回りなら (dy, -dx) が外向き法線）
            #   法線ベクトル（進行方向に対して右側を外向きとする定義が多いが、
            #   ここでは反時計回りなら「右側＝内側」「左側＝外側」。
            #   一般的なBEMでは外向き法線を使う。
            #   反時計回りの場合、(dx, dy) -> (dy, -dx) が外向き法線
            self._normals[i] = np.array([dy, -dx]) / length

class BEMPotentialFlow:
    """境界要素法による定常ポテンシャル流れ解析

    以下の手順で計算を行う:
    1. 境界条件の設定 (BoundaryConditionsオブジェクトの作成)
    2. 1で作成したオブジェクトを使用し、本クラスを初期化
    3. `calculate_velocity()` または `calculate_velocity_mesh()` メソッドで
       内部点の速度ベクトルを計算
    """

    def __init__(self, bcs: BoundaryConditions):
        self.bcs: BoundaryConditions
        """境界条件オブジェクト"""
        self._path: Path
        """ポリゴンのPathオブジェクト（境界判定用）"""
        self._is_solved: bool
        """境界積分方程式が解かれたかどうかのフラグ"""

        self.phi_boundary: np.ndarray
        """境界上のポテンシャル値 φ"""
        self.q_boundary: np.ndarray
        """境界上の法線方向速度 q = dφ/dn"""

        # φとqの方程式の係数H,G
        self._H: np.ndarray
        """境界積分方程式の係数行列 H"""
        self._G: np.ndarray
        """境界積分方程式の係数行列 G"""

        # 境界条件を設定して境界積分方程式を解く
        self.setup(bcs)

    def setup(self, bcs: BoundaryConditions) -> None:
        """境界条件を設定して境界積分方程式を解く

        Parameters
        ----------
        bcs : BoundaryConditions
            境界条件オブジェクト
        """
        # 境界条件を設定
        self.bcs = bcs
        self._path = Path(self.bcs.original_vertices())
        self._is_solved = False

        n = self.bcs.count()

        self.phi_boundary = np.zeros(n)
        self.q_boundary = np.zeros(n)

        # 境界積分方程式を解く
        self._init_boundary_integrals()
        self._solve()

    def _init_boundary_integrals(self) -> None:
        """境界積分方程式の係数行列 H, G を初期化する"""
        n = self.bcs.count()
        self._H = np.zeros((n, n))
        self._G = np.zeros((n, n))

        # 対角成分 (i==j) の計算
        # --> H_ii = c(xi) = 1/2
        np.fill_diagonal(self._H, 0.5)
        # --> G_ii: (l_j/(2*pi)) * (1 - ln(l_j/2))
        g_diag = self.bcs.lengths() / (2*np.pi) * (1 - np.log(self.bcs.lengths()/2))
        self._G[np.arange(n), np.arange(n)] = g_diag

        for i in range(n): # 観測点 i
            xi = self.bcs.midpoint(i)

            # 対角要素を取得
            h_ii = self._H[i, i]
            g_ii = self._G[i, i]

            # 各要素の端点
            p1s, p2s = self.bcs.all_endpoints()

            # H_ij: ∂G/∂n = -1/(2π r^2) * (r・n) の積分
            alpha1s = np.arctan2(p1s[:,1] - xi[1], p1s[:,0] - xi[0])
            alpha2s = np.arctan2(p2s[:,1] - xi[1], p2s[:,0] - xi[0])
            angle_diffs = _normalize_angles(alpha2s - alpha1s)
            self._H[i, :] = - _INV_2PI * angle_diffs
            self._H[i, i] = h_ii  # 対角成分を元に戻す

            # 局所座標 t1s, t2s の計算
            tangents = self.bcs.tangents()
            # r・t の計算
            r1s = p1s - xi
            r2s = p2s - xi
            t1s = np.sum(r1s * tangents, axis=1, keepdims=False)
            t2s = np.sum(r2s * tangents, axis=1, keepdims=False)

            # G_ij: G = -1/(2π)*ln(r) の積分
            norm_r1s = np.linalg.norm(r1s, axis=1, keepdims=False)
            norm_r2s = np.linalg.norm(r2s, axis=1, keepdims=False)
            ds = np.sum(r1s * self.bcs.normals(), axis=1, keepdims=False)
            self._G[i, :] = _INV_2PI * (
                + t1s * np.log(norm_r1s) - t2s * np.log(norm_r2s) + (t2s - t1s)
                - ds * angle_diffs)
            self._G[i, i] = g_ii  # 対角成分を元に戻す

    def _solve(self) -> None:
        """境界積分方程式を解く (境界におけるポテンシャル値 φ と法線方向速度 q を求める)"""
        # 連立方程式 Ax = b の構築
        A = np.zeros((self.bcs.count(), self.bcs.count()))
        b = np.zeros(self.bcs.count())
        for j, cond in enumerate(self.bcs.conditions):
            if cond.boundary_type == BoundaryType.NEUMANN: # q (流速) が既知, φ が未知
                # 左辺: H[i,j]*φ_j -> A[i,j] = H[i,j]
                # 右辺: G[i,j]*q_known
                A[:, j] = self._H[:, j]
                b += self._G[:, j] * cond.value

            elif cond.boundary_type == BoundaryType.DIRICHLET: # φ が既知, q が未知
                # 左辺: -G[i,j]*q_j -> A[i,j] = -G[i,j]
                # 右辺: -H[i,j]*φ_known (右辺に移項するため符号反転)
                A[:, j] = -self._G[:, j]
                b -= self._H[:, j] * cond.value

        # Dirichlet条件が1つでもあれば行列は正則になるため、厳密解法 (solve) を試みる。
        try:
            x = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            # 特異行列になる場合は最小二乗解を求める
            print("Warning: Singular matrix encountered. Using least squares solution.")
            x, _, _, _ = np.linalg.lstsq(A, b, rcond=None)

        # 解の割り当て
        for j, cond in enumerate(self.bcs.conditions):
            if cond.boundary_type == BoundaryType.NEUMANN:
                self.phi_boundary[j] = x[j]
                self.q_boundary[j] = cond.value
            elif cond.boundary_type == BoundaryType.DIRICHLET:
                self.phi_boundary[j] = cond.value
                self.q_boundary[j] = x[j]

        self._is_solved = True

    def _calculate_potential(self, points: np.ndarray) -> np.ndarray:
        """
        複数の内部点におけるポテンシャル値 φ(x, y) を計算する

        Parameters
        ----------
        points : np.ndarray
            内部点の座標リスト (M x 2)

        Returns
        -------
        potentials : np.ndarray
            内部点のポテンシャル値リスト (M,)

        Notes
        -----
        - 点が境界外、または境界に極端に近い場合は np.nan を返す
        - pointsは境界内部にあることを前提とする
        """
        # 次元を拡張する
        # points: (M, 2) -> (M, 1, 2)
        # segments (境界要素): (N, 2) -> (1, N, 2)
        # (M, N, 2) の形状で全点 x 境界要素の組み合わせによる演算を行う
        ps = points[:, np.newaxis, :]  # (M, 1, 2)

        # 端点p1,p2を取得: (1, N, 2)
        p1s_2d, p2s_2d = self.bcs.all_endpoints()
        p1s = p1s_2d[np.newaxis, :, :]
        p2s = p2s_2d[np.newaxis, :, :]

        # 端点へのベクトルr1,r2を計算: (M, N, 2)
        r1s = p1s - ps
        r2s = p2s - ps
        # 極端に境界に近い点は後でnanにする
        sq_tol = (1e-3)**2  # 1単位を1mmとして1μmの距離
        sq_r1s = np.sum(r1s**2, axis=2)  # |r1|^2 (M, N)
        sq_r2s = np.sum(r2s**2, axis=2)  # |r2|^2 (M, N)
        close_mask = (sq_r1s < sq_tol) | (sq_r2s < sq_tol)

        # 接線ベクトルと法線ベクトルを取得: (1, N, 2)
        ts = self.bcs.tangents()[np.newaxis, :, :]
        ns = self.bcs.normals()[np.newaxis, :, :]

        # 内積r・t, r・nを計算: (M, N)
        s1s = np.sum(r1s * ts, axis=2)
        s2s = np.sum(r2s * ts, axis=2)
        ds = np.sum(r1s * ns, axis=2)

        # ∠(p2s-point-p1s) を計算: (M, N)
        theta1s = np.arctan2(r1s[:, :, 1], r1s[:, :, 0])
        theta2s = np.arctan2(r2s[:, :, 1], r2s[:, :, 0])
        angle_diffs = _normalize_angles(theta2s - theta1s)

        # 境界値 q, φ を取得: (1, N)
        qs = self.q_boundary[np.newaxis, :]
        phis = self.phi_boundary[np.newaxis, :]

        # ポテンシャル値の計算: (M, N)
        with np.errstate(divide='ignore', invalid='ignore'):
            phi_spec = - qs * (s2s/2 * np.log(sq_r2s) - s1s/2 * np.log(sq_r1s) - (s2s - s1s)) \
                     + (phis - qs * ds) * angle_diffs
        phi = _INV_2PI * np.sum(phi_spec, axis=1)  # (M,)

        # 極端に境界に近い点はnanにする
        phi[close_mask.any(axis=1)] = np.nan
        return phi

    def calculate_potential(self, x: float, y: float) -> float:
        """単一点のポテンシャル値 φ(x, y) を計算

        Parameters
        ----------
        x : float
            内部点のx座標
        y : float
            内部点のy座標

        Returns
        -------
        phi : float
            内部点のポテンシャル値

        Notes
        -----
        - 点が境界外、または境界に極端に近い場合は np.nan を返す
        """
        if not self._is_solved:
            raise RuntimeError("Boundary integral equation is not solved yet.")

        # 内外判定
        if not self._path.contains_point((x, y)):
            return np.nan

        point = np.array([[x, y]])  # (1, 2)
        phi = self._calculate_potential(point)[0]
        return phi

    def calculate_potential_mesh(self, xx: np.ndarray, yy: np.ndarray) -> np.ndarray:
        """格子点群におけるポテンシャル値 φ(x, y) を計算する

        Parameters
        ----------
        xx : np.ndarray
            格子点のx座標の2次元配列
        yy : np.ndarray
            格子点のy座標の2次元配列

        Returns
        -------
        phi_mesh : np.ndarray
            格子点のポテンシャル値の2次元配列
        """
        if not self._is_solved:
            raise RuntimeError("Boundary integral equation is not solved yet.")

        # 格子点を内部点リストに変換
        points = np.column_stack((xx.ravel(), yy.ravel()))  # (M, 2)

        # 内外判定マスクを作成
        inside_mask = self._path.contains_points(points)  # (M,)
        if not np.any(inside_mask):
            # 全点が外部の場合、全てnanを返す
            return np.full(xx.shape, np.nan)

        # 内部点のみ抽出して計算
        internal_points = points[inside_mask]  # (K, 2)
        potentials_internal = self._calculate_potential(internal_points)  # (K,)

        # 結果を元の格子形状に戻す
        potentials = np.full(points.shape[0], np.nan)  # (M,)
        potentials[inside_mask] = potentials_internal
        phi_mesh = potentials.reshape(xx.shape)  # (grid_y, grid_x)

        return phi_mesh

    def _calculate_velocity(self, points: np.ndarray) -> np.ndarray:
        """
        複数の内部点における速度ベクトルを計算する

        Parameters
        ----------
        points : np.ndarray
            内部点の座標リスト (M x 2)

        Returns
        -------
        velocities : np.ndarray
            内部点の速度ベクトルリスト (M x 2)

        Notes
        -----
        - 点が境界外、または境界に極端に近い場合は np.nan を返す
        - pointsは境界内部にあることを前提とする
        """
        # 次元を拡張する
        # points: (M, 2) -> (M, 1, 2)
        # segments (境界要素): (N, 2) -> (1, N, 2)
        # (M, N, 2) の形状で全点 x 境界要素の組み合わせによる演算を行う
        ps = points[:, np.newaxis, :]  # (M, 1, 2)

        # 端点p1,p2を取得: (1, N, 2)
        p1s_2d, p2s_2d = self.bcs.all_endpoints()
        p1s = p1s_2d[np.newaxis, :, :]
        p2s = p2s_2d[np.newaxis, :, :]

        # 端点へのベクトルr1,r2を計算: (M, N, 2)
        r1s = p1s - ps
        r2s = p2s - ps
        # 極端に境界に近い点は後でnanにする
        sq_tol = (1e-3)**2  # 1単位を1mmとして1μmの距離
        sq_r1s = np.sum(r1s**2, axis=2)  # |r1|^2 (M, N)
        sq_r2s = np.sum(r2s**2, axis=2)  # |r2|^2 (M, N)
        close_mask = (sq_r1s < sq_tol) | (sq_r2s < sq_tol)

        # 接線ベクトルと法線ベクトルを取得: (1, N, 2)
        ts = self.bcs.tangents()[np.newaxis, :, :]
        ns = self.bcs.normals()[np.newaxis, :, :]

        # 内積r・t, r・nを計算: (M, N)
        s1s = np.sum(r1s * ts, axis=2)
        s2s = np.sum(r2s * ts, axis=2)
        ds = np.sum(r1s * ns, axis=2)

        # ∠(p2s-point-p1s) を計算: (M, N)
        theta1s = np.arctan2(r1s[:, :, 1], r1s[:, :, 0])
        theta2s = np.arctan2(r2s[:, :, 1], r2s[:, :, 0])
        angle_diffs = _normalize_angles(theta2s - theta1s)

        # 境界値 q, φ を取得: (1, N)
        qs = self.q_boundary[np.newaxis, :]
        phis = self.phi_boundary[np.newaxis, :]

        # 係数計算: (M, N)
        # 境界近傍での警告は無視する
        with np.errstate(divide='ignore', invalid='ignore'):
            t_coef = (qs * np.log(sq_r2s / sq_r1s) / 2 - phis * (ds / sq_r2s - ds / sq_r1s))
            n_coef = (qs * angle_diffs - phis * (s1s / sq_r1s - s2s / sq_r2s))

        # 速度成分の計算: (M, N, 2)
        v_spec = t_coef[:, :, np.newaxis] * ts + n_coef[:, :, np.newaxis] * ns
        v = _INV_2PI * np.sum(v_spec, axis=1)  # (M, 2)

        # 極端に境界に近い点はnanにする
        v[close_mask.any(axis=1), :] = np.nan

        return v

    def calculate_velocity(self, x: float, y: float) -> tuple[float, float]:
        """単一点の速度ベクトル (u, v) を計算

        Parameters
        ----------
        x : float
            内部点のx座標
        y : float
            内部点のy座標

        Returns
        -------
        u : float
            内部点のx方向速度成分
        v : float
            内部点のy方向速度成分

        Notes
        -----
        - 点が境界外、または境界に極端に近い場合は (np.nan, np.nan) を返す
        """
        if not self._is_solved:
            raise RuntimeError("Boundary integral equation is not solved yet. "
                               "Call the 'solve()' method first.")

        point = np.array([x, y])

        # ポリゴン外判定
        if not self._path.contains_point((x, y)):
            return np.nan, np.nan

        velocity = self._calculate_velocity(point[np.newaxis, :])  # (1, 2)
        return velocity[0, 0], velocity[0, 1]

    def calculate_velocity_mesh(self, xx: np.ndarray, yy: np.ndarray
                                ) -> tuple[np.ndarray, np.ndarray]:
        """
        内部点の速度ベクトル (u, v) を計算

        Parameters
        ----------
        xx : np.ndarray
            内部点のx座標グリッド (2D配列)
        yy : np.ndarray
            内部点のy座標グリッド (2D配列)

        Returns
        -------
        u_grid : np.ndarray
            内部点のx方向速度成分グリッド (2D配列)
            サイズはX,Yと同じ
        v_grid : np.ndarray
            内部点のy方向速度成分グリッド (2D配列)
            サイズはX,Yと同じ

        Notes
        -----
        - 点が境界外、または境界に極端に近い場合は np.nan を設定
        """
        u_grid = np.zeros_like(xx)
        v_grid = np.zeros_like(yy)

        rows, cols = xx.shape

        # グリット状の点を (N, 2) 配列に変換
        points = np.column_stack((xx.ravel(), yy.ravel()))

        # ポリゴン内部の点のみ抽出
        is_inside = self._path.contains_points(points)

        if not np.any(is_inside):
            # 全て外部点の場合はnanを返す
            return np.full(xx.shape, np.nan), np.full(yy.shape, np.nan)

        # 速度計算
        u_flat = np.full(points.shape[0], np.nan)
        v_flat = np.full(points.shape[0], np.nan)
        velocities = self._calculate_velocity(points[is_inside])
        u_flat[is_inside] = velocities[:, 0]
        v_flat[is_inside] = velocities[:, 1]

        # 元のグリッド形状に戻す
        u_grid = u_flat.reshape((rows, cols))
        v_grid = v_flat.reshape((rows, cols))
        return u_grid, v_grid




# --- 描画関数 ---

def _plot_result(xx: np.ndarray, yy: np.ndarray,
                 poly: np.ndarray,
                 pp: Optional[np.ndarray] = None,
                 uu_vv: Optional[tuple[np.ndarray, np.ndarray]] = None,
                 fluid_color: str = 'jet',
                 skip: int = 2) -> None:
    """計算結果のプロット (ポテンシャル場と流線)

    Parameters
    ----------
    xx : np.ndarray
        x座標グリッド (2D配列)
    yy : np.ndarray
        y座標グリッド (2D配列)
    poly : np.ndarray
        境界の頂点リスト ((N+1) x 2 配列)
        終点を含むこと (N+1番目の要素は最初の頂点と同じ点であること)
    pp : Optional[np.ndarray], optional
        ポテンシャル場グリッド (2D配列), by default None
    uu_vv : Optional[tuple[np.ndarray, np.ndarray]], optional
        x,y方向速度成分グリッドのタプル (uu, vv), by default None
    fluid_color : str, optional
        流体の色マップ (デフォルトは 'jet')
        'jet'以外が指定された場合は単色として扱う
    skip : int, optional
        ベクトル場の表示間引き率 (デフォルトは2)
        例えばxxが100x100の場合、skip=2なら50x50のベクトルが表示される

    Notes
    -----
    - `xx`, `yy` と `pp` または `uuvv` の中の配列は同じサイズの2D配列であること
    """
    if pp is None and uu_vv is None:
        print("Warning: No data to plot.")
        return

    num_plots = (1 if pp is not None else 0) + (1 if uu_vv is not None else 0)
    plot_index = 1

    x_width = xx.max() - xx.min()
    y_width = yy.max() - yy.min()

    fig_width = 10
    if x_width > y_width:
        fig_size = (fig_width, fig_width * num_plots * y_width / x_width)
    else:
        fig_size = (fig_width * x_width / y_width, fig_width * num_plots)

    plt.figure(figsize=fig_size)

    # ポテンシャル場の描画
    if pp is not None:
        plt.subplot(num_plots, 1, plot_index)
        plt.contourf(xx, yy, pp, cmap='jet', levels=50)
        plt.plot(poly[:,0], poly[:,1], 'k-', linewidth=2, label='Boundary')
        plt.title("Potential Field")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.axis('equal')
        plt.colorbar(label='Potential')
        plt.grid(True, linestyle=':', alpha=0.6)
        plot_index += 1

    # 流線の描画
    if uu_vv is not None:
        uu, vv = uu_vv
        plt.subplot(num_plots, 1, plot_index)
        # nanを0に置換
        uu = np.nan_to_num(uu, nan=0.0)
        vv = np.nan_to_num(vv, nan=0.0)
        color_map = None
        color = fluid_color
        if fluid_color == 'jet':
            # 色づけのために速度の大きさを計算
            speed = np.log(np.sqrt(uu**2 + vv**2) + 1)
            color = speed
            color_map = fluid_color
        plt.streamplot(xx, yy, uu, vv, color=color,
                       cmap=color_map, density=2, linewidth=2)

        # 境界の描画
        plt.plot(poly[:,0], poly[:,1], 'k-', linewidth=2, label='Boundary')

        # ベクトル場（間引いて表示）
        plt.quiver(xx[::skip, ::skip], yy[::skip, ::skip],
                uu[::skip, ::skip], vv[::skip, ::skip],
                scale=100, color='black', alpha=0.2)

        plt.title("Streamlines")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.axis('equal')
        if fluid_color == 'jet':
            plt.colorbar(label='Log Speed')
        plt.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.show()



# --- メイン実行部分 ---
# エイリアス定義
_BT = BoundaryType
_BC = BoundaryCondition
_VERTEX = tuple[float, float]
_VERTICES = list[_VERTEX]

def _get_boundary_conditions(outer_vertices: _VERTICES,
                             io_values: dict[tuple[_VERTEX, _VERTEX], float],
                             max_length: Optional[float] = None
                             ) -> BoundaryConditions:
    """境界条件オブジェクトを作成するユーティリティ関数

    Parameters
    ----------
    vertices : _VERTICES
        境界の頂点リスト (反時計回り)
    io_values : dict[tuple[_VERTEX, _VERTEX], float]
        Dirichlet境界条件を指定する辞書.
        キーは (x1, y1), (x2, y2) のタプルで境界要素の端点座標を表す.
        値はその境界要素に対応するDirichlet条件の値を表す.
        ここで指定した要素以外はNeumann条件 (q=0) として扱う.
    max_length : Optional[float], optional
        境界要素の最大長さ (デフォルトは None, 指定しない場合は分割しない)

    Returns
    -------
    bc_obj : BoundaryConditions
        境界条件オブジェクト

    Raises
    ------
    ValueError
        指定された境界要素が頂点リストに存在しない場合
    """
    n = len(outer_vertices)
    conditions = []
    o_verts = np.array(outer_vertices)

    # io_valuesのキーを検索しやすい形式に変換
    # (x1, y1, x2, y2) -> frozenset({(x1, y1), (x2, y2)})
    # これにより、(p1, p2) と (p2, p1) を同じキーとして扱える
    io_map = {frozenset({k[0], k[1]}): v for k, v in io_values.items()}

    # 各辺の境界条件を決定
    for i in range(n):
        p1 = tuple(o_verts[i])
        p2 = tuple(o_verts[(i + 1) % n])
        segment_key = frozenset({p1, p2})

        if segment_key in io_map:
            # io_valuesに指定された辺はDirichlet条件
            conditions.append(_BC(_BT.DIRICHLET, io_map[segment_key]))
        else:
            # それ以外はNeumann条件 (壁)
            conditions.append(_BC(_BT.NEUMANN, 0.0))

    # 整合性チェック: io_valuesで指定された辺がすべて使われたか
    used_keys = {frozenset({tuple(o_verts[i]), tuple(o_verts[(i + 1) % n])}) for i in range(n)}
    for key in io_map:
        if key not in used_keys:
            p1, p2 = tuple(key)
            raise ValueError(f"Boundary segment {p1} -> {p2} from io_values not found in vertices list.")

    return BoundaryConditions(o_verts, conditions, max_length)

def _boundary_example_1():
    """コの字型パイプの境界条件例"""
    # 形状定義（矩形水槽の中に障害物があるような単純な形状、あるいはパイプ）
    # ここでは単純な「コ」の字型パイプを定義（反時計回り）
    vertices = [
        (0, 0), (5, 0), (5, 2), (3, 2),
        (3, 1), (2, 1), (2, 2), (0, 2)
    ]

    # 境界条件の設定
    # NOTE: 質量保存則 (積分値の合計=0) が厳密には必要
    # NOTE: 電磁気学等のポテンシャル場と異なり、ポテンシャル流れv=∇φでは
    #       ポテンシャルの低い方から高い方へ流れることに注意. したがって、
    #       流入部は低い値、流出部は高い値を指定する.
    bc_list = [
        _BC(_BT.NEUMANN,    0.0),  # 0: (0,0)->(5,0) 底面 (壁)
        _BC(_BT.NEUMANN,    0.0),  # 1: (5,0)->(5,2) 右面 (壁)
        _BC(_BT.DIRICHLET, -1.0),  # 2: (5,2)->(3,2) 上面右 (流入)
        _BC(_BT.NEUMANN,    0.0),  # 3: (3,2)->(3,1) 凹み右 (壁)
        _BC(_BT.NEUMANN,    0.0),  # 4: (3,1)->(2,1) 凹み底 (壁)
        _BC(_BT.NEUMANN,    0.0),  # 5: (2,1)->(2,2) 凹み左 (壁)
        _BC(_BT.DIRICHLET,  1.0),  # 6: (2,2)->(0,2) 上面左 (流出)
        _BC(_BT.NEUMANN,    0.0),  # 7: (0,2)->(0,0) 左面 (壁)
    ]
    return BoundaryConditions(np.array(vertices), bc_list, max_length=0.1)

def _boundary_example_2():
    """複雑な形状の境界条件例"""
    vertices = [
        (0, 2), (1, 2), (1, 1), (1, 0), (2, 0), (2, 1), (3, 1), (4, 0),
        (6, 0), (6, 1), (5, 1), (4, 3), (3, 3), (2.5, 2), (2, 2), (1, 3), (0, 3)
    ]
    return _get_boundary_conditions(
        vertices,
        io_values = {
            ((1, 0), (2, 0)): -1.0,  # (1,0)->(2,0) 流入部
            ((6, 0), (6, 1)): -1.0,  # (6,0)->(6,1) 流入部
            ((0, 3), (0, 2)):  2.0   # (0,3)->(0,2) 流出部
        },
        max_length=0.1
    )

def _parse_args():
    """コマンドライン引数を解析する. 計算を行う必要がない場合にはquit()"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--ex", type=int, choices=[1, 2],
                        default=1,
                        help="境界条件の例番号 (1または2)")
    parser.add_argument("--pot", action="store_true",
                        help="ポテンシャル場を表示")
    parser.add_argument("--vel", action="store_true",
                        help="速度場を表示")
    args = parser.parse_args()
    if not args.pot and not args.vel:
        print("Warning: At least one of --pot or --vel must be specified.")
        quit()
    return args

def main():
    """メイン実行関数: 境界条件設定、計算、描画"""
    args = _parse_args()

    print("Setting up boundary conditions and solving BEM...")
    start_time = perf_counter()
    if args.ex == 1:
        bcs = _boundary_example_1()
    else:
        bcs = _boundary_example_2()
    print(f"  -> Boundary conditions set up. ({perf_counter() - start_time:.2f} sec)")

    # 1. 境界積分方程式を解く
    print("Solving BEM...")
    start_time = perf_counter()
    bem = BEMPotentialFlow(bcs)
    print(f"  -> BEM solved. ({perf_counter() - start_time:.2f} sec)")

    # 2.1 グリッド作成
    n_grids = 400
    min_x, max_x, min_y, max_y = bcs.aabb()
    x_range = np.linspace(min_x - 0.5, max_x + 0.5, n_grids)
    y_range = np.linspace(min_y - 0.5, max_y + 0.5, n_grids)
    xx, yy = np.meshgrid(x_range, y_range)
    total_points = xx.size

    # 2.2 各点でのポテンシャル値 φ を計算
    pp = None
    if args.pot:
        print("Calculating internal potential field...")
        start_time = perf_counter()
        pp = bem.calculate_potential_mesh(xx, yy)
        total_time = perf_counter() - start_time
        print(f"  -> Internal potential field calculated. ({total_time:.2f} sec; "
            f"{total_time / total_points * 1e3:.2f} ms/point)")

    # 2.3 各点での速度ベクトル (u, v) を計算
    uu_vv = None
    if args.vel:
        print("Calculating internal velocity field...")
        start_time = perf_counter()
        uu_vv = bem.calculate_velocity_mesh(xx, yy)
        total_time = perf_counter() - start_time
        print(f"  -> Internal velocity field calculated. ({total_time:.2f} sec; "
          f"{total_time / total_points * 1e3:.2f} ms/point)")

    # 3. 結果のプロット
    print("Plotting results...")
    skip = max(1, n_grids // 50)  # 最大50x50のベクトル表示
    _plot_result(xx, yy, np.vstack([bcs.vertices, bcs.vertices[0]]),
                 pp = pp, uu_vv=uu_vv, skip=skip)

if __name__ == "__main__":
    main()
