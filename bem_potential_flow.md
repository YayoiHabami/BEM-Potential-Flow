# ポテンシャル流れ

> WARN: GitHubのMarkdownレンダラーでは数式が正しく表示されない場合があります。特に太字 (`\bm`や`\boldsymbol`) が太字として表示されないことがあります。 数式を正しく表示するには、VSCodeのMarkdownプレビューなどを利用してください。

## 目次

- [目次](#目次)
- [1. シミュレーションの概要](#1-シミュレーションの概要)
  - [1.1 本シミュレーションの目的と解析フロー](#11-本シミュレーションの目的と解析フロー)
  - [1.2 ポテンシャル流れとは](#12-ポテンシャル流れとは)
  - [1.3 BEMとは](#13-bemとは)
- [2. 境界積分方程式の定式化](#2-境界積分方程式の定式化)
  - [2.1 支配方程式](#21-支配方程式)
  - [2.2 境界条件の設定](#22-境界条件の設定)
  - [2.3 境界積分方程式による解析](#23-境界積分方程式による解析)
- [3. BEMによる計算](#3-bemによる計算)
  - [3.1 基本解とその微分の導出](#31-基本解とその微分の導出)
  - [3.2 境界積分方程式の離散化](#32-境界積分方程式の離散化)
  - [3.3 内部点での速度ベクトルの計算](#33-内部点での速度ベクトルの計算)
- [参考文献](#参考文献)
- [Appendix](#appendix)
  - [A. 計算](#a-計算)
  - [A-Eq2.6](#a-eq26)
    - [A-Eq3.9](#a-eq39)
    - [A-Eq3.14](#a-eq314)
  - [B. 理論的な背景](#b-理論的な背景)
    - [B.1 ポテンシャル流れにおけるナビエ–ストークス方程式の簡略化](#b1-ポテンシャル流れにおけるナビエストークス方程式の簡略化)

## 1. シミュレーションの概要

### 1.1 本シミュレーションの目的と解析フロー

　本シミュレーションの目的は、障害物が存在する複雑な環境において、始点から終点へ至る流れ場を計算することである。具体的には、迷路のような環境における入口から出口への流れを、[ポテンシャル流れ](#12-ポテンシャル流れとは)と[境界要素法（BEM: Boundary Element Method）](#13-bemとは)を用いて解析する。

　解析は以下の手順で行う。

1. **形状定義**：解析領域の境界（壁、入口、出口）を閉曲線として定義し、線分要素に分割する。
2. **境界条件の設定**：入口と出口にはポテンシャル値を与え、壁面には法線流速ゼロ（壁を突き抜けない）という条件を設定する。
3. **境界値の決定**：BEMを用いて境界積分方程式を解き、境界上の未知量（ポテンシャルまたは流速）を決定する。
4. **内部流速の計算**：得られた境界値を用いて、領域内部の任意点における速度ベクトルを計算する。

> **背景メモ**：本シミュレーションの最終的な目標は、3次元のポリゴン領域において、z座標が最大となる面へ向かう流れ場を計算することである。今回目的とする領域では、z方向に進めない行き詰まり（局所的な頂点）が多数存在する（可能性がある）。単純にz方向への流れを計算すると、これらの行き詰まりに捕捉される問題があるため、ポテンシャル流れの概念を利用してz方向に連続な経路を見つけることを目指している。ただし、いきなり3次元の解析を行うのは困難であるため、まずは2次元での実装と検証を行う。

### 1.2 ポテンシャル流れとは

　**ポテンシャル流れ**とは、流体の速度場 $\boldsymbol{v}$ がスカラー関数である速度ポテンシャル $\phi$ の勾配として表される流れである。すなわち、

$$\boldsymbol{v} = \nabla \phi \qquad\text{(1.1)}$$

を満たす流れを指す。

　ポテンシャル流れが成立するためには、以下の条件を満たす必要がある。

1. **非回転性** (渦なし)：流体の回転成分である渦度 $\boldsymbol{\omega} = \nabla \times \boldsymbol{v}$ がゼロであること。渦なし流れであれば、必ず速度ポテンシャルが存在し、ポテンシャル流れとなる。
2. **非圧縮性**：流体の密度 $\rho$ が一定であること。連続の式 $\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \boldsymbol{v}) = 0$ により、非圧縮性の場合は $\nabla \cdot \boldsymbol{v} = 0$ が成り立つ。
3. **非粘性**：流体の粘性が無視できること。

　また、本シミュレーションでは、さらに以下の条件を仮定する。

4. **定常性**：時間に依存しない流れであること ($\frac{\partial \boldsymbol{v}}{\partial t} = 0$).
5. **外力の不存在**：重力などの外力が存在しないこと。

> **注**：条件3の非粘性は、式 $(1.1)$ を満たすための必要条件ではない。しかし、流体が存在する領域全体で非回転性を保つためには必要となる。詳細については[補遺B.1](#b1-ポテンシャル流れにおけるナビエストークス方程式の簡略化)を参照のこと。

　ポテンシャル流れは、計算コストが低く理論的な解析が容易であるという利点がある。一方で、粘性の影響が支配的な壁面近傍の流れや、渦が発生する流れなど、多くの実際の流体現象を正確に表現することはできない。それでも、航空力学や水理学などの分野では、ポテンシャル流れを用いた近似解析が広く利用されている。特に、次節で説明する境界要素法（BEM）を用いた数値解析が盛んに行われている。

### 1.3 BEMとは

　**境界要素法**（BEM: Boundary Element Method）は、境界上の情報のみを用いて領域内部の物理量を計算する数値解析手法である。有限要素法（FEM）や有限差分法（FDM）といった他の主要な数値解析手法が領域全体を離散化（メッシュ化）するのに対し、BEMは**境界のみ**を離散化する点に特徴がある。

　BEMの主な利点は以下の通りである。

1. **計算コストの削減**：解析対象の次元を1つ減らすことができる。例えば、2次元領域の解析では境界が1次元となるため、計算コストが大幅に削減される。
2. **無限領域への適用**：領域が無限に広がる問題であっても、境界上の情報のみを扱うため、容易に解析できる。これは、無限遠まで広がる流れ場の解析に特に適している。

　一方、BEMには欠点もある。離散化によって得られる連立方程式の係数行列が密行列かつ非対称となるため、要素数が増加すると計算コストが急激に増大する。この問題は、高速多重極法（FMM: Fast Multipole Method）などの手法を組み合わせることである程度緩和できる。

## 2. 境界積分方程式の定式化

### 2.1 支配方程式

　前節で述べたBEMを用いて、領域内部の任意の点における速度ベクトルを計算する方法について説明する。まず、ポテンシャル流れの支配方程式を確認する。

　点 $\boldsymbol{x}$ における速度ベクトル $\boldsymbol{v} = (u, v)$ は、流れの非回転性により、式 $(1.1)$ を満たすスカラー関数 $\phi(\boldsymbol{x})$（速度ポテンシャル）を用いて表される。さらに、非圧縮性の条件 $\nabla \cdot \boldsymbol{v} = 0$ より, $\phi$ はラプラス方程式

$$\nabla^2 \phi = \frac{\partial^2 \phi}{\partial x^2} + \frac{\partial^2 \phi}{\partial y^2} = 0 \qquad\text{(2.1)}$$

を満たす。したがって、ポテンシャル流れの解析は、与えられた領域内でラプラス方程式 $(2.1)$ と境界条件を満たす速度ポテンシャル $\phi$ を求める問題に帰着される。

### 2.2 境界条件の設定

　閉曲線境界 $\Gamma$ 上での条件を設定する。閉曲線によって区切られた、流体が存在する側を内部、その補集合を外部として、境界上の外向き法線ベクトルを $\boldsymbol{n}$ とする。境界の一部を**開口部**とし、そこでは流体の出入りがあるものとする。その他の部分は**壁面**であり、流体は壁を貫通しない。

**開口部**：流体の出入りがある辺では、法線方向の速度成分を既知とする。

$$\frac{\partial \phi}{\partial n} = \boldsymbol{v} \cdot \boldsymbol{n} = V_{\text{in/out}} \qquad\text{(2.2)}$$

ここで, $V_{\text{in/out}}$ は既知の値である（流入なら負、流出なら正など、法線の向きの定義に依存する）。

**壁面**：壁面では流体が壁を突き抜けないため、法線方向の速度成分はゼロである。

$$\frac{\partial \phi}{\partial n} = \boldsymbol{v} \cdot \boldsymbol{n} = 0 \qquad\text{(2.3)}$$

　非圧縮性流体では、質量保存則により、閉じた領域への総流入量と総流出量は一致しなければならない。すなわち、

$$\int_{\Gamma} \frac{\partial \phi}{\partial n} ds = 0 \qquad\text{(2.4)}$$

が成り立つ。

> **注**：本プログラムでは、開口部では速度ポテンシャル $\phi$ の値を指定し（ディリクレ条件）、壁面では法線速度成分 $\partial \phi / \partial n$ を指定する（ノイマン条件）形で境界条件を設定する。

### 2.3 境界積分方程式による解析

　任意形状の境界を持つ問題において、ラプラス方程式の解を内部点 $\boldsymbol{x}$ で求めるため、グリーンの第2恒等式を利用する（参考文献[1]ではGreen-Gaussの定理から導出）。観測点を $\boldsymbol{x}$、境界上の積分点を $\boldsymbol{\xi}$ とし、相対位置ベクトルを $\boldsymbol{r} = \boldsymbol{\xi} - \boldsymbol{x}$、その大きさを $r = \|\boldsymbol{r}\|$ と表す。2次元ラプラス方程式の基本解（点源によるポテンシャル）は次式で与えられる。

$$G(\boldsymbol{x}, \boldsymbol{\xi}) = -\frac{1}{2\pi} \ln r \qquad\text{(2.5)}$$

　文献[1]の第2章より、グリーンの第2恒等式を適用すると、領域内部の点 $\boldsymbol{x}$ における速度ポテンシャルは、境界 $\Gamma$ 上の積分で次のように表される。

$$c(\boldsymbol{x}) \phi(\boldsymbol{x}) = \oint_{\Gamma} \left(G(\boldsymbol{x}, \boldsymbol{\xi}) \frac{\partial \phi(\boldsymbol{\xi})}{\partial n} - \phi(\boldsymbol{\xi}) \frac{\partial G(\boldsymbol{x}, \boldsymbol{\xi})}{\partial n} \right) ds(\boldsymbol{\xi}) \qquad\text{(2.6)}$$

ここで, $\phi(\boldsymbol{\xi})$ は境界上の点 $\boldsymbol{\xi}$ における速度ポテンシャル, $\frac{\partial \phi(\boldsymbol{\xi})}{\partial n}$ はその法線方向微分（すなわち法線方向流速）である。また, $c(\boldsymbol{x})$ は観測点 $\boldsymbol{x}$ の位置に依存する幾何学的係数であり、式 $(2.7)$ のように定義される<sup>[1, Eq. (29)]</sup>。本シミュレーションでは、境界を線分要素で近似しているため, $\boldsymbol{x} \in \Gamma$ の場合は $c(\boldsymbol{x}) = 1/2$ とする。

$$c(\boldsymbol{x}) = \begin{cases}
    1 & (\boldsymbol{x} \text{ が領域内部にある場合}) \\\
    \frac{1}{2} & (\boldsymbol{x} \in \Gamma \text{ かつ境界が滑らかな場合}) \\\
    1 - \frac{\theta}{2\pi} & (\boldsymbol{x} \in \Gamma \text{ かつ境界が角を持つ場合}) \\\
    0 & (\boldsymbol{x} \text{ が領域外にある場合})
\end{cases} \qquad\text{(2.7)}$$

　境界条件から, $\frac{\partial \phi(\boldsymbol{\xi})}{\partial n}$ は開口部では既知の値 $V_{\text{in/out}}(\boldsymbol{\xi})$、壁面ではゼロである。一方, $\phi(\boldsymbol{\xi})$ は未知量である。したがって、式 $(2.7)$ は境界上の未知量 $\phi(\boldsymbol{\xi})$ に関する積分方程式となる。この積分方程式を離散化して連立方程式を解くことで、境界上の全ての点における速度ポテンシャル $\phi(\boldsymbol{\xi})$ を求めることができる。

**内部点での速度ベクトルの計算**

　境界上の $\phi$ と $\frac{\partial \phi}{\partial n}$ が求まれば、内部点 $\boldsymbol{x}$ における速度ベクトル $\boldsymbol{v}(\boldsymbol{x}) = \nabla \phi(\boldsymbol{x})$ を計算できる。内部点 $\boldsymbol{x}$ では $c(\boldsymbol{x}) = 1$ であるから、式 $(2.7)$ の両辺を $\boldsymbol{x}$ で微分すると、

$$\begin{aligned}
    \boldsymbol{v}(\boldsymbol{x}) &= \nabla_{\boldsymbol{x}} \phi(\boldsymbol{x}) \\\
    &= \oint_{\Gamma} \left( \nabla_{\boldsymbol{x}} G(\boldsymbol{x}, \boldsymbol{\xi}) \frac{\partial \phi(\boldsymbol{\xi})}{\partial n} - \phi(\boldsymbol{\xi}) \nabla_{\boldsymbol{x}} \frac{\partial G(\boldsymbol{x}, \boldsymbol{\xi})}{\partial n} \right) ds(\boldsymbol{\xi}) \\\
    &= \sum_{j=1}^{N} \int_{\Gamma_j} \left( \nabla_{\boldsymbol{x}} G(\boldsymbol{x}, \boldsymbol{\xi}) \frac{\partial \phi(\boldsymbol{\xi})}{\partial n} - \phi(\boldsymbol{\xi}) \nabla_{\boldsymbol{x}} \frac{\partial G(\boldsymbol{x}, \boldsymbol{\xi})}{\partial n} \right) ds(\boldsymbol{\xi})
\end{aligned} \qquad\text{(2.8)}$$

が得られる。最後の式は境界を $N$ 個の要素 $\Gamma_j$ に分割した場合の離散化表現である。

　以降では、式 $(2.8)$ に現れる4つの量、すなわち $\nabla_{\boldsymbol{x}} G(\boldsymbol{x}, \boldsymbol{\xi})$, $\nabla_{\boldsymbol{x}} \frac{\partial G(\boldsymbol{x}, \boldsymbol{\xi})}{\partial n}$, $\phi(\boldsymbol{\xi})$, $\frac{\partial \phi(\boldsymbol{\xi})}{\partial n}$ の計算方法を説明する。前者2つ（基本解の微分）は[節3.1](#31-基本解とその微分の導出)で数学的に導出し、後者2つ（境界上の物理量）は[節3.2](#32-境界積分方程式の離散化)でBEMによる数値的な求め方を説明する。

<!-- 疑問点: 上の計算は3次元においても成り立つか. 少なくとも、式(2.5)は異なるようである -->

## 3. BEMによる計算

### 3.1 基本解とその微分の導出

　前述の通り、2次元ラプラス方程式の基本解は式 $(2.5)$ で与えられる。式 $(2.8)$ から、速度ベクトル $\boldsymbol{v}(\boldsymbol{x})$ を計算するためには、基本解の $\boldsymbol{x}$ に関する勾配 $\nabla_{\boldsymbol{x}} G(\boldsymbol{x}, \boldsymbol{\xi})$ と、基本解の法線方向微分の $\boldsymbol{x}$ に関する勾配 $\nabla_{\boldsymbol{x}} \frac{\partial G(\boldsymbol{x}, \boldsymbol{\xi})}{\partial n}$ を求める必要がある。

　まず、基本解の勾配を求める. $\boldsymbol{\xi} = (\xi_x, \xi_y)$ とすれば $\boldsymbol{r} = (\xi_x - x, \xi_y - y)$, $r = \sqrt{(\xi_x - x)^2 + (\xi_y - y)^2}$ となるから, $r$ の $\boldsymbol{x}$ に関する勾配は以下のようになる。

$$\begin{aligned}
    \nabla_{\boldsymbol{x}} r &= \left( \frac{\partial r}{\partial x}, \frac{\partial r}{\partial y} \right) \\\
    &= \left( \frac{-(\xi_x - x)}{r}, \frac{-(\xi_y - y)}{r} \right) \\\
    &= -\frac{\boldsymbol{r}}{r}
\end{aligned} \qquad\text{(3.1)}$$

　したがって, $\nabla_{\boldsymbol{x}} G(\boldsymbol{x}, \boldsymbol{\xi})$ は以下のように求められる。

$$\begin{aligned}
    \nabla_{\boldsymbol{x}} G(\boldsymbol{x}, \boldsymbol{\xi}) &= -\frac{1}{2\pi} \nabla_{\boldsymbol{x}} \ln r \\\
    &= -\frac{1}{2\pi} \cdot \frac{1}{r} \nabla_{\boldsymbol{x}} r \\\
    &= -\frac{1}{2\pi} \cdot \frac{1}{r} \left( -\frac{\boldsymbol{r}}{r} \right) \\\
    &= \frac{1}{2\pi r^2} \boldsymbol{r}
\end{aligned} \qquad\text{(3.2)}$$

　次に、基本解の法線方向微分の $\boldsymbol{x}$ に関する勾配を求める。まず, $r$ の $\boldsymbol{\xi}$ に関する勾配は以下のようになる。

$$\begin{aligned}
    \nabla_{\boldsymbol{\xi}} r &= \left( \frac{\partial r}{\partial \xi_x}, \frac{\partial r}{\partial \xi_y} \right) \\\
    &= \left( \frac{\xi_x - x}{r}, \frac{\xi_y - y}{r} \right) \\\
    &= \frac{\boldsymbol{r}}{r}
\end{aligned} \qquad\text{(3.3)}$$

したがって、点 $\boldsymbol{\xi}$ における法線ベクトルを $\boldsymbol{n} = (n_x, n_y)$ としたときの法線微分 $\frac{\partial G}{\partial n}$ は、次のように表される。

$$\begin{aligned}
    \frac{\partial G(\boldsymbol{x}, \boldsymbol{\xi})}{\partial n} &= \nabla_{\boldsymbol{\xi}} G(\boldsymbol{x}, \boldsymbol{\xi}) \cdot \boldsymbol{n} \\\
    &= \left(-\frac{1}{2\pi} \cdot \frac{1}{r} \nabla_{\boldsymbol{\xi}} r\right) \cdot \boldsymbol{n} \\\
    &= \left(-\frac{1}{2\pi} \cdot \frac{1}{r} \cdot \frac{\boldsymbol{r}}{r}\right) \cdot \boldsymbol{n} \\\
    &= -\frac{\boldsymbol{r} \cdot \boldsymbol{n}}{2\pi r^2}
\end{aligned} \qquad\text{(3.4)}$$

　次に、このスカラー量を $\boldsymbol{x}$ で微分する。商の微分公式を用いると、

$$\nabla_{\boldsymbol{x}} \frac{\partial G(\boldsymbol{x}, \boldsymbol{\xi})}{\partial n} = -\frac{1}{2\pi} \frac{r^2 \nabla_{\boldsymbol{x}} (\boldsymbol{r} \cdot \boldsymbol{n}) - (\boldsymbol{r} \cdot \boldsymbol{n}) \nabla_{\boldsymbol{x}} (r^2)}{r^4} \qquad\text{(3.5)}$$
となる。ここで, $\nabla_{\boldsymbol{x}} (\boldsymbol{r} \cdot \boldsymbol{n}) = \nabla_{\boldsymbol{x}} ((\xi_x - x) n_x + (\xi_y - y) n_y) = -\boldsymbol{n}$ であり, $\nabla_{\boldsymbol{x}} (r^2) = \nabla_{\boldsymbol{x}} ((\xi_x - x)^2 + (\xi_y - y)^2) = -2\boldsymbol{r}$ である。これらを式 $(3.5)$ に代入すると、式 $(3.6)$ が得られる。

$$\begin{aligned}
    \nabla_{\boldsymbol{x}} \frac{\partial G(\boldsymbol{x}, \boldsymbol{\xi})}{\partial n} &= -\frac{1}{2\pi} \frac{r^2 (-\boldsymbol{n}) - (\boldsymbol{r} \cdot \boldsymbol{n}) (-2\boldsymbol{r})}{r^4} \\\
    &= \frac{1}{2\pi r^2} \left( \boldsymbol{n} - \frac{2 (\boldsymbol{r} \cdot \boldsymbol{n})}{r^2} \boldsymbol{r} \right)
\end{aligned} \qquad\text{(3.6)}$$

### 3.2 境界積分方程式の離散化

　本節では、境界積分方程式 $(2.6)$ を離散化し、連立一次方程式 $A\mathbf{x} = \mathbf{b}$ の形に変換する手順を説明する。本稿では、閉曲線 $\Gamma$ を $N$ 個の線分要素 $\Gamma_j \ (j = 1, 2, \ldots, N)$ に分割し、各要素 $\Gamma_j$ 上でポテンシャル $\phi$ と流速 $q = \partial \phi / \partial n$ を一定と仮定する。これにより、境界上の未知量を有限個の値で近似できる。

$$\phi(\boldsymbol{\xi}) \approx \phi_j, \quad q(\boldsymbol{\xi}) \approx q_j, \quad \left( \boldsymbol{\xi} \in \Gamma_j \right) \qquad\text{(3.7)}$$

　以降、各要素の中点を $\boldsymbol{\xi}_j$、各要素の長さを $l_j$ と表す。中点 $\boldsymbol{\xi}_j$ における境界積分方程式 $(2.6)$ は, $\phi$ の項を左側に, $q$ の項を右側に移項すると以下のように書ける。

$$\begin{aligned}
    &\quad c(\boldsymbol{\xi}_i) \phi(\boldsymbol{\xi}_i) + \oint_{\Gamma} \phi(\boldsymbol{\xi}) \frac{\partial G(\boldsymbol{\xi}_i, \boldsymbol{\xi})}{\partial n} ds(\boldsymbol{\xi}) = \oint_{\Gamma} G(\boldsymbol{\xi}_i, \boldsymbol{\xi}) q(\boldsymbol{\xi}) ds(\boldsymbol{\xi}) \\\
    &\Rightarrow c(\boldsymbol{\xi}_i) \phi_i + \sum_{j=1}^{N} \phi_j \int_{\Gamma_j} \frac{\partial G(\boldsymbol{\xi}_i, \boldsymbol{\xi})}{\partial n} ds(\boldsymbol{\xi}) = \sum_{j=1}^{N} q_j \int_{\Gamma_j} G(\boldsymbol{\xi}_i, \boldsymbol{\xi}) ds(\boldsymbol{\xi})
    &\Rightarrow \sum_{j=1}^{N} H_{ij} \phi_j = \sum_{j=1}^{N} G_{ij} q_j
\end{aligned} \qquad\text{(3.8)}$$

ここで、行列要素 $H_{ij}$ と $G_{ij}$ は、クロネッカーのデルタ $\delta_{ij}$ を用いて以下のように定義される。この計算の詳細は、[付録A](#a-eq39)を参照されたい。

$$\begin{aligned}
    H_{ij} &= \int_{\Gamma_j} \frac{\partial G(\boldsymbol{\xi}_i, \boldsymbol{\xi})}{\partial n} ds(\boldsymbol{\xi}) + \delta_{ij} c(\boldsymbol{\xi}_i) \\\
    &= \int_{\Gamma_j} \left( -\frac{\boldsymbol{r} \cdot \boldsymbol{n}}{2\pi r^2} \right) ds(\boldsymbol{\xi}) + \delta_{ij} c(\boldsymbol{\xi}_i) \\\
    G_{ij} &= \int_{\Gamma_j} G(\boldsymbol{\xi}_i, \boldsymbol{\xi}) ds(\boldsymbol{\xi}) \\\
    &= \int_{\Gamma_j} \left( -\frac{1}{2\pi} \ln r \right) ds(\boldsymbol{\xi})
\end{aligned} \qquad\text{(3.9)}$$

> ただし、前述のように線分として境界を近似しているため、境界上の点 $\boldsymbol{\xi} _j$ について位置ベクトル $\boldsymbol{r}$ と法線ベクトル $\boldsymbol{n}$ は直交する。したがって, $\boldsymbol{r} \cdot \boldsymbol{n} = 0$ となるから, $H _{ii} = c(\boldsymbol{\xi}_i) = 1/2$ となる。

　式 $(3.8)$ は, $i = 1, 2, \ldots, N$ について連立させることで、以下の行列方程式として表される。

$$\begin{bmatrix}
    H_{11} & H_{12} & \cdots & H_{1N} \\\
    H_{21} & H_{22} & \cdots & H_{2N} \\\
    \vdots & \vdots & \ddots & \vdots \\\
    H_{N1} & H_{N2} & \cdots & H_{NN}
  \end{bmatrix}
  \begin{bmatrix}
    \phi_1 \\\ \phi_2 \\\ \vdots \\\ \phi_N
  \end{bmatrix} =
  \begin{bmatrix}
    G_{11} & G_{12} & \cdots & G_{1N} \\\
    G_{21} & G_{22} & \cdots & G_{2N} \\\
    \vdots & \vdots & \ddots & \vdots \\\
    G_{N1} & G_{N2} & \cdots & G_{NN}
  \end{bmatrix}
  \begin{bmatrix}
    q_1 \\\ q_2 \\\ \vdots \\\ q_N
\end{bmatrix} \qquad\text{(3.10)}$$

　[節2.2](#22-境界条件の設定)で説明したように、各要素 $\Gamma_j$ において, $\phi_j$ または $q_j$ のどちらかが指定されるため、式 $(3.10)$ は未知量に関する連立一次方程式として解くことができる。このとき、既知の値を右辺にベクトル $\mathbf{b}$ としてまとめ、未知量を左辺にベクトル $\mathbf{x}$ としてまとめることで、以下の形に変換できる。

$$A \mathbf{x} = \mathbf{b} \qquad\text{(3.11)}$$

**具体例**：境界 $\Gamma$ が4つの要素 $\Gamma_1, \Gamma_2, \Gamma_3, \Gamma_4$ に分割されており, $\Gamma_1$ と $\Gamma_3$ でポテンシャル $\phi$ が指定され, $\Gamma_2$ と $\Gamma_4$ で流速 $q$ が指定されている場合を考える。このとき、式 $(3.10)$ から、以下のように連立方程式を構築できる。ただし、区別のため、既知の変数には添え字 "k" を付ける。また, $H_{ij}$ と $G_{ij}$ については、いずれも既に計算されているものとする。

$$\begin{bmatrix}
  H_{11} & H_{12} & H_{13} & H_{14} \\\
  H_{21} & H_{22} & H_{23} & H_{24} \\\
  H_{31} & H_{32} & H_{33} & H_{34} \\\
  H_{41} & H_{42} & H_{43} & H_{44}
\end{bmatrix}
\begin{bmatrix}
  \phi_{1,k} \\\ \phi_2 \\\ \phi_{3,k} \\\ \phi_4
\end{bmatrix} =
\begin{bmatrix}
  G_{11} & G_{12} & G_{13} & G_{14} \\\
  G_{21} & G_{22} & G_{23} & G_{24} \\\
  G_{31} & G_{32} & G_{33} & G_{34} \\\
  G_{41} & G_{42} & G_{43} & G_{44}
\end{bmatrix}
\begin{bmatrix}
  q_1 \\\ q_{2,k} \\\ q_3 \\\ q_{4,k}
\end{bmatrix}$$

$$\Rightarrow \begin{bmatrix}
    -G_{11} & H_{12} & -G_{13} & H_{14} \\\
    -G_{21} & H_{22} & -G_{23} & H_{24} \\\
    -G_{31} & H_{32} & -G_{33} & H_{34} \\\
    -G_{41} & H_{42} & -G_{43} & H_{44}
\end{bmatrix} \begin{bmatrix}
  q_1 \\\ \phi_2 \\\ q_3 \\\ \phi_4
\end{bmatrix} = \begin{bmatrix}
    - H_{11} \phi_{1,k} + G_{12} q_{2,k} - H_{13} \phi_{3,k} + G_{14} q_{4,k} \\\
    - H_{21} \phi_{1,k} + G_{22} q_{2,k} - H_{23} \phi_{3,k} + G_{24} q_{4,k} \\\
    - H_{31} \phi_{1,k} + G_{32} q_{2,k} - H_{33} \phi_{3,k} + G_{34} q_{4,k} \\\
    - H_{41} \phi_{1,k} + G_{42} q_{2,k} - H_{43} \phi_{3,k} + G_{44} q_{4,k}
\end{bmatrix} \qquad\text{(3.12)}$$

### 3.3 内部点での速度ベクトルの計算

　境界上の未知量 $\phi_j$ と $q_j$ が求まれば、式 $(2.8)$ を用いて、任意の内部点 $\boldsymbol{x}$ における速度ベクトル $\boldsymbol{v}(\boldsymbol{x})$ を計算できる。具体的には、式 $(3.2)$ と式 $(3.6)$ を式 $(2.8)$ に代入することで、以下のように表される。

$$\begin{aligned}
    \boldsymbol{v}(\boldsymbol{x}) &= \sum_{j=1}^{N} \int_{\Gamma_j} \left( \nabla_{\boldsymbol{x}} G(\boldsymbol{x}, \boldsymbol{\xi}) \frac{\partial \phi(\boldsymbol{\xi})}{\partial n} - \phi(\boldsymbol{\xi}) \nabla_{\boldsymbol{x}} \frac{\partial G(\boldsymbol{x}, \boldsymbol{\xi})}{\partial n} \right) ds(\boldsymbol{\xi}) \\\
    &= \sum_{j=1}^{N} \int_{\Gamma_j} \left( \frac{1}{2\pi r^2} \boldsymbol{r} \cdot q_j - \phi_j \cdot \frac{1}{2\pi r^2} \left( \boldsymbol{n} - \frac{2 (\boldsymbol{r} \cdot \boldsymbol{n})}{r^2} \boldsymbol{r} \right) \right) ds(\boldsymbol{\xi}) \\\
    &= \frac{1}{2\pi} \sum_{j=1}^{N} \left( q_j \int_{\Gamma_j} \frac{\boldsymbol{r}}{r^2} ds(\boldsymbol{\xi}) - \phi_j \int_{\Gamma_j} \frac{1}{r^2} \left( \boldsymbol{n} - \frac{2 (\boldsymbol{r} \cdot \boldsymbol{n})}{r^2} \boldsymbol{r} \right) ds(\boldsymbol{\xi}) \right)
\end{aligned} \qquad\text{(3.13)}$$

　したがって、速度ベクトル $\boldsymbol{v}(\boldsymbol{x})$ は、上の式に現れる2つの積分を計算することで求められる。線分要素 $\Gamma_j$ の単位接線ベクトルを $\boldsymbol{t}$, 単位法線ベクトルを $\boldsymbol{n}$, 要素の両端点へのベクトルを $\boldsymbol{r}_1, \boldsymbol{r}_2$ と表す。これらを用いて、上の2つの積分は以下のように計算できる。導出の詳細は、[付録A](#a-eq314)を参照されたい。

$$\boldsymbol{v}(\boldsymbol{x}) = \frac{1}{2\pi} \sum_{j=1}^{N} \left( \left( q_j \ln \frac{r_2}{r_1} - \phi_j \left( \frac{d}{r_2^2} - \frac{d}{r_1^2} \right) \right) \boldsymbol{t} + \left( q_j (\theta_2 - \theta_1) - \phi_j \left( -\frac{s_2}{r_2^2} + \frac{s_1}{r_1^2} \right) \right) \boldsymbol{n} \right) \qquad\text{(3.14)}$$

ここで, $s_1 = \boldsymbol{r}_1 \cdot \boldsymbol{t}$, $s_2 = \boldsymbol{r}_2 \cdot \boldsymbol{t}$, および $d = \boldsymbol{r}_1 \cdot \boldsymbol{n} = \boldsymbol{r}_2 \cdot \boldsymbol{n}$ である。

## 参考文献

[1] Tara LaForce, "PE281 Boundary Element Method Course Notes", Stanford University, 2006, https://web.stanford.edu/class/energy281/BoundaryElementMethod.pdf

## Appendix

### A. 計算

　本節では、本文で述べた一部の式について、これを整理・導出する過程を示す。

### A-Eq2.6

　式 $(2.6)$ を展開する。式 $(2.5)$, 式 $(3.4)$ を代入すると、以下のようになる。

$$\begin{aligned}
    c(\boldsymbol{x}) \phi(\boldsymbol{x}) &= \oint_{\Gamma} \left(G(\boldsymbol{x}, \boldsymbol{\xi}) \frac{\partial \phi(\boldsymbol{\xi})}{\partial n} - \phi(\boldsymbol{\xi}) \frac{\partial G(\boldsymbol{x}, \boldsymbol{\xi})}{\partial n} \right) ds(\boldsymbol{\xi}) \\\
    &= \sum_{j=1}^{N} \int_{\Gamma_j} \left( G(\boldsymbol{x}, \boldsymbol{\xi}) q_j - \phi_j \frac{\partial G(\boldsymbol{x}, \boldsymbol{\xi})}{\partial n} \right) ds(\boldsymbol{\xi}) \\\
    &= \sum_{j=1}^{N} \int_{\Gamma_j} \left( -\frac{q_j}{2\pi} \ln r + \phi_j \frac{\boldsymbol{r} \cdot \boldsymbol{n}}{2\pi r^2} \right) ds(\boldsymbol{\xi})
\end{aligned} \qquad\text{(A26-1)}$$

　この計算のため、相対ベクトル $\boldsymbol{r}$ を分解し、局所座標系 $s$ を導入する。線分要素 $\Gamma_j$ の単位接線ベクトルを $\boldsymbol{t}$、単位法線ベクトルを $\boldsymbol{n}$ とすれば, $\boldsymbol{r}$ は以下のように分解できる。

$$\begin{aligned}
    &\boldsymbol{r} = s \boldsymbol{t} + d \boldsymbol{n} \\\
    & \text{where} \quad \left\lbrace \begin{aligned}
        s &= \boldsymbol{r} \cdot \boldsymbol{t} \\\
        d &= \boldsymbol{r} \cdot \boldsymbol{n}
    \end{aligned} \right.
\end{aligned} \qquad\text{(A26-2)}$$

ここで, $d$ は観測点 $\boldsymbol{x}$ から線分要素 $\Gamma_j$ までの符号付最短距離であり, $s$ は線分要素上の局所座標である。したがって, $r = \sqrt{s^2 + d^2}$ となる。以降、積分変数を $ds(\boldsymbol{\xi}) = ds$ とし、線分の始点を $s = s_1$、終点を $s = s_2$ と表す。

　これらを式 $\text{(A26-1)}$ に代入すると、要素 $\Gamma_j$ に関する積分は以下のように書き換えられる。

$$\begin{aligned}
    I_j &= \int_{\Gamma_j} \left( -\frac{q_j}{2\pi} \ln r + \phi_j \frac{\boldsymbol{r} \cdot \boldsymbol{n}}{2\pi r^2} \right) ds(\boldsymbol{\xi}) \\\
    &= -\frac{q_j}{2\pi} \int_{s_1}^{s_2} \ln \sqrt{s^2 + d^2} ds + \frac{\phi_j}{2\pi} \int_{s_1}^{s_2} \frac{d}{s^2 + d^2} ds \\\
\end{aligned} \qquad\text{(A26-3)}$$

　これらの積分は、以下のように計算できる。

$$\begin{aligned}
    \int_{s_1}^{s_2} \ln \sqrt{s^2 + d^2} ds &= \frac{1}{2} \int_{s_1}^{s_2} \ln (s^2 + d^2) ds \\\
    &= \frac{1}{2} \left[ s \ln (s^2 + d^2) - 2s + 2d \arctan \left( \frac{s}{d} \right) \right]_{s_1}^{s_2} \\\
    &= \frac{s_2}{2} \ln r_2^2 - \frac{s_1}{2} \ln r_1^2 - (s_2 - s_1) + d \left( \arctan \frac{s_2}{d} - \arctan \frac{s_1}{d} \right) \\\
    \int_{s_1}^{s_2} \frac{d}{s^2 + d^2} ds &= \left[ \arctan \left( \frac{s}{d} \right) \right]_{s_1}^{s_2} \\\
    &= \arctan \frac{s_2}{d} - \arctan \frac{s_1}{d}
\end{aligned}$$

　以降, $\Delta\theta_j = \arctan \frac{s_2}{d} - \arctan \frac{s_1}{d}$ と表す。これらを式 $\text{(A26-3)}$ に代入すると、以下のようになる。この $I_j$ を全要素について合計することで、式 $(2.6)$ の右辺が得られる。

$$\begin{aligned}
    &I_j = -\frac{q_j}{2\pi} \left( \frac{s_2}{2} \ln r_2^2 - \frac{s_1}{2} \ln r_1^2 - (s_2 - s_1) \right) + \frac{\phi_j - q_j d}{2\pi} \Delta\theta_j \\\
\end{aligned} \qquad\text{(A26-4)}$$

> 実際のプログラムでは、`arctan2`関数を用いて $\Delta \theta_j$ を計算している。

#### A-Eq3.9

　式 $(3.9)$ の積分を計算する。

$$\begin{aligned}
    H_{ij} &= \int_{\Gamma_j} \left( -\frac{\boldsymbol{r} \cdot \boldsymbol{n}}{2\pi r^2} \right) ds(\boldsymbol{\xi}) + \delta_{ij} c(\boldsymbol{\xi}_i) \\\
    G_{ij} &= \int_{\Gamma_j} \left( -\frac{1}{2\pi} \ln r \right) ds(\boldsymbol{\xi})
\end{aligned} \qquad\text{(3.9)}$$

$H_{ij}$ **の計算：**

　まず、対角項 ($i = j$) について考える。このとき、前述の通り $\boldsymbol{r} \cdot \boldsymbol{n} = 0$ となるため、

$$H_{ii} = c(\boldsymbol{\xi}_i) = \frac{1}{2} \qquad\text{(A39-1)}$$

とできる。次に、非対角項 ($i \neq j$) について考える。このとき, $\boldsymbol{r}$ と $\boldsymbol{n}$ のなす角を $\theta$ とすると、以下のように変形できる。

$$\begin{aligned}
    H_{ij} &= \int_{\Gamma_j} \left( -\frac{\boldsymbol{r} \cdot \boldsymbol{n}}{2\pi r^2} \right) ds(\boldsymbol{\xi}) \\\
    &= -\frac{1}{2\pi} \int_{\Gamma_j} \frac{r \cos \theta}{r^2} ds(\boldsymbol{\xi}) \\\
    &= -\frac{1}{2\pi} \int_{\Gamma_j} \frac{\cos \theta}{r} ds(\boldsymbol{\xi})
\end{aligned}$$

ここで、観測点 $i$ から見て、微小要素 $ds$ が張る角度を $d\alpha$ とすると、演習の接線方向に対する要素 $ds$ の射影成分は $ds \cos \theta$ であり、これが弧長 $r d\alpha$ に等しいことから, $ds \cos \theta = r d\alpha$ が成り立つ。したがって、上の式は以下のように変形できる。

$$\begin{aligned}
    H_{ij} &= -\frac{1}{2\pi} \int_{\Gamma_j} \frac{\cos \theta}{r} ds(\boldsymbol{\xi}) \\\
    &= -\frac{1}{2\pi} \int_{\alpha_1}^{\alpha_2} d\alpha \\\
    &= -\frac{1}{2\pi} (\alpha_2 - \alpha_1)
\end{aligned} \qquad\text{(A39-2)}$$

すなわち、線分要素 $\Gamma_j$ の両端点 $\boldsymbol{\xi} _{j1}, \boldsymbol{\xi} _{j2}$ が観測点 $\boldsymbol{\xi} _i$ となす角度 $(\alpha_2 - \alpha_1)$ を用いて, $H _{ij}$ を表すことができる。

> プログラムでは、要素の両端点の座標を $(x_{j1}, y_{j1}), (x_{j2}, y_{j2})$、観測点の座標を $(x_i, y_i)$ として、以下のように計算できる。この差分 $(\alpha_2 - \alpha_1)$ を $[-\pi, \pi]$ の範囲に正規化し、上の式に代入することで $H_{ij}$ を求める。
>
> $$\begin{aligned}\alpha_1 &= \text{atan2}(y_{j1} - y_i, x_{j1} - x_i) \\\ \alpha_2 &= \text{atan2}(y_{j2} - y_i, x_{j2} - x_i)\end{aligned}$$

$G_{ij}$ **の計算：**

　まず、対角項 ($i = j$) について考える。線分要素 $\Gamma_j$ の長さを $l_j$ とすると、中点 $\boldsymbol{\xi}_j$ を原点に取った局所座標 $s \in \left[-\frac{l_j}{2}, \frac{l_j}{2}\right]$ を導入できる。このとき $r = |s|$ となるから、以下のように計算できる。

$$\begin{aligned}
    G_{ii} &= -\frac{1}{2\pi} \int_{-\frac{l_j}{2}}^{\frac{l_j}{2}} \ln |s| ds \\\
    &= -\frac{1}{\pi} \int_{0}^{\frac{l_j}{2}} \ln s \ ds \\\
    &= -\frac{1}{\pi} \left[ s \ln s - s \right]_{0}^{\frac{l_j}{2}} \\\
    &= -\frac{l_j}{2\pi} \left( 1 - \ln\frac{l_j}{2} \right)
\end{aligned} \qquad\text{(A39-3)}$$

　次に、非対角項 ($i \neq j$) について考える。観測点 $\boldsymbol{\xi}_i$ から線分 $\Gamma_j$ までの最短距離を $d$ とすると, $r = \sqrt{t^2 + d^2}$ となる。両端点における局所座標 $t$ の値をそれぞれ $t_1, t_2$ とすると、以下のように計算できる。

$$\begin{aligned}
    G_{ij} &= -\frac{1}{2\pi} \int_{t_1}^{t_2} \ln \sqrt{t^2 + d^2} \ dt \\\
    &= -\frac{1}{4\pi} \int_{t_1}^{t_2} \ln (t^2 + d^2) \ dt
\end{aligned}$$

　ここで、積分部分を $I$ と置いて部分積分を行うと、

$$\begin{aligned}
    I &= \int \ln (t^2 + d^2) \ dt \\\
    &= [t \ln (t^2 + d^2)] - \int t \cdot \frac{2t}{t^2 + d^2} \ dt \\\
    &= t\ln(t^2 + d^2) - 2 \int \frac{t^2}{t^2 + d^2} \ dt \\\
    &= t\ln(t^2 + d^2) - 2 \left( t - d^2 \int \frac{1}{t^2 + d^2} \ dt \right) \\\
    &= t\ln(t^2 + d^2) - 2t + 2d\ \text{arctan} \frac{t}{d} \\\
    &= 2\left(t \ln r - t + d\ \text{arctan} \frac{t}{d} \right)
\end{aligned}$$

　したがって, $G_{ij}$ は以下のように求められる。

$$\begin{aligned}
    G_{ij} &= -\frac{1}{4\pi} [I]_{t_1}^{t_2} \\\
    &= -\frac{1}{2\pi} \left[  t \ln r - t + d\ \text{arctan} \frac{t}{d} \right]_{t_1}^{t_2} \\\
    &= \frac{1}{2\pi} \left( t_1 \ln r_1 - t_2 \ln r_2 + (t_2 - t_1) - d(\theta_2 - \theta_1) \right)
\end{aligned} \qquad\text{(A39-4)}$$

#### A-Eq3.14

　式 $(3.13)$ の積分を計算する。

$$\boldsymbol{v}(\boldsymbol{x}) = \frac{1}{2\pi} \sum_{j=1}^{N} \left( q_j \int_{\Gamma_j} \frac{\boldsymbol{r}}{r^2} ds(\boldsymbol{\xi}) - \phi_j \int_{\Gamma_j} \frac{1}{r^2} \left( \boldsymbol{n} - \frac{2 (\boldsymbol{r} \cdot \boldsymbol{n})}{r^2} \boldsymbol{r} \right) ds(\boldsymbol{\xi}) \right) \qquad\text{(3.13)}$$

　この計算のため、[A-Eq2.6](#a-eq26)と同様に、相対ベクトル $\boldsymbol{r}$ を分解し、局所座標系 $s$ を導入する。以降、積分変数を $ds(\boldsymbol{\xi}) = ds$ とし、線分の始点を $s = s_1$、終点を $s = s_2$ と表す。

(1) $\int_{\Gamma_j} \frac{\boldsymbol{r}}{r^2} ds(\boldsymbol{\xi})$ の計算：

$$\begin{aligned}
    \int_{\Gamma_j} \frac{\boldsymbol{r}}{r^2} ds &= \int_{s_1}^{s_2} \frac{s \boldsymbol{t} + d \boldsymbol{n}}{s^2 + d^2} ds \\\
    &= \boldsymbol{t} \int_{s_1}^{s_2} \frac{s}{s^2 + d^2} ds + \boldsymbol{n} \int_{s_1}^{s_2} \frac{d}{s^2 + d^2} ds \\\
    &= \boldsymbol{t} \left[ \ln(r) \right]_{s_1}^{s_2} + \boldsymbol{n} \left[ \text{arctan} \frac{s}{d} \right]_{s_1}^{s_2} \\\
    &= \left( \ln \frac{r_2}{r_1} \right) \boldsymbol{t} + \left( \theta_2 - \theta_1 \right) \boldsymbol{n}
\end{aligned} \qquad\text{(A313-2)}$$

(2) $\int_{\Gamma_j} \frac{1}{r^2} \left( \boldsymbol{n} - \frac{2 (\boldsymbol{r} \cdot \boldsymbol{n})}{r^2} \boldsymbol{r} \right) ds(\boldsymbol{\xi})$ の計算：

$$\begin{aligned}
    &\quad \int_{\Gamma_j} \frac{1}{r^2} \left( \boldsymbol{n} - \frac{2 (\boldsymbol{r} \cdot \boldsymbol{n})}{r^2} \boldsymbol{r} \right) ds \\\
    &= \int_{s_1}^{s_2} \left( \frac{\boldsymbol{n}}{s^2 + d^2} - \frac{2d(s \boldsymbol{t} + d^2 \boldsymbol{n})}{(s^2 + d^2)^2} \right) ds \\\
    &= \boldsymbol{t} \int_{s_1}^{s_2} \left( -\frac{2ds}{(s^2 + d^2)^2} \right) ds + \boldsymbol{n} \int_{s_1}^{s_2} \left( \frac{1}{s^2 + d^2} - \frac{2d^2}{(s^2 + d^2)^2} \right) ds \\\
    &= \left[ \frac{d}{r^2} \right]_{s_1}^{s_2} \boldsymbol{t} + \left[ -\frac{s}{r^2} \right]_{s_1}^{s_2} \boldsymbol{n} \\\
    &= \left( \frac{d}{r_2^2} - \frac{d}{r_1^2} \right) \boldsymbol{t} + \left( -\frac{s_2}{r_2^2} + \frac{s_1}{r_1^2} \right) \boldsymbol{n}
\end{aligned} \qquad\text{(A313-3)}$$

　以上より、式 $(3.13)$ は以下のように表される。

$$\begin{aligned}
    \boldsymbol{v}(\boldsymbol{x}) &= \frac{1}{2\pi} \sum_{j=1}^{N} \left( q_j \left( \left( \ln \frac{r_2}{r_1} \right) \boldsymbol{t} + \left( \theta_2 - \theta_1 \right) \boldsymbol{n} \right) - \phi_j \left( \left( \frac{d}{r_2^2} - \frac{d}{r_1^2} \right) \boldsymbol{t} + \left( -\frac{s_2}{r_2^2} + \frac{s_1}{r_1^2} \right) \boldsymbol{n} \right) \right) \\\
    &= \frac{1}{2\pi} \sum_{j=1}^{N} \left( \left( q_j \ln \frac{r_2}{r_1} - \phi_j \left( \frac{d}{r_2^2} - \frac{d}{r_1^2} \right) \right) \boldsymbol{t} + \left( q_j (\theta_2 - \theta_1) - \phi_j \left( -\frac{s_2}{r_2^2} + \frac{s_1}{r_1^2} \right) \right) \boldsymbol{n} \right)
\end{aligned} \qquad\text{(A313-4)}$$

ここで, $s_1 = \boldsymbol{r}_1 \cdot \boldsymbol{t}$, $s_2 = \boldsymbol{r}_2 \cdot \boldsymbol{t}$, および $d = \boldsymbol{r}_1 \cdot \boldsymbol{n} = \boldsymbol{r}_2 \cdot \boldsymbol{n}$ である。

### B. 理論的な背景

#### B.1 ポテンシャル流れにおけるナビエ–ストークス方程式の簡略化

　流体の速度を $\boldsymbol{v}(\boldsymbol{x}, t)$, 流体にかかる圧力を $p(\boldsymbol{x}, t)$ とする。また、流体の密度を $\rho$, 粘性係数を $\mu$, 単位体積あたりの外力を $\boldsymbol{f}(\boldsymbol{x}, t)$ とする。このとき、ナビエ–ストークス方程式は以下のように表される。

$$\rho \left( \frac{\partial \boldsymbol{v}}{\partial t} + (\boldsymbol{v} \cdot \nabla) \boldsymbol{v} \right) = -\nabla p + \mu \nabla^2 \boldsymbol{v} + \rho \boldsymbol{f} \qquad\text{(A11-1)}$$

　また、連続の式は次のように表される。

$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \boldsymbol{v}) = 0 \qquad\text{(A11-2)}$$

　非圧縮性を仮定すると, $\rho = \text{const.}$ となるため、式 $(A11-2)$ は以下のように簡略化される。

$$\nabla \cdot \boldsymbol{v} = 0 \qquad\text{(A11-3)}$$

　さらに、非回転性 ($\nabla \times \boldsymbol{v} = 0$) を仮定すると、式 $(A11-1)$ は式 $(A11-4)$ のように簡略化される。この時点で、粘性項 $\mu \nabla^2 \boldsymbol{v}$ はゼロになり、粘性の影響は無視される。

$$\begin{aligned}
    \text{LHS} &= \rho \left( \frac{\partial \boldsymbol{v}}{\partial t} + (\boldsymbol{v} \cdot \nabla) \boldsymbol{v} \right) \\\
    &= \rho \left( \frac{\partial \boldsymbol{v}}{\partial t} + \nabla \left( \frac{1}{2} |\boldsymbol{v}|^2 \right) - \boldsymbol{v} \times (\nabla \times \boldsymbol{v}) \right) \\\
    &= \rho \left( \frac{\partial \boldsymbol{v}}{\partial t} + \nabla \left( \frac{1}{2} |\boldsymbol{v}|^2 \right) \right) \\\
    \text{RHS} &= -\nabla p + \mu \nabla^2 \boldsymbol{v} + \rho \boldsymbol{f} \\\
    &= -\nabla p + (\nabla (\nabla \cdot \boldsymbol{v}) - \nabla \times (\nabla \times \boldsymbol{v})) + \rho \boldsymbol{f} \\\
    &= -\nabla p + \rho \boldsymbol{f} \\\
\end{aligned}$$

$$\therefore \quad \rho \left( \frac{\partial \boldsymbol{v}}{\partial t} + \nabla \left( \frac{1}{2} |\boldsymbol{v}|^2 \right) \right) = -\nabla p + \rho \boldsymbol{f} \qquad\text{(A11-4)}$$

> 本来、粘性を考慮する場合、壁面では滑りなし ($\boldsymbol{v} = 0$) の境界条件を課す必要がある。この場合、壁面近傍では流体に速度差が生じるため、非回転性の仮定が破れる。したがって、ポテンシャル流れでは、非粘性流体を仮定し、壁面では流体が壁面に対して滑る（すなわち、壁面に垂直な速度成分のみゼロ）と考えることで、非回転性を保っている。
>
> 本シミュレーションは、迷路を解くような（終点領域が存在する場合に、そこへ向かう流れが発生する）流れを計算することを目的としている。この場合、壁面で滑りなしの境界条件を課す必要はなく、ポテンシャル流れの仮定が妥当と考えた。

　さらに、流れの定常性 ($\frac{\partial \boldsymbol{v}}{\partial t} = 0$) と外力の不存在 ($\boldsymbol{f} = 0$) を仮定すると、式 $(A11-4)$ は以下のように簡略化される（ベルヌーイの定理の微分形）。

$$\frac{\rho}{2} \nabla |\boldsymbol{v}|^2 = -\nabla p \qquad\text{(A11-5)}$$
