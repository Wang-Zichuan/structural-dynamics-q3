# 结构动力学 Q3 有限元计算报告

## 1. 有限元模型

采用二节点 Euler-Bernoulli Hermite 梁单元。每个节点包含两个自由度：横向位移 $w$ 与转角 $\theta$。计算采用无量纲形式，梁长为 $4$，弯曲刚度与分布质量密度取为 $1$。

单元刚度矩阵为

$$
\mathbf k_e=\frac{1}{h^3}\begin{bmatrix}12&6h&-12&6h\\6h&4h^2&-6h&2h^2\\-12&-6h&12&-6h\\6h&2h^2&-6h&4h^2\end{bmatrix}
$$

一致质量矩阵为

$$
\mathbf m_e=\frac{h}{420}\begin{bmatrix}156&22h&54&-13h\\22h&4h^2&13h&-3h^2\\54&13h&156&-22h\\-13h&-3h^2&-22h&4h^2\end{bmatrix}
$$

边界和附加构件处理如下：

- $B$ 点约束位移自由度 $w_B=0$，保留转角自由度。
- $C$ 点集中质量加入平动质量项。
- $D$ 点竖向弹簧加入平动刚度项。

代表算例参数为：

- 单元数：80
- 节点数：81
- 单元长度：$h=0.05$
- 集中质量参数：$\alpha=0.5$
- 弹簧刚度参数：$\kappa=20$
- 初速度尺度：$\bar v_0=1$

## 2. Fortran 输出文件

Fortran 程序 `fem_q3_solver.f90` 组装并求解广义特征值问题

$$
\mathbf K\boldsymbol\phi_n=\bar{\omega}_n^2\mathbf M\boldsymbol\phi_n
$$

并导出如下 CSV 文件：

- `outputs/fem/fem_parameters.csv`
- `outputs/fem/fem_frequencies.csv`
- `outputs/fem/fem_modes.csv`
- `outputs/fem/fem_c_response.csv`
- `outputs/fem/fem_snapshots.csv`

## 3. 固有频率结果

| Mode | $\bar\omega_n$ | $\beta_n=\sqrt{\bar\omega_n}$ | $Y_C$ | $Y_D$ |
|---:|---:|---:|---:|---:|
| 1 | 1.29291121 | 1.13706253 | -0.47061026 | 0.00682808 |
| 2 | 1.95408264 | 1.39788506 | 0.18938252 | 0.25452487 |
| 3 | 3.17948233 | 1.78311030 | -0.58440807 | -0.36727120 |
| 4 | 6.93522170 | 2.63348091 | 0.35593038 | -0.91072216 |
| 5 | 12.23652472 | 3.49807443 | 0.05501317 | -0.35799071 |
| 6 | 16.98171295 | 4.12088740 | 0.36805753 | -0.13155424 |
| 7 | 21.44393578 | 4.63075974 | -0.28544241 | 0.35790956 |
| 8 | 30.19935982 | 5.49539442 | 0.33860352 | 0.89620846 |

表中列出前 8 阶模态；动力响应由 Fortran 输出的前 24 阶模态叠加得到。

频率谱如下图所示。

![FEM frequency spectrum](D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_fem/fem_frequency_spectrum.png)

前六阶频率与解析传递矩阵结果的相对误差如下。误差很小，说明有限元离散对低阶模态已经收敛良好。

![FEM frequency error](D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_fem/fem_frequency_error.png)

## 4. 振型与响应

前四阶归一化振型如下。

![FEM mode shapes](D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_fem/fem_mode_shapes.png)

$C$ 点位移和速度响应如下。

![FEM C point response](D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_fem/fem_c_point_response.png)

全梁时空响应如下。

![FEM spacetime response](D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_fem/fem_spacetime_response.png)

## 5. 图像文件

- frequency: `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_fem/fem_frequency_spectrum.png`, `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_fem/fem_frequency_spectrum.pdf`
- error: `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_fem/fem_frequency_error.png`, `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_fem/fem_frequency_error.pdf`
- modes: `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_fem/fem_mode_shapes.png`, `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_fem/fem_mode_shapes.pdf`
- c_response: `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_fem/fem_c_point_response.png`, `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_fem/fem_c_point_response.pdf`
- spacetime: `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_fem/fem_spacetime_response.png`, `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_fem/fem_spacetime_response.pdf`

## 6. 备注

本报告中的图像由 `postprocess_fem_q3.py` 从 Fortran 导出的 CSV 文件读取生成。若修改单元数、$\alpha$、$\kappa$ 或响应时间，只需重新编译运行 Fortran 程序，再运行 Python 后处理脚本即可更新全部图表与报告。
