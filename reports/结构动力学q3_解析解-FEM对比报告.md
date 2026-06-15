# 结构动力学 Q3 解析解与有限元法对比报告

## 1. 对比设置

本报告比较传递矩阵解析解与 Euler-Bernoulli 梁单元有限元法。两者均采用相同无量纲参数：

$$
\alpha=\frac{m}{\rho Sl}=0.5,\qquad \kappa=\frac{kl^3}{EJ}=20.0
$$

频率采用无量纲形式

$$
\bar{\omega}_n=\omega_n l^2\sqrt{\frac{\rho S}{EJ}}=\beta_n^2
$$

FEM 结果来自 Fortran 程序 `fem_q3_solver.f90` 导出的 CSV；解析结果由传递矩阵特征方程重新计算得到。响应比较均使用前 24 阶模态截断。

## 2. 固有频率对比

| Mode | $\beta_n$ analytical | $\bar\omega_n$ analytical | $\bar\omega_n$ FEM | Error (%) |
|---:|---:|---:|---:|---:|
| 1 | 1.13706253 | 1.29291120 | 1.29291121 | 6.561e-07 |
| 2 | 1.39788505 | 1.95408261 | 1.95408264 | 1.614e-06 |
| 3 | 1.78311026 | 3.17948221 | 3.17948233 | 3.652e-06 |
| 4 | 2.63348065 | 6.93522035 | 6.93522170 | 1.953e-05 |
| 5 | 3.49807330 | 12.23651680 | 12.23652472 | 6.477e-05 |
| 6 | 4.12088500 | 16.98169320 | 16.98171295 | 1.163e-04 |
| 7 | 4.63075532 | 21.44389486 | 21.44393578 | 1.908e-04 |
| 8 | 5.49538421 | 30.19924764 | 30.19935982 | 3.715e-04 |
| 9 | 6.62342308 | 43.86973326 | 43.87010941 | 8.574e-04 |
| 10 | 7.22338705 | 52.17732044 | 52.17760933 | 5.537e-04 |
| 11 | 7.75075013 | 60.07412758 | 60.07546580 | 2.228e-03 |
| 12 | 8.58187870 | 73.64864210 | 73.67590498 | 3.702e-02 |

前 8 阶频率柱状对比如下。

![Frequency comparison](D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_frequency_spectrum.png)

前 12 阶频率相对误差如下。

![Frequency error](D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_frequency_error.png)

## 3. 振型对比

下图叠加前四阶归一化振型。由于振型符号可任意取反，绘图时按与 FEM 振型内积为正进行符号对齐。

![Mode overlay](D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_mode_overlay.png)

前四阶振型差值如下。

![Mode difference](D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_mode_difference.png)

前四阶振型节点 $L_2$ 相对误差为：

- Mode 1: $2.455e-09$
- Mode 2: $1.727e-09$
- Mode 3: $2.404e-08$
- Mode 4: $5.903e-08$

## 4. 动力响应对比

$C$ 点位移与速度响应对比如下。

![C response comparison](D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_c_response.png)

$C$ 点响应差值如下。

![C response error](D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_c_response_error.png)

全梁时空响应差值云图如下。

![Spacetime difference](D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_spacetime_difference.png)

响应误差指标：

- $C$ 点位移 RMS 误差：$1.270e-03$
- $C$ 点速度 RMS 误差：$1.008e-01$
- 全梁时空位移 RMS 误差：$3.086e-03$
- 全梁时空位移最大绝对误差：$1.529e-02$

## 5. 图像文件

- frequency: `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_frequency_spectrum.png`, `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_frequency_spectrum.pdf`
- frequency_error: `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_frequency_error.png`, `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_frequency_error.pdf`
- mode_overlay: `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_mode_overlay.png`, `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_mode_overlay.pdf`
- mode_difference: `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_mode_difference.png`, `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_mode_difference.pdf`
- c_response: `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_c_response.png`, `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_c_response.pdf`
- c_response_error: `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_c_response_error.png`, `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_c_response_error.pdf`
- field_difference: `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_spacetime_difference.png`, `D:/DeskTopProjects/jiegoudonglixueq3/assets/q3_comparison/comparison_spacetime_difference.pdf`

## 6. 结论

在当前离散参数下，FEM 前几阶固有频率与解析传递矩阵解高度一致，振型曲线几乎重合。动力响应的差异主要来自 FEM 空间离散和模态截断；随着单元数与参与模态数增加，该差异可继续降低。
