from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scienceplots  # noqa: F401


ROOT = Path(__file__).resolve().parent
CSV_DIR = ROOT / "outputs" / "fem"
FIG_DIR = ROOT / "assets" / "q3_fem"
REPORT = ROOT / "结构动力学q3_FEM计算报告.md"


def set_style() -> None:
    plt.style.use(["science", "ieee", "no-latex"])
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.unicode_minus": False,
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "lines.linewidth": 1.25,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
        }
    )


def save(fig: plt.Figure, name: str) -> tuple[Path, Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    png = FIG_DIR / f"{name}.png"
    pdf = FIG_DIR / f"{name}.pdf"
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def read_parameters() -> dict[str, float]:
    df = pd.read_csv(CSV_DIR / "fem_parameters.csv")
    params: dict[str, float] = {}
    for _, row in df.iterrows():
        key = str(row["parameter"])
        val = float(row["value"])
        params[key] = val
    return params


def plot_frequency(freq: pd.DataFrame) -> tuple[Path, Path]:
    freq_plot = freq.head(8)
    fig, ax = plt.subplots(figsize=(3.5, 2.35))
    ax.bar(freq_plot["mode"], freq_plot["omega_bar"], width=0.66, color="#1f6f6b", edgecolor="black", linewidth=0.45)
    ax.set_xlabel("Mode number")
    ax.set_ylabel(r"Dimensionless frequency $\bar{\omega}_n$")
    ax.set_xticks(freq_plot["mode"])
    ymax = float(freq_plot["omega_bar"].max())
    ax.set_ylim(0, ymax * 1.14)
    for _, row in freq_plot.iterrows():
        ax.text(row["mode"], row["omega_bar"] + ymax * 0.018, f"{row['omega_bar']:.2f}", ha="center", va="bottom", fontsize=6)
    return save(fig, "fem_frequency_spectrum")


def plot_modes(modes: pd.DataFrame) -> tuple[Path, Path]:
    fig, ax = plt.subplots(figsize=(7.0, 3.1))
    styles = [
        ("#1f4e79", "-"),
        ("#b23a48", "--"),
        ("#2f6f3e", ":"),
        ("#6c4c99", "-."),
    ]
    for idx, (color, ls) in enumerate(styles, start=1):
        ax.plot(modes["xi"], modes[f"mode_{idx}"], color=color, ls=ls, label=f"Mode {idx}")
    for x, label in zip([0, 1, 2, 3, 4], ["A", "B", "C", "D", "E"]):
        ax.axvline(x, color="0.72", lw=0.65, zorder=0)
        ax.text(x, 1.13, label, ha="center", va="bottom", fontsize=7)
    ax.scatter([1], [0], marker="^", s=38, color="black", zorder=4, label="pin support")
    ax.scatter([2], [0], marker="o", s=34, facecolor="white", edgecolor="black", zorder=4, label="lumped mass")
    ax.scatter([3], [0], marker="s", s=30, facecolor="white", edgecolor="black", zorder=4, label="spring")
    ax.set_xlabel(r"Dimensionless coordinate $\xi=x/l$")
    ax.set_ylabel(r"Normalized displacement $w/\max|w|$")
    ax.set_xlim(0, 4)
    ax.set_ylim(-1.15, 1.25)
    ax.legend(ncols=2, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.04), fontsize=7)
    return save(fig, "fem_mode_shapes")


def plot_c_response(resp: pd.DataFrame) -> tuple[Path, Path]:
    fig, ax = plt.subplots(figsize=(7.0, 2.6))
    ax.plot(resp["tau"], resp["w_c"], color="#1f4e79", label=r"Displacement $\bar{w}_C$")
    ax.plot(resp["tau"], resp["v_c"], color="#b23a48", label=r"Velocity $\dot{\bar{w}}_C$")
    ax.set_xlabel(r"Dimensionless time $\tau=t/t_0$")
    ax.set_ylabel("Dimensionless response")
    ax.set_xlim(float(resp["tau"].min()), float(resp["tau"].max()))
    ax.legend(frameon=False, ncols=2, loc="upper right")
    return save(fig, "fem_c_point_response")


def plot_spacetime(snap: pd.DataFrame) -> tuple[Path, Path]:
    pivot = snap.pivot(index="tau", columns="xi", values="w")
    tau = pivot.index.to_numpy(dtype=float)
    xi = pivot.columns.to_numpy(dtype=float)
    w = pivot.to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(7.0, 3.0))
    lim = float(np.max(np.abs(w)))
    levels = np.linspace(-lim, lim, 81)
    cf = ax.contourf(xi, tau, w, levels=levels, cmap="RdBu_r", extend="both")
    for x in [1, 2, 3]:
        ax.axvline(x, color="black", lw=0.45, alpha=0.45)
    ax.set_xlabel(r"Dimensionless coordinate $\xi=x/l$")
    ax.set_ylabel(r"Dimensionless time $\tau=t/t_0$")
    cbar = fig.colorbar(cf, ax=ax, pad=0.015, shrink=0.95)
    cbar.set_label(r"Dimensionless displacement $\bar{w}$")
    return save(fig, "fem_spacetime_response")


def plot_frequency_error(freq: pd.DataFrame) -> tuple[Path, Path]:
    analytic = np.array(
        [
            1.29291120,
            1.95408261,
            3.17948221,
            6.93522035,
            12.23651680,
            16.98169320,
        ]
    )
    fem = freq["omega_bar"].to_numpy(dtype=float)[: analytic.size]
    err = np.abs((fem - analytic) / analytic) * 100.0
    fig, ax = plt.subplots(figsize=(3.5, 2.35))
    ax.semilogy(np.arange(1, analytic.size + 1), err, marker="o", color="#6c4c99", markersize=3.5)
    ax.set_xlabel("Mode number")
    ax.set_ylabel("Relative error (%)")
    ax.set_xticks(np.arange(1, analytic.size + 1))
    ax.grid(True, which="both", axis="y", lw=0.35, alpha=0.45)
    return save(fig, "fem_frequency_error")


def write_report(
    params: dict[str, float],
    freq: pd.DataFrame,
    figures: dict[str, tuple[Path, Path]],
) -> None:
    lines = [
        "# 结构动力学 Q3 有限元计算报告",
        "",
        "## 1. 有限元模型",
        "",
        "采用二节点 Euler-Bernoulli Hermite 梁单元。每个节点包含两个自由度：横向位移 $w$ 与转角 $\\theta$。计算采用无量纲形式，梁长为 $4$，弯曲刚度与分布质量密度取为 $1$。",
        "",
        "单元刚度矩阵为",
        "",
        "$$",
        r"\mathbf k_e=\frac{1}{h^3}\begin{bmatrix}12&6h&-12&6h\\6h&4h^2&-6h&2h^2\\-12&-6h&12&-6h\\6h&2h^2&-6h&4h^2\end{bmatrix}",
        "$$",
        "",
        "一致质量矩阵为",
        "",
        "$$",
        r"\mathbf m_e=\frac{h}{420}\begin{bmatrix}156&22h&54&-13h\\22h&4h^2&13h&-3h^2\\54&13h&156&-22h\\-13h&-3h^2&-22h&4h^2\end{bmatrix}",
        "$$",
        "",
        "边界和附加构件处理如下：",
        "",
        "- $B$ 点约束位移自由度 $w_B=0$，保留转角自由度。",
        "- $C$ 点集中质量加入平动质量项。",
        "- $D$ 点竖向弹簧加入平动刚度项。",
        "",
        "代表算例参数为：",
        "",
        f"- 单元数：{int(params['nel'])}",
        f"- 节点数：{int(params['nodes'])}",
        f"- 单元长度：$h={params['h']:.6g}$",
        f"- 集中质量参数：$\\alpha={params['alpha']:.6g}$",
        f"- 弹簧刚度参数：$\\kappa={params['kappa']:.6g}$",
        f"- 初速度尺度：$\\bar v_0={params['v0_bar']:.6g}$",
        "",
        "## 2. Fortran 输出文件",
        "",
        "Fortran 程序 `fem_q3_solver.f90` 组装并求解广义特征值问题",
        "",
        "$$",
        r"\mathbf K\boldsymbol\phi_n=\bar{\omega}_n^2\mathbf M\boldsymbol\phi_n",
        "$$",
        "",
        "并导出如下 CSV 文件：",
        "",
        "- `outputs/fem/fem_parameters.csv`",
        "- `outputs/fem/fem_frequencies.csv`",
        "- `outputs/fem/fem_modes.csv`",
        "- `outputs/fem/fem_c_response.csv`",
        "- `outputs/fem/fem_snapshots.csv`",
        "",
        "## 3. 固有频率结果",
        "",
        "| Mode | $\\bar\\omega_n$ | $\\beta_n=\\sqrt{\\bar\\omega_n}$ | $Y_C$ | $Y_D$ |",
        "|---:|---:|---:|---:|---:|",
    ]
    for _, row in freq.head(8).iterrows():
        lines.append(
            f"| {int(row['mode'])} | {row['omega_bar']:.8f} | {row['beta']:.8f} | {row['Y_C']:.8f} | {row['Y_D']:.8f} |"
        )
    lines += [
        "",
        "表中列出前 8 阶模态；动力响应由 Fortran 输出的前 24 阶模态叠加得到。",
        "",
        "频率谱如下图所示。",
        "",
        f"![FEM frequency spectrum]({figures['frequency'][0].as_posix()})",
        "",
        "前六阶频率与解析传递矩阵结果的相对误差如下。误差很小，说明有限元离散对低阶模态已经收敛良好。",
        "",
        f"![FEM frequency error]({figures['error'][0].as_posix()})",
        "",
        "## 4. 振型与响应",
        "",
        "前四阶归一化振型如下。",
        "",
        f"![FEM mode shapes]({figures['modes'][0].as_posix()})",
        "",
        "$C$ 点位移和速度响应如下。",
        "",
        f"![FEM C point response]({figures['c_response'][0].as_posix()})",
        "",
        "全梁时空响应如下。",
        "",
        f"![FEM spacetime response]({figures['spacetime'][0].as_posix()})",
        "",
        "## 5. 图像文件",
        "",
    ]
    for key, (png, pdf) in figures.items():
        lines.append(f"- {key}: `{png.as_posix()}`, `{pdf.as_posix()}`")
    lines += [
        "",
        "## 6. 备注",
        "",
        "本报告中的图像由 `postprocess_fem_q3.py` 从 Fortran 导出的 CSV 文件读取生成。若修改单元数、$\\alpha$、$\\kappa$ 或响应时间，只需重新编译运行 Fortran 程序，再运行 Python 后处理脚本即可更新全部图表与报告。",
        "",
    ]
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    set_style()
    params = read_parameters()
    freq = pd.read_csv(CSV_DIR / "fem_frequencies.csv")
    modes = pd.read_csv(CSV_DIR / "fem_modes.csv")
    resp = pd.read_csv(CSV_DIR / "fem_c_response.csv")
    snap = pd.read_csv(CSV_DIR / "fem_snapshots.csv")

    figures = {
        "frequency": plot_frequency(freq),
        "error": plot_frequency_error(freq),
        "modes": plot_modes(modes),
        "c_response": plot_c_response(resp),
        "spacetime": plot_spacetime(snap),
    }
    write_report(params, freq, figures)
    print(f"Report written: {REPORT}")
    for key, paths in figures.items():
        print(key, paths[0], paths[1])


if __name__ == "__main__":
    main()
