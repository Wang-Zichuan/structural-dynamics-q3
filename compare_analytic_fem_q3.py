from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scienceplots  # noqa: F401

from visualize_q3 import BeamParams, find_roots, raw_modes, state_at


ROOT = Path(__file__).resolve().parent
FEM_DIR = ROOT / "outputs" / "fem"
FIG_DIR = ROOT / "assets" / "q3_comparison"
REPORT = ROOT / "结构动力学q3_解析解-FEM对比报告.md"


def set_style() -> None:
    plt.style.use(["science", "ieee", "no-latex"])
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.unicode_minus": False,
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "lines.linewidth": 1.2,
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


def read_fem() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    return (
        pd.read_csv(FEM_DIR / "fem_frequencies.csv"),
        pd.read_csv(FEM_DIR / "fem_modes.csv"),
        pd.read_csv(FEM_DIR / "fem_c_response.csv"),
        pd.read_csv(FEM_DIR / "fem_snapshots.csv"),
    )


def analytic_mode_on_nodes(mode: dict, xi: np.ndarray, params: BeamParams) -> np.ndarray:
    beta = float(mode["beta"])
    cvec = mode["cvec"]
    y = np.array([state_at(beta, cvec, float(x), params)[0] for x in xi])
    y /= np.max(np.abs(y))
    return y


def align_to_fem(y_ana: np.ndarray, y_fem: np.ndarray) -> np.ndarray:
    return -y_ana if np.dot(y_ana, y_fem) < 0.0 else y_ana


def analytic_response(modes: list[dict], xi: np.ndarray, tau: np.ndarray, params: BeamParams) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    c_disp = np.zeros_like(tau)
    c_vel = np.zeros_like(tau)
    field = np.zeros((tau.size, xi.size))
    for mode in modes:
        beta = float(mode["beta"])
        omega = float(mode["omega_bar"])
        cvec = mode["cvec"]
        y_c = float(mode["y_c"])
        mbar = float(mode["mass_bar"])
        phi_x = np.array([state_at(beta, cvec, float(x), params)[0] for x in xi])
        coef = params.alpha * params.v0_bar * y_c / (mbar * omega)
        sinv = np.sin(omega * tau)
        cosv = np.cos(omega * tau)
        field += np.outer(sinv, coef * phi_x)
        c_disp += coef * y_c * sinv
        c_vel += coef * y_c * omega * cosv
    return c_disp, c_vel, field


def fem_field_from_snapshots(snap: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    pivot = snap.pivot(index="tau", columns="xi", values="w")
    return pivot.columns.to_numpy(dtype=float), pivot.index.to_numpy(dtype=float), pivot.to_numpy(dtype=float)


def plot_frequency_comparison(fem_freq: pd.DataFrame, ana_modes: list[dict]) -> tuple[Path, Path]:
    n = 8
    idx = np.arange(1, n + 1)
    ana = np.array([float(m["omega_bar"]) for m in ana_modes[:n]])
    fem = fem_freq["omega_bar"].to_numpy(dtype=float)[:n]
    width = 0.36
    fig, ax = plt.subplots(figsize=(3.5, 2.35))
    ax.bar(idx - width / 2, ana, width=width, color="#1f4e79", edgecolor="black", linewidth=0.35, label="Analytical")
    ax.bar(idx + width / 2, fem, width=width, color="#b23a48", edgecolor="black", linewidth=0.35, label="FEM")
    ax.set_xlabel("Mode number")
    ax.set_ylabel(r"Dimensionless frequency $\bar{\omega}_n$")
    ax.set_xticks(idx)
    ax.set_ylim(0, max(ana.max(), fem.max()) * 1.14)
    ax.legend(frameon=False, loc="upper left")
    return save(fig, "comparison_frequency_spectrum")


def plot_frequency_error(fem_freq: pd.DataFrame, ana_modes: list[dict]) -> tuple[Path, Path]:
    n = 12
    idx = np.arange(1, n + 1)
    ana = np.array([float(m["omega_bar"]) for m in ana_modes[:n]])
    fem = fem_freq["omega_bar"].to_numpy(dtype=float)[:n]
    err = np.abs((fem - ana) / ana) * 100.0
    fig, ax = plt.subplots(figsize=(3.5, 2.35))
    ax.semilogy(idx, err, marker="o", markersize=3.2, color="#6c4c99")
    ax.set_xlabel("Mode number")
    ax.set_ylabel("Relative error (%)")
    ax.set_xticks(idx)
    ax.grid(True, which="both", axis="y", lw=0.35, alpha=0.45)
    return save(fig, "comparison_frequency_error")


def plot_mode_overlay(fem_modes: pd.DataFrame, ana_modes: list[dict], params: BeamParams) -> tuple[Path, Path]:
    xi = fem_modes["xi"].to_numpy(dtype=float)
    colors = ["#1f4e79", "#b23a48", "#2f6f3e", "#6c4c99"]
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 4.6), sharex=True, sharey=True)
    axes = axes.ravel()
    for i, ax in enumerate(axes, start=1):
        fem = fem_modes[f"mode_{i}"].to_numpy(dtype=float)
        ana = analytic_mode_on_nodes(ana_modes[i - 1], xi, params)
        ana = align_to_fem(ana, fem)
        ax.plot(xi, ana, color=colors[i - 1], label="Analytical")
        ax.plot(xi, fem, color="black", ls="--", lw=0.95, label="FEM")
        for x in [1, 2, 3]:
            ax.axvline(x, color="0.75", lw=0.55, zorder=0)
        ax.set_title(f"Mode {i}", fontsize=8)
        ax.set_xlim(0, 4)
        ax.set_ylim(-1.12, 1.12)
        if i in [3, 4]:
            ax.set_xlabel(r"$\xi=x/l$")
        if i in [1, 3]:
            ax.set_ylabel(r"Normalized displacement")
    axes[0].legend(frameon=False, fontsize=7, loc="lower left")
    return save(fig, "comparison_mode_overlay")


def plot_mode_difference(fem_modes: pd.DataFrame, ana_modes: list[dict], params: BeamParams) -> tuple[Path, Path]:
    xi = fem_modes["xi"].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(7.0, 2.7))
    for i, color in zip(range(1, 5), ["#1f4e79", "#b23a48", "#2f6f3e", "#6c4c99"]):
        fem = fem_modes[f"mode_{i}"].to_numpy(dtype=float)
        ana = align_to_fem(analytic_mode_on_nodes(ana_modes[i - 1], xi, params), fem)
        ax.plot(xi, fem - ana, color=color, label=f"Mode {i}")
    for x in [1, 2, 3]:
        ax.axvline(x, color="0.75", lw=0.55, zorder=0)
    ax.set_xlabel(r"Dimensionless coordinate $\xi=x/l$")
    ax.set_ylabel("FEM - analytical")
    ax.set_xlim(0, 4)
    ax.legend(frameon=False, ncols=4, loc="upper center")
    return save(fig, "comparison_mode_difference")


def plot_c_response_comparison(fem_resp: pd.DataFrame, ana_disp: np.ndarray, ana_vel: np.ndarray) -> tuple[Path, Path]:
    tau = fem_resp["tau"].to_numpy(dtype=float)
    fig, axes = plt.subplots(2, 1, figsize=(7.0, 4.0), sharex=True)
    axes[0].plot(tau, ana_disp, color="#1f4e79", label="Analytical")
    axes[0].plot(tau, fem_resp["w_c"], color="black", ls="--", lw=0.95, label="FEM")
    axes[0].set_ylabel(r"$\bar{w}_C$")
    axes[0].legend(frameon=False, ncols=2, loc="upper right")
    axes[1].plot(tau, ana_vel, color="#b23a48", label="Analytical")
    axes[1].plot(tau, fem_resp["v_c"], color="black", ls="--", lw=0.95, label="FEM")
    axes[1].set_xlabel(r"Dimensionless time $\tau=t/t_0$")
    axes[1].set_ylabel(r"$\dot{\bar{w}}_C$")
    axes[1].legend(frameon=False, ncols=2, loc="upper right")
    axes[1].set_xlim(tau.min(), tau.max())
    return save(fig, "comparison_c_response")


def plot_c_response_error(fem_resp: pd.DataFrame, ana_disp: np.ndarray, ana_vel: np.ndarray) -> tuple[Path, Path]:
    tau = fem_resp["tau"].to_numpy(dtype=float)
    dw = fem_resp["w_c"].to_numpy(dtype=float) - ana_disp
    dv = fem_resp["v_c"].to_numpy(dtype=float) - ana_vel
    fig, ax = plt.subplots(figsize=(7.0, 2.6))
    ax.plot(tau, dw, color="#1f4e79", label=r"$\Delta \bar{w}_C$")
    ax.plot(tau, dv, color="#b23a48", label=r"$\Delta \dot{\bar{w}}_C$")
    ax.set_xlabel(r"Dimensionless time $\tau=t/t_0$")
    ax.set_ylabel("FEM - analytical")
    ax.set_xlim(tau.min(), tau.max())
    ax.legend(frameon=False, ncols=2, loc="upper right")
    return save(fig, "comparison_c_response_error")


def plot_field_difference(xi: np.ndarray, tau: np.ndarray, diff: np.ndarray) -> tuple[Path, Path]:
    fig, ax = plt.subplots(figsize=(7.0, 3.0))
    lim = max(float(np.max(np.abs(diff))), 1e-14)
    levels = np.linspace(-lim, lim, 81)
    cf = ax.contourf(xi, tau, diff, levels=levels, cmap="RdBu_r", extend="both")
    for x in [1, 2, 3]:
        ax.axvline(x, color="black", lw=0.45, alpha=0.45)
    ax.set_xlabel(r"Dimensionless coordinate $\xi=x/l$")
    ax.set_ylabel(r"Dimensionless time $\tau=t/t_0$")
    cbar = fig.colorbar(cf, ax=ax, pad=0.015, shrink=0.95)
    cbar.set_label(r"Difference $\Delta\bar{w}$")
    return save(fig, "comparison_spacetime_difference")


def comparison_table(fem_freq: pd.DataFrame, ana_modes: list[dict], n: int = 12) -> pd.DataFrame:
    rows = []
    for i in range(n):
        ana = float(ana_modes[i]["omega_bar"])
        fem = float(fem_freq.loc[i, "omega_bar"])
        rows.append(
            {
                "mode": i + 1,
                "beta_analytical": float(ana_modes[i]["beta"]),
                "omega_analytical": ana,
                "omega_fem": fem,
                "rel_error_percent": abs((fem - ana) / ana) * 100.0,
            }
        )
    return pd.DataFrame(rows)


def write_report(
    params: BeamParams,
    table: pd.DataFrame,
    mode_l2: list[float],
    c_metrics: dict[str, float],
    field_metrics: dict[str, float],
    figures: dict[str, tuple[Path, Path]],
) -> None:
    lines = [
        "# 结构动力学 Q3 解析解与有限元法对比报告",
        "",
        "## 1. 对比设置",
        "",
        "本报告比较传递矩阵解析解与 Euler-Bernoulli 梁单元有限元法。两者均采用相同无量纲参数：",
        "",
        "$$",
        rf"\alpha=\frac{{m}}{{\rho Sl}}={params.alpha},\qquad \kappa=\frac{{kl^3}}{{EJ}}={params.kappa}",
        "$$",
        "",
        "频率采用无量纲形式",
        "",
        "$$",
        r"\bar{\omega}_n=\omega_n l^2\sqrt{\frac{\rho S}{EJ}}=\beta_n^2",
        "$$",
        "",
        "FEM 结果来自 Fortran 程序 `fem_q3_solver.f90` 导出的 CSV；解析结果由传递矩阵特征方程重新计算得到。响应比较均使用前 24 阶模态截断。",
        "",
        "## 2. 固有频率对比",
        "",
        "| Mode | $\\beta_n$ analytical | $\\bar\\omega_n$ analytical | $\\bar\\omega_n$ FEM | Error (%) |",
        "|---:|---:|---:|---:|---:|",
    ]
    for _, row in table.iterrows():
        lines.append(
            f"| {int(row['mode'])} | {row['beta_analytical']:.8f} | {row['omega_analytical']:.8f} | "
            f"{row['omega_fem']:.8f} | {row['rel_error_percent']:.3e} |"
        )
    lines += [
        "",
        "前 8 阶频率柱状对比如下。",
        "",
        f"![Frequency comparison]({figures['frequency'][0].as_posix()})",
        "",
        "前 12 阶频率相对误差如下。",
        "",
        f"![Frequency error]({figures['frequency_error'][0].as_posix()})",
        "",
        "## 3. 振型对比",
        "",
        "下图叠加前四阶归一化振型。由于振型符号可任意取反，绘图时按与 FEM 振型内积为正进行符号对齐。",
        "",
        f"![Mode overlay]({figures['mode_overlay'][0].as_posix()})",
        "",
        "前四阶振型差值如下。",
        "",
        f"![Mode difference]({figures['mode_difference'][0].as_posix()})",
        "",
        "前四阶振型节点 $L_2$ 相对误差为：",
        "",
    ]
    for i, err in enumerate(mode_l2, start=1):
        lines.append(f"- Mode {i}: ${err:.3e}$")
    lines += [
        "",
        "## 4. 动力响应对比",
        "",
        "$C$ 点位移与速度响应对比如下。",
        "",
        f"![C response comparison]({figures['c_response'][0].as_posix()})",
        "",
        "$C$ 点响应差值如下。",
        "",
        f"![C response error]({figures['c_response_error'][0].as_posix()})",
        "",
        "全梁时空响应差值云图如下。",
        "",
        f"![Spacetime difference]({figures['field_difference'][0].as_posix()})",
        "",
        "响应误差指标：",
        "",
        f"- $C$ 点位移 RMS 误差：${c_metrics['disp_rms']:.3e}$",
        f"- $C$ 点速度 RMS 误差：${c_metrics['vel_rms']:.3e}$",
        f"- 全梁时空位移 RMS 误差：${field_metrics['rms']:.3e}$",
        f"- 全梁时空位移最大绝对误差：${field_metrics['max_abs']:.3e}$",
        "",
        "## 5. 图像文件",
        "",
    ]
    for key, (png, pdf) in figures.items():
        lines.append(f"- {key}: `{png.as_posix()}`, `{pdf.as_posix()}`")
    lines += [
        "",
        "## 6. 结论",
        "",
        "在当前离散参数下，FEM 前几阶固有频率与解析传递矩阵解高度一致，振型曲线几乎重合。动力响应的差异主要来自 FEM 空间离散和模态截断；随着单元数与参与模态数增加，该差异可继续降低。",
        "",
    ]
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    set_style()
    fem_freq, fem_modes, fem_resp, fem_snap = read_fem()

    params = BeamParams(alpha=0.50, kappa=20.0, n_modes=24, beta_max=40.0, scan_points=36000, response_modes=24)
    betas = find_roots(params)
    ana_modes = raw_modes(betas, params)

    xi_nodes = fem_modes["xi"].to_numpy(dtype=float)
    mode_l2 = []
    for i in range(1, 5):
        fem = fem_modes[f"mode_{i}"].to_numpy(dtype=float)
        ana = align_to_fem(analytic_mode_on_nodes(ana_modes[i - 1], xi_nodes, params), fem)
        mode_l2.append(float(np.linalg.norm(fem - ana) / np.linalg.norm(ana)))

    xi_field, tau_field, fem_field = fem_field_from_snapshots(fem_snap)
    tau_resp = fem_resp["tau"].to_numpy(dtype=float)
    ana_disp, ana_vel, _ = analytic_response(ana_modes, xi_nodes, tau_resp, params)
    _, _, ana_field = analytic_response(ana_modes, xi_field, tau_field, params)
    field_diff = fem_field - ana_field

    c_metrics = {
        "disp_rms": float(np.sqrt(np.mean((fem_resp["w_c"].to_numpy(dtype=float) - ana_disp) ** 2))),
        "vel_rms": float(np.sqrt(np.mean((fem_resp["v_c"].to_numpy(dtype=float) - ana_vel) ** 2))),
    }
    field_metrics = {
        "rms": float(np.sqrt(np.mean(field_diff**2))),
        "max_abs": float(np.max(np.abs(field_diff))),
    }

    table = comparison_table(fem_freq, ana_modes, n=12)
    (ROOT / "outputs" / "comparison").mkdir(parents=True, exist_ok=True)
    table.to_csv(ROOT / "outputs" / "comparison" / "frequency_comparison.csv", index=False)

    figures = {
        "frequency": plot_frequency_comparison(fem_freq, ana_modes),
        "frequency_error": plot_frequency_error(fem_freq, ana_modes),
        "mode_overlay": plot_mode_overlay(fem_modes, ana_modes, params),
        "mode_difference": plot_mode_difference(fem_modes, ana_modes, params),
        "c_response": plot_c_response_comparison(fem_resp, ana_disp, ana_vel),
        "c_response_error": plot_c_response_error(fem_resp, ana_disp, ana_vel),
        "field_difference": plot_field_difference(xi_field, tau_field, field_diff),
    }
    write_report(params, table, mode_l2, c_metrics, field_metrics, figures)
    print(f"Report written: {REPORT}")
    print(table.head(12).to_string(index=False))
    print("C metrics:", c_metrics)
    print("Field metrics:", field_metrics)


if __name__ == "__main__":
    main()
