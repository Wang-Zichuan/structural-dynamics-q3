from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scienceplots  # noqa: F401
from scipy.linalg import null_space
from scipy.optimize import brentq


@dataclass(frozen=True)
class BeamParams:
    alpha: float = 0.50
    kappa: float = 20.0
    n_modes: int = 6
    beta_max: float = 16.0
    scan_points: int = 12000
    response_modes: int = 6
    tau_max: float = 16.0
    v0_bar: float = 1.0


E1 = np.array([1.0, 0.0, 0.0, 0.0])
E3 = np.array([0.0, 0.0, 1.0, 0.0])
E4 = np.array([0.0, 0.0, 0.0, 1.0])
G1 = E1.copy()
G2 = np.array([0.0, 1.0, 0.0, 0.0])


def transfer_matrix(beta: float, s: float) -> np.ndarray:
    ch = np.cosh(beta * s)
    sh = np.sinh(beta * s)
    c = np.cos(beta * s)
    sn = np.sin(beta * s)
    b = beta
    return np.array(
        [
            [(ch + c) / 2, (sh + sn) / (2 * b), (ch - c) / (2 * b**2), (sh - sn) / (2 * b**3)],
            [b * (sh - sn) / 2, (ch + c) / 2, (sh + sn) / (2 * b), (ch - c) / (2 * b**2)],
            [b**2 * (ch - c) / 2, b * (sh - sn) / 2, (ch + c) / 2, (sh + sn) / (2 * b)],
            [b**3 * (sh + sn) / 2, b**2 * (ch - c) / 2, b * (sh - sn) / 2, (ch + c) / 2],
        ],
        dtype=float,
    )


def mass_jump(beta: float, alpha: float) -> np.ndarray:
    mat = np.eye(4)
    mat[3, 0] = alpha * beta**4
    return mat


def spring_jump(kappa: float) -> np.ndarray:
    mat = np.eye(4)
    mat[3, 0] = -kappa
    return mat


def characteristic_matrix(beta: float, params: BeamParams) -> np.ndarray:
    p1 = transfer_matrix(beta, 1.0)
    jm = mass_jump(beta, params.alpha)
    jk = spring_jump(params.kappa)
    t = p1 @ jk @ p1 @ jm @ p1
    return np.array(
        [
            [E1 @ p1 @ G1, E1 @ p1 @ G2, 0.0],
            [E3 @ t @ p1 @ G1, E3 @ t @ p1 @ G2, E3 @ t @ E4],
            [E4 @ t @ p1 @ G1, E4 @ t @ p1 @ G2, E4 @ t @ E4],
        ],
        dtype=float,
    )


def characteristic_det(beta: float, params: BeamParams) -> float:
    return float(np.linalg.det(characteristic_matrix(beta, params)))


def signed_log_det(beta: np.ndarray, params: BeamParams) -> np.ndarray:
    values = np.array([characteristic_det(float(b), params) for b in beta])
    return np.sign(values) * np.log10(1.0 + np.abs(values))


def find_roots(params: BeamParams) -> np.ndarray:
    beta_grid = np.linspace(0.05, params.beta_max, params.scan_points)
    det_values = np.array([characteristic_det(float(b), params) for b in beta_grid])
    roots: list[float] = []

    for i in range(len(beta_grid) - 1):
        f0, f1 = det_values[i], det_values[i + 1]
        if not np.isfinite(f0) or not np.isfinite(f1):
            continue
        if f0 == 0.0:
            candidate = beta_grid[i]
        elif f0 * f1 < 0.0:
            try:
                candidate = brentq(lambda x: characteristic_det(x, params), beta_grid[i], beta_grid[i + 1], maxiter=200)
            except ValueError:
                continue
        else:
            continue
        if candidate > 0.0 and all(abs(candidate - r) > 1e-3 for r in roots):
            roots.append(candidate)
        if len(roots) >= params.n_modes:
            break

    if len(roots) < params.n_modes:
        raise RuntimeError(f"Only found {len(roots)} roots below beta={params.beta_max}. Increase beta_max.")
    return np.array(roots[: params.n_modes])


def modal_vector(beta: float, params: BeamParams) -> np.ndarray:
    ns = null_space(characteristic_matrix(beta, params))
    if ns.size == 0:
        raise RuntimeError(f"No null vector found at beta={beta:.8g}.")
    c = ns[:, 0]
    c = c / np.max(np.abs(c[:2]))
    return np.real_if_close(c).astype(float)


def state_at(beta: float, cvec: np.ndarray, xi: float, params: BeamParams) -> np.ndarray:
    a, b, r = cvec
    z0 = np.array([a, b, 0.0, 0.0])
    p1 = transfer_matrix(beta, 1.0)

    if xi <= 1.0:
        return transfer_matrix(beta, xi) @ z0

    z1_minus = p1 @ z0
    z1_plus = z1_minus + r * E4
    if xi <= 2.0:
        return transfer_matrix(beta, xi - 1.0) @ z1_plus

    z2_minus = p1 @ z1_plus
    z2_plus = mass_jump(beta, params.alpha) @ z2_minus
    if xi <= 3.0:
        return transfer_matrix(beta, xi - 2.0) @ z2_plus

    z3_minus = p1 @ z2_plus
    z3_plus = spring_jump(params.kappa) @ z3_minus
    return transfer_matrix(beta, xi - 3.0) @ z3_plus


def mode_shape(beta: float, cvec: np.ndarray, xi: np.ndarray, params: BeamParams) -> np.ndarray:
    y = np.array([state_at(beta, cvec, float(x), params)[0] for x in xi])
    max_abs = np.max(np.abs(y))
    return y / max_abs


def modal_mass(beta: float, cvec: np.ndarray, params: BeamParams) -> float:
    xi = np.linspace(0.0, 4.0, 2001)
    y = np.array([state_at(beta, cvec, float(x), params)[0] for x in xi])
    distributed = np.trapz(y * y, xi)
    concentrated = params.alpha * state_at(beta, cvec, 2.0, params)[0] ** 2
    return float(distributed + concentrated)


def raw_modes(betas: np.ndarray, params: BeamParams) -> list[dict[str, float | np.ndarray]]:
    modes = []
    for i, beta in enumerate(betas, start=1):
        cvec = modal_vector(float(beta), params)
        mbar = modal_mass(float(beta), cvec, params)
        y_c = state_at(float(beta), cvec, 2.0, params)[0]
        y_d = state_at(float(beta), cvec, 3.0, params)[0]
        modes.append(
            {
                "index": i,
                "beta": float(beta),
                "omega_bar": float(beta**2),
                "cvec": cvec,
                "mass_bar": mbar,
                "y_c": float(y_c),
                "y_d": float(y_d),
            }
        )
    return modes


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


def save_figure(fig: plt.Figure, out_dir: Path, name: str) -> tuple[Path, Path]:
    png = out_dir / f"{name}.png"
    pdf = out_dir / f"{name}.pdf"
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def plot_characteristic(params: BeamParams, betas: np.ndarray, out_dir: Path) -> tuple[Path, Path]:
    beta_grid = np.linspace(0.05, params.beta_max, 3000)
    det_log = signed_log_det(beta_grid, params)
    fig, ax = plt.subplots(figsize=(3.5, 2.35))
    ax.plot(beta_grid, det_log, color="#1f4e79")
    ax.axhline(0, color="0.25", lw=0.7)
    for beta in betas:
        ax.axvline(beta, color="#b23a48", lw=0.65, alpha=0.55)
    ax.scatter(betas, np.zeros_like(betas), s=12, color="#b23a48", zorder=3, label="roots")
    ax.set_xlabel(r"Frequency parameter $\beta$")
    ax.set_ylabel(r"$\operatorname{sgn}(D)\log_{10}(1+|D|)$")
    ax.set_xlim(0, params.beta_max)
    ax.legend(frameon=False, loc="best")
    return save_figure(fig, out_dir, "q3_characteristic_roots")


def plot_frequency_spectrum(modes: list[dict[str, float | np.ndarray]], out_dir: Path) -> tuple[Path, Path]:
    idx = np.array([m["index"] for m in modes], dtype=int)
    omega = np.array([m["omega_bar"] for m in modes], dtype=float)
    fig, ax = plt.subplots(figsize=(3.5, 2.35))
    ax.bar(idx, omega, width=0.64, color="#246b69", edgecolor="black", linewidth=0.45)
    ax.set_xlabel("Mode number")
    ax.set_ylabel(r"Dimensionless frequency $\bar{\omega}_n=\beta_n^2$")
    ax.set_xticks(idx)
    ax.set_ylim(0, omega.max() * 1.16)
    for x, y in zip(idx, omega):
        ax.text(x, y + omega.max() * 0.025, f"{y:.2f}", ha="center", va="bottom", fontsize=6)
    return save_figure(fig, out_dir, "q3_frequency_spectrum")


def plot_mode_shapes(modes: list[dict[str, float | np.ndarray]], params: BeamParams, out_dir: Path) -> tuple[Path, Path]:
    xi = np.linspace(0.0, 4.0, 1201)
    fig, ax = plt.subplots(figsize=(7.0, 3.1))
    colors = ["#1f4e79", "#b23a48", "#2f6f3e", "#6c4c99", "#a85f00", "#2b7a78"]
    for mode, color in zip(modes[:4], colors):
        y = mode_shape(float(mode["beta"]), mode["cvec"], xi, params)
        ax.plot(xi, y, color=color, label=fr"Mode {int(mode['index'])}, $\beta$={float(mode['beta']):.3f}")
    for x, label in zip([0, 1, 2, 3, 4], ["A", "B", "C", "D", "E"]):
        ax.axvline(x, color="0.72", lw=0.65, zorder=0)
        ax.text(x, 1.13, label, ha="center", va="bottom", fontsize=7)
    ax.scatter([1], [0], marker="^", s=38, color="black", zorder=4, label="pin support")
    ax.scatter([2], [0], marker="o", s=34, facecolor="white", edgecolor="black", zorder=4, label="lumped mass")
    ax.scatter([3], [0], marker="s", s=30, facecolor="white", edgecolor="black", zorder=4, label="spring")
    ax.set_xlabel(r"Dimensionless coordinate $\xi=x/l$")
    ax.set_ylabel(r"Normalized mode shape $Y_n/\max|Y_n|$")
    ax.set_xlim(0, 4)
    ax.set_ylim(-1.15, 1.25)
    ax.legend(ncols=2, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.04), fontsize=7)
    return save_figure(fig, out_dir, "q3_mode_shapes")


def response_fields(modes: list[dict[str, float | np.ndarray]], params: BeamParams) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    xi = np.linspace(0.0, 4.0, 401)
    tau = np.linspace(0.0, params.tau_max, 1601)
    response = np.zeros((tau.size, xi.size))
    c_response = np.zeros_like(tau)
    c_velocity = np.zeros_like(tau)

    for mode in modes[: params.response_modes]:
        beta = float(mode["beta"])
        omega = float(mode["omega_bar"])
        cvec = mode["cvec"]
        mbar = float(mode["mass_bar"])
        yc = float(mode["y_c"])
        phi_x = np.array([state_at(beta, cvec, float(x), params)[0] for x in xi])
        coef = params.alpha * params.v0_bar * yc / (mbar * omega)
        response += np.outer(np.sin(omega * tau), coef * phi_x)
        c_response += coef * yc * np.sin(omega * tau)
        c_velocity += coef * yc * omega * np.cos(omega * tau)
    return xi, tau, response, np.vstack([c_response, c_velocity])


def plot_c_response(tau: np.ndarray, c_data: np.ndarray, out_dir: Path) -> tuple[Path, Path]:
    fig, ax = plt.subplots(figsize=(7.0, 2.6))
    ax.plot(tau, c_data[0], color="#1f4e79", label=r"Displacement $\bar{w}_C$")
    ax.plot(tau, c_data[1], color="#b23a48", label=r"Velocity $\dot{\bar{w}}_C$")
    ax.set_xlabel(r"Dimensionless time $\tau=t/t_0$")
    ax.set_ylabel("Dimensionless response")
    ax.set_xlim(tau.min(), tau.max())
    ax.legend(frameon=False, ncols=2, loc="upper right")
    return save_figure(fig, out_dir, "q3_c_point_response")


def plot_spacetime(xi: np.ndarray, tau: np.ndarray, response: np.ndarray, out_dir: Path) -> tuple[Path, Path]:
    fig, ax = plt.subplots(figsize=(7.0, 3.0))
    lim = np.max(np.abs(response))
    levels = np.linspace(-lim, lim, 81)
    contour = ax.contourf(xi, tau, response, levels=levels, cmap="RdBu_r", extend="both")
    for x in [1, 2, 3]:
        ax.axvline(x, color="black", lw=0.45, alpha=0.45)
    ax.set_xlabel(r"Dimensionless coordinate $\xi=x/l$")
    ax.set_ylabel(r"Dimensionless time $\tau=t/t_0$")
    cbar = fig.colorbar(contour, ax=ax, pad=0.015, shrink=0.95)
    cbar.set_label(r"Dimensionless displacement $\bar{w}$")
    return save_figure(fig, out_dir, "q3_spacetime_response")


def write_report(
    out_path: Path,
    params: BeamParams,
    modes: list[dict[str, float | np.ndarray]],
    image_paths: dict[str, tuple[Path, Path]],
) -> None:
    lines = [
        "# 结构动力学 Q3 可视化计算报告",
        "",
        "## 1. 计算模型",
        "",
        "本报告对解析解进行无量纲数值可视化。采用的无量纲参数为",
        "",
        "$$",
        r"\alpha=\frac{m}{\rho Sl},\qquad \kappa=\frac{kl^3}{EJ}",
        "$$",
        "",
        "代表算例参数如下：",
        "",
        f"- 集中质量参数：$\\alpha={params.alpha:.3g}$",
        f"- 弹簧刚度参数：$\\kappa={params.kappa:.3g}$",
        f"- 初速度尺度：$\\bar v_0={params.v0_bar:.3g}$",
        f"- 响应叠加模态数：{params.response_modes}",
        "",
        "无量纲固有频率定义为",
        "",
        "$$",
        r"\bar{\omega}_n=\omega_n l^2\sqrt{\frac{\rho S}{EJ}}=\beta_n^2",
        "$$",
        "",
        "实际圆频率可由",
        "",
        "$$",
        r"\omega_n=\frac{\bar{\omega}_n}{l^2}\sqrt{\frac{EJ}{\rho S}}",
        "$$",
        "",
        "恢复。",
        "",
        "## 2. 特征根与固有频率",
        "",
        "| Mode | $\\beta_n$ | $\\bar{\\omega}_n=\\beta_n^2$ | $\\bar M_n$ | $Y_n(2)$ | $Y_n(3)$ |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for mode in modes:
        lines.append(
            f"| {int(mode['index'])} | {float(mode['beta']):.6f} | {float(mode['omega_bar']):.6f} | "
            f"{float(mode['mass_bar']):.6f} | {float(mode['y_c']):.6f} | {float(mode['y_d']):.6f} |"
        )

    lines += [
        "",
        "特征方程扫描与根的位置如下图所示。",
        "",
        f"![Characteristic roots]({image_paths['characteristic'][0].as_posix()})",
        "",
        "前几阶无量纲频率谱如下。",
        "",
        f"![Frequency spectrum]({image_paths['spectrum'][0].as_posix()})",
        "",
        "## 3. 振型函数",
        "",
        "下图给出前四阶振型，所有振型均按各自最大绝对位移归一化。图中 $B$ 为外部简单支承点，$C$ 为集中质量位置，$D$ 为弹簧位置。",
        "",
        f"![Mode shapes]({image_paths['modes'][0].as_posix()})",
        "",
        "## 4. 动力响应",
        "",
        "响应采用如下无量纲时间",
        "",
        "$$",
        r"\tau=\frac{t}{t_0},\qquad t_0=l^2\sqrt{\frac{\rho S}{EJ}}",
        "$$",
        "",
        "并取初位移为零，初始动量只由 $C$ 点集中质量的速度贡献。$C$ 点响应如下。",
        "",
        f"![C point response]({image_paths['c_response'][0].as_posix()})",
        "",
        "全梁时空响应云图如下。",
        "",
        f"![Spacetime response]({image_paths['spacetime'][0].as_posix()})",
        "",
        "## 5. 输出文件",
        "",
    ]
    for label, (png, pdf) in image_paths.items():
        lines.append(f"- {label}: `{png.as_posix()}`, `{pdf.as_posix()}`")

    lines += [
        "",
        "## 6. 说明",
        "",
        "由于题目没有给出 $m$、$k$、$EJ$、$\\rho S$ 和 $l$ 的具体数值，本报告采用无量纲代表算例展示解析解的计算流程。若需要换成具体工程参数，只需修改 Python 脚本中的 `BeamParams(alpha=..., kappa=...)`，并通过频率恢复公式换算为实际单位。",
        "",
    ]
    out_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    params = BeamParams()
    root = Path(__file__).resolve().parent
    out_dir = root / "assets" / "q3_visualization"
    out_dir.mkdir(parents=True, exist_ok=True)

    set_style()
    betas = find_roots(params)
    modes = raw_modes(betas, params)

    image_paths: dict[str, tuple[Path, Path]] = {}
    image_paths["characteristic"] = plot_characteristic(params, betas, out_dir)
    image_paths["spectrum"] = plot_frequency_spectrum(modes, out_dir)
    image_paths["modes"] = plot_mode_shapes(modes, params, out_dir)
    xi, tau, response, c_data = response_fields(modes, params)
    image_paths["c_response"] = plot_c_response(tau, c_data, out_dir)
    image_paths["spacetime"] = plot_spacetime(xi, tau, response, out_dir)

    write_report(root / "结构动力学q3_可视化计算报告.md", params, modes, image_paths)
    print("Computed roots:")
    for mode in modes:
        print(f"mode {int(mode['index'])}: beta={float(mode['beta']):.8f}, omega_bar={float(mode['omega_bar']):.8f}")
    print(f"Figures and report written to {out_dir}")


if __name__ == "__main__":
    main()
