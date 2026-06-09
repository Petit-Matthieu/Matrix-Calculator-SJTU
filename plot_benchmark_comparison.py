#!/usr/bin/env python3

import csv
import html
import math
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "benchmark_comparison.svg"
WIDTH = 1800
HEIGHT = 1800


def run_benchmark(name):
    completed = subprocess.run(
        [str(ROOT / name)],
        cwd=ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return completed.stdout


def extract_csv(output, header_prefix):
    lines = output.splitlines()
    start = next(index for index, line in enumerate(lines) if line.startswith(header_prefix))
    rows = []
    for line in lines[start + 1 :]:
        if not line.strip() or line.startswith("validation_guard"):
            break
        rows.append(line)
    return list(csv.DictReader([lines[start], *rows]))


def escape(text):
    return html.escape(str(text), quote=True)


class Svg:
    def __init__(self):
        self.parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" '
            f'viewBox="0 0 {WIDTH} {HEIGHT}">',
            "<style>",
            "text { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; }",
            ".title { font-size: 34px; font-weight: 700; fill: #172033; }",
            ".subtitle { font-size: 17px; fill: #5b6472; }",
            ".panel-title { font-size: 22px; font-weight: 700; fill: #172033; }",
            ".panel-subtitle { font-size: 15px; fill: #667085; }",
            ".axis { font-size: 14px; fill: #667085; }",
            ".legend { font-size: 14px; fill: #344054; }",
            ".note { font-size: 15px; fill: #475467; }",
            "</style>",
            f'<rect width="{WIDTH}" height="{HEIGHT}" fill="#f7f9fc"/>',
        ]

    def add(self, markup):
        self.parts.append(markup)

    def rect(self, x, y, width, height, fill, stroke="none", radius=0, stroke_width=1):
        self.add(
            f'<rect x="{x}" y="{y}" width="{width}" height="{height}" rx="{radius}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}"/>'
        )

    def line(self, x1, y1, x2, y2, stroke, stroke_width=1, dash=None):
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(
            f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
            f'stroke="{stroke}" stroke-width="{stroke_width}"{dash_attr}/>'
        )

    def circle(self, cx, cy, radius, fill, stroke="#ffffff", stroke_width=2):
        self.add(
            f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{radius}" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{stroke_width}"/>'
        )

    def text(self, x, y, value, css_class, anchor="start", weight=None):
        weight_attr = f' font-weight="{weight}"' if weight else ""
        self.add(
            f'<text x="{x:.2f}" y="{y:.2f}" class="{css_class}" '
            f'text-anchor="{anchor}"{weight_attr}>{escape(value)}</text>'
        )

    def rotated_text(self, x, y, value, css_class):
        self.add(
            f'<text x="{x:.2f}" y="{y:.2f}" class="{css_class}" text-anchor="middle" '
            f'transform="rotate(-90 {x:.2f} {y:.2f})">{escape(value)}</text>'
        )

    def path(self, points, stroke, stroke_width=3):
        coordinates = " ".join(
            ("M" if index == 0 else "L") + f" {x:.2f} {y:.2f}"
            for index, (x, y) in enumerate(points)
        )
        self.add(
            f'<path d="{coordinates}" fill="none" stroke="{stroke}" '
            f'stroke-width="{stroke_width}" stroke-linejoin="round" stroke-linecap="round"/>'
        )

    def finish(self):
        self.parts.append("</svg>")
        OUTPUT.write_text("\n".join(self.parts) + "\n", encoding="utf-8")


def nice_linear_max(value):
    if value <= 0:
        return 1
    magnitude = 10 ** math.floor(math.log10(value))
    normalized = value / magnitude
    if normalized <= 1:
        top = 1
    elif normalized <= 2:
        top = 2
    elif normalized <= 5:
        top = 5
    else:
        top = 10
    return top * magnitude


def draw_panel(svg, x, y, width, height, title, subtitle):
    svg.rect(x, y, width, height, "#ffffff", "#d9e0ea", 8)
    svg.text(x + 26, y + 38, title, "panel-title")
    svg.text(x + 26, y + 62, subtitle, "panel-subtitle")
    return {
        "x": x + 82,
        "y": y + 91,
        "width": width - 112,
        "height": height - 150,
    }


def evenly_spaced_x(plot, labels):
    if len(labels) == 1:
        return [plot["x"] + plot["width"] / 2]
    return [
        plot["x"] + plot["width"] * index / (len(labels) - 1)
        for index in range(len(labels))
    ]


def draw_x_axis(svg, plot, labels, axis_label):
    xs = evenly_spaced_x(plot, labels)
    bottom = plot["y"] + plot["height"]
    svg.line(plot["x"], bottom, plot["x"] + plot["width"], bottom, "#98a2b3", 1.2)
    for x, label in zip(xs, labels):
        svg.line(x, bottom, x, bottom + 6, "#98a2b3")
        svg.text(x, bottom + 24, label, "axis", "middle")
    svg.text(plot["x"] + plot["width"] / 2, bottom + 48, axis_label, "axis", "middle")
    return xs


def draw_log_y_axis(svg, plot, values, axis_label):
    minimum = 10 ** math.floor(math.log10(min(values)))
    maximum = 10 ** math.ceil(math.log10(max(values)))
    log_minimum = math.log10(minimum)
    log_maximum = math.log10(maximum)

    def map_value(value):
        ratio = (math.log10(value) - log_minimum) / (log_maximum - log_minimum)
        return plot["y"] + plot["height"] * (1 - ratio)

    exponent = int(log_minimum)
    while exponent <= int(log_maximum):
        value = 10**exponent
        y = map_value(value)
        svg.line(plot["x"], y, plot["x"] + plot["width"], y, "#e8ecf2")
        svg.text(plot["x"] - 12, y + 5, f"{value:g}", "axis", "end")
        exponent += 1
    svg.rotated_text(plot["x"] - 58, plot["y"] + plot["height"] / 2, axis_label, "axis")
    return map_value


def draw_linear_y_axis(svg, plot, values, axis_label):
    maximum = nice_linear_max(max(values) * 1.1)
    steps = 5

    def map_value(value):
        return plot["y"] + plot["height"] * (1 - value / maximum)

    for step in range(steps + 1):
        value = maximum * step / steps
        y = map_value(value)
        svg.line(plot["x"], y, plot["x"] + plot["width"], y, "#e8ecf2")
        svg.text(plot["x"] - 12, y + 5, f"{value:g}", "axis", "end")
    svg.rotated_text(plot["x"] - 58, plot["y"] + plot["height"] / 2, axis_label, "axis")
    return map_value


def draw_series(svg, xs, values, map_y, color):
    points = [(x, map_y(value)) for x, value in zip(xs, values)]
    svg.path(points, color)
    for x, y in points:
        svg.circle(x, y, 5, color)


def draw_legend(svg, x, y, items, columns=1, column_width=170):
    for index, (label, color) in enumerate(items):
        column = index % columns
        row = index // columns
        item_x = x + column * column_width
        item_y = y + row * 24
        svg.line(item_x, item_y - 4, item_x + 24, item_y - 4, color, 3)
        svg.circle(item_x + 12, item_y - 4, 4, color)
        svg.text(item_x + 33, item_y, label, "legend")


def float_column(rows, name):
    return [float(row[name]) for row in rows]


def draw_multiplication_chart(svg, rows):
    plot = draw_panel(
        svg,
        48,
        118,
        1704,
        570,
        "Matrix multiplication runtime",
        "Lower is better. Accelerate uses BLAS cblas_dgemm; repository versions use matrix_ops.c.",
    )
    labels = [row["n"] for row in rows]
    xs = draw_x_axis(svg, plot, labels, "square matrix size n")
    series = [
        ("Original i-j-k", "#e5484d", float_column(rows, "ijk_ms")),
        ("Reordered i-k-j", "#f59e0b", float_column(rows, "ikj_ms")),
        ("Blocked", "#0f9f8f", float_column(rows, "blocked_ms")),
        ("Optimized portable C", "#7c3aed", float_column(rows, "optimized_ms")),
        ("Accelerate BLAS", "#2563eb", float_column(rows, "blas_ms")),
    ]
    map_y = draw_log_y_axis(
        svg, plot, [value for _, _, values in series for value in values], "runtime (ms, log scale)"
    )
    for _, color, values in series:
        draw_series(svg, xs, values, map_y, color)
    draw_legend(svg, plot["x"] + 18, plot["y"] + 22, [(label, color) for label, color, _ in series], 5, 185)


def draw_operation_speedup_chart(svg, rows):
    plot = draw_panel(
        svg,
        48,
        730,
        828,
        1000,
        "New operations: Accelerate speedup",
        "Custom runtime divided by Accelerate runtime. Above 1x means Accelerate is faster.",
    )
    labels = [row["n"] for row in rows]
    xs = draw_x_axis(svg, plot, labels, "square matrix size n")
    series = [
        ("Add", "#2563eb", float_column(rows, "add_speedup")),
        ("Subtract", "#7c3aed", float_column(rows, "sub_speedup")),
        ("Scale", "#059669", float_column(rows, "scale_speedup")),
        ("Transpose", "#f97316", float_column(rows, "transpose_speedup")),
        ("Frobenius norm", "#e5484d", float_column(rows, "norm_speedup")),
    ]
    map_y = draw_linear_y_axis(
        svg, plot, [value for _, _, values in series for value in values], "Accelerate speedup (x)"
    )
    svg.line(plot["x"], map_y(1.0), plot["x"] + plot["width"], map_y(1.0), "#98a2b3", 1.5, "6 5")
    svg.text(plot["x"] + plot["width"] - 2, map_y(1.0) - 8, "equal speed", "axis", "end")
    for _, color, values in series:
        draw_series(svg, xs, values, map_y, color)
    draw_legend(svg, plot["x"] + 14, plot["y"] + 20, [(label, color) for label, color, _ in series], 2, 178)


def draw_lu_chart(svg, rows):
    plot = draw_panel(
        svg,
        924,
        730,
        828,
        1000,
        "LU factorization runtime",
        "Lower is better. Repository no-pivot LU versus Accelerate LAPACK dgetrf.",
    )
    labels = [row["n"] for row in rows]
    xs = draw_x_axis(svg, plot, labels, "square matrix size n")
    series = [
        ("Repository LU", "#e5484d", float_column(rows, "custom_ms")),
        ("Accelerate LAPACK", "#2563eb", float_column(rows, "lapack_ms")),
    ]
    map_y = draw_log_y_axis(
        svg, plot, [value for _, _, values in series for value in values], "runtime (ms, log scale)"
    )
    for _, color, values in series:
        draw_series(svg, xs, values, map_y, color)
    draw_legend(svg, plot["x"] + 15, plot["y"] + 20, [(label, color) for label, color, _ in series], 2, 205)
    maximum_speedup = max(float_column(rows, "lapack_speedup"))
    svg.text(
        plot["x"] + plot["width"],
        plot["y"] + plot["height"] - 14,
        f"Peak measured LAPACK speedup: {maximum_speedup:.1f}x",
        "note",
        "end",
        "600",
    )


def main():
    operations_output = run_benchmark("benchmark_ops_vs_accelerate")
    lu_output = run_benchmark("benchmark_lu_vs_accelerate")
    operation_rows = extract_csv(operations_output, "n,add_custom_ms")
    multiplication_rows = extract_csv(operations_output, "n,ijk_ms")
    lu_rows = extract_csv(lu_output, "n,runs,custom_ms")

    svg = Svg()
    svg.text(48, 58, "Matrix Calculator Performance Comparison", "title")
    svg.text(
        48,
        88,
        "Repository implementation compared with Apple Accelerate on this Mac. Timed regions exclude initialization.",
        "subtitle",
    )
    draw_multiplication_chart(svg, multiplication_rows)
    draw_operation_speedup_chart(svg, operation_rows)
    draw_lu_chart(svg, lu_rows)
    svg.finish()
    print(OUTPUT)


if __name__ == "__main__":
    main()
