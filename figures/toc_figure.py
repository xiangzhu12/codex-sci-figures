#!/usr/bin/env python3
"""Vector TOC schematic for the uploaded FeN4–COF manuscript.

Run: python figures/toc_figure.py --font 'Liberation Sans'
For Arial: python figures/toc_figure.py --font Arial [--font-dir /path/to/fonts]
Only matplotlib and numpy are required. No source images are embedded.
All molecular drawings are schematic, not atomistic geometry models.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Circle, Polygon, FancyArrowPatch, FancyBboxPatch, Ellipse

OUT = Path(__file__).resolve().parent
C = dict(ink='#243447', grey='#687887', blue='#307EA3', pale='#EFF6FA',
         purple='#80549D', lilac='#F6F1F9', orange='#D28A21', green='#197B70',
         red='#C45651', carbon='#536473')
RATE = 77.33  # μmol g−1 h−1; average over the reported 3 h run
SELECTIVITY = 90  # % CO
FE_N_SHORT, FE_N_LONG = 1.939, 2.035  # Å; EXAFS fit, not drawing scale
TA_BEFORE, TA_AFTER = 52.71, 849.20  # ps; slow fs-TA decay components


def line(ax, a, b, color=C['carbon'], lw=1.15, **kw):
    ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=lw,
            solid_capstyle='round', **kw)


def text(ax, x, y, s, size=11, color=C['ink'], weight='normal', **kw):
    # Some Arial-compatible fonts lack Unicode superscript +/- glyphs.
    # MathText draws superscripts using the selected family, avoiding tofu.
    s = s.replace('e⁻', r'$\mathrm{e}^{-}$').replace('H⁺', r'$\mathrm{H}^{+}$')
    s = s.replace('g⁻¹', r'$\mathrm{g}^{-1}$').replace('h⁻¹', r'$\mathrm{h}^{-1}$')
    return ax.text(x, y, s, ha='center', va='center', fontsize=size,
                   color=color, fontweight=weight, **kw)


def arrow(ax, a, b, color=C['purple'], lw=2, scale=15, rad=0, **kw):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle='-|>', mutation_scale=scale,
                               lw=lw, color=color,
                               connectionstyle=f'arc3,rad={rad}', **kw))


def benzene(ax, center, radius, rotation=0, color=C['carbon'], lw=1):
    center = np.asarray(center)
    angles = rotation + np.arange(6) * np.pi / 3
    pts = center + radius * np.c_[np.cos(angles), np.sin(angles)]
    ax.add_patch(Polygon(pts, closed=True, fill=False, ec=color, lw=lw,
                         joinstyle='round'))
    for i in (0, 2, 4):
        a, b = pts[i], pts[(i + 1) % 6]
        aa, bb = center + .73 * (a - center), center + .73 * (b - center)
        line(ax, aa, bb, color, lw * .75)
    return pts


def porphyrin(ax, center, r=.35, color=C['blue'], enlarged=False):
    """Four pyrrole rings joined through meso carbons; N atoms face inward.

    Phenyl substituents emerge from the four meso positions. Alternating
    bonds are schematic, not a unique Kekulé assignment of a metalloporphyrin.
    """
    c = np.asarray(center)
    nitrogens = []
    for theta in np.pi / 4 + np.arange(4) * np.pi / 2:
        u = np.array([np.cos(theta), np.sin(theta)])
        v = np.array([-np.sin(theta), np.cos(theta)])
        local = np.array([[.52, 0], [.87, .31], [1.23, .20],
                          [1.23, -.20], [.87, -.31]])
        pts = c + r * (local[:, :1] * u + local[:, 1:] * v)
        ax.add_patch(Polygon(pts, closed=True, fc='white', ec=color,
                             lw=1.25 if enlarged else .9, zorder=3))
        for i, j in ((1, 2), (3, 4)):
            mid = pts.mean(axis=0)
            line(ax, mid + .78 * (pts[i] - mid),
                 mid + .78 * (pts[j] - mid), color,
                 .9 if enlarged else .6, zorder=4)
        nitrogens.append(pts[0])
        if enlarged:
            text(ax, *pts[0], 'N', 12, color, 'bold',
                 bbox=dict(fc='white', ec='none', pad=.2), zorder=8)
        else:
            text(ax, *pts[0], 'N', 5.2, color,
                 bbox=dict(fc='white', ec='none', pad=0), zorder=8)
    # Join neighboring pyrroles with meso carbons and phenyl arms.
    for theta in np.arange(4) * np.pi / 2:
        u = np.array([np.cos(theta), np.sin(theta)])
        v = np.array([-np.sin(theta), np.cos(theta)])
        meso = c + .95 * r * u
        line(ax, meso, c + r * (.835 * u + .395 * v), color)
        line(ax, meso, c + r * (.835 * u - .395 * v), color)
        phc = c + 1.58 * r * u
        phr = .32 * r
        line(ax, meso, phc - phr * u, color)
        benzene(ax, phc, phr, theta, color, 1.05 if enlarged else .8)
        line(ax, phc + phr * u, c + 2.08 * r * u, color)
    return np.asarray(nitrogens)


def linker(ax, a, b):
    """Simplified phenyl–N=CH–(2,5-dihydroxyphenyl)–CH=N–phenyl linker."""
    a, b = np.asarray(a), np.asarray(b)
    u = (b - a) / np.linalg.norm(b - a)
    v = np.array([-u[1], u[0]])
    mid = (a + b) / 2
    r = .15
    theta = np.arctan2(u[1], u[0])
    ring = benzene(ax, mid, r, theta, C['blue'], .9)
    for p, q in ((a, mid - r * u), (mid + r * u, b)):
        line(ax, p, q, C['blue'], .95)
        # A second short stroke conveys the imine double bond.
        p1, q1 = p + .22 * (q - p) + .035 * v, p + .65 * (q - p) + .035 * v
        line(ax, p1, q1, C['blue'], .65)
        nc = p + .12 * (q - p)
        text(ax, *nc, 'N', 5.3, C['blue'],
             bbox=dict(fc='white', ec='none', pad=.1), zorder=7)
    # Two para-related hydroxy substituents on the dialdehyde-derived ring.
    for idx, sign in ((1, 1), (4, -1)):
        p = ring[idx]
        q = p + sign * .10 * v
        line(ax, p, q, C['blue'], .7)
        text(ax, *(q + sign * .06 * v), 'OH', 5.2, C['blue'])


def cof(ax):
    centers = [(1.03, 2.05), (2.77, 2.05), (1.03, 3.79), (2.77, 3.79)]
    r = .285
    for c in centers:
        ns = porphyrin(ax, c, r)
        for n in ns:
            line(ax, c, n, C['blue'], .65, zorder=5)
        ax.add_patch(Circle(c, .078, fc=C['orange'], ec='white', lw=.6, zorder=9))
        text(ax, *c, 'Fe', 5.5, 'white', 'bold', zorder=10)
    edge = 2.08 * r
    for a, b in ((centers[0], centers[1]), (centers[2], centers[3]),
                 (centers[0], centers[2]), (centers[1], centers[3])):
        a, b = np.array(a), np.array(b)
        u = (b - a) / np.linalg.norm(b - a)
        linker(ax, a + edge * u, b - edge * u)
    text(ax, 1.90, 1.18, 'Porphyrin–imine COF', 12, C['blue'], 'bold')
    text(ax, 1.90, .86, 'Near-equivalent Fe–N sites', 9.8, C['grey'])


def molecule(ax, x, y, co2=True):
    positions = [(x - .44, y), (x, y), (x + .44, y)] if co2 else [(x, y), (x + .44, y)]
    if co2:
        for d in (-1, 1):
            for dy in (-.034, .034):
                line(ax, (x + d * .13, y + dy), (x + d * .30, y + dy), C['carbon'], 1.7)
        kinds = [('O', C['red']), ('C', C['carbon']), ('O', C['red'])]
    else:
        for dy in (-.052, 0, .052):
            line(ax, (x + .13, y + dy), (x + .30, y + dy), C['carbon'], 1.3)
        kinds = [('C', C['carbon']), ('O', C['red'])]
    for p, (label, color) in zip(positions, kinds):
        ax.add_patch(Circle(p, .15, fc=color, ec='white', lw=.8, zorder=4))
        text(ax, *p, label, 10, 'white', 'bold', zorder=5)


def lightning_path(a, b, rng, roughness=.12, levels=5):
    """Seeded midpoint displacement; true vector polylines, no image filter."""
    pts = np.array([a, b], dtype=float)
    for depth in range(levels):
        result = [pts[0]]
        for p, q in zip(pts[:-1], pts[1:]):
            d = q - p
            normal = np.array([-d[1], d[0]])
            normal /= max(np.linalg.norm(normal), 1e-9)
            offset = rng.normal(0, roughness / (2 ** (.72 * depth)))
            result.extend([(p + q) / 2 + offset * normal, q])
        pts = np.array(result)
    return pts


def plasma(ax, center, width, height, electrodes=True, seed=41):
    """Purple DBD discharge illustration with branching luminous filaments.

    This is a conceptual plasma icon, not a reconstruction of the apparatus.
    Transparent ellipses and nested line strokes remain editable SVG vectors.
    No metal particles are introduced into the COF scientific story.
    """
    rng = np.random.default_rng(seed)
    cx, cy = center
    scale = min(width, height)
    # Smooth visual glow built from translucent nested vector ellipses.
    for i in range(30):
        f = 1 - .69 * i / 29
        ax.add_patch(Ellipse((cx, cy), width * f, height * f,
                            fc='#9E42DD', ec='none', alpha=.028, zorder=1))
    for _ in range(17):
        px = cx + rng.uniform(-.30, .30) * width
        py = cy + rng.uniform(-.32, .32) * height
        ax.add_patch(Ellipse((px, py), width * .25, height * .22,
                            angle=rng.uniform(-35, 35), fc='#CF78F1',
                            ec='none', alpha=.07, zorder=1.1))

    def filament(points, strength=1.):
        for lw, color, alpha in ((6.4, '#7D29BA', .09),
                                 (3.3, '#B747EC', .19),
                                 (1.45, '#ECA4FF', .70),
                                 (.55, '#FFFFFF', .97)):
            ax.plot(points[:, 0], points[:, 1], lw=lw * strength,
                    color=color, alpha=alpha, solid_capstyle='round',
                    solid_joinstyle='round', zorder=2)

    # Filaments originate at opposing electrode edges and branch inward.
    for side in (-1, 1):
        for j in range(6):
            a = np.array([cx + side * width * .39,
                          cy + height * (.13 * (j - 2.5))])
            b = np.array([cx + rng.uniform(-.11, .11) * width,
                          cy + rng.uniform(-.28, .28) * height])
            pts = lightning_path(a, b, rng, roughness=scale * .073)
            filament(pts)
            k = int(rng.integers(8, len(pts) - 4))
            endpoint = pts[k] + np.array([rng.uniform(-.18, .18) * width,
                                         rng.choice([-1, 1]) * height * .16])
            branch = lightning_path(pts[k], endpoint, rng,
                                    roughness=scale * .035, levels=3)
            filament(branch, .62)
    # Small ionization glints contribute detail without photographic rasters.
    for _ in range(20):
        x = cx + rng.uniform(-.33, .33) * width
        y = cy + rng.uniform(-.31, .31) * height
        r = scale * rng.uniform(.003, .009)
        ax.add_patch(Circle((x, y), r, fc='white', ec='none', alpha=.8, zorder=2.1))
    if electrodes:
        for side in (-1, 1):
            x = cx + side * width * .45
            # Stacked rectangles suggest a metallic electrode without raster shading.
            ax.add_patch(FancyBboxPatch((x - width * .060, cy - height * .23),
                                       width * .12, height * .46,
                                       boxstyle=f'round,pad=0,rounding_size={scale * .025}',
                                       fc='#292737', ec='#5D556D', lw=.8, zorder=3))
            for j in range(12):
                yy = cy - height * .205 + height * .41 * j / 12
                line(ax, (x - width * .048, yy),
                     (x + width * .048, yy), '#97909F', .8,
                     alpha=.12 + .24 * np.sin(np.pi * j / 12), zorder=3.1)
            lip = x - side * width * .054
            line(ax, (lip, cy - height * .21), (lip, cy + height * .21),
                 '#E777E8', 2.2, zorder=3.2)
            line(ax, (lip, cy - height * .17), (lip, cy + height * .17),
                 '#FFE9FF', .6, zorder=3.3)


def export_plasma_assets():
    for name, electrodes in (('plasma_background', False), ('plasma_discharge', True)):
        fig = plt.figure(figsize=(8, 4), facecolor='none')
        ax = fig.add_axes([0, 0, 1, 1], xlim=(0, 8), ylim=(0, 4))
        ax.set_aspect('equal'); ax.axis('off')
        plasma(ax, (4, 2), 7.2, 3.65, electrodes=electrodes)
        for ext in ('png', 'svg'):
            fig.savefig(OUT / f'{name}.{ext}', dpi=300, transparent=True)
        plt.close(fig)


def build(font):
    plt.rcParams.update({'font.family': font, 'font.size': 11,
                         'svg.fonttype': 'path', 'pdf.fonttype': 42,
                         'mathtext.fontset': 'custom', 'mathtext.rm': font,
                         'mathtext.it': f'{font}:italic',
                         'mathtext.bf': f'{font}:bold',
                         'mathtext.cal': font, 'mathtext.sf': font,
                         'mathtext.tt': font,
                         'axes.unicode_minus': False})
    fig = plt.figure(figsize=(12, 4.8), facecolor='white')
    ax = fig.add_axes([0, 0, 1, 1], xlim=(0, 12), ylim=(0, 4.8))
    ax.set_aspect('equal'); ax.axis('off')
    text(ax, 6, 4.52, 'Plasma-distorted FeN₄ sites for CO₂-to-CO photoreduction',
         17, C['ink'], 'bold')
    cof(ax)
    # Plasma treatment is separate from the light-driven catalytic operation.
    text(ax, 4.10, 4.00, 'DBD plasma', 11, C['purple'], 'bold')
    plasma(ax, (4.10, 3.24), 1.13, 1.21)
    arrow(ax, (3.49, 2.36), (4.72, 2.36), C['purple'], 2.4, 18)
    text(ax, 4.10, 2.06, 'Local distortion', 9.6, C['purple'])
    # Expanded site: four N donors remain; two unequal Fe–N bond sets.
    center = (6.13, 2.83)
    ax.add_patch(Circle(center, 1.28, fc=C['lilac'], ec='none', zorder=0))
    ns = porphyrin(ax, center, .59, C['purple'], enlarged=True)
    fe = np.array(center) + np.array([.052, 0])
    for n in ns:
        short = n[0] > center[0]
        line(ax, fe, n, C['orange'] if short else C['purple'], 2.7, zorder=6)
    ax.add_patch(Circle(fe, .145, fc=C['orange'], ec='white', lw=1.1, zorder=9))
    text(ax, *fe, 'Fe', 11, 'white', 'bold', zorder=10)
    text(ax, 6.13, 1.43, 'Asym-FeN₄-COF', 13, C['purple'], 'bold')
    line(ax, (5.11, 1.06), (5.41, 1.06), C['orange'], 2.7)
    text(ax, 5.87, 1.06, f'{FE_N_SHORT:.3f} Å', 10, C['grey'])
    line(ax, (6.49, 1.06), (6.79, 1.06), C['purple'], 2.7)
    text(ax, 7.22, 1.06, f'{FE_N_LONG:.3f} Å', 10, C['grey'])
    text(ax, 6.13, .66, 'FeN₄ motif retained', 10, C['grey'])
    # Visible photons absorbed by the conjugated framework, electron transfer
    # toward the catalytic site. Direction is conceptual, not a measured trajectory.
    sx, sy = 7.85, 3.92
    ax.add_patch(Circle((sx, sy), .20, fc='#F6C45D', ec=C['orange'], lw=.8))
    for ang in np.arange(8) * np.pi / 4:
        u = np.array([np.cos(ang), np.sin(ang)])
        line(ax, np.array([sx, sy]) + .27 * u,
             np.array([sx, sy]) + .37 * u, C['orange'], 1)
    text(ax, 9.12, 3.96, 'Visible light', 10.8, C['orange'], 'bold')
    t = np.linspace(0, 1, 90)
    a, b = np.array([7.65, 3.70]), np.array([7.02, 3.23])
    u = b - a; normal = np.array([-u[1], u[0]]) / np.linalg.norm(u)
    pts = a + t[:, None] * u + .05 * np.sin(7 * np.pi * t)[:, None] * normal
    ax.plot(pts[:, 0], pts[:, 1], color=C['orange'], lw=1.5)
    arrow(ax, tuple(pts[-4]), tuple(b - .10 * normal), C['orange'], 1, 10)
    arrow(ax, (7.38, 2.79), (6.38, 2.87), C['blue'], 1.9, 14, rad=.25)
    text(ax, 7.31, 2.38, 'e⁻ transfer', 10, C['blue'], 'bold')
    arrow(ax, (7.43, 3.14), (8.15, 3.14), C['green'], 1.8, 14)
    # Photocatalytic CO2-to-CO conversion and a conventional proposed pathway.
    molecule(ax, 8.94, 3.16, True)
    molecule(ax, 10.83, 3.16, False)
    arrow(ax, (9.56, 3.16), (10.59, 3.16), C['green'], 2.2, 16)
    text(ax, 10.07, 3.53, '2e⁻ + 2H⁺', 9.5, C['green'])
    text(ax, 8.94, 2.76, 'CO₂', 14, C['ink'], 'bold')
    text(ax, 11.05, 2.76, 'CO', 14, C['green'], 'bold')
    text(ax, 10.05, 2.35, 'CO₂  →  *COOH  →  *CO  →  CO', 10.1, C['grey'])
    text(ax, 10.05, 2.09, 'Proposed surface pathway', 8.6, C['grey'])
    # Results are mass-normalized production rate and CO selectivity.
    ax.add_patch(FancyBboxPatch((8.17, .55), 3.48, 1.25,
                               boxstyle='round,pad=0.08,rounding_size=0.12',
                               fc='#F0F7F5', ec='none'))
    text(ax, 9.18, 1.37, f'{RATE:.2f}', 27, C['green'], 'bold')
    text(ax, 9.18, .90, 'μmol g⁻¹ h⁻¹', 10.5, C['green'])
    line(ax, (10.12, .78), (10.12, 1.56), '#C7DED8', 1)
    text(ax, 10.93, 1.37, f'{SELECTIVITY}%', 27, C['green'], 'bold')
    text(ax, 10.93, .90, 'CO selectivity', 10.5, C['green'])
    # Slow TA component is not equated to the average TRPL lifetime.
    text(ax, 4.78, .22, f'Slow carrier decay: {TA_BEFORE:.2f} → {TA_AFTER:.2f} ps (≈16×)',
         9.6, C['blue'])
    # No raster artist is used; svg text paths preserve the chosen glyphs.
    for ext in ('png', 'svg'):
        fig.savefig(OUT / f'toc_figure.{ext}', dpi=300, facecolor='white',
                    metadata={'Title': 'Plasma-distorted FeN4 COF for CO2 reduction',
                              'Description': 'Schematic based on the uploaded manuscript; surface pathway proposed.'}
                    if ext == 'svg' else None)
    plt.close(fig)
    export_plasma_assets()
    # Keep the earlier filenames synchronized for existing download links.
    for ext in ('png', 'svg'):
        (OUT / f'TOC_FeN4_COF.{ext}').write_bytes((OUT / f'toc_figure.{ext}').read_bytes())
    return {'font': font, 'png_pixels': [3600, 1440], 'dpi': 300,
            'co_rate_umol_g_h': RATE, 'co_selectivity_percent': SELECTIVITY,
            'fe_n_distances_A': [FE_N_SHORT, FE_N_LONG],
            'slow_TA_decay_ps': [TA_BEFORE, TA_AFTER],
            'svg_text': 'glyph paths', 'surface_pathway': 'proposed schematic',
            'plasma_style': 'seeded branching vector discharge with translucent glow',
            'plasma_assets': ['plasma_background.png', 'plasma_background.svg',
                              'plasma_discharge.png', 'plasma_discharge.svg']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--font', default='Arial', help='Exact installed family; no silent substitution')
    parser.add_argument('--font-dir', type=Path, help='Directory containing legally obtained TTF/OTF fonts')
    args = parser.parse_args()
    if args.font_dir:
        for p in args.font_dir.rglob('*'):
            if p.suffix.lower() in ('.ttf', '.otf'):
                fm.fontManager.addfont(str(p))
    try:
        font_file = fm.findfont(fm.FontProperties(family=args.font), fallback_to_default=False)
    except ValueError:
        parser.error(f'{args.font} is unavailable. Supply --font-dir or explicitly select --font "Liberation Sans".')
    report = build(args.font)
    report['font_file'] = font_file
    (OUT / 'TOC_metadata.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
