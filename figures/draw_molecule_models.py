#!/usr/bin/env python3
"""Render separate vector ball-and-stick CO2 and CO illustrations.
Run with Python, numpy and matplotlib: python figures/draw_molecule_models.py
Carbon: charcoal; oxygen: red. No labels, transparent canvas, PNG at 300 dpi.
CO2 is linear (180 degrees); CO is diatomic. Bond lengths set center spacing
in angstrom units (CO2: 1.16 Å; CO: 1.13 Å); sphere radii are illustrative,
not atomic radii. One cylindrical stick represents each chemical bond;
its bond order is not encoded by the number of sticks.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, Ellipse
from matplotlib.colors import to_rgb

OUT = Path(__file__).resolve().parent


def sphere(ax, x, y, radius, color):
    base = np.array(to_rgb(color))
    # Overlapping vector disks create highlights without embedded raster images.
    ax.add_patch(Circle((x, y), radius, fc=base * .48, ec='none', zorder=3))
    for i in range(85):
        t = (i + 1) / 86
        r = radius * (1 - .84 * t)
        shift = radius * t * .25
        col = base * (.54 + .44 * t) + (1 - base) * .38 * t ** 2
        ax.add_patch(Circle((x - shift, y + shift), r,
                            fc=np.clip(col, 0, 1), ec='none', zorder=3 + t))
    ax.add_patch(Ellipse((x - radius * .31, y + radius * .39),
                         radius * .31, radius * .17, angle=28,
                         fc='white', alpha=.35, ec='none', zorder=5))


def stick(ax, x1, x2, radius=.095):
    # Vertical brightness profile suggests a cylinder; all bands are vectors.
    bands = 64
    for j in range(bands):
        t = -1 + 2 * j / bands
        brightness = .35 + .50 * max(0, 1 - (t - .22) ** 2)
        ax.add_patch(Rectangle((x1, t * radius), x2 - x1,
                               2 * radius / bands + .0003,
                               fc=(brightness * .96, brightness * .98, brightness),
                               ec='none', zorder=1))


def render(name, atoms):
    fig = plt.figure(figsize=(6, 2.5), facecolor='none')
    ax = fig.add_axes([0, 0, 1, 1], xlim=(-1.95, 1.95), ylim=(-.8125, .8125))
    ax.set_aspect('equal'); ax.axis('off')
    for a, b in zip(atoms[:-1], atoms[1:]):
        stick(ax, a[1], b[1])
    for element, x in atoms:
        sphere(ax, x, 0, .43 if element == 'O' else .39,
               '#D33332' if element == 'O' else '#4E5660')
    for ext in ('png', 'svg'):
        fig.savefig(OUT / f'{name}_ball_and_stick.{ext}', dpi=300,
                    transparent=True,
                    metadata={'Title': f'{name} ball-and-stick schematic'} if ext == 'svg' else None)
    plt.close(fig)


def main():
    render('CO2', [('O', -1.16), ('C', 0), ('O', 1.16)])
    render('CO', [('C', -.565), ('O', .565)])
    print('Generated CO2 and CO: transparent 1800 × 750 PNG (300 dpi) and vector SVG.')


if __name__ == '__main__':
    main()
