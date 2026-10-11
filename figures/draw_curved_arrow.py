#!/usr/bin/env python3
"""S-shaped tapered arrow, inspired by the user-provided reference.
Run: python figures/draw_curved_arrow.py
Requires numpy and matplotlib; outputs transparent PNG (300 dpi) and SVG.
The broad tail narrows at the inflection and widens toward the arrowhead.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.colors import to_rgb

OUT = Path(__file__).resolve().parent


def bezier(points, t):
    p = np.asarray(points, dtype=float)
    q = t[:, None]
    xy = (1-q)**3*p[0] + 3*(1-q)**2*q*p[1] + 3*(1-q)*q**2*p[2] + q**3*p[3]
    d = 3*(1-q)**2*(p[1]-p[0]) + 6*(1-q)*q*(p[2]-p[1]) + 3*q**2*(p[3]-p[2])
    return xy, d


def draw():
    t = np.linspace(0, 1, 601)
    # Two cubic curves reproduce the left reference's nearly horizontal
    # endpoints and steep central bend; equal derivatives keep the join smooth.
    upper, du = bezier([(0.65,3.32),(1.75,3.32),(2.92,2.74),(2.92,2.12)], t)
    lower, dl = bezier([(2.92,2.12),(2.92,1.50),(4.10,.66),(5.22,.66)], t)
    xy = np.vstack([upper,lower[1:]])
    d = np.vstack([du,dl[1:]])
    tangents = d / np.linalg.norm(d,axis=1)[:,None]
    normals = np.c_[-tangents[:,1],tangents[:,0]]
    # The right reference has a broad end, a narrow waist, and a broad
    # arrow end. Preserve that width trend along the S-shaped centerline.
    s = np.linspace(0,1,len(xy))
    width = np.interp(s,[0,.30,.50,.70,1],[.28,.18,.025,.15,.27])
    left = xy + normals*width[:,None]/2
    right = xy - normals*width[:,None]/2
    fig = plt.figure(figsize=(6.4,4), facecolor='none')
    ax = fig.add_axes([0,0,1,1],xlim=(0,6.4),ylim=(0,4))
    ax.set_aspect('equal'); ax.axis('off')
    start = np.array(to_rgb('#BD8564'))
    end = np.array(to_rgb('#FFB185'))
    # Overlapping narrow quadrilaterals implement a vector color gradient.
    for i in range(len(xy)-1):
        color = start*(1-s[i]) + end*s[i]
        polygon = np.array([left[i],left[min(i+2,len(xy)-1)],
                            right[min(i+2,len(xy)-1)],right[i]])
        ax.add_patch(Polygon(polygon,closed=True,fc=color,ec='none',antialiased=False))
    # Horizontal arrowhead aligned to the end tangent, broad like the reference.
    tip = xy[-1] + np.array([.47,0])
    head = [tip, xy[-1]+[-.055,.41], xy[-1]+[-.055,.135],
            xy[-1]+[-.12,.135], xy[-1]+[-.12,-.135],
            xy[-1]+[-.055,-.135], xy[-1]+[-.055,-.41]]
    ax.add_patch(Polygon(head,closed=True,fc=end,ec='none',joinstyle='miter'))
    for ext in ('png','svg'):
        fig.savefig(OUT/f'curved_tapered_arrow.{ext}',dpi=300,transparent=True)
    plt.close(fig)
    print('Saved curved_tapered_arrow.png (1920 × 1200, 300 dpi) and pure-vector SVG.')


if __name__ == '__main__':
    draw()
