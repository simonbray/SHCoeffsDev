#!/usr/bin/env python
import sys

import numpy as np

from math import sin, cos
from mayavi import mlab


def rtp2xyz(vec):
    x = vec[0] * sin(vec[1]) * cos(vec[2])
    y = vec[0] * sin(vec[1]) * sin(vec[2])
    z = vec[0] * cos(vec[1])

    return [x, y, z]


def main():
    # Create and open an application window.
    #mlab.figure(1, fgcolor=(1, 1, 1), bgcolor=(0, 0, 0), size=(1600, 1600))
    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(1600, 1600))


    x = []
    y = []
    z = []
    t = []
    with open("reconstruction.csv") as in_file:
        for line in in_file:
            coords = np.asarray(line.split(","), dtype=np.float64, order='C')
            xyz = rtp2xyz([coords[3], coords[0], coords[1]])
            x.append(xyz[0])
            y.append(xyz[1])
            z.append(xyz[2])
            t.append(coords[5])

    max_abs_charge = np.max(np.abs(t))


    surface = mlab.points3d(x, y, z, t, mode='sphere',
                            colormap='RdYlBu',
                            scale_factor=0.05,
                            scale_mode='none',
                            vmin=-max_abs_charge,
                            vmax=max_abs_charge)

    mlab.quiver3d([0.0, 0.0, 0.0], [0.0, 0.0, 0.0],[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0])

    surface.module_manager.scalar_lut_manager.reverse_lut = False

    mlab.scalarbar(surface, orientation='vertical', nb_labels=3, title='kcal/mol')

    mlab.draw()
    mlab.show()


if __name__ == "__main__":
    main()
