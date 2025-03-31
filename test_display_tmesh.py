#!/usr/bin/env python
import sys

import numpy as np

from mayavi import mlab


def main():
    # Create and open an application window.
    #mlab.figure(1, fgcolor=(1, 1, 1), bgcolor=(0, 0, 0), size=(1600, 1600))
    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(1600, 1600))

    file_name = "ZINC67007472_origin.tmesh"

    with open(file_name) as in_file:
        x = []
        y = []
        z = []
        t = []
        lines = in_file.readlines()
        for line in lines:
            tokens = line.split()
            if len(tokens) == 9:
                coords = np.asarray(tokens[:3], dtype=np.float64, order='C')
                x.append(coords[0])
                y.append(coords[1])
                z.append(coords[2])
                t.append(float(tokens[7]))
        triangles = []
        for i in range(len(lines)):
            line = lines[i]
            if line[:3] == '3 0':
                triangles.append((int(lines[i+1]), int(lines[i + 2]), int(lines[i + 3])))

    max_abs_charge = np.max(np.abs(t))

    #surface = mlab.triangular_mesh(x, y, z, triangles, scalars=t, colormap='RdYlBu', opacity=0.8, representation='fancymesh')
    surface = mlab.triangular_mesh(x, y, z, triangles, scalars=t, colormap='RdYlBu',
                                   representation='surface',
                                   opacity=0.5,
                                   resolution=64, 
                                   vmin=-max_abs_charge,
                                   vmax=max_abs_charge)

    surface.module_manager.scalar_lut_manager.reverse_lut = False

    mlab.scalarbar(surface, orientation='vertical', nb_labels=3, title='kcal/mol')



    x = []
    y = []
    z = []
    t = []
    with open("unit_intersect.csv") as in_file:
        for line in in_file:
            coords = np.asarray(line.split(","), dtype=np.float64, order='C')
            
            x.append(coords[0])
            y.append(coords[1])
            z.append(coords[2])
            t.append(coords[3])

    mlab.points3d(x, y, z, t, mode='sphere',
                  colormap='RdYlBu',
                  scale_factor=0.05,
                  scale_mode='none',
                  vmin=-max_abs_charge,
                  vmax=max_abs_charge)

    mlab.quiver3d([0.0, 0.0, 0.0], [0.0, 0.0, 0.0],[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0])

    mlab.draw()
    mlab.show()


if __name__ == "__main__":
    main()
