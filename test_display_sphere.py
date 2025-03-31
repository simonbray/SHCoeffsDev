#!/usr/bin/env python
import ast
import sys
import math

import numpy as np

from mayavi import mlab


def xyz2rtp(x, y, z):
    r = math.sqrt(x**2+y**2+z**2)
    theta = math.acos(z/r)
    phi = math.atan2(y, x)
    return (r, theta, phi)


def rtp2xyz(r, theta, phi):
    x = r*math.sin(theta)*math.cos(phi)
    y = r*math.sin(theta)*math.sin(phi)
    z = r*math.cos(theta)
    return (x, y, z)


def main():
    file_name = "unit_sphere.csv"

    # Create and open an application window.
    mlab.figure(1, bgcolor=(0, 0, 0), fgcolor=(1, 1, 1), size=(1600, 1600))

    x = []
    y = []
    z = []
    with open(file_name, "r") as in_file:
        for line in in_file:
            tokens = list(map(float, line.split(",")))

            a, b, c = rtp2xyz(tokens[0], tokens[1], tokens[2])
            print(a, b, c)
            x.append(a)
            y.append(b)
            z.append(c)
    
    surface = mlab.points3d(x, y, z, mode='point')

    mlab.draw()
    mlab.show()


if __name__ == "__main__":
    main()