# -*- coding: utf-8 -*-
# Copyright (c) 2011 Oliver Beckstein
# Placed into the Public Domain
# NO WARRANTY

"""
Load HOLLOW constraints file into PyMOL and show the box

In PyMOL:

 PyMOL> run /path/to/drawBrick.py
 PyMOL> drawHollowConstraints hollow.constraints

This will draw the constraints region. Modify the constraints file,
call drawHollowConstraints again and cycle through the frames of the
object to see how the box changes.

At the moment, only the *brick* constraints are implemented.
"""


# 1. Load data from constraints file
# 2. Translate into PyMOL selections
# 3. Have PyMOL get the coordinates
# 4. Set up geometry
# 5. Draw in PyMOL.


from pymol.cgo import *
from pymol import cmd
from random import randint

import numpy

def read_parameters(fname):
    # hollow.util.read_parameters
    class DataHolder:
        pass
    f = open(fname, 'r')
    try:
        result = DataHolder()
        result.__dict__ = eval(f.read())
    finally:
        f.close()
    return result

def str2array(s):
    seq = str(s).replace('[','').replace(']','').replace('(','').replace(')','').split(",")
    return numpy.array([float(x.strip()) for x in seq])


# based on http://www.pymolwiki.org/index.php/DrawBoundingBox

def box(lower_left, upper_right, linewidth=2, rgb=None):
    minX, minY, minZ = lower_left
    maxX, maxY, maxZ = upper_right
    if rgb is None:
        r,g,b = (1., 1., 1.)
    else:
        r,g,b = str2array(rgb)
    return [
        LINEWIDTH, float(linewidth),

        BEGIN, LINES,
        COLOR, r, g, b,

        VERTEX, minX, minY, minZ,       #1
        VERTEX, minX, minY, maxZ,       #2

        VERTEX, minX, maxY, minZ,       #3
        VERTEX, minX, maxY, maxZ,       #4

        VERTEX, maxX, minY, minZ,       #5
        VERTEX, maxX, minY, maxZ,       #6

        VERTEX, maxX, maxY, minZ,       #7
        VERTEX, maxX, maxY, maxZ,       #8


        VERTEX, minX, minY, minZ,       #1
        VERTEX, maxX, minY, minZ,       #5

        VERTEX, minX, maxY, minZ,       #3
        VERTEX, maxX, maxY, minZ,       #7

        VERTEX, minX, maxY, maxZ,       #4
        VERTEX, maxX, maxY, maxZ,       #8

        VERTEX, minX, minY, maxZ,       #2
        VERTEX, maxX, minY, maxZ,       #6


        VERTEX, minX, minY, minZ,       #1
        VERTEX, minX, maxY, minZ,       #3

        VERTEX, maxX, minY, minZ,       #5
        VERTEX, maxX, maxY, minZ,       #7

        VERTEX, minX, minY, maxZ,       #2
        VERTEX, minX, maxY, maxZ,       #4

        VERTEX, maxX, minY, maxZ,       #6
        VERTEX, maxX, maxY, maxZ,       #8

        END
        ]

def get_selection_first_coord(sel, name="(sel)"):
    """Return *first* coordinate of selection"""

    #print "--> selection %(name)s:  %(sel)s" % vars()
    cmd.select(name, sel)
    try:
        return numpy.array(cmd.get_model(name).atom[0].coord)
    except IndexError:
        print "--> empty selection: %(sel)s" % vars()
        raise

def drawBrick(res_num1, atom1, res_num2, atom2, **kwargs):
    """Draw a box with loer left corner and upper right corner specified by atoms.

    drawBrick res_num1,atom1, res_num2,atom2[, chain1[,chain2[,offset1[,offset2[,boxname]]]]]

    ARGUMENTS

    res_num1, atom1 - lower left
    res_num2, atom2 - upper right

    OPTIONAL

    chain1, chain2  - defaults to "A"

    offset1         - add to lower left front (0,0,0)
    offset2         - add to upper right back; default (0,0,0)

    rgb             - line colour, e.g (1,1,1) for white
    linewidth       - width of box lines
    """

    chain1 = str(kwargs.pop("chain1", "A")).strip()
    chain2 = str(kwargs.pop("chain2", "A")).strip()
    res_num1 = int(res_num1)
    res_num2 = int(res_num2)
    offset1 = str2array(kwargs.pop("offset1", "(0,0,0)"))
    offset2 = str2array(kwargs.pop("offset2", "(0,0,0)"))

    boxname = kwargs.pop("boxname", "brick")
    selname = kwargs.pop("name", "(sel)")
    p1 = get_selection_first_coord("///%(chain1)s/`%(res_num1)d/%(atom1)s" % vars(), selname)
    p2 = get_selection_first_coord("///%(chain2)s/`%(res_num2)d/%(atom2)s" % vars(), selname)

    #print "lower_left: %r" % lower_left
    #print "offset1:    %r" % offset1

    p1 += offset1
    p2 += offset2

    lower_left = numpy.where(p1 < p2, p1, p2)
    upper_right = numpy.where(p1 < p2, p2, p1)

    print "--> lower left front corner: %(lower_left)r" % vars()
    print "--> upper right back corner: %(upper_right)r" % vars()

    boxCGO = box(lower_left, upper_right, **kwargs)
    cmd.load_cgo(boxCGO, boxname)
    return boxname

def drawHollowConstraints(fname, linewidth=2, rgb="(1,1,1)"):
    """Draw the constraints volume from a Hollow input file.

    Only brick/box constraints implemented at the moment.

    rgb             - line colour, e.g (1,1,1) for white
    linewidth       - width of box lines
    """
    constraints = read_parameters(fname)
    if not constraints.type in ("box", "brick"):
        raise NotImplementedError("sorry, constraint %(type)r not implemented" % vars(constraints))
    c = constraints
    return drawBrick(c.res_num1, c.atom1, c.res_num2, c.atom2,
                     chain1=c.chain1, chain2=c.chain2,
                     offset1=c.offset1, offset2=c.offset2,
                     boxname=c.type,
                     linewidth=linewidth, rgb=rgb)

cmd.extend("drawBrick",  drawBrick)
cmd.extend("drawHollowConstraints", drawHollowConstraints)
