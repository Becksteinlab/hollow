# Hollow --- http://hollow.sourceforge.net
# Copyright (c) 2009 Bosco Ho & Franz Gruswitz
# Note: This used to be hollow.py in the original Hollow 1.1 distribution.
#       In order to Hollow conform to standard Python distributions, the core
#       functionality was moved into the hollow.core module and the hollow
#       script is only a thin wrapper around make_hollow_spheres(). -- OB 2010
"""
hollow.core module
==================

This module containts code to generate the lattics of spheres that
fills voids in protein cavities. Users really only need
:func:`make_hollow_spheres`, which can imported in user code as ::

  from hollow.core import make_hollow_spheres

(see the ``hollow`` script as an example).

"""

import math
import os
import array
import sys

# should really make it use numpy throughout... [orbeckst]
import numpy

# hollow modules
import util, pdbtext, vector3d, pdbstruct, asa

from pkg_resources import resource_filename

# load default values from "hollow.txt"
default_fname = resource_filename(__name__, "hollow.txt")
defaults = util.read_parameters(default_fname)


class BaseGrid:
  def __init__(self, n):
    self.n = n
    self.n_sq = self.n * self.n
    self.array = array.array('b')
    self.n_cube = self.n_sq * self.n
    for i in xrange(self.n_cube):
      self.array.append(0)

  def is_set(self, i, j, k):
    l = i*self.n_sq + j*self.n + k
    return self.array[l] > 0

  def set(self, i, j , k, is_state):
    l = i*self.n_sq + j*self.n + k
    if is_state:
      self.array[l] = 1
    else:
      self.array[l] = 0


class Grid:
  def __init__(self, grid_spacing, width, center):
    self.center = center
    self.width = width
    self.spacing = grid_spacing
    self.inv_spacing = 1.0 / self.spacing

    self.n = int(math.ceil(self.width * self.inv_spacing))
    self.half_n = self.n // 2

    self.excluded = BaseGrid(self.n)
    self.is_excluded = self.excluded.is_set
    self.set_excluded = self.excluded.set

    self.drilled = BaseGrid(self.n)
    self.is_drilled = self.drilled.is_set
    self.set_drilled = self.drilled.set

    self.x = [self.center.x + (i-self.half_n)*self.spacing
              for i in xrange(self.n)]
    self.y = [self.center.y + (i-self.half_n)*self.spacing
              for i in xrange(self.n)]
    self.z = [self.center.z + (i-self.half_n)*self.spacing
              for i in xrange(self.n)]

  def is_excluded_or_drilled(self, i, j, k):
    return self.is_excluded(i, j, k) or self.is_drilled(i, j, k)

  def indices(self, pos):
    return \
       ((pos.x-self.center.x)*self.inv_spacing + self.half_n,
        (pos.y-self.center.y)*self.inv_spacing + self.half_n,
        (pos.z-self.center.z)*self.inv_spacing + self.half_n)

  def pos(self, i, j, k):
    return vector3d.Vector3d(self.x[i], self.y[j], self.z[k])

  def is_grid_point_near_sphere(self, i, j, k, pos, r_sq):
    d_x = self.x[i] - pos.x
    d_y = self.y[j] - pos.y
    d_z = self.z[k] - pos.z
    return d_x*d_x + d_y*d_y + d_z*d_z < r_sq

  def int_range(self, low_f, high_f):
    low = max(0, int(math.floor(low_f)))
    high = min(self.n, int(math.ceil(high_f)) + 1)
    return xrange(low, high)

  def exclude_sphere(self, pos, r):
    r_sq = r*r
    low = vector3d.Vector3d(pos.x-r, pos.y-r, pos.z-r)
    high = vector3d.Vector3d(pos.x + r, pos.y + r, pos.z + r)
    low_i, low_j, low_k = self.indices(low)
    high_i, high_j, high_k = self.indices(high)
    for i in self.int_range(low_i, high_i):
      for j in self.int_range(low_j, high_j):
        for k in self.int_range(low_k, high_k):
          if not self.is_excluded(i, j, k):
            if self.is_grid_point_near_sphere(i, j, k, pos, r_sq):
              self.set_excluded(i, j, k, True)

  def permutation(self, i, j, k, dim):
    if dim == 0:
      return i, j, k
    if dim == 1:
      return j, k, i
    if dim == 2:
      return k, i, j

  def drill_in_dim(self, is_reversed, i, j, dim):
    drill_range = xrange(self.n)
    if is_reversed:
      drill_range = reversed(drill_range)
    for k in drill_range:
      a, b, c = self.permutation(i, j, k, dim)
      if self.is_excluded(a, b, c):
        return
      self.set_drilled(a, b, c, True)

  def exclude_edge_to_interior(self):
    for i in xrange(self.n):
      for j in xrange(self.n):
        self.drill_in_dim(True,  i, j, 0)
        self.drill_in_dim(False, i, j, 0)
        self.drill_in_dim(True,  i, j, 1)
        self.drill_in_dim(False, i, j, 1)
        self.drill_in_dim(True,  i, j, 2)
        self.drill_in_dim(False, i, j, 2)

  def is_surrounded(self, i, j, k):
    indices_list = [ \
        (i,j,k),
        (i+1,j,k), (i-1,j,k),
        (i,j+1,k), (i,j-1,k),
        (i,j,k-1), (i,j,k+1)]
    for (a, b, c) in indices_list:
      if 0 <= a < self.n and 0 <= b < self.n and 0 <= c < self.n:
        if self.is_excluded_or_drilled(a, b, c):
          return False
    return True

  def exclude_surrounded(self, skip):
    surrounded_grid_points = []
    for i in xrange(self.n):
      for j in xrange(self.n):
        for k in xrange(self.n):
          if self.is_surrounded(i,j,k):
            surrounded_grid_points.append([i,j,k])
    for (i,j,k) in surrounded_grid_points:
      if skip > 0:
        if i % skip == 0 and j % skip == 0 and k % skip == 0:
          continue
      self.set_excluded(i, j, k, True)

  def exclude_points_in_constraint(self, constraint_fn):
    for i in xrange(self.n):
      for j in xrange(self.n):
        for k in xrange(self.n):
          if not self.is_excluded_or_drilled(i, j, k):
            if not constraint_fn(self.pos(i,j,k)):
              self.set_excluded(i, j, k, True)

  def make_mol(self, res_type, atom_type):
    mol = pdbstruct.Molecule()
    element = ""
    for c in atom_type[:2]:
      if not c.isdigit() and c != " ":
        element += c
    i_res = 1
    for i in xrange(self.n):
      for j in xrange(self.n):
        for k in xrange(self.n):
          if not (self.is_excluded_or_drilled(i,j,k)):
            hollow_atom = pdbstruct.Atom()
            hollow_atom.pos = self.pos(i,j,k)
            hollow_atom.type = atom_type
            hollow_atom.element = element
            hollow_atom.res_type = res_type
            hollow_atom.res_num = i_res
            hollow_atom.num = i_res
            mol.insert_atom(hollow_atom)
            i_res += 1
    return mol


def exclude_atoms_from_grid(grid, atoms, probe):
  for atom in atoms:
    grid.exclude_sphere(atom.pos, atom.radius + probe)


def exclude_surface(grid, atoms, probe):
  test_point = vector3d.Vector3d()
  sphere_points = asa.generate_sphere_points(960)
  for i, atom_i in enumerate(atoms):
    is_interior_atom = "asa" in atom_i.__dict__ and atom_i.asa < 9
    if not is_interior_atom:
      neighbor_indices = asa.find_neighbor_indices(atoms, probe, i)
      n_neighbor = len(neighbor_indices)
      j_closest_neighbor = 0
      radius = probe + atom_i.radius
      for point in sphere_points:
        is_point_accessible = True
        test_point.x = point[0]*radius + atom_i.pos.x
        test_point.y = point[1]*radius + atom_i.pos.y
        test_point.z = point[2]*radius + atom_i.pos.z
        cycled_indices = range(j_closest_neighbor, n_neighbor)
        cycled_indices.extend(range(j_closest_neighbor))
        for j in cycled_indices:
          atom_j = atoms[neighbor_indices[j]]
          r = atom_j.radius + probe
          d_x = test_point.x - atom_j.pos.x
          d_y = test_point.y - atom_j.pos.y
          d_z = test_point.z - atom_j.pos.z
          if d_x*d_x + d_y*d_y + d_z*d_z < r*r:
              j_closest_neighbor = j
              is_point_accessible = False
              break
        if is_point_accessible:
          grid.exclude_sphere(test_point, probe)


def calculate_average_bfactor(
    grid_chain, protein_atoms, bfactor_probe):
  max_bfactor = 0.0
  for atom in protein_atoms:
    if atom.bfactor > max_bfactor:
      max_bfactor = atom.bfactor
  for grid_atom in grid_chain.atoms():
    bfactors = []
    for protein_atom in protein_atoms:
      if protein_atom.element != "H":
        radius = bfactor_probe
        dist = vector3d.pos_distance(protein_atom.pos, grid_atom.pos)
        if dist < radius:
          bfactors.append(protein_atom.bfactor)
    n_bfactor = len(bfactors)
    if n_bfactor == 0:
      grid_atom.bfactor = max_bfactor
    else:
      grid_atom.bfactor = sum(bfactors)/float(n_bfactor)


def add_asa(atoms, probe):
  areas = asa.calculate_asa(atoms, probe)
  for area, atom in zip(areas, atoms):
    atom.asa = area


def get_molecule(pdb):
  txt = open(pdb, 'r').read()
  txt = pdbtext.strip_other_nmr_models(txt)
  txt = pdbtext.strip_lines(txt, lambda l: l.startswith('ANISOU'))
  txt = pdbtext.strip_lines(txt, lambda l: l.startswith('CONECT'))
  txt = pdbtext.strip_lines(txt, lambda l: l.startswith('MASTER'))
  txt = pdbtext.strip_alternative_atoms(txt)
  bare_fname = util.temp_fname('.pdb')
  open(bare_fname, 'w').write(txt)
  mol = pdbstruct.Molecule(bare_fname)
  if os.path.isfile(bare_fname):
    os.remove(bare_fname)
  return mol


def find_atom(atoms, chain, res_num, atom_type):
  if chain == "":
    chain = " "
  for atom in atoms:
    if atom.chain_id == chain \
      and atom.res_num == res_num \
      and atom.type == atom_type:
      return atom
  raise ValueError("Can't find '%s' atom %s of res %d" \
          % (chain, atom_type, res_num))


def get_sphere_constraint_fn(center, radius):
  return lambda pos: vector3d.pos_distance(center, pos) <= radius


def get_cylinder_constraint_fn(center1, center2, radius):
  axis12 = center2 - center1

  def cylinder_constraint_fn(pos):
    pos1 = pos - center1
    if vector3d.dot(pos1, axis12) < 0:
      return False
    pos1_perp = pos1.perpendicular_vec(axis12)
    if pos1_perp.length() > radius:
      return False
    pos2 = pos - center2
    if vector3d.dot(pos2, axis12) > 0:
      return False
    return True

  return cylinder_constraint_fn

def get_brick_constraint_fn(lower_left_front, upper_right_back):
  """Constraint to be within in a orthorhombic box."""
  lo = lower_left_front
  hi = upper_right_back
  def brick_constraint_fn(pos):
    return lo.x < pos.x < hi.x and lo.y < pos.y < hi.y and lo.z < pos.z < hi.z
  return brick_constraint_fn

def print_time(timer):
  print "   >>>>>", timer.str().split()[1]


def make_hollow_spheres(
    pdb,
    out_pdb="",
    grid_spacing=defaults.grid_spacing,
    size_interior_probe=defaults.interior_probe,
    is_skip_waters=defaults.is_skip_waters,
    size_surface_probe=defaults.surface_probe,
    constraint_file="",
    is_hollow=defaults.is_hollow
    size_bfactor_probe=defaults.bfactor_probe):
  """Generate spheres that fill the void in structure *pdb*.

  make_hollow_spheres(pdb[,out_pdb,grid_spacing,size_interior_probe,is_skip_waters,size_surface_probe,constraint_file,size_bfactor_probe])

  Read the structure from PDB file *pdb* and write the fake atoms that show the voids
  to file *out_pdb*. If a *constraint_file* is provided, only regions within the
  defined spherical or cylindrical region are analyzed. All optional arguments have
  default values that are shown below in square brackets.

  :Arguments:
     pdb
        a PDB structure file (remove ligands)

  :Keywords:
     out_pdb
        name of the PDB file containing the spheres generated by Hollow;
        if left empty, then the file is named <pdb>-hollow.pdb
     grid_spacing
        spacing of the grid on which spheres are placed, in Angstroem. It determines
        how detailed the resultant fake atoms will be in the output PDB file. Values
        smaller than 0.2 A are typically not necessary and can lead to memory
        problems (especially when not used with a *constraint_file*). Values of 0.5 -
        1.0 A tend to produce good results in automatic mode (the default) when the
        whole protein is scanned [%(grid_spacing)f]
     size_interior_probe
        The definition of the accessible surface involves rolling a surface probe
        over the atomic surface defined by the van der Waals radius of the atoms. The
        *size_interior probe* defines the size of the surface probe. [%(interior_probe)f]
     is_skip_waters
         [%(is_skip_waters)s]
     size_surface_probe
         In the automatic analysis mode, we want to generate fake atoms in clefts. We
         start the calculation with a cubic grid of fake atoms and system atically
         eliminate fake atoms that like within the protein structure. However, we
         need to also eliminate fake atoms outside the accessible surface of the
         protein as well, whilst allowing fake atoms to fill in surface clefts. To do
         that, we need to define a wrapping surface around the protein structure,
         that wraps around the accessible surface but contains spaces corresponding
         to surface clefts. We do this by defining a large exterior surface probe
         through *size_surface_probe*   [%(surface_probe)f]
     constraint_file
        read constraints from file; this enables **constraint mode** (Without this
        file, Hollow operates in **automatic analysis mode**.) See the README file and
        the examples (channel and binding_site) for example files.
     size_bfactor_probe
        We also calculate appropriate B-factors of every fake atom, by averaging over
        the heavy protein atoms around each fake atom. [%(bfactor_probe)f]

  For further help, type ``hollow.help()`` and read the README file in the
  distribution or see the online documentation at http://hollow.sourceforge.net/.
  """

  # load protein and get protein parameters
  print "Loading", pdb
  mol = get_molecule(pdb)
  atoms = mol.atoms()
  if is_skip_waters:
    print "Skipping water molecules"
    atoms = [a for a in atoms if a.res_type != "HOH"]
  pdbstruct.add_radii(atoms)

  # setup constraints and grid size in width
  timer = util.Timer()
  if not constraint_file:
    center = pdbstruct.get_center(atoms)
    width = pdbstruct.get_width(atoms, center)
  else:
    print "Loading constraints from %s" % constraint_file
    constraints = util.read_parameters(constraint_file)
    if constraints.type == 'sphere':
      atom1 = find_atom(
          atoms, constraints.chain1, constraints.res_num1,
          constraints.atom1)
      radius = constraints.radius
      constraint_fn = get_sphere_constraint_fn(atom1.pos, radius)
      center = atom1.pos
      width = 2.0*constraints.radius + 2.0*grid_spacing
      radius = radius - 3.0*grid_spacing
      inner_constraint_fn = get_sphere_constraint_fn(atom1.pos, radius)
      print "%s constraints: center = %r, radius = %g" % \
          (constraints.type, center, constraints.radius)
    elif constraints.type == 'coordsphere':
      try:
        coord = vector3d.Vector3d(*constraints.coord)
      except AttributeError:
        coord = vector3d.Vector3d(0,0,0)
      radius = constraints.radius
      constraint_fn = get_sphere_constraint_fn(coord, radius)
      center = coord
      width = 2.0*constraints.radius + 2.0*grid_spacing
      radius = radius - 3.0*grid_spacing
      inner_constraint_fn = get_sphere_constraint_fn(coord, radius)
      print "%s constraints: center = %r, radius = %g" % \
          (constraints.type, center, constraints.radius)
    elif constraints.type == 'cylinder':
      atom1 = find_atom(
          atoms, constraints.chain1, constraints.res_num1,
          constraints.atom1)
      atom2 = find_atom(
          atoms, constraints.chain2, constraints.res_num2,
          constraints.atom2)
      axis12 = atom2.pos - atom1.pos

      offset1 = -axis12.normal_vec()
      offset1.scale(constraints.axis_offset1)
      center1 = atom1.pos + offset1

      offset2 = axis12.normal_vec()
      offset2.scale(constraints.axis_offset2)
      center2 = atom2.pos + offset2

      center = center1 + center2
      center.scale(0.5)
      radius = constraints.radius
      constraint_fn = get_cylinder_constraint_fn(
          center1, center2, radius)

      half_length = vector3d.pos_distance(center, center1)
      width = 2.0 * grid_spacing + 2.0 * \
                math.sqrt(half_length*half_length \
                  + constraints.radius*constraints.radius)
      border_length = 3.0*grid_spacing
      center1 = center1 + axis12.normal_vec().scaled_vec(border_length)
      center2 = center2 - axis12.normal_vec().scaled_vec(border_length)
      radius = radius - border_length
      inner_constraint_fn = get_cylinder_constraint_fn(
          center1, center2, radius)
      print "%s constraints: center = %r, axis = %r, radius = %g" % \
          (constraints.type, center, axis12, constraints.radius)
    elif constraints.type == 'coordcylinder':
      try:
        coord1 = vector3d.Vector3d(*constraints.coord1)
      except AttributeError:
        coord1 = vector3d.Vector3d(0,0,0)
      try:
        coord2 = vector3d.Vector3d(*constraints.coord2)
      except AttributeError:
        coord2 = vector3d.Vector3d(0,0,0)
      axis12 = coord2 - coord1

      offset1 = -axis12.normal_vec()
      offset1.scale(constraints.axis_offset1)
      center1 = coord1 + offset1

      offset2 = axis12.normal_vec()
      offset2.scale(constraints.axis_offset2)
      center2 = coord2 + offset2

      center = center1 + center2
      center.scale(0.5)
      radius = constraints.radius
      constraint_fn = get_cylinder_constraint_fn(
          center1, center2, radius)

      half_length = vector3d.pos_distance(center, center1)
      width = 2.0 * grid_spacing + 2.0 * \
                math.sqrt(half_length*half_length \
                  + constraints.radius*constraints.radius)
      border_length = 3.0*grid_spacing
      center1 = center1 + axis12.normal_vec().scaled_vec(border_length)
      center2 = center2 - axis12.normal_vec().scaled_vec(border_length)
      radius = radius - border_length
      inner_constraint_fn = get_cylinder_constraint_fn(
          center1, center2, radius)
      print "%s constraints: center = %r, axis = %r, radius = %g" % \
          (constraints.type, center, axis12, constraints.radius)
    elif constraints.type in ("brick", "box"):
      atom1 = find_atom(
          atoms, constraints.chain1, constraints.res_num1,
          constraints.atom1)
      atom2 = find_atom(
          atoms, constraints.chain2, constraints.res_num2,
          constraints.atom2)
      try:
        offset1 = vector3d.Vector3d(*constraints.offset1)
      except AttributeError:
        offset1 = vector3d.Vector3d(0,0,0)
      try:
        offset2 = vector3d.Vector3d(*constraints.offset2)
      except AttributeError:
        offset2 = vector3d.Vector3d(0,0,0)

      p1 = atom1.pos + offset1
      p2 = atom2.pos + offset2

      # true lower-left-front ("min") and upper-right-back ("max")
      def toarray(v):
        return numpy.array([v.x, v.y, v.z])

      p1 = toarray(p1)
      p2 = toarray(p2)

      pmin = numpy.where(p1 < p2, p1, p2)
      pmax = numpy.where(p1 < p2, p2, p1)
      lower_left_front_corner = vector3d.Vector3d(*pmin)
      upper_right_back_corner = vector3d.Vector3d(*pmax)

      longest_side = max(pmax - pmin) 
      width = longest_side + 2.0*grid_spacing
      center = lower_left_front_corner + upper_right_back_corner
      center.scale(0.5)
      constraint_fn = get_brick_constraint_fn(lower_left_front_corner, upper_right_back_corner)
      inner_constraint_fn = constraint_fn   # XXX should be reduced in size
      print "%s constraints: loleftfr = %r uprightbk = %r" % \
          (constraints.type, lower_left_front_corner, upper_right_back_corner)
    elif constraints.type in ('coordbrick', 'coordbox'):
      try:
        coord1 = vector3d.Vector3d(*constraints.coord1)
      except AttributeError:
        coord1 = vector3d.Vector3d(0,0,0)
      try:
        coord2 = vector3d.Vector3d(*constraints.coord2)
      except AttributeError:
        coord2 = vector3d.Vector3d(0,0,0)
      p1 = coord1
      p2 = coord2
      # true lower-left-front ("min") and upper-right-back ("max")
      def toarray(v):
        return numpy.array([v.x, v.y, v.z])

      p1 = toarray(p1)
      p2 = toarray(p2)

      pmin = numpy.where(p1 < p2, p1, p2)
      pmax = numpy.where(p1 < p2, p2, p1)
      lower_left_front_corner = vector3d.Vector3d(*pmin)
      upper_right_back_corner = vector3d.Vector3d(*pmax)

      longest_side = max(pmax - pmin)
      width = longest_side + 2.0*grid_spacing
      center = lower_left_front_corner + upper_right_back_corner
      center.scale(0.5)
      constraint_fn = get_brick_constraint_fn(lower_left_front_corner, upper_right_back_corner)
      inner_constraint_fn = constraint_fn   # XXX should be reduced in size
      print "%s constraints: loleftfr = %r uprightbk = %r" % \
          (constraints.type, lower_left_front_corner, upper_right_back_corner)

    else:
      raise TypeError("Don't understand constraint type %r" % constraint.type)

  # Make the grid
  n_point = width / grid_spacing
  print "Setting up grid: %d x %d x %d, spacing %.3f, width %.1f" \
         % (n_point, n_point, n_point, grid_spacing, width)
  grid = Grid(grid_spacing, width, center)
  print_time(timer)

  print "Excluding protein bulk from grid with %.1f angstrom probe" \
          % size_interior_probe
  timer.start()
  exclude_atoms_from_grid(grid, atoms, size_interior_probe)
  print_time(timer)

  if constraint_file:
    # Eliminate all grid points outside constraints
    print "Excluding grid points outside constraint"
    timer.start()
    grid.exclude_points_in_constraint(constraint_fn)
    print_time(timer)


  is_calculate_asa_shell = True
  if constraint_file:
    if not constraints.remove_asa_shell:
      is_calculate_asa_shell = False

  if is_calculate_asa_shell:
    print "Calculating asa of atoms in protein"
    # Roll large ball over surface residues to eliminate
    # grid points over surface, then drill in from the
    # edge to eliminate the rest
    timer.start()
    add_asa(atoms, 1.4)
    print_time(timer)

    print "Excluding surface shell from grid with %.1f angstrom probe" \
            % size_surface_probe
    timer.start()
    exclude_surface(grid, atoms, size_surface_probe)
    print_time(timer)

    print "Excluding edges"
    grid.exclude_edge_to_interior()

  print "Excluding surrounded points"
  timer.start()
  hole_size = int(1.5 * 1.4 / grid_spacing)

  if is_hollow == 'yes':
    grid.exclude_surrounded(hole_size)
  elif is_hollow == 'no':
    continue
  else:
    raise ValueError('Please enter yes or no for the is_hollow option.')
  print_time(timer)

  # Make hollow spheres from grid-points
  if not out_pdb:
    out_pdb = pdb.replace('.pdb', '-hollow.pdb')
  print "Saving hollow spheres to", out_pdb
  grid_chain = grid.make_mol(defaults.res_type, defaults.atom_type)
  if size_bfactor_probe:
    print "Averaging nearby protein b-factors for each hollow atom"
    timer.start()
    if constraint_file:
      atoms = [a for a in atoms if constraint_fn(a.pos)]
    calculate_average_bfactor(grid_chain, atoms, size_bfactor_probe)
    print_time(timer)
  if constraint_file:
    for atom in grid_chain.atoms():
      if inner_constraint_fn(atom.pos):
        atom.occupancy = 1.0
      else:
        atom.occupancy = 0.0
  is_hetatm = not defaults.atom_field.startswith("ATOM")
  pdbstruct.save_atoms(grid_chain.atoms(), out_pdb, is_hetatm)

  return grid

# monkey-patch the defaults into the doc string
make_hollow_spheres.__doc__ = make_hollow_spheres.__doc__ % defaults.__dict__

idle_help_txt = """
  ----------------------

  If you see this message, then you are probably trying to run
  hollow.py using the IDLE interpreter.

  Load the hollow.py module in the interactive IDLE interpreter:

    >>> import hollow

  At any point after this, you can retrieve this message by typing:

    >>> hollow.help()

  Let's say for arguments sake that your pdb file is 3hbs.pdb. Then
  to generate the hollow spheres in the automated mode, where the
  hollow spheres will be written to the file '3hbs-hollow.pdb':

    >>> hollow.make_hollow_spheres('3hbs.pdb')

  If you want to choose your own output hollow filename:

    >>> hollow.make_hollow_spheres('3hbs.pdb','hollow.pdb')

  If you want to use a cylinder constraint file or a spherical
  constraint file, then run the command using the following parameters:

    >>> hollow.make_hollow_spheres('3hbs.pdb', constraint_file='constraint')

  To change the grid spacing:

    >>> hollow.make_hollow_spheres('3hbs.pdb', grid_spacing=0.5)

  Indeed, all the options recognized by hollow.make_hollow_spheres are:
"""

definition_str = """
   make_hollow_spheres(
      pdb,
      out_pdb="",
      grid_spacing=%f,
      size_interior_probe=%f,
      is_skip_waters=%s,
      size_surface_probe=%f,
      constraint_file="",
      size_bfactor_probe=%f)
  """ % \
  (defaults.grid_spacing, defaults.interior_probe,
   str(defaults.is_skip_waters), defaults.surface_probe,
   defaults.bfactor_probe)


function_txt = """
  Only the first parameter is required, all other parameters have
  default values. To override the default values, just replace them
  in order up to the ones you want to replace, such as:

    >>> hollow.make_hollow_spheres('3hbs.pdb', 'hollow.pdb',
            0.25, 1.4, False, 3.3, '', 0)

  Another method is to use keyword replacement:

    >>> hollow.make_hollow_spheres('3hbs.pdb', is_skip_waters=False)
"""


def help():
  print idle_help_txt, definition_str, function_txt

