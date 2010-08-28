#!/usr/bin/env python
# coding: utf-8


import math
import vector3d
import array
import pdbstruct
import util


class Grid:
  def __init__(self, grid_spacing, width, center):
    self.width = float(width)
    half_width = self.width / 2.0
    self.center = center.copy()
    self.spacing = float(grid_spacing)
    self.inv_spacing = 1.0 / self.spacing

    self.n = 1
    cover = 0
    self.low = vector3d.Vector3d()
    while cover < half_width:
      self.n += 1
      half_n_point = int(self.n / 2)
      self.low.x = self.center.x - half_n_point*self.spacing
      self.low.y = self.center.y - half_n_point*self.spacing
      self.low.z = self.center.z - half_n_point*self.spacing
      width_1 = abs(self.center.x - self.low.x)
      high_x = self.low.x + self.n*self.spacing
      width_2 = abs(high_x - self.center.x)
      cover = min(width_1, width_2)
      
    self.actual_width = self.n*self.spacing
    self.total_volume = self.actual_width**3
    self.total_point = self.n**3
    
    self.array = array.array('b')
    for i in xrange(self.total_point):
      self.array.append(0)
    self.n_sq = self.n**2
      
    self.x = [self.low.x + i*self.spacing for i in xrange(self.n)]
    self.y = [self.low.y + j*self.spacing for j in xrange(self.n)]
    self.z = [self.low.z + k*self.spacing for k in xrange(self.n)]
      
  def reset(self):
    for i in xrange(self.total_point):
      self.array[i] = 0
    
  def indices(self, pos):
    return ((pos.x-self.low.x)*self.inv_spacing,
            (pos.y-self.low.y)*self.inv_spacing,
            (pos.z-self.low.z)*self.inv_spacing)

  def pos(self, i, j, k):
    return vector3d.Vector3d(self.x[i], self.y[j], self.z[k])

  def is_grid_point_near_sphere(self, i, j, k, pos, r_sq):
    d_x = self.x[i] - pos.x
    d_y = self.y[j] - pos.y 
    d_z = self.z[k] - pos.z
    return (d_x*d_x + d_y*d_y + d_z*d_z) < r_sq
    
  def int_range(self, low_f, high_f):
    low = max(0, int(math.floor(low_f-1)))
    high = min(self.n, int(math.ceil(high_f) + 2))
    return range(low, high)

  def exclude_sphere(self, pos, r):
    low = vector3d.Vector3d(pos.x - r, pos.y - r, pos.z - r)
    low_i, low_j, low_k = self.indices(low)
    high = vector3d.Vector3d(pos.x + r, pos.y + r, pos.z + r)
    high_i, high_j, high_k = self.indices(high)
    r_sq = r*r
    for i in self.int_range(low_i, high_i):
      for j in self.int_range(low_j, high_j):
        for k in self.int_range(low_k, high_k):
          l = i*self.n_sq + j*self.n + k 
          if self.array[l] == 0:
            if self.is_grid_point_near_sphere(i, j, k, pos, r_sq):
              self.array[l] = 1

  def n_excluded(self):
    return sum(self.array)
    
  def make_mol(self, res_type, atom_type):
    mol = pdbstruct.Molecule()
    element = ""
    for c in atom_type[:2]:
      if not c.isdigit() and c != " ":
        element += c
    i_res = 1
    for i in range(self.n):
      for j in range(self.n):
        for k in range(self.n):
          l = i*self.n_sq + j*self.n + k 
          if self.array[l]:
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

        
def volume(atoms, grid_spacing, out_fname=""):
  center = pdbstruct.get_center(atoms)
  width = pdbstruct.get_width(atoms, center) + 4.0
  grid = Grid(grid_spacing, width, center)
  print "Grid %d x %d x %d; Width %.2f Å" % \
          (grid.n, grid.n, grid.n, grid.actual_width)
  for atom in atoms:
    grid.exclude_sphere(atom.pos, atom.radius)
  d_volume = float(grid_spacing)**3
  volume = grid.n_excluded()*d_volume
  print u"Volume %.1f Å^3 (%d x %.3f Å^3)" \
           % (volume, grid.n_excluded(), d_volume)
  if out_fname:
    atoms = grid.make_mol("HOH", "O").atoms()
    pdbstruct.save_atoms(atoms, out_fname)


def load_psyco_if_in_system():
  try:
    import psyco
    psyco.full()
    print "Loaded optional psyco JIT compiler"
  except:
    pass


if __name__ == "__main__":
  usage = \
  """
  Calculates the volume of the atoms in a PDB file.
  Copyright (c) 2007 Bosco Ho.  

  Usage: volume.py in_pdb [spacing]
  
  -spacing    the spacing in the grid used to calculate the 
              volume grid (default=0.5). The smaller value,
              the more accurate, but slower.
  """
  import sys
  if len(sys.argv) < 2:
    print usage
    sys.exit(1)

  spacing = 0.5
  if len(sys.argv) > 2:
    spacing = float(sys.argv[2])
    
  load_psyco_if_in_system()
  
  mol = pdbstruct.Molecule(sys.argv[1])
  atoms = mol.atoms()
  pdbstruct.add_radii(atoms)

  print "%d atoms in input pdb file" % len(atoms)
  volume(atoms, spacing, sys.argv[1].replace('.pdb', '.grid.pdb'))
