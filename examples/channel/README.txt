==========================================
 Visualizing the OprP channel with hollow
==========================================

Generate surface
================

See http://hollow.sourceforge.net/channel.html for details.

Commandline::

  hollow -c constraint -o hollow.pdb 2o4v-a.pdb 


Visualization
=============

In PyMOL, load the Pdb file and the hollow.pdb::

  pymol hollow.pdb 2o4v-a.pdb 

Once in pymol, we show the channel surface by showing the surface of
hollow spheres inside the cylinder (q=1) using two-sided lighting::

  pymol> show surface, hollow and q>0
  pymol> hide nonbonded
  pymol> set two_sided_lighting, on 


Using the hollow spheres, choosing the channel residues is trivial::

  pymol> select lining, byres hollow around 5
  pymol> show sticks, lining
  pymol> cartoon tube





