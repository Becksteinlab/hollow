

def strip_lines(pdb_txt, tag_func):
  new_lines = []
  for line in pdb_txt.splitlines():
    if tag_func(line):
      continue
    new_lines.append(line)
  return '\n'.join(new_lines)


def sort_file_by_line(fname, key_func):
  lines = open(fname, 'r').read().splitlines()
  pairs = [(key_func(l), l) for l in lines]
  new_lines = [l for (v, l) in sorted(pairs)]
  open(fname, 'w').write('\n'.join(new_lines))


def replace_txt_in_file(fname, tag, replacement):
  f = open(fname, 'r')
  txt = f.read().replace(tag, replacement)
  f.close()
  f = open(fname, 'w')
  f.write(txt)
  f.close()


def strip_hydrogens(pdb_txt):
  new_lines = []
  for line in pdb_txt.splitlines():
    if line.startswith("ATOM"):
      atom_type = line[12:16]
      if atom_type[0] == "H" or atom_type[1] == "H":
        continue
    new_lines.append(line)
  return '\n'.join(new_lines)


def renumber_residues(pdb_txt):

  get_res_tag = lambda line: line[17:27]

  sorted_res_tags = []
  lines = pdb_txt.splitlines()
  for line in lines:
    if line.startswith('ATOM') or line.startswith('HETATM'):
      tag = get_res_tag(line)
      if tag not in sorted_res_tags:
        sorted_res_tags.append(tag)
  
  res_tag_to_new_resnum = {}
  for i, tag in enumerate(sorted_res_tags):
    res_tag_to_new_resnum[tag] = "%3d " % (i+1)

  new_lines = []
  for line in lines:
    new_line = line
    for start_tag in ['ATOM', 'ANISOU', 'HETATM']:
      if line.startswith(start_tag):
        tag = get_res_tag(line)
        resnum = res_tag_to_new_resnum[tag]
        new_line = line[:23] + resnum + line[27:]
    new_lines.append(new_line)
  return '\n'.join(new_lines)


def strip_other_nmr_models(pdb_txt):
  new_lines = []
  for line in pdb_txt.splitlines():
    new_lines.append(line)
    if line.startswith("ENDMDL"):
      break
  return '\n'.join(new_lines)


def strip_alternative_atoms(pdb_txt):
  new_lines = []
  for line in pdb_txt.splitlines():
    new_line = line
    if line.startswith('ATOM'):
      alt_loc = line[16]
      if not alt_loc in [' ']:
        if alt_loc in ['A', 'a']:
          new_line = line[:16] + ' ' + line[17:]
        else:
          continue
    new_lines.append(new_line)
  return '\n'.join(new_lines)


def clean_pdb(in_pdb, out_pdb):
  txt = open(in_pdb, 'r').read()
  txt = strip_other_nmr_models(txt)
  txt = strip_lines(txt, lambda l: l.startswith('HETATM'))
  txt = strip_lines(txt, lambda l: l.startswith('ANISOU'))
  txt = strip_lines(txt, lambda l: l.startswith('CONECT'))
  txt = strip_lines(txt, lambda l: l.startswith('MASTER'))
  txt = strip_lines(txt, lambda l: l.startswith("ATOM") and \
                                   l[17:20] in ['WAT', 'TIP', 'HOH'])
  txt = strip_alternative_atoms(txt)
  txt = strip_hydrogens(txt)
  txt = renumber_residues(txt)
  open(out_pdb, 'w').write(txt)
