import os, copy, tempfile, re, glob, subprocess, time


def re_glob(dir_tag, reg_exp=""):
  fnames = glob.glob(dir_tag)
  return [f for f in fnames if re.search(reg_exp, f)]


def goto_dir(new_dir):
  if not os.path.isdir(new_dir):
    os.makedirs(new_dir)
  os.chdir(new_dir)


def insert_path(path, insert):
  if path.startswith('/'):
    return path
  else:
    return os.path.join(insert, path)


def temp_fname(suffix=''):
  fd, fname = tempfile.mkstemp(suffix, 'tmp-', '.')
  f = os.fdopen(fd, 'w')
  f.close()
  os.unlink(fname)
  return fname


def fname_variant(fname):
  root, ext = os.path.splitext(fname)
  i = 1
  new_fname = "%s-%d%s" % (root, i, ext)
  while os.path.isfile(new_fname):
    i += 1
    new_fname = "%s-%d%s" % (root, i, ext)
  return new_fname


def clean_fname(fname):
  try:
    os.remove(fname)
  except:
    pass


def run_with_output(cmd):
  p = subprocess.Popen(cmd, shell=True,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT)
  return p.stdout.read()


def is_same_dict_in_file(parms, fname):
  try:
    saved_parms = eval(open(fname).read())
    result = saved_parms == parms
  except:
    result = False
  return result


def replace_dict(t, d, max_len = 0):
  """Replace $token in string t using a dictionary

  By default, converts floats to %.3f format:
    - t: string whose tokens will be replace_dict
    - d: dictionary of tokens (string keys) and their 
         replacements (values, may be of any type)
    - max_len: maximum length (in characters) of 
               replacement values.
  """

  s = copy.deepcopy(t)
  for (k, v) in d.iteritems():
    token = "$%s" % k
    if max_len == 0:
      if type(v) is float:
        s = s.replace(token, "%.3f" % v)
      else:
        s = s.replace(token, str(v))
    else:
      if type(v) is float:
        s = s.replace(token, ("%.3f" % v).rjust(max_len))
      else:
        s = s.replace(token, str(v).rjust(max_len))
  return s


def write_dict(d, fname):
  """Outputs the dictionary to a file in repr format."""
  open(fname, "w").write(format_dict(d))


def format_dict(d):
  """Makes repr() comaptible string of a dictionary"""
  s = "{ "
  n = len(d)
  items = sorted(d.items())
  max_len = max([len(repr(k)) for k in d.keys()])
  for i in range(n):
    if i > 0:
      s += "  "
    key, value = items[i]
    s += "%s : %s" % (repr(key).ljust(max_len), 
                      repr(value))
    if i < n-1:
      s += ", \n"
  s += " }"
  return s


def words_in_file(fname):
  result = []
  for line in open(fname).readlines():
    result.extend(line.split())
  return result


def read_parameters(fname):
  class DataHolder: pass
  f = open(fname, 'r')
  result = DataHolder()
  result.__dict__ = eval(f.read())
  f.close()
  return result


def elapsed_time_str(time):
  s = str(time) + ' '
  minute = time / 60.0
  if minute > 60:
    s += "%.f:%02.f:" % (minute / 60, minute % 60)
  elif minute >= 1:
    s += "%.f:" % minute
  sec = time % 60.0
  if sec < 0.01:
    s += "%07.4fs" % sec
  else:
    s += "%05.2fs" % sec
  return s
  
  
class Timer:
  def __init__(self):
    self._elapsed = 0;
    self._start = time.time()

  def start(self):
    self._start = time.time()
    self._elapsed = 0

  def stop(self):
    self._elapsed = time.time() - self._start

  def elapsed(self):
    if self._elapsed == 0:
      return time.time() - self._start
    else:
      return self._elapsed

  def str(self):
    elapsed_time = self.elapsed()
    return elapsed_time_str(elapsed_time)
    
  def __str__(self):
    return self.str()


