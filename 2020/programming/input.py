Settings = dict()

Settings["basis"] = "STO-3G"
Settings["molecule"] = """
  0 1
  O
  H 1 R
  H 1 R 2 A
  R = 0.9
  A = 104.5
  symmetry c1
"""
Settings["nalpha"] = 5
Settings["nbeta"] = 5
Settings["scf_max_iter"] = 50

# 'o' frozen doubly occupied orbitals
# 'a' active orbitals
Settings["active_space"] = 'ooooaaa'  # Any remaining that remain refer to frozen unoccupied orbitals