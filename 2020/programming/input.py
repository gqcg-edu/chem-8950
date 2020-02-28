Settings = dict()

Settings["basis"] = "6-31G"
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
# 'u' frozen unoccupied orbitals
#       If insufficient u's are provided, you're code will append additional 
#       u's to get the length to the number of molecular orbitals
# 'full' automatically set the active space to 'a' * nmo
Settings["active_space"] = 'oaaaaaa'

# If 'full' perform a full CI
# If an integer then is the excitation level
#    (1 = CIS, 2 = CISD, 3 = CISDT, etc.)
#Settings["excitation_level"] = 2
Settings["excitation_level"] = 'full' 
