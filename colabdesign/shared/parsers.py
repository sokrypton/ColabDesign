####################################
# some code borrowed from AlphaFold
####################################
import string
from colabdesign.af.alphafold.common import residue_constants
import numpy as np

aa_order = {a:n for n,a in enumerate(residue_constants.restypes + ["X","-"])}
order_aa = {n:a for a,n in aa_order.items()}
deletion_table = str.maketrans('', '', string.ascii_lowercase)

def aa2num(seq):
  return [aa_order.get(res,20) for res in seq.translate(deletion_table)]

def num2aa(seq):
  return "".join([order_aa.get(res,"X") for res in seq])

def parse_fasta(fasta_string):
  """Parses FASTA string and returns list of strings with amino-acid sequences.

  Arguments:
    fasta_string: The string contents of a FASTA file.

  Returns:
    A tuple of two lists:
    * A list of sequences.
    * A list of sequence descriptions taken from the comment lines. In the
      same order as the sequences.
  """
  sequences = []
  descriptions = []
  index = -1
  for line in fasta_string.splitlines():
    line = line.strip()
    if line.startswith('>'):
      index += 1
      descriptions.append(line[1:])  # Remove the '>' at the beginning.
      sequences.append('')
      continue
    elif not line:
      continue  # Skip blank lines.
    sequences[index] += line
  return sequences, descriptions

def parse_a3m(a3m_filename=None, a3m_string=None):
  """Parses sequences and deletion matrix from a3m format alignment.

  Args:
    a3m_string: The string contents of a a3m file. The first sequence in the
      file should be the query sequence.

  Returns:
    A tuple of:
      * A list of sequences that have been aligned to the query. These
        might contain duplicates.
      * The deletion matrix for the alignment as a list of lists. The element
        at `deletion_matrix[i][j]` is the number of residues deleted from
        the aligned sequence i at residue position j.
      * A list of descriptions, one per sequence, from the a3m file.
  """
  if a3m_filename is not None:
    a3m_string = open(a3m_filename,"r").read()
  sequences, descriptions = parse_fasta(a3m_string)
  deletion_matrix = []
  for msa_sequence in sequences:
    deletion_vec = []
    deletion_count = 0
    for j in msa_sequence:
      if j.islower():
        deletion_count += 1
      else:
        deletion_vec.append(deletion_count)
        deletion_count = 0
    deletion_matrix.append(deletion_vec)

  # Make the MSA matrix out of aligned (deletion-free) sequences.
  msa = [aa2num(s) for s in sequences]
  return np.array(msa), np.array(deletion_matrix)