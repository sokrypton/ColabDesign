import numpy as np
from colabdesign.af.alphafold.common import residue_constants

def make_atom14_positions(batch):
  """Constructs denser atom positions (14 dimensions instead of 37)."""
  restype_atom14_to_atom37 = []  # mapping (restype, atom14) --> atom37
  restype_atom37_to_atom14 = []  # mapping (restype, atom37) --> atom14
  restype_atom14_mask = []

  for rt in residue_constants.restypes:
    atom_names = residue_constants.restype_name_to_atom14_names[
        residue_constants.restype_1to3[rt]]

    restype_atom14_to_atom37.append([
        (residue_constants.atom_order[name] if name else 0)
        for name in atom_names
    ])

    atom_name_to_idx14 = {name: i for i, name in enumerate(atom_names)}
    restype_atom37_to_atom14.append([
        (atom_name_to_idx14[name] if name in atom_name_to_idx14 else 0)
        for name in residue_constants.atom_types
    ])

    restype_atom14_mask.append([(1. if name else 0.) for name in atom_names])

  # Add dummy mapping for restype 'UNK'.
  restype_atom14_to_atom37.append([0] * 14)
  restype_atom37_to_atom14.append([0] * 37)
  restype_atom14_mask.append([0.] * 14)

  restype_atom14_to_atom37 = np.array(restype_atom14_to_atom37, dtype=np.int32)
  restype_atom37_to_atom14 = np.array(restype_atom37_to_atom14, dtype=np.int32)
  restype_atom14_mask = np.array(restype_atom14_mask, dtype=np.float32)

  # Create the mapping for (residx, atom14) --> atom37, i.e. an array
  # with shape (num_res, 14) containing the atom37 indices for this protein.
  residx_atom14_to_atom37 = restype_atom14_to_atom37[batch["aatype"]]
  residx_atom14_mask = restype_atom14_mask[batch["aatype"]]

  # Create a mask for known ground truth positions.
  residx_atom14_gt_mask = residx_atom14_mask * np.take_along_axis(
      batch["all_atom_mask"], residx_atom14_to_atom37, axis=1).astype(np.float32)

  # Gather the ground truth positions.
  residx_atom14_gt_positions = residx_atom14_gt_mask[:, :, None] * (
      np.take_along_axis(batch["all_atom_positions"],
                         residx_atom14_to_atom37[..., None],
                         axis=1))

  prot = {}
  prot["atom14_atom_exists"] = residx_atom14_mask
  prot["atom14_gt_exists"] = residx_atom14_gt_mask
  prot["atom14_gt_positions"] = residx_atom14_gt_positions

  prot["residx_atom14_to_atom37"] = residx_atom14_to_atom37

  # Create the gather indices for mapping back.
  residx_atom37_to_atom14 = restype_atom37_to_atom14[batch["aatype"]]
  prot["residx_atom37_to_atom14"] = residx_atom37_to_atom14

  # Create the corresponding mask.
  restype_atom37_mask = np.zeros([21, 37], dtype=np.float32)
  for restype, restype_letter in enumerate(residue_constants.restypes):
    restype_name = residue_constants.restype_1to3[restype_letter]
    atom_names = residue_constants.residue_atoms[restype_name]
    for atom_name in atom_names:
      atom_type = residue_constants.atom_order[atom_name]
      restype_atom37_mask[restype, atom_type] = 1

  residx_atom37_mask = restype_atom37_mask[batch["aatype"]]
  prot["atom37_atom_exists"] = residx_atom37_mask

  # As the atom naming is ambiguous for 7 of the 20 amino acids, provide
  # alternative ground truth coordinates where the naming is swapped
  restype_3 = [
      residue_constants.restype_1to3[res] for res in residue_constants.restypes
  ]
  restype_3 += ["UNK"]

  # Matrices for renaming ambiguous atoms.
  all_matrices = {res: np.eye(14, dtype=np.float32) for res in restype_3}
  for resname, swap in residue_constants.residue_atom_renaming_swaps.items():
    correspondences = np.arange(14)
    for source_atom_swap, target_atom_swap in swap.items():
      source_index = residue_constants.restype_name_to_atom14_names[
          resname].index(source_atom_swap)
      target_index = residue_constants.restype_name_to_atom14_names[
          resname].index(target_atom_swap)
      correspondences[source_index] = target_index
      correspondences[target_index] = source_index
      renaming_matrix = np.zeros((14, 14), dtype=np.float32)
      for index, correspondence in enumerate(correspondences):
        renaming_matrix[index, correspondence] = 1.
    all_matrices[resname] = renaming_matrix.astype(np.float32)
  renaming_matrices = np.stack([all_matrices[restype] for restype in restype_3])

  # Pick the transformation matrices for the given residue sequence
  # shape (num_res, 14, 14).
  renaming_transform = renaming_matrices[batch["aatype"]]

  # Apply it to the ground truth positions. shape (num_res, 14, 3).
  alternative_gt_positions = np.einsum("rac,rab->rbc",
                                       residx_atom14_gt_positions,
                                       renaming_transform)
  prot["atom14_alt_gt_positions"] = alternative_gt_positions

  # Create the mask for the alternative ground truth (differs from the
  # ground truth mask, if only one of the atoms in an ambiguous pair has a
  # ground truth position).
  alternative_gt_mask = np.einsum("ra,rab->rb",
                                  residx_atom14_gt_mask,
                                  renaming_transform)

  prot["atom14_alt_gt_exists"] = alternative_gt_mask

  # Create an ambiguous atoms mask.  shape: (21, 14).
  restype_atom14_is_ambiguous = np.zeros((21, 14), dtype=np.float32)
  for resname, swap in residue_constants.residue_atom_renaming_swaps.items():
    for atom_name1, atom_name2 in swap.items():
      restype = residue_constants.restype_order[
          residue_constants.restype_3to1[resname]]
      atom_idx1 = residue_constants.restype_name_to_atom14_names[resname].index(
          atom_name1)
      atom_idx2 = residue_constants.restype_name_to_atom14_names[resname].index(
          atom_name2)
      restype_atom14_is_ambiguous[restype, atom_idx1] = 1
      restype_atom14_is_ambiguous[restype, atom_idx2] = 1

  # From this create an ambiguous_mask for the given sequence.
  prot["atom14_atom_is_ambiguous"] = (restype_atom14_is_ambiguous[batch["aatype"]])
  return prot