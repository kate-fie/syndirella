from config import config


def getNeigIdxs(mol, atomIdxs, n_hops):
    visitedAtomsIx = set([])
    for idx in atomIdxs:
        current_atom = mol.GetAtomWithIdx(idx)
        visitedAtomsIx = _getNeigIdxs(current_atom, n_hops, visitedAtomsIdx=visitedAtomsIx)
    return visitedAtomsIx


def _getNeigIdxs(current_atom, remaining_hops, visitedAtomsIdx= None):
  if visitedAtomsIdx is None:
      visitedAtomsIdx = set([])
  visitedAtomsIdx.add(current_atom.GetIdx())
  if remaining_hops > 0:
    for neig in current_atom.GetNeighbors():
        if neig.GetIdx() not in visitedAtomsIdx:
            _getNeigIdxs(neig, remaining_hops - 1, visitedAtomsIdx)

    return visitedAtomsIdx


def expand2DIdxsToNeigs(mol, attachmentIdxs, n_hops_from_attachment=config.SMARTS_N_HOPS_FROM_ATTACHMENT):
    assert n_hops_from_attachment >= 0, "Error, n_hops_from_attachment must be 0 or larger"
    if isinstance(attachmentIdxs, int):
        attachmentIdxs = [attachmentIdxs]
    attachment_region_idxs = set(attachmentIdxs)
    if n_hops_from_attachment:
        neigs = getNeigIdxs(mol, attachmentIdxs, n_hops_from_attachment)
        attachment_region_idxs = attachment_region_idxs.union(neigs)
    return attachment_region_idxs