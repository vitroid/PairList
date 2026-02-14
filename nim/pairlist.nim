## Nim bindings for pairlist (pairs_fine equivalent)
## Returns pairs (i, j) within rc under PBC. Cartesian coordinates.
## Compile from project root: nim c nim/pairlist.nim

{.compile: "csource/pairlist.c".}
{.passC: "-Icsource".}

proc c_free(p: pointer) {.cdecl, importc: "free".}

type
  Vec3* = array[3, float64]
  Pair* = tuple[i, j: int]

proc pairlist*(nAtoms0: cint, atoms0: ptr cdouble, nAtoms1: cint, atoms1: ptr cdouble,
               lower: cdouble, higher: cdouble, cell: ptr cdouble,
               pairs: ptr ptr cint): cint {.cdecl, importc: "pairlist".}

proc pairlistSeq*(pos: seq[Vec3], rc: float64, cell: Vec3,
                  pos2: seq[Vec3] = @[]): seq[Pair] =
  ## pairs_fine equivalent. Returns (i, j) pairs within rc (Cartesian, PBC).
  ## pos2 empty = homo, otherwise hetero.
  var pairs: ptr cint
  let
    n0 = cint pos.len
    n1 = cint pos2.len
    lower = 0.0
  if n0 == 0:
    return @[]
  var atoms0 = newSeq[cdouble](3 * pos.len)
  for i, p in pos:
    atoms0[i*3+0] = p[0]
    atoms0[i*3+1] = p[1]
    atoms0[i*3+2] = p[2]
  var cellArr: array[3, cdouble] = [cell[0].cdouble, cell[1].cdouble, cell[2].cdouble]
  let nPairs = if pos2.len == 0:
    pairlist(n0, addr atoms0[0], 0, nil, lower.cdouble, rc.cdouble,
             addr cellArr[0], addr pairs)
  else:
    var atoms1 = newSeq[cdouble](3 * pos2.len)
    for i, p in pos2:
      atoms1[i*3+0] = p[0]
      atoms1[i*3+1] = p[1]
      atoms1[i*3+2] = p[2]
    pairlist(n0, addr atoms0[0], n1, addr atoms1[0], lower.cdouble, rc.cdouble,
             addr cellArr[0], addr pairs)
  result = newSeq[Pair](nPairs)
  let arr = cast[ptr UncheckedArray[cint]](pairs)
  for k in 0 ..< nPairs:
    result[k] = (arr[k*2+0].int, arr[k*2+1].int)
  c_free(pairs)

iterator pairlistIter*(pos: seq[Vec3], rc: float64, cell: Vec3,
                       pos2: seq[Vec3] = @[]): Pair =
  ## Iterator version of pairlistSeq (pairs_iter equivalent).
  for p in pairlistSeq(pos, rc, cell, pos2):
    yield p

when isMainModule:
  # Quick test: 2x2x2 simple cubic
  var pos: seq[Vec3]
  for x in 0..1:
    for y in 0..1:
      for z in 0..1:
        pos.add [float64(x), float64(y), float64(z)]
  let cell: Vec3 = [2.0, 2.0, 2.0]
  let pairs = pairlistSeq(pos, 1.5, cell)
  echo "Pairs (rc=1.5): ", pairs.len
  for (i, j) in pairlistIter(pos, 1.5, cell):
    echo "  ", i, " - ", j
