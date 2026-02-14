## Nim bindings for pairlist (pairs_fine equivalent)
## Uses Pairs + distance filter for triclinic cell support.
import std/math
## Compile from project root: nim c nim/pairlist.nim

{.compile: "csource/pairlist.c".}
{.passC: "-Icsource".}

proc c_free(p: pointer) {.cdecl, importc: "free".}

type
  Vec3* = array[3, float64]
  Mat3* = array[3, Vec3]  ## Cell matrix: rows = a, b, c
  Pair* = tuple[i, j: int]

proc Pairs*(npos: cint, rpos: ptr cdouble, ngrid: ptr cint,
            pairs: ptr ptr cint): cint {.cdecl, importc: "Pairs".}

# --- Math helpers ---
func dot*(a, b: Vec3): float64 = a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
func cross*(a, b: Vec3): Vec3 =
  [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
func norm*(a: Vec3): float64 = sqrt(a.dot(a))

## d @ cell: fractional delta -> Cartesian
func matVec*(cell: Mat3, d: Vec3): Vec3 =
  [d[0]*cell[0][0] + d[1]*cell[1][0] + d[2]*cell[2][0],
   d[0]*cell[0][1] + d[1]*cell[1][1] + d[2]*cell[2][1],
   d[0]*cell[0][2] + d[1]*cell[1][2] + d[2]*cell[2][2]]

func inv*(m: Mat3): Mat3 =
  let det = m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) -
            m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) +
            m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0])
  let idet = 1.0 / det
  result[0][0] = (m[1][1]*m[2][2]-m[1][2]*m[2][1]) * idet
  result[0][1] = (m[0][2]*m[2][1]-m[0][1]*m[2][2]) * idet
  result[0][2] = (m[0][1]*m[1][2]-m[0][2]*m[1][1]) * idet
  result[1][0] = (m[1][2]*m[2][0]-m[1][0]*m[2][2]) * idet
  result[1][1] = (m[0][0]*m[2][2]-m[0][2]*m[2][0]) * idet
  result[1][2] = (m[0][2]*m[1][0]-m[0][0]*m[1][2]) * idet
  result[2][0] = (m[1][0]*m[2][1]-m[1][1]*m[2][0]) * idet
  result[2][1] = (m[0][1]*m[2][0]-m[0][0]*m[2][1]) * idet
  result[2][2] = (m[0][0]*m[1][1]-m[0][1]*m[1][0]) * idet

## v @ inv(m): Cartesian -> fractional
func toFractional*(pos: Vec3, cell: Mat3): Vec3 = matVec(inv(cell), pos)

func determineGrid*(cell: Mat3, radius: float64): tuple[gx, gy, gz: int] =
  ## Grid divisions for triclinic cell (Python determine_grid)
  let
    a = cell[0]
    b = cell[1]
    c = cell[2]
    al = a.norm
    bl = b.norm
    cl = c.norm
    ae: Vec3 = [a[0]/al, a[1]/al, a[2]/al]
    be: Vec3 = [b[0]/bl, b[1]/bl, b[2]/bl]
    ce: Vec3 = [c[0]/cl, c[1]/cl, c[2]/cl]
    an = dot(a, cross(be, ce))
    bn = dot(b, cross(ce, ae))
    cn = dot(c, cross(ae, be))
  var
    gf0 = an / radius
    gf1 = bn / radius
    gf2 = cn / radius
  if gf0 < 1.0: gf0 = 1.0
  if gf1 < 1.0: gf1 = 1.0
  if gf2 < 1.0: gf2 = 1.0
  (floor(gf0).int, floor(gf1).int, floor(gf2).int)

proc pairsFine*(pos: seq[Vec3], rc: float64, cell: Mat3,
                fractional: bool = true): seq[Pair] =
  ## pairs_fine equivalent. Triclinic cell supported.
  ## pos: fractional coords (or Cartesian if fractional=false)
  if pos.len == 0:
    return @[]
  var rpos: seq[cdouble]
  var frac: seq[Vec3]  ## Fractional coords for distance calc
  if fractional:
    rpos = newSeq[cdouble](3 * pos.len)
    frac = pos
    for i, p in pos:
      rpos[i*3+0] = p[0].cdouble
      rpos[i*3+1] = p[1].cdouble
      rpos[i*3+2] = p[2].cdouble
  else:
    let celli = inv(cell)
    rpos = newSeq[cdouble](3 * pos.len)
    frac = newSeq[Vec3](pos.len)
    for i, p in pos:
      let f = matVec(celli, p)
      frac[i] = f
      rpos[i*3+0] = f[0].cdouble
      rpos[i*3+1] = f[1].cdouble
      rpos[i*3+2] = f[2].cdouble
  let grid = determineGrid(cell, rc)
  var ngrid: array[3, cint] = [grid.gx.cint, grid.gy.cint, grid.gz.cint]
  var pairs: ptr cint
  let nRough = Pairs(pos.len.cint, addr rpos[0], addr ngrid[0], addr pairs)
  let arr = cast[ptr UncheckedArray[cint]](pairs)
  result = newSeqOfCap[Pair](nRough)
  let rc2 = rc * rc
  for k in 0 ..< nRough:
    let i = arr[k*2+0].int
    let j = arr[k*2+1].int
    var d: Vec3 = [frac[i][0] - frac[j][0], frac[i][1] - frac[j][1], frac[i][2] - frac[j][2]]
    for dd in 0..2:
      d[dd] -= floor(d[dd] + 0.5)
    let cart = matVec(cell, d)
    let L2 = cart.dot(cart)
    if L2 < rc2:
      result.add (i, j)
  c_free(pairs)

iterator pairsIter*(pos: seq[Vec3], rc: float64, cell: Mat3,
                    fractional: bool = true): Pair =
  for p in pairsFine(pos, rc, cell, fractional):
    yield p

when isMainModule:
  # Test: 2x2x2 simple cubic, fractional coords, orthogonal cell
  var pos: seq[Vec3]
  for x in 0..1:
    for y in 0..1:
      for z in 0..1:
        pos.add [float64(x)/2.0, float64(y)/2.0, float64(z)/2.0]  # fractional [0, 0.5]
  # Orthogonal cell: a=(2,0,0), b=(0,2,0), c=(0,0,2)
  let cell: Mat3 = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]
  let pairs = pairsFine(pos, 1.5, cell)
  echo "Pairs (rc=1.5, fractional): ", pairs.len
  for (i, j) in pairsIter(pos, 1.5, cell):
    echo "  ", i, " - ", j
