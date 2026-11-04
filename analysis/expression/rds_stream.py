"""Reduce a dense R matrix stored in a .rds file without ever holding it in memory.

R writes a matrix column by column, so a single pass over the stream can add each
nucleus into the group it belongs to. Peak memory is one column plus the result.
"""
import gzip, struct, numpy as np

class Reader:
    def __init__(self, path):
        self.f = gzip.open(path, 'rb')
    def i4(self):
        return struct.unpack('>i', self.f.read(4))[0]
    def header(self):
        assert self.f.read(2) == b'X\n', 'not an XDR rds'
        ver, wr, mn = self.i4(), self.i4(), self.i4()
        if ver >= 3:
            n = self.i4()
            self.f.read(n)                      # native encoding name
        return ver
    def flags(self):
        v = self.i4()
        return v & 0xFF, bool(v & (1 << 9))     # type, has attributes
    def length(self):
        n = self.i4()
        if n == -1:                             # long vector
            hi, lo = self.i4(), self.i4()
            n = (hi << 32) | (lo & 0xFFFFFFFF)
        return n


def _item(r, refs):
    """enough of R's serialization to walk the attributes that trail the data"""
    flags = r.i4()
    t = flags & 0xFF
    has_attr, has_tag = bool(flags & (1 << 9)), bool(flags & (1 << 10))
    if t == 255:                                   # a back-reference
        idx = flags >> 8
        if idx == 0: idx = r.i4()
        return refs[idx - 1]
    if t in (254, 253, 242): return None           # nil, global env, empty env
    if t == 1:                                     # symbol
        name = _item(r, refs); refs.append(name); return name
    if t in (2, 17, 18):                           # pairlist
        if has_attr: _item(r, refs)
        tag = _item(r, refs) if has_tag else None
        car = _item(r, refs)
        cdr = _item(r, refs)
        out = [(tag, car)]
        if isinstance(cdr, list): out += cdr
        return out
    if t == 9:                                     # a single string
        n = r.i4()
        return None if n < 0 else r.f.read(n).decode("utf-8", "replace")
    if t in (13, 14, 16, 19):                      # int, double, character, list
        n = r.length()
        if t == 13:   v = np.frombuffer(r.f.read(4 * n), dtype=">i4")
        elif t == 14: v = np.frombuffer(r.f.read(8 * n), dtype=">f8")
        else:         v = [_item(r, refs) for _ in range(n)]
        if has_attr: _item(r, refs)
        return v
    raise ValueError("unhandled SEXP type %d" % t)

def read_attributes(r):
    """call straight after the data; returns {tag: value}"""
    return dict((k, v) for k, v in _item(r, []) if k)

def reduce_matrix(path, group=None, ngroup=1, nrow=None, cols_per_read=256):
    """Stream a dense R matrix, summing each column into its group.
    group[j] is the group of column j, -1 to skip. Returns a dict."""
    r = Reader(path); r.header()
    t, has_attr = r.flags()
    assert t == 14, "expected a double matrix, got SEXP type %d" % t
    n = r.length()
    if nrow is None:
        raise ValueError("nrow required")
    assert n % nrow == 0, "length %d not divisible by nrow %d" % (n, nrow)
    ncol = n // nrow
    if group is None: group = np.zeros(ncol, dtype=np.int64)
    assert len(group) == ncol, "group has %d entries, matrix has %d columns" % (len(group), ncol)
    num = np.zeros((nrow, ngroup)); colsum = np.zeros(ncol)
    j = 0
    while j < ncol:
        k = min(cols_per_read, ncol - j)
        blk = np.frombuffer(r.f.read(nrow * 8 * k), dtype=">f8").reshape(k, nrow).T
        colsum[j:j+k] = blk.sum(axis=0)
        g = group[j:j+k]
        for gi in np.unique(g):
            if gi < 0: continue
            num[:, gi] += blk[:, g == gi].sum(axis=1)
        j += k
    attrs = read_attributes(r) if has_attr else {}
    dim = attrs.get("dim"); dn = attrs.get("dimnames")
    genes = dn[0] if dn else None
    cells = dn[1] if dn and len(dn) > 1 else None
    return dict(num=num, colsum=colsum, ncol=ncol, nrow=nrow, genes=genes, cells=cells,
                dim=(list(dim) if dim is not None else None))

def peek_shape(path, nrow_hint=None):
    """length only, without reading the data"""
    r = Reader(path); r.header(); t, _ = r.flags()
    assert t == 14
    return r.length()
