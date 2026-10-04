"""Immutable file-only cumulative FIMA and deposited-power histories.

No OpenMC runtime import; existing FIMA-only tables remain supported.
"""
import csv
import hashlib
from pathlib import Path

import numpy as np


class DepletionHistory:
    """Physical annular bins in cm and saved times in seconds.

    FIMA is fissions / initial U+Zr atoms. Interpolate cumulative quantities
    linearly between recorded times; never extrapolate or infer FIMA from BU.
    """

    def __init__(self, path):
        self.path = Path(path).resolve()
        self.sha256 = hashlib.sha256(self.path.read_bytes()).hexdigest()
        with self.path.open(newline='') as stream:
            rows = list(csv.DictReader(stream))
        if not rows:
            raise ValueError('Empty depletion history')
        required = ('time_s', 'radial_index', 'axial_index', 'r_in_cm',
                    'r_out_cm', 'z_low_cm', 'z_high_cm', 'volume_cm3',
                    'initial_U_Zr_atoms', 'cumulative_fissions', 'FIMA')
        if any(k not in rows[0] for k in required):
            raise ValueError(f'Depletion table requires {required}')
        data = np.array([[float(row[k]) for k in required] for row in rows])
        if not np.isfinite(data).all():
            raise ValueError('Non-finite depletion data')
        if not np.equal(data[:, 1:3], np.floor(data[:, 1:3])).all():
            raise ValueError('Domain indices must be integers')
        self.times = np.unique(data[:, 0])
        self.keys = sorted({tuple(map(int, row)) for row in data[:, 1:3]})
        nr = max(k[0] for k in self.keys)
        nz = max(k[1] for k in self.keys)
        if self.keys != [(r, z) for r in range(1, nr+1) for z in range(1, nz+1)]:
            raise ValueError('Domain indices must form a complete positive grid')
        self.shape = (nr, nz)
        if self.times[0] != 0 or len(self.times) < 2:
            raise ValueError('History must start at zero and contain >=2 times')
        ordered = []
        for t in self.times:
            block = data[data[:, 0] == t]
            block = block[np.lexsort((block[:, 2], block[:, 1]))]
            if len(block) != len(self.keys) or list(map(tuple, block[:, 1:3].astype(int))) != self.keys:
                raise ValueError('Duplicate or missing time/domain key')
            ordered.append(block)
        ordered = np.array(ordered)
        if not np.allclose(ordered[:, :, 3:9], ordered[0, :, 3:9], rtol=1e-13, atol=0):
            raise ValueError('Geometry or initial inventories change with time')
        self.bounds_cm = ordered[0, :, 3:7].copy()
        self.volumes_cm3 = ordered[0, :, 7].copy()
        self.initial_atoms = ordered[0, :, 8].copy()
        self.cumulative_fissions = ordered[:, :, 9].copy()
        self.fima = ordered[:, :, 10].copy()
        ri, ro, zl, zh = self.bounds_cm.T
        volume = np.pi*(ro**2-ri**2)*(zh-zl)
        if (ri < 0).any() or (ro <= ri).any() or (zh <= zl).any() or (self.initial_atoms <= 0).any():
            raise ValueError('Invalid physical bins or initial atom inventory')
        if not np.allclose(volume, self.volumes_cm3, rtol=1e-11, atol=1e-12):
            raise ValueError('Bin volume disagrees with physical coordinates')
        grid = self.bounds_cm.reshape(nr, nz, 4)
        if not np.allclose(grid[:, :, :2], grid[:, :1, :2], rtol=0, atol=1e-12) or not np.allclose(grid[:, :, 2:], grid[:1, :, 2:], rtol=0, atol=1e-12):
            raise ValueError('Annular/axial boundaries must be consistent')
        if not np.allclose(grid[1:, 0, 0], grid[:-1, 0, 1], rtol=0, atol=1e-12) or not np.allclose(grid[0, 1:, 2], grid[0, :-1, 3], rtol=0, atol=1e-12):
            raise ValueError('Depletion grid has gaps or overlaps')
        if (self.fima < 0).any() or (self.cumulative_fissions < 0).any():
            raise ValueError('Negative cumulative exposure')
        if np.any(self.fima[0] != 0) or np.any(self.cumulative_fissions[0] != 0):
            raise ValueError('Initial exposure must be zero')
        if (np.diff(self.fima, axis=0) < 0).any() or (np.diff(self.cumulative_fissions, axis=0) < 0).any():
            raise ValueError('Cumulative exposure must be monotone')
        if not np.allclose(self.fima*self.initial_atoms, self.cumulative_fissions, rtol=1e-11, atol=0):
            raise ValueError('FIMA disagrees with cumulative fissions / initial atoms')
        for array in (self.times, self.bounds_cm, self.volumes_cm3, self.initial_atoms, self.cumulative_fissions, self.fima):
            array.setflags(write=False)

    def fissions_at(self, time_s):
        t = float(time_s)
        if not np.isfinite(t) or t < self.times[0] or t > self.times[-1]:
            raise ValueError(f'Time {t} outside saved history [{self.times[0]}, {self.times[-1]}]')
        hi = int(np.searchsorted(self.times, t, side='left'))
        if self.times[hi] == t:
            return self.cumulative_fissions[hi].copy()
        lo = hi-1
        alpha = (t-self.times[lo])/(self.times[hi]-self.times[lo])
        return (1-alpha)*self.cumulative_fissions[lo]+alpha*self.cumulative_fissions[hi]

    def fima_at(self, time_s):
        return self.fissions_at(time_s)/self.initial_atoms


class HeatingHistory(DepletionHistory):
    """Saved normalized deposited_power_W, never inferred from exposure.

    Linear interpolation of instantaneous power at the saved transport times;
    preserves nonnegativity and the interpolated total power. No extrapolation.
    This is distinct from BOS integration of cumulative energy in depletion.
    """

    def __init__(self, path):
        super().__init__(path)
        with self.path.open(newline='') as stream:
            rows=list(csv.DictReader(stream))
        if 'deposited_power_W' not in rows[0]:
            raise ValueError('Native heating requires deposited_power_W in watts')
        by_key={(float(r['time_s']),int(r['radial_index']),int(r['axial_index'])):float(r['deposited_power_W']) for r in rows}
        self.power_W=np.array([[by_key[(t,*key)] for key in self.keys] for t in self.times])
        if not np.isfinite(self.power_W).all() or (self.power_W<0).any():
            raise ValueError('Power must be finite and nonnegative')
        self.power_W.setflags(write=False)

    def power_at(self, time_s):
        # Reuse strict range validation; no clamping at the end of history.
        self.fissions_at(time_s)
        t=float(time_s);hi=int(np.searchsorted(self.times,t))
        if self.times[hi]==t:
            return self.power_W[hi].copy()
        lo=hi-1;a=(t-self.times[lo])/(self.times[hi]-self.times[lo])
        return (1-a)*self.power_W[lo]+a*self.power_W[hi]

    def qdot_at(self, time_s):
        return self.power_at(time_s)/(self.volumes_cm3*1e-6)
