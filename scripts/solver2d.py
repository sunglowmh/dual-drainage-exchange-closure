"""Explicit diffusive-wave 2D surface solver on the campus domain.

Scheme (LISFLOOD-FP style, "flow over the higher bed"):

    eta     = z + h
    dh_face = max(eta_i, eta_j) - max(z_i, z_j)        depth over the higher bed
    q       = sign(eta_i - eta_j) * dh_face^(5/3) * sqrt(|eta_i - eta_j| / dx) / n

    dh/dt   = (q_in_from_W - q_out_to_E + q_in_from_S - q_out_to_N) / dx

Two numerical points matter here:

* **donor volume limiter** -- dem_rev1 lowers the river bed by a uniform 2 m,
  which leaves ~2.4 m near-vertical steps along the channel; and the raised
  buildings give similar steps.  Without limiting the transferred volume to a
  quarter of the donor cell's stored water the scheme overflows to inf within
  ~10 steps.  With it, depths stay non-negative and the run is stable.
* static neighbour/index arrays and a `dt` update every few steps -- that takes
  the cost from ~29 ms to a few ms per step.

Boundary treatment
  * perimeter wall -> faces to inactive cells do not conduct: a closed basin;
  * gates -> one-way faces to a ghost cell holding a FIXED water surface (the
    municipal road level): a flap gate that only lets water out;
  * the river/lake inside the campus is ordinary storage (bed already lowered
    by dem_rev1).
"""
import numpy as np
import os
ROOT = os.environ.get("EMS_SWMM_BASE",
                   os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class Surface2D:
    def __init__(self, z, active, cell=2.0, manning=0.05, gates=(),
                 dry=1e-4, cfl=0.7, dt_max=2.0, dt_min=0.02,
                 dt_update_every=8, smax_cap=0.1, dtype=np.float64):
        self.cell = cell
        self.n_man = manning
        self.dry = dry
        self.cfl = cfl
        self.dt_max = dt_max
        self.dt_min = dt_min
        self.dt_update_every = dt_update_every
        self.smax_cap = smax_cap
        self.H, self.W = z.shape

        self.act = active
        self.k = np.full(z.shape, -1, np.int64)
        rr, cc = np.nonzero(active)
        self.k[rr, cc] = np.arange(len(rr))
        self.n = len(rr)
        self.rr, self.cc = rr, cc
        self.z = z[rr, cc].astype(dtype)
        self.h = np.zeros(self.n, dtype=dtype)

        def nb(dr, dc):
            r2 = np.clip(rr + dr, 0, self.H - 1)
            c2 = np.clip(cc + dc, 0, self.W - 1)
            v = self.k[r2, c2]
            valid = ((rr + dr >= 0) & (rr + dr < self.H) &
                     (cc + dc >= 0) & (cc + dc < self.W))
            return np.where(valid, v, -1)

        self.nE = nb(0, 1).astype(np.int64)
        self.nW = nb(0, -1).astype(np.int64)
        self.nN = nb(-1, 0).astype(np.int64)
        self.nS = nb(1, 0).astype(np.int64)

        # ---- static face bookkeeping (computed once) ----
        self._prep(which="E", nb=self.nE)
        iN, jN, zmaxN = self._prep(which="N", nb=self.nN, ret=True)
        self.iN, self.jN, self.zmaxN = iN, jN, zmaxN
        # inflow sources: cell k receives qE[nW[k]] and qN[nS[k]]
        self.iW = np.nonzero(self.nW >= 0)[0]
        self.srcW = self.nW[self.iW]
        self.iS = np.nonzero(self.nS >= 0)[0]
        self.srcS = self.nS[self.iS]

        # gate ghost cells
        self.gates = []
        for (gr, gc, stage, verdict) in gates:
            if 0 <= gr < self.H and 0 <= gc < self.W:
                kk = self.k[gr, gc]
                if kk >= 0:
                    self.gates.append((int(kk), float(stage), verdict))
        self.out_vol = 0.0
        self.clip_vol = 0.0
        self.lake = None                 # indices of the lake-hollow cells
        self.lake_open_t = None          # manual gate: only open after this time (s)
        # ---- surface -> sewer return inlets (second half of a two-way coupling)
        # Each manhole acts as an inlet: Q = coef * A * sqrt(2 g h_surface),
        # capped by the downstream conduit's full-flow capacity.  Water removed
        # here is accumulated per node and written back into SWMM as an
        # external inflow on the next coupling iteration.
        self.inlets = []                 # [(cell_index, node_id, q_cap, area)]
        self.inlet_coef = 0.6
        self.in_k = None                 # np.array of cell indices
        self.in_cap = None               # np.array of capacities (m3/s)
        self.in_area = None              # np.array of inlet areas (m2)
        self.in_ok = None                # (n_inlets, n_steps) bool: sewer can accept
        self.in_head = None              # (n_inlets, n_steps) Z_node (m, datum)
        self.in_zr = None                # (n_inlets,) manhole rim elevation
        self.in_step = 60.0              # time resolution of `in_ok`
        self.ret_ids = []                # node ids aligned with in_k
        self.ret_vol = 0.0               # total volume returned to the sewers
        self.src_vol = 0.0               # volume the solver actually INJECTED
        self.dry_loss = 0.0              # volume removed by dry-cell clipping
        self._ret_acc = {}
        self.ret_hist = []               # [(t, {node_id: m3})]
        self._t_sim = 0.0                # accumulated simulated time (s)
        self.lake_out_vol = 0.0
        self._dt = dt_max
        self._step_i = 0

    # ------------------------------------------------------- lake gate outlet
    def set_lake_outlet(self, cells, sill_stage, capacity, ramp=0.20):
        """Discharge of the lake through the sluice gate to the pump chamber.

        The gate is opened in flood; the pump station (0.722 m3/s on the main
        unit) is the binding constraint, so the outlet is modelled as a
        capacity-limited drain that switches on once the lake rises above the
        normal landscape water level (`sill_stage`) and ramps in over `ramp` m.
        Water is drawn from the hollow in proportion to depth, which mimics the
        slow draw-down of a flat lake surface.
        """
        self.lake = np.asarray(sorted(cells), dtype=np.int64)
        self.lake_sill = float(sill_stage)
        self.lake_cap = float(capacity)
        self.lake_ramp = float(ramp)

    def set_return_inlets(self, items, coef=0.6, ok=None, step=60.0, cap=None,
                          head=None, rim=None):
        """Register surface->sewer inlets.

        items : iterable of (node_id, cell_index, q_cap_m3s, area_m2).
                `q_cap` = full-flow capacity of the node's downstream conduit
                (used when `cap` is not supplied).
        cap   : optional (n_inlets, n_steps) array of the SPARE capacity per
                coupling interval = full-flow capacity - current conduit flow.
                Falls back to the scalar `q_cap` when omitted.
        ok    : optional (n_inlets, n_steps) bool array -- whether the SEWER can
                accept water during that coupling interval (False while the node
                is surcharged/flooding).
        """
        self.inlet_coef = coef
        self.in_step = step
        self.inlets = [(int(k), nid, float(qc), float(a))
                       for (nid, k, qc, a) in items if k is not None and k >= 0]
        self.ret_ids = [nid for (_, nid, _, _) in self.inlets]
        self.in_k = np.array([k for (k, _, _, _) in self.inlets], np.int64)
        qc = np.array([q for (_, _, q, _) in self.inlets], float)
        self.in_area = np.array([a for (_, _, _, a) in self.inlets], float)
        if cap is not None:
            cap = np.asarray(cap, float)
            if cap.shape[0] != self.in_k.size:
                raise ValueError("cap has %d rows, expected %d"
                                 % (cap.shape[0], self.in_k.size))
        self.in_cap = cap if cap is not None else qc
        if ok is not None:
            ok = np.asarray(ok, bool)
            if ok.shape[0] != self.in_k.size:
                raise ValueError("ok has %d rows, expected %d" % (ok.shape[0], self.in_k.size))
        self.in_ok = ok
        # HEAD-DRIVEN mode: `head` = Z_node(t) for every inlet x coupling step,
        # `rim` = manhole rim elevation.  When supplied the flow-ratio mask and
        # the spare-capacity cap are both bypassed.
        if head is not None:
            head = np.asarray(head, np.float32)
            rim = np.asarray(rim, float)
            if head.shape[0] != self.in_k.size or rim.size != self.in_k.size:
                raise ValueError("head/rim shape %s / %d, expected %d rows"
                                 % (head.shape, rim.size, self.in_k.size))
        self.in_head, self.in_zr = head, rim
        return self.in_k.size

    def set_lake_init(self, stage):
        """Fill the lake hollow up to `stage` (m) as the initial condition.

        The landscape pond normally holds ~1 m of water; the 2D domain is
        otherwise built dry, which silently grants the model the whole hollow
        as free storage.  Pre-storm draw-down is expressed by lowering `stage`.
        """
        if self.lake is None:
            return 0.0
        z = self.z[self.lake]
        h = np.maximum(stage - z, 0.0)
        self.h[self.lake] = h
        return float(h.sum()) * self.cell ** 2

    def lake_stage(self):
        if self.lake is None or self.lake.size == 0:
            return None
        eta = self.z[self.lake] + self.h[self.lake]
        wet = self.h[self.lake] > 0.01
        if not wet.any():
            return None
        return float(eta[wet].max())

    def _lake_drain(self, dt):
        if self.lake is None or self.lake.size == 0:
            return 0.0
        # A MANUAL gate stays shut until an operator opens it: with
        # lake_open_t set the drain is disabled before that time.
        # (`_t_sim` is accumulated in `step()`.)
        if self.lake_open_t is not None and self._t_sim < self.lake_open_t:
            return 0.0
        stage = self.lake_stage()
        if stage is None:
            return 0.0
        head = stage - self.lake_sill
        if head <= 0.0:
            return 0.0
        # Only the water ABOVE the gate sill may leave, and a DRY cell must not
        # register any drainable water: a cell whose bed sits above the sill
        # would otherwise contribute (bed - sill) of phantom storage and the
        # gate would over-draw (it removed 10 400 m3 of a 12 374 m3 injection in
        # a test).  For a cell below the sill the drainable column is the water
        # above the sill; for one above the sill it is the whole water column.
        z = self.z[self.lake]
        h = self.h[self.lake]
        eta = z + h
        drainable = np.where(z < self.lake_sill,
                             np.maximum(eta - self.lake_sill, 0.0), h)
        tot = float(drainable.sum()) * self.cell ** 2
        if tot <= 0.0:
            return 0.0
        q = self.lake_cap * min(1.0, head / self.lake_ramp)
        v = min(q * dt, tot)
        if v <= 0.0:
            return 0.0
        self.h[self.lake] -= (v / self.cell ** 2) * (drainable / drainable.sum())
        self.lake_out_vol += v
        return v

    def _prep(self, which, nb, ret=False):
        i = np.nonzero(nb >= 0)[0]
        j = nb[i]
        zmax = np.maximum(self.z[i], self.z[j])
        if which == "E":
            self.iE, self.jE, self.zmaxE = i, j, zmax
        if ret:
            return i, j, zmax

    # ------------------------------------------------------------------ flux
    def _flux(self, h, dt, i, j, zmax):
        """Flux from each cell to its east/north neighbour (m2/s), donor-limited."""
        ei = self.z[i] + h[i]
        ej = self.z[j] + h[j]
        dh = np.maximum(ei, ej) - zmax
        d = ei - ej
        q = np.zeros(self.n, dtype=h.dtype)
        pos = np.nonzero((dh > 0) & (d != 0))[0]
        if pos.size:
            ii = i[pos]
            val = (np.sign(d[pos]) * np.power(dh[pos], 5.0 / 3.0) *
                   np.sqrt(np.abs(d[pos]) / self.cell) / self.n_man)
            cap = self.cell / (4.0 * dt)
            hd = np.where(val > 0, h[ii], h[j[pos]])
            q[ii] = np.clip(val, -hd * cap, hd * cap)
        return q

    def step(self, dt, src=None):
        self._t_sim += dt          # single time accumulator for all sub-models
        z, h = self.z, self.h
        qE = self._flux(h, dt, self.iE, self.jE, self.zmaxE)
        qN = self._flux(h, dt, self.iN, self.jN, self.zmaxN)

        inflow = -qE
        inflow -= qN
        # each cell has at most one west / south neighbour -> plain fancy indexing
        inflow[self.iW] += qE[self.srcW]
        inflow[self.iS] += qN[self.srcS]

        h += dt / self.cell * inflow
        neg = h < 0
        if neg.any():
            self.clip_vol += float(-h[neg].sum()) * self.cell ** 2
            h[neg] = 0.0
        if src is not None:
            _sv = float(src.sum()) * dt
            self.src_vol += _sv
            h += dt * src / (self.cell ** 2)

        # gates: one-way outflow to a fixed outside stage
        for kk, gz, _v in self.gates:
            eta = z[kk] + h[kk]
            if eta > gz + 1e-6:
                dh = eta - max(z[kk], gz)
                if dh > 0:
                    q = (dh ** (5.0 / 3.0) * np.sqrt((eta - gz) / self.cell)
                         / self.n_man)
                    q = min(q, h[kk] * self.cell / (4.0 * dt))
                    v = min(q * dt / self.cell, h[kk] * self.cell)
                    h[kk] -= v / self.cell
                    self.out_vol += v * self.cell

        # lake -> sluice gate -> pump chamber -> 国康路 municipal sewer
        self._lake_drain(dt)

        # surface -> sewer return inlets (two-way coupling); the removed water
        # leaves the 2D domain, so it is booked in `ret_vol` like the gates.
        # A manhole can only accept water while its SEWER has SPARE capacity:
        # `in_ok` says whether the node was overflowing (a surcharged sewer
        # cannot take water back) and `in_cap` may vary per coupling interval,
        # carrying (full-flow capacity - current conduit flow).
        if self.in_k is not None and self.in_k.size:
            j = int(self._t_sim / self.in_step)
            hh = h[self.in_k]
            if self.in_head is not None:
                # HEAD-DRIVEN (C3): the manhole accepts water only while the
                # SURFACE water surface stands above the manhole water surface.
                #   dH = (rim + h_surface) - Z_node
                # State-only criterion -> no self-cancelling feedback; and a
                # sewer that fills up raises Z_node, so dH self-limits.
                if j < self.in_head.shape[1]:
                    dH = (self.in_zr + hh) - self.in_head[:, j]
                    act = (hh > 1e-3) & (dH > 1e-3)
                    if act.any():
                        idx = np.nonzero(act)[0]
                        q = (self.inlet_coef * self.in_area[idx]
                             * np.sqrt(2.0 * 9.81 * dH[idx]))
                        # the pipe must still have spare conveyance: a manhole
                        # whose sewer is full cannot pass surface water on, even
                        # though the head difference looks favourable
                        cap = (self.in_cap[idx, j] if self.in_cap.ndim == 2
                               else self.in_cap[idx])
                        q = np.minimum(q, cap)
                        v_ok = q > 0.0
                        idx, q = idx[v_ok], q[v_ok]
                        ha = hh[idx]
                        v = np.minimum(q * dt, ha * self.cell ** 2)
                        # np.add.at: duplicate inlets may share a cell; plain fancy-index
                        # subtraction keeps only the last value, under-removing
                        # water while still counting it in ret_vol
                        np.add.at(h, self.in_k[idx], -(v / (self.cell ** 2)))
                        self.ret_vol += float(v.sum())
                        for m, vv in zip(idx, v):
                            if vv > 0.0:
                                nid = self.ret_ids[m]
                                self._ret_acc[nid] = self._ret_acc.get(nid, 0.0) + float(vv)
            else:
                if self.in_ok is None:
                    ok = np.ones(self.in_k.size, bool)
                elif j < self.in_ok.shape[1]:
                    ok = self.in_ok[:, j]
                else:
                    ok = np.zeros(self.in_k.size, bool)
                act = (hh > 1e-3) & ok
                if act.any():
                    # operate on the active subset only: sqrt of a negative/NaN
                    # depth otherwise raises an "invalid value" warning每步
                    idx = np.nonzero(act)[0]
                    ha = hh[idx]
                    q = self.inlet_coef * self.in_area[idx] * np.sqrt(2.0 * 9.81 * ha)
                    cap = (self.in_cap[idx, j] if self.in_cap.ndim == 2
                           else self.in_cap[idx])
                    q = np.minimum(q, cap)
                    v = np.minimum(q * dt, ha * self.cell ** 2)
                    # np.add.at: duplicate inlets may share a cell; plain fancy-index
                    # subtraction keeps only the last value, under-removing
                    # water while still counting it in ret_vol
                    np.add.at(h, self.in_k[idx], -(v / (self.cell ** 2)))
                    self.ret_vol += float(v.sum())
                    for m, vv in zip(idx, v):
                        if vv > 0.0:
                            nid = self.ret_ids[m]
                            self._ret_acc[nid] = self._ret_acc.get(nid, 0.0) + float(vv)

        # NOTE: no depth thresholding here.  Zeroing h below `dry` looks
        # harmless but it silently destroys mass -- with dt = 0.25 s over
        # 7 200 steps it removed 32 % of the injected volume, and the loss
        # grows as dt shrinks.  Huge numbers of 1e-4 m cells are re-wetted and
        # re-zeroed every step.  Threshold only when reporting.
        return h

    # -------------------------------------------------------------- stepping
    def volume(self):
        return float(self.h.sum() * self.cell ** 2)

    def safe_dt(self):
        h = self.h
        if h.max() <= self.dry:
            return self.dt_max
        smax = 0.0
        for i, j in ((self.iE, self.jE), (self.iN, self.jN)):
            dh = h[i]
            if dh.size == 0:
                continue
            de = np.abs((self.z[i] + h[i]) - (self.z[j] + h[j]))
            smax = max(smax, float((de / self.cell).max()))
        if smax <= 0:
            return self.dt_max
        smax = min(smax, self.smax_cap)
        v = (1.0 / self.n_man) * float(h.max()) ** (2.0 / 3.0) * np.sqrt(smax)
        return float(np.clip(self.cfl * self.cell / max(v, 0.02),
                             self.dt_min, self.dt_max))

    def run(self, total_time, src_of=None, dt=None, record=None,
            record_every=60.0, verbose=False):
        t = 0.0
        nxt = record_every
        k = 0
        while t < total_time - 1e-9:
            if dt is None:
                if k % self.dt_update_every == 0:
                    self._dt = self.safe_dt()
                dt_ = self._dt
            else:
                dt_ = dt
            dt_ = min(dt_, total_time - t)
            src = src_of(t + 0.5 * dt_) if src_of is not None else None
            self.step(dt_, src)
            t += dt_
            k += 1
            if record is not None and t >= nxt - 1e-9:
                record(t, self)
                if self.inlets:
                    self.ret_hist.append([t, dict(self._ret_acc)])
                nxt += record_every
                if verbose:
                    print("  t=%7.1f  V=%10.1f m3  out=%10.1f  hmax=%.3f" % (
                        t, self.volume(), self.out_vol, self.h.max()), flush=True)
        return t
