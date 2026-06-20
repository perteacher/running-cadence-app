// Image y grows DOWNWARD. A foot at its LOWEST point (ground contact) is therefore a LOCAL MAXIMUM of raw y.

(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = factory();
  } else {
    root.CadenceCore = factory();
  }
}(typeof globalThis !== 'undefined' ? globalThis : this, function () {

  // --- Signal helpers ---

  function median(arr) {
    const s = arr.slice().sort((a, b) => a - b);
    const m = Math.floor(s.length / 2);
    return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
  }

  function std(arr) {
    const n = arr.length;
    if (n < 2) return 0;
    const mu = arr.reduce((a, b) => a + b, 0) / n;
    return Math.sqrt(arr.reduce((s, x) => s + (x - mu) ** 2, 0) / (n - 1));
  }

  // 2nd-order Butterworth biquad low-pass filter (RBJ cookbook).
  // Cutoff=7 Hz kills jitter while preserving 2.5–3.2 Hz step band.
  function lowpass(sig, fps, cutoffHz) {
    const f0 = cutoffHz / fps;
    const w0 = 2 * Math.PI * f0;
    const cosW = Math.cos(w0);
    const sinW = Math.sin(w0);
    const alpha = sinW / (2 * Math.SQRT2); // Q = sqrt(2)/2 => Butterworth
    const b0 = (1 - cosW) / 2, b1 = 1 - cosW, b2 = (1 - cosW) / 2;
    const a0 = 1 + alpha, a1 = -2 * cosW, a2 = 1 - alpha;
    const B0 = b0 / a0, B1 = b1 / a0, B2 = b2 / a0;
    const A1 = a1 / a0, A2 = a2 / a0;
    const out = new Array(sig.length);
    let x1 = 0, x2 = 0, y1 = 0, y2 = 0;
    for (let i = 0; i < sig.length; i++) {
      const x0 = sig[i];
      const y0 = B0 * x0 + B1 * x1 + B2 * x2 - A1 * y1 - A2 * y2;
      out[i] = y0;
      x2 = x1; x1 = x0; y2 = y1; y1 = y0;
    }
    return out;
  }

  // Centered moving-average detrend; removes slow vertical drift/panning.
  function detrend(sig, fps, winSec) {
    const half = Math.round(fps * winSec / 2);
    const n = sig.length;
    return sig.map((v, i) => {
      const lo = Math.max(0, i - half), hi = Math.min(n - 1, i + half);
      let sum = 0;
      for (let j = lo; j <= hi; j++) sum += sig[j];
      return v - sum / (hi - lo + 1);
    });
  }

  // Find local maxima with refractory minDistance and prominence threshold.
  function findPeaks(sig, { minDistance = 5, prominence = 0 } = {}) {
    const n = sig.length;
    const peaks = [];
    for (let i = 1; i < n - 1; i++) {
      if (sig[i] > sig[i - 1] && sig[i] > sig[i + 1]) {
        // prominence: min difference to the lowest valley on either side up to next higher peak
        const lMin = Math.min(...sig.slice(Math.max(0, i - minDistance * 3), i));
        const rMin = Math.min(...sig.slice(i + 1, Math.min(n, i + minDistance * 3 + 1)));
        const prom = sig[i] - Math.max(lMin, rMin);
        if (prom >= prominence) peaks.push(i);
      }
    }
    // Enforce refractory: keep only the highest peak within each minDistance window
    const kept = [];
    let lastKept = -Infinity;
    for (const p of peaks) {
      if (p - lastKept >= minDistance) {
        kept.push(p);
        lastKept = p;
      } else if (sig[p] > sig[kept[kept.length - 1]]) {
        kept[kept.length - 1] = p;
        lastKept = p;
      }
    }
    return kept;
  }

  // Normalized autocorrelation of a signal.
  function autocorr(sig) {
    const n = sig.length;
    const mu = sig.reduce((a, b) => a + b, 0) / n;
    const s = sig.map(x => x - mu);
    const denom = s.reduce((a, x) => a + x * x, 0);
    if (denom === 0) return new Array(n).fill(0);
    const ac = new Array(n).fill(0);
    for (let lag = 0; lag < n; lag++) {
      let sum = 0;
      for (let i = 0; i < n - lag; i++) sum += s[i] * s[i + lag];
      ac[lag] = sum / denom;
    }
    return ac;
  }

  // --- Gap filling by linear interpolation ---
  function interpolateGaps(sig, vis, visThreshold) {
    const out = sig.slice();
    const n = out.length;
    for (let i = 0; i < n; i++) {
      if (isNaN(out[i]) || (vis && vis[i] < visThreshold)) out[i] = NaN;
    }
    // forward fill edges, then linear interpolate interior gaps
    let first = 0;
    while (first < n && isNaN(out[first])) first++;
    if (first === n) return out;
    for (let i = 0; i < first; i++) out[i] = out[first];
    let last = n - 1;
    while (last >= 0 && isNaN(out[last])) last--;
    for (let i = last + 1; i < n; i++) out[i] = out[last];
    for (let i = 0; i < n; i++) {
      if (isNaN(out[i])) {
        let j = i + 1;
        while (j < n && isNaN(out[j])) j++;
        const vStart = out[i - 1], vEnd = out[j];
        for (let k = i; k < j; k++) {
          out[k] = vStart + (vEnd - vStart) * (k - i + 1) / (j - i + 1);
        }
        i = j - 1;
      }
    }
    return out;
  }

  // --- Main API ---

  function computeSPM(samples, options) {
    const {
      cutoffHz = 7,
      detrendWinSec = 1.0,
      minDistanceSec = 0.30,
      prominenceK = 0.3,
      visThreshold = 0.5
    } = options || {};

    const { leftAnkleY, rightAnkleY, hipMidY, leftVis, rightVis, tMs } = samples;
    const n = tMs.length;

    // 1. fps and duration
    const diffs = [];
    for (let i = 1; i < n; i++) diffs.push(tMs[i] - tMs[i - 1]);
    const fps = 1000 / median(diffs);
    const durationSec = (tMs[n - 1] - tMs[0]) / 1000;

    const minDist = Math.max(2, Math.round(fps * minDistanceSec));

    // 2. Interpolate missing frames
    function prepFoot(ankleY, vis) {
      const missing = ankleY.filter((v, i) => isNaN(v) || (vis && vis[i] < visThreshold)).length;
      const missingFrac = missing / n;
      const interp = interpolateGaps(ankleY, vis, visThreshold);
      return { interp, missingFrac };
    }

    const leftPrep  = prepFoot(leftAnkleY,  leftVis);
    const rightPrep = prepFoot(rightAnkleY, rightVis);

    // 3 & 4. Lowpass + detrend
    function process(sig) {
      const lp = lowpass(sig, fps, cutoffHz);
      return detrend(lp, fps, detrendWinSec);
    }

    const leftSig  = process(leftPrep.interp);
    const rightSig = process(rightPrep.interp);

    // 5. Find peaks (ground contact = local max in y-down convention).
    // Edge guard: skip first/last 5 frames to avoid filter-startup transients.
    const EDGE_GUARD = 5;
    function getPeaks(sig) {
      const trimmed = sig.slice(EDGE_GUARD, sig.length - EDGE_GUARD);
      const prom = prominenceK * std(trimmed);
      return findPeaks(trimmed, { minDistance: minDist, prominence: prom })
        .map(i => i + EDGE_GUARD); // restore original indices
    }

    const leftPeaks  = getPeaks(leftSig);
    const rightPeaks = getPeaks(rightSig);

    // 6. Peak-based SPM
    const totalSteps = leftPeaks.length + rightPeaks.length;
    const spmPeak = durationSec > 0 ? totalSteps / durationSec * 60 : 0;

    // Inter-peak CV for quality
    function intervalCV(peaks) {
      if (peaks.length < 3) return Infinity;
      const intervals = [];
      for (let i = 1; i < peaks.length; i++) intervals.push(peaks[i] - peaks[i - 1]);
      const m = intervals.reduce((a, b) => a + b, 0) / intervals.length;
      return m > 0 ? std(intervals) / m : Infinity;
    }

    const leftCV  = intervalCV(leftPeaks);
    const rightCV = intervalCV(rightPeaks);
    const combinedPeaks = [...leftPeaks, ...rightPeaks].sort((a, b) => a - b);
    const combinedCV = intervalCV(combinedPeaks);

    // 7. Autocorrelation fallback
    const useSig = leftPrep.missingFrac <= rightPrep.missingFrac ? leftSig : rightSig;
    const ac = autocorr(useSig);
    const minLag = Math.round(fps * 0.25); // shortest plausible stride ~240 spm
    const maxLag = Math.round(fps * 1.5);  // longest plausible stride ~40 spm
    const hiLag = Math.min(maxLag, ac.length - 1);
    let bestLag = minLag, bestAC = -Infinity;
    for (let lag = minLag; lag <= hiLag; lag++) {
      if (ac[lag] > bestAC) { bestAC = ac[lag]; bestLag = lag; }
    }
    // Parabolic (3-point) interpolation of the autocorr peak -> sub-sample lag.
    // Integer lags quantize cadence (e.g. ~189.5 vs 200 SPM across the lag boundary),
    // so the same clip can jump between runs; the refined lag smooths that out.
    // Skip at the search-window edges where a 3-point fit isn't defined.
    let refinedLag = bestLag;
    if (bestLag > minLag && bestLag < hiLag) {
      const am1 = ac[bestLag - 1], a0 = ac[bestLag], ap1 = ac[bestLag + 1];
      const denom = am1 - 2 * a0 + ap1;
      let delta = denom !== 0 ? 0.5 * (am1 - ap1) / denom : 0;
      delta = Math.max(-1, Math.min(1, delta));
      refinedLag = bestLag + delta;
    }
    const spmAuto = refinedLag > 0 ? (fps / refinedLag) * 60 * 2 : 0;

    // 8. Choose method — cross-validate peak detection against autocorrelation
    //    to catch double-counting (spurious peaks => spmPeak ~2x the true cadence).
    const PLAUSIBLE_LO = 140, PLAUSIBLE_HI = 205;   // valid running cadence band
    const peakInBand = spmPeak >= PLAUSIBLE_LO && spmPeak <= PLAUSIBLE_HI;
    const autoInBand = spmAuto >= PLAUSIBLE_LO && spmAuto <= PLAUSIBLE_HI;
    const enoughPeaks = totalSteps >= 6;
    const regularSpacing = combinedCV < 0.30;

    // Autocorr locks onto the fundamental period, so it resists double-counting.
    // If spmPeak disagrees with it by >12% (≈2x => classic double-count), distrust peaks.
    const peakAutoDev = spmAuto > 0 ? Math.abs(spmPeak - spmAuto) / spmAuto : Infinity;
    const agreesAuto = peakAutoDev <= 0.12;

    // Heavy left/right peak-count asymmetry (>1.5x) also signals spurious peaks.
    const lpc = leftPeaks.length, rpc = rightPeaks.length;
    const peakRatio = (lpc > 0 && rpc > 0) ? Math.max(lpc, rpc) / Math.min(lpc, rpc) : Infinity;
    const balanced = peakRatio <= 1.5;

    let method, spmRaw, confidence, recheck = false;

    // Trust peaks only when ALL hold: enough + regular + in-band + agrees-autocorr + balanced.
    if (enoughPeaks && regularSpacing && peakInBand && agreesAuto && balanced) {
      method = 'peak';
      spmRaw = spmPeak;
      confidence = combinedCV < 0.15 ? 'high' : 'medium';
    } else if (autoInBand) {
      // Peaks untrustworthy but autocorr is plausible -> adopt autocorr.
      method = 'autocorr';
      spmRaw = spmAuto;
      const heavilyInterp = leftPrep.missingFrac > 0.4 && rightPrep.missingFrac > 0.4;
      confidence = heavilyInterp ? 'low' : 'medium';
    } else if (peakInBand && enoughPeaks && regularSpacing) {
      // Autocorr implausible (wrong harmonic) but peaks are sane -> fall back to peaks, low conf.
      method = 'peak';
      spmRaw = spmPeak;
      confidence = 'low';
    } else {
      // Neither estimate lands in the plausible band -> flag for re-measure.
      method = 'autocorr';
      spmRaw = spmAuto;
      confidence = 'low';
      recheck = true;
    }

    // 9. Round
    const spm = Math.round(spmRaw);

    return {
      spm,
      confidence,
      method,
      recheck,
      stepCount: totalSteps,
      durationSec,
      fps,
      peaks: { left: leftPeaks, right: rightPeaks },
      debug: {
        fps, spmPeak, spmAuto, combinedCV, leftCV, rightCV,
        peakAutoDev, peakRatio,
        leftPeakCount: leftPeaks.length, rightPeakCount: rightPeaks.length,
        leftMissingFrac: leftPrep.missingFrac, rightMissingFrac: rightPrep.missingFrac,
        bestLag, bestAC
      }
    };
  }

  // --- Overstride analysis ---

  function computeOverstride(samples, peaks, dims) {
    const { w, h } = dims;
    const {
      leftAnkleX, rightAnkleX,
      leftAnkleY, rightAnkleY,
      leftKneeX, leftKneeY,
      rightKneeX, rightKneeY,
      hipMidX,
      leftVis, rightVis,
      leftKneeVis, rightKneeVis
    } = samples;

    // Running direction from the overall hipX trend (mean per-frame change).
    // forward = +x for zero/rightward travel, -x for leftward. (research §1-4)
    function meanDiff(arr) {
      let sum = 0, cnt = 0, prev = null;
      if (!arr) return 0;
      for (const v of arr) {
        if (isFinite(v)) { if (prev !== null) { sum += v - prev; cnt++; } prev = v; }
      }
      return cnt > 0 ? sum / cnt : 0;
    }
    const dir = meanDiff(hipMidX) < 0 ? -1 : 1;

    // Per-side stride length (frames) -> search window for the initial-contact edge.
    function strideWin(idxArr) {
      if (!idxArr || idxArr.length < 2) return 10;
      const gaps = [];
      for (let i = 1; i < idxArr.length; i++) gaps.push(idxArr[i] - idxArr[i - 1]);
      const g = median(gaps);
      return Math.max(4, Math.min(20, Math.round(g * 0.6)));
    }

    // Touchdown via velocity zero-crossing: follow the descent (foot falling, y rising) and
    // stop where it halts -> ground contact. More precise than a fixed % threshold.
    function touchdownIdx(yArr, p, win) {
      if (!Number.isFinite(yArr && yArr[p])) return p;
      const lo = Math.max(1, p - win);
      let vmaxIdx = lo, vmax = -Infinity;
      for (let i = lo; i <= p; i++) {
        if (Number.isFinite(yArr[i]) && Number.isFinite(yArr[i - 1])) {
          const v = yArr[i] - yArr[i - 1];     // y-down: descending => v>0
          if (v > vmax) { vmax = v; vmaxIdx = i; }
        }
      }
      if (vmax <= 0) return p;                  // no clear descent -> fall back to peak
      const HALT = 0.25 * vmax;                 // velocity fallen to 25% of peak = contact settled
      let i = vmaxIdx;
      while (i < p && Number.isFinite(yArr[i + 1]) && (yArr[i + 1] - yArr[i]) > HALT) i++;
      return i;
    }

    // --- Per-frame signed shin angle at one index (or null if coords missing) ---
    function shinAt(k, ankleXArr, ankleYArr, kneeXArr, kneeYArr, kneeVisArr) {
      const ax = ankleXArr && ankleXArr[k], ay = ankleYArr && ankleYArr[k];
      const kx = kneeXArr  && kneeXArr[k],  ky = kneeYArr  && kneeYArr[k];
      const hx = hipMidX   && hipMidX[k],   vis = kneeVisArr && kneeVisArr[k];
      if (!Number.isFinite(ax) || !Number.isFinite(ay) || !Number.isFinite(kx) ||
          !Number.isFinite(ky) || !Number.isFinite(hx)) return null;
      const axpx = ax * w, aypx = ay * h, kxpx = kx * w, kypx = ky * h, hxpx = hx * w;
      // Signed: only the forward component (foot ahead of knee in travel direction) counts.
      // Propulsion posture (knee ahead of foot) clamps to 0 -- not overstride.
      const fwd    = Math.max(0, dir * (axpx - kxpx));
      const shin   = Math.atan2(fwd, Math.abs(aypx - kypx)) * 180 / Math.PI;
      const offset = dir * (axpx - hxpx);
      return { shin, offset, vis: Number.isFinite(vis) ? vis : 0 };
    }

    // Overstride is EXPERIMENTAL on 30fps side-view 2D, so gate hard against noise:
    //  - measure only the NEAR leg (higher ankle+knee visibility); the occluded far leg is dropped
    //  - per stride: median of a short touchdown window, high-visibility frames only
    //  - discard shin > 35deg (impossible at contact => swing frame / broken keypoint)
    const MAX_SHIN = 35, VIS_MIN = 0.5, WIN_AHEAD = 3;

    function sideVis(idxArr, ankleVisArr, kneeVisArr) {
      let sum = 0, n = 0;
      for (const p of idxArr) {
        const av = ankleVisArr && ankleVisArr[p], kv = kneeVisArr && kneeVisArr[p];
        const parts = [];
        if (Number.isFinite(av)) parts.push(av);
        if (Number.isFinite(kv)) parts.push(kv);
        if (parts.length) { sum += parts.reduce((a, b) => a + b, 0) / parts.length; n++; }
      }
      return n ? sum / n : 0;
    }

    const strides = [];

    function addSide(foot, idxArr, ankleXArr, ankleYArr, kneeXArr, kneeYArr, kneeVisArr) {
      const win = strideWin(idxArr);
      for (const p of idxArr) {
        const td = touchdownIdx(ankleYArr, p, win);
        const windowAngles = [];
        let offSum = 0, offN = 0, visSum = 0, visN = 0;
        for (let j = td; j <= td + WIN_AHEAD; j++) {
          const r = shinAt(j, ankleXArr, ankleYArr, kneeXArr, kneeYArr, kneeVisArr);
          if (r && r.vis >= VIS_MIN && r.shin <= MAX_SHIN) {  // high-vis + physically-possible only
            windowAngles.push(r.shin);
            offSum += r.offset; offN++;
            visSum += r.vis;    visN++;
          }
        }
        if (windowAngles.length === 0) continue;
        const aytd = (ankleYArr && Number.isFinite(ankleYArr[td])) ? Math.round(ankleYArr[td] * 1000) / 1000 : null;
        strides.push({
          foot,
          peakFrame: p,
          touchdownFrame: td,
          ankleY_at_td: aytd,
          windowAngles: windowAngles.map(a => Math.round(a * 10) / 10),
          strideMedian: Math.round(median(windowAngles) * 10) / 10,
          offset: offN ? offSum / offN : 0,
          vis: visN ? visSum / visN : 0
        });
      }
    }

    // Near leg = the side with higher mean (ankle+knee) visibility across its strides.
    const leftSideVis  = sideVis(peaks.left,  leftVis,  leftKneeVis);
    const rightSideVis = sideVis(peaks.right, rightVis, rightKneeVis);
    const nearFoot = (rightSideVis > leftSideVis) ? 'R' : 'L';
    if (nearFoot === 'L') addSide('L', peaks.left,  leftAnkleX,  leftAnkleY,  leftKneeX,  leftKneeY,  leftKneeVis);
    else                  addSide('R', peaks.right, rightAnkleX, rightAnkleY, rightKneeX, rightKneeY, rightKneeVis);

    if (strides.length === 0) {
      return { shinAngleDeg: null, horizOffset: null, classification: 'unknown',
               confidence: 'low', strikeCount: 0, strideCount: 0, nearFoot,
               debug: { strides: [], finalMedian: null, nearFoot } };
    }

    const strideMedians = strides.map(s => s.strideMedian);
    const offsets       = strides.map(s => s.offset);
    const visArr        = strides.map(s => s.vis);

    const shinAngleDeg = Math.round(median(strideMedians) * 10) / 10;
    const horizOffset  = Math.round(median(offsets));
    const strikeCount  = strides.length;
    const avgVis       = visArr.reduce((a, b) => a + b, 0) / visArr.length;

    let classification;
    if (shinAngleDeg < 10)       classification = 'good';
    else if (shinAngleDeg < 15)  classification = 'borderline';
    else                         classification = 'overstride';

    // Experimental metric: trust only with enough clean strides; too few (<6) -> low.
    let confidence;
    if (strikeCount >= 6 && avgVis >= 0.6)      confidence = 'high';
    else if (strikeCount >= 4 && avgVis >= 0.4) confidence = 'medium';
    else                                         confidence = 'low';
    if (strikeCount < 6) confidence = 'low';

    return {
      shinAngleDeg, horizOffset, classification, confidence,
      strikeCount, strideCount: strikeCount, nearFoot,
      debug: { strides, finalMedian: shinAngleDeg, nearFoot }
    };
  }

  // --- Trunk forward lean (shoulder-mid -> hip-mid vs vertical) ---
  // Reliable on 30fps side-view 2D (large sagittal angle, no contact-frame dependency):
  // take the MEDIAN over every frame with good shoulder+hip visibility. Deterministic.
  function computeTrunkLean(samples, dims) {
    const { w, h } = dims;
    const { shoulderMidX, shoulderMidY, hipMidX, hipMidY, shoulderVis, hipVis } = samples;
    const n = (hipMidY && hipMidY.length) || (shoulderMidY && shoulderMidY.length) || 0;

    function meanDiff(arr) {
      let sum = 0, cnt = 0, prev = null;
      if (!arr) return 0;
      for (const v of arr) { if (Number.isFinite(v)) { if (prev !== null) { sum += v - prev; cnt++; } prev = v; } }
      return cnt > 0 ? sum / cnt : 0;
    }
    const dir = meanDiff(hipMidX) < 0 ? -1 : 1;   // forward = travel direction

    const VIS = 0.5;
    const angles = [];
    for (let i = 0; i < n; i++) {
      const sx = shoulderMidX && shoulderMidX[i], sy = shoulderMidY && shoulderMidY[i];
      const hx = hipMidX && hipMidX[i],           hy = hipMidY && hipMidY[i];
      const sv = shoulderVis && shoulderVis[i],   hv = hipVis && hipVis[i];
      if (!Number.isFinite(sx) || !Number.isFinite(sy) || !Number.isFinite(hx) || !Number.isFinite(hy)) continue;
      if (Number.isFinite(sv) && sv < VIS) continue;
      if (Number.isFinite(hv) && hv < VIS) continue;
      // pixel-space shoulder->hip vector; angle from vertical. Sign: + leaning forward (travel dir).
      const fwd  = dir * (sx - hx) * w;
      const vert = Math.abs((sy - hy) * h);
      if (vert < 1) continue;   // skip degenerate frames (shoulder/hip same height) => no ±90° outlier
      angles.push(Math.atan2(fwd, vert) * 180 / Math.PI);
    }

    if (angles.length < 3) {
      return { leanDeg: null, classification: 'unknown', confidence: 'low', frameCount: angles.length };
    }
    const leanDeg = Math.round(median(angles) * 10) / 10;

    // Dead-zone: an upright torso (small +/- angle) is NORMAL. Only a clear backward
    // lean is flagged. This is torso lean (shoulder->hip), NOT whole-body lean-from-ankle.
    let classification;
    if (leanDeg < -3)       classification = 'back';        // clearly leaning back -> nudge forward
    else if (leanDeg <= 12) classification = 'good';        // -3..+12 = upright/natural torso, fine
    else                    classification = 'excessive';   // >12 over-lean / hip-flexion

    let confidence;
    if (angles.length >= 30)      confidence = 'high';
    else if (angles.length >= 12) confidence = 'medium';
    else                          confidence = 'low';

    return { leanDeg, classification, confidence, frameCount: angles.length };
  }

  // --- Vertical oscillation (hip bounce), scale-free ---
  // amplitude = robust peak-to-peak (p95-p5) of the DETRENDED hipMidY, in pixels,
  // normalized by torso length (|shoulder-hip|) => % independent of camera distance/zoom.
  // Whole-clip based => deterministic & robust.
  function computeVerticalOscillation(samples, dims) {
    const { w, h } = dims;
    const { hipMidX, hipMidY, shoulderMidX, shoulderMidY, hipVis, shoulderVis, tMs } = samples;
    const n = (hipMidY && hipMidY.length) || 0;

    // fps from timestamps (fallback 30)
    let fps = 30;
    if (tMs && tMs.length > 2) {
      const diffs = [];
      for (let i = 1; i < tMs.length; i++) { const d = tMs[i] - tMs[i - 1]; if (d > 0) diffs.push(d); }
      if (diffs.length) fps = 1000 / median(diffs);
    }

    const VIS = 0.5;
    const hipSeries = [], torso = [];
    for (let i = 0; i < n; i++) {
      const hy = hipMidY && hipMidY[i], hx = hipMidX && hipMidX[i];
      const hv = hipVis && hipVis[i];
      if (!Number.isFinite(hy)) continue;
      if (Number.isFinite(hv) && hv < VIS) continue;
      hipSeries.push(hy);
      const sy = shoulderMidY && shoulderMidY[i], sx = shoulderMidX && shoulderMidX[i];
      const sv = shoulderVis && shoulderVis[i];
      if (Number.isFinite(sy) && Number.isFinite(sx) && Number.isFinite(hx) &&
          (!Number.isFinite(sv) || sv >= VIS)) {
        const dxpx = (sx - hx) * w, dypx = (sy - hy) * h;
        torso.push(Math.sqrt(dxpx * dxpx + dypx * dypx));
      }
    }

    if (hipSeries.length < 10 || torso.length < 3) {
      return { ratioPct: null, ampPx: null, torsoPx: null, classification: 'unknown', confidence: 'low', frameCount: hipSeries.length };
    }

    // detrend out slow vertical drift / panning, then robust peak-to-peak amplitude.
    const detr = detrend(hipSeries, fps, 1.0);
    function pct(arr, p) {
      const s = arr.slice().sort((a, b) => a - b);
      const idx = Math.min(s.length - 1, Math.max(0, Math.round((p / 100) * (s.length - 1))));
      return s[idx];
    }
    const ampPx   = (pct(detr, 95) - pct(detr, 5)) * h;
    const torsoPx = median(torso);
    const ratioPct = torsoPx > 0 ? Math.round((ampPx / torsoPx) * 1000) / 10 : null;

    if (ratioPct == null) {
      return { ratioPct: null, ampPx: null, torsoPx: null, classification: 'unknown', confidence: 'low', frameCount: hipSeries.length };
    }

    // Provisional bands (to be calibrated on real data, like the cadence intercept).
    let classification;
    if (ratioPct < 16)       classification = 'low';     // economical
    else if (ratioPct <= 22) classification = 'normal';
    else                     classification = 'high';    // bouncy

    let confidence;
    if (hipSeries.length >= 60)      confidence = 'high';
    else if (hipSeries.length >= 20) confidence = 'medium';
    else                             confidence = 'low';

    return {
      ratioPct, ampPx: Math.round(ampPx), torsoPx: Math.round(torsoPx),
      classification, confidence, frameCount: hipSeries.length
    };
  }

  const CadenceCore = { computeSPM, computeOverstride, computeTrunkLean, computeVerticalOscillation, lowpass, detrend, findPeaks, autocorr };

  // --- Self-test (Node only) ---
  if (typeof require !== 'undefined' && require.main === module) {
    // Seeded LCG PRNG
    function makePRNG(seed) {
      let s = seed >>> 0;
      return function () {
        s = (Math.imul(1664525, s) + 1013904223) >>> 0;
        return s / 4294967296;
      };
    }
    function gaussNoise(rng, sigma) {
      // Box-Muller
      const u = 1 - rng(), v = rng();
      return sigma * Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
    }
    function buildSamples(n, fps, strideHz, noiseS, driftAmp, leftNaNFrac, rng) {
      const PI2 = 2 * Math.PI;
      const leftAnkleY = [], rightAnkleY = [], hipMidY = [];
      const leftVis = [], rightVis = [], tMs = [];
      for (let i = 0; i < n; i++) {
        const t = i / fps;
        tMs.push(t * 1000);
        const noise = noiseS > 0 ? gaussNoise(rng, noiseS) : 0;
        const drift = driftAmp * (t / (n / fps));
        leftAnkleY.push(0.8 + 0.05 * Math.cos(PI2 * strideHz * t) + noise + drift);
        rightAnkleY.push(0.8 + 0.05 * Math.cos(PI2 * strideHz * t + Math.PI) + noise + drift);
        hipMidY.push(0.5);
        leftVis.push(1); rightVis.push(1);
      }
      if (leftNaNFrac > 0) {
        const nBad = Math.floor(n * leftNaNFrac);
        const step = Math.floor(n / nBad);
        for (let i = 0; i < nBad; i++) {
          const idx = i * step;
          leftVis[idx] = 0;
        }
      }
      return { leftAnkleY, rightAnkleY, hipMidY, leftVis, rightVis, tMs };
    }

    const FPS = 30, DUR = 10, N = FPS * DUR;
    const rng = makePRNG(42);

    const tests = [
      { label: '1', samples: buildSamples(N, FPS, 1.5, 0,    0,   0,   rng), expect: 180, tol: 5  },
      { label: '2', samples: buildSamples(N, FPS, 1.4, 0,    0,   0,   rng), expect: 168, tol: 5  },
      { label: '3', samples: buildSamples(N, FPS, 1.5, 0.01, 0.1, 0,   rng), expect: 180, tol: 6  },
      { label: '4', samples: buildSamples(N, FPS, 1.5, 0,    0,   0.3, rng), expect: 180, tol: 8  },
    ];

    for (const t of tests) {
      const r = CadenceCore.computeSPM(t.samples);
      const ok = Math.abs(r.spm - t.expect) <= t.tol;
      console.log(`${ok ? 'PASS' : 'FAIL'} Test${t.label}: got ${r.spm} expected ${t.expect}±${t.tol} method=${r.method} conf=${r.confidence}`);
    }

    // --- Autocorr sub-sample (parabolic) interpolation ---
    // A cadence whose true lag falls BETWEEN two integers must not snap to either.
    // period=18.5 frames -> 194.6 SPM (between lag18=200 and lag19=189.5); the refined
    // lag should land on the true value (integer-only autocorr is ~5 SPM off here).
    // A genuine ~190 case must also stay ~190 -> adjacent lags both read correctly.
    {
      function pureAnkle(n, fps, periodFrames) {
        const lA = [], rA = [], lV = [], rV = [], ts = [];
        for (let i = 0; i < n; i++) {
          const ph = 2 * Math.PI * i / periodFrames;
          lA.push(0.8 + 0.05 * Math.cos(ph));
          rA.push(0.8 + 0.05 * Math.cos(ph + Math.PI));
          lV.push(1); rV.push(1); ts.push(i / fps * 1000);
        }
        return { leftAnkleY: lA, rightAnkleY: rA, hipMidY: new Array(n).fill(0.5), leftVis: lV, rightVis: rV, tMs: ts };
      }
      const fps = 30;
      const trueMid = (fps / 18.5) * 120;                                   // 194.6
      const rMid = CadenceCore.computeSPM(pureAnkle(300, fps, 18.5));
      const r190 = CadenceCore.computeSPM(pureAnkle(300, fps, 3600 / 190)); // period 18.947 -> 190
      const okMid = Math.abs(rMid.debug.spmAuto - trueMid) <= 2;
      const ok190 = Math.abs(r190.debug.spmAuto - 190) <= 2;
      const ok = okMid && ok190;
      console.log(`${ok ? 'PASS' : 'FAIL'} AC-interp: between-lag=${Math.round(rMid.debug.spmAuto * 10) / 10}(exp≈${Math.round(trueMid * 10) / 10}) ~190=${Math.round(r190.debug.spmAuto * 10) / 10}`);
    }

    // --- Overstride tests ---
    function makeOSsamples(n, leftAnkleX, rightAnkleX, leftKneeX, leftKneeY, rightKneeX, rightKneeY, hipMidX, leftAnkleY, rightAnkleY, leftKneeVis, rightKneeVis) {
      // Build flat arrays of length n filled with constant values
      function fill(v) { return new Array(n).fill(v); }
      return {
        leftAnkleX:   fill(leftAnkleX),
        rightAnkleX:  fill(rightAnkleX),
        leftAnkleY:   fill(leftAnkleY),
        rightAnkleY:  fill(rightAnkleY),
        leftKneeX:    fill(leftKneeX),
        leftKneeY:    fill(leftKneeY),
        rightKneeX:   fill(rightKneeX),
        rightKneeY:   fill(rightKneeY),
        hipMidX:      fill(hipMidX),
        leftKneeVis:  fill(leftKneeVis),
        rightKneeVis: fill(rightKneeVis),
      };
    }

    // (a) Vertical shin: ankle directly below knee (same X), kneeVis 0.9
    {
      const idxs = [10, 20, 30, 40];
      const s = makeOSsamples(50, 0.5, 0.5, 0.5, 0.4, 0.5, 0.4, 0.5, 0.6, 0.6, 0.9, 0.9);
      const peaks = { left: idxs, right: [] };
      const r = CadenceCore.computeOverstride(s, peaks, { w: 1000, h: 1000 });
      const ok = r.shinAngleDeg !== null && Math.abs(r.shinAngleDeg) <= 1 && r.classification === 'good';
      console.log(`${ok ? 'PASS' : 'FAIL'} OSTest-a: vertical shin shinAngleDeg=${r.shinAngleDeg} class=${r.classification}`);
    }

    // (b) Known 30 deg + enough strides -> high confidence. knee(0.5,0.4), ankle(0.61547,0.6):
    //     dx=115.47px, dy=200px -> atan2(115.47,200)=30deg. (45deg would now be discarded as >35.)
    {
      const idxs = [8, 14, 20, 26, 32, 38];   // 6 strides => 'high' eligible
      const s = makeOSsamples(44, 0.61547, 0.61547, 0.5, 0.4, 0.5, 0.4, 0.5, 0.6, 0.6, 0.9, 0.9);
      const peaks = { left: idxs, right: [] };
      const r = CadenceCore.computeOverstride(s, peaks, { w: 1000, h: 1000 });
      const ok = r.shinAngleDeg !== null && Math.abs(r.shinAngleDeg - 30) <= 1 && r.confidence === 'high';
      console.log(`${ok ? 'PASS' : 'FAIL'} OSTest-b: 30deg shinAngleDeg=${r.shinAngleDeg} conf=${r.confidence} strides=${r.strideCount}`);
    }

    // (c) Aspect-ratio correctness: same normalized coords, portrait vs landscape -> DIFFERENT pixel angles
    {
      // knee(0.5,0.4) ankle(0.6,0.6) with portrait w=720 h=1280:
      //   kxpx=360,kypx=512; axpx=432,aypx=768 -> dx=72,dy=256 -> atan2(72,256)
      // with landscape w=1280 h=720:
      //   kxpx=640,kypx=288; axpx=768,aypx=432 -> dx=128,dy=144 -> atan2(128,144)
      const expectedPortrait  = Math.atan2(Math.abs(0.55*720 - 0.5*720), Math.abs(0.6*1280 - 0.4*1280)) * 180 / Math.PI;
      const expectedLandscape = Math.atan2(Math.abs(0.55*1280 - 0.5*1280), Math.abs(0.6*720 - 0.4*720)) * 180 / Math.PI;
      const idxs = [5, 10, 15, 20];
      const s = makeOSsamples(30, 0.55, 0.55, 0.5, 0.4, 0.5, 0.4, 0.5, 0.6, 0.6, 0.9, 0.9);
      const peaks = { left: idxs, right: [] };
      const rP = CadenceCore.computeOverstride(s, peaks, { w: 720,  h: 1280 });
      const rL = CadenceCore.computeOverstride(s, peaks, { w: 1280, h: 720  });
      const okP = Math.abs(rP.shinAngleDeg - Math.round(expectedPortrait  * 10) / 10) <= 0.15;
      const okL = Math.abs(rL.shinAngleDeg - Math.round(expectedLandscape * 10) / 10) <= 0.15;
      const okDiff = Math.abs(rP.shinAngleDeg - rL.shinAngleDeg) > 1; // must differ
      const ok = okP && okL && okDiff;
      console.log(`${ok ? 'PASS' : 'FAIL'} OSTest-c: portrait=${rP.shinAngleDeg}(exp≈${Math.round(expectedPortrait*10)/10}) landscape=${rL.shinAngleDeg}(exp≈${Math.round(expectedLandscape*10)/10}) differ=${okDiff}`);
    }

    // (d) Classification bands: craft normalized coords to yield median ~12 and ~20 degrees pixel-space
    {
      // dims 1000x1000 (square, so normalized == pixel proportions).
      // knee at (0.5, 0.4) -> kxpx=500, kypx=400.
      // Want ankle shin = angle: atan2(|axpx-kxpx|, |aypx-kypx|) = targetDeg
      // Fix aypx = 600 -> dypx = 200. Then dxpx = dypx * tan(targetDeg).
      // ankleX = (kxpx + dxpx) / 1000, ankleY = 600/1000 = 0.6
      const W = 1000, H = 1000;
      const kxpx = 500, kypx = 400, aypx = 600;
      const dypx = aypx - kypx; // 200
      function ankleXNorm(deg) { return (kxpx + dypx * Math.tan(deg * Math.PI / 180)) / W; }
      const ax12 = ankleXNorm(12);
      const ax20 = ankleXNorm(20);
      const idxs = [5, 10, 15, 20];

      const s12 = makeOSsamples(30, ax12, ax12, 0.5, 0.4, 0.5, 0.4, 0.5, 0.6, 0.6, 0.9, 0.9);
      const s20 = makeOSsamples(30, ax20, ax20, 0.5, 0.4, 0.5, 0.4, 0.5, 0.6, 0.6, 0.9, 0.9);
      const peaks = { left: idxs, right: [] };
      const r12 = CadenceCore.computeOverstride(s12, peaks, { w: W, h: H });
      const r20 = CadenceCore.computeOverstride(s20, peaks, { w: W, h: H });
      const ok12 = r12.classification === 'borderline';
      const ok20 = r20.classification === 'overstride';
      const ok = ok12 && ok20;
      console.log(`${ok ? 'PASS' : 'FAIL'} OSTest-d: shin12=${r12.shinAngleDeg}=>${r12.classification} shin20=${r20.shinAngleDeg}=>${r20.classification}`);
    }

    // --- Double-count guard: spurious peaks at 2x must NOT inflate cadence ---
    {
      // True stride 1.5 Hz (=>180 SPM) plus a strong 2nd harmonic, so the peak detector
      // finds ~2 peaks per stride. The guard should distrust peaks and adopt autocorr.
      const PI2 = 2 * Math.PI, fps = 30, nn = 300, f = 1.5;
      const lA = [], rA = [], lV = [], rV = [], ts = [];
      for (let i = 0; i < nn; i++) {
        const t = i / fps;
        ts.push(t * 1000);
        lA.push(0.8 + 0.05 * Math.cos(PI2 * f * t)            + 0.045 * Math.cos(PI2 * 2 * f * t));
        rA.push(0.8 + 0.05 * Math.cos(PI2 * f * t + Math.PI)  + 0.045 * Math.cos(PI2 * 2 * f * t + Math.PI));
        lV.push(1); rV.push(1);
      }
      const s = { leftAnkleY: lA, rightAnkleY: rA, hipMidY: new Array(nn).fill(0.5), leftVis: lV, rightVis: rV, tMs: ts };
      const r = CadenceCore.computeSPM(s);
      const guarded = r.spm >= 165 && r.spm <= 195 && r.method === 'autocorr';
      console.log(`${guarded ? 'PASS' : 'FAIL'} DoubleCountGuard: spm=${r.spm} method=${r.method} spmPeak≈${Math.round(r.debug.spmPeak)} spmAuto≈${Math.round(r.debug.spmAuto)} conf=${r.confidence}`);
    }

    // --- Overstride: initial-contact + signed shin angle ---
    function osArrays(n, obj) {
      const o = {};
      for (const k in obj) o[k] = new Array(n).fill(obj[k]);
      return o;
    }

    // (e) Initial contact, foot AHEAD of knee in travel direction -> positive angle.
    {
      const n = 30, idxs = [5, 12, 19, 26];
      const s = osArrays(n, {
        leftAnkleX: 0.62, rightAnkleX: 0.62, leftAnkleY: 0.62, rightAnkleY: 0.62,
        leftKneeX: 0.50, leftKneeY: 0.42, rightKneeX: 0.50, rightKneeY: 0.42,
        leftKneeVis: 0.9, rightKneeVis: 0.9
      });
      s.hipMidX = []; for (let i = 0; i < n; i++) s.hipMidX.push(0.30 + 0.01 * i); // rightward trend -> dir=+1
      const r = CadenceCore.computeOverstride(s, { left: idxs, right: [] }, { w: 1000, h: 1000 });
      const ok = r.shinAngleDeg > 5 && r.classification !== 'unknown';
      console.log(`${ok ? 'PASS' : 'FAIL'} OSTest-e: IC foot-ahead shin=${r.shinAngleDeg} class=${r.classification}`);
    }

    // (f) Propulsion posture: knee AHEAD of foot -> clamps to 0 (not overstride).
    {
      const n = 30, idxs = [5, 12, 19, 26];
      const s = osArrays(n, {
        leftAnkleX: 0.40, rightAnkleX: 0.40, leftAnkleY: 0.62, rightAnkleY: 0.62,  // ankle BEHIND knee
        leftKneeX: 0.55, leftKneeY: 0.42, rightKneeX: 0.55, rightKneeY: 0.42,
        leftKneeVis: 0.9, rightKneeVis: 0.9
      });
      s.hipMidX = []; for (let i = 0; i < n; i++) s.hipMidX.push(0.30 + 0.01 * i); // rightward -> dir=+1
      const r = CadenceCore.computeOverstride(s, { left: idxs, right: [] }, { w: 1000, h: 1000 });
      const ok = r.shinAngleDeg === 0 && r.classification === 'good';
      console.log(`${ok ? 'PASS' : 'FAIL'} OSTest-f: propulsion shin=${r.shinAngleDeg} class=${r.classification}`);
    }

    // (g) Direction-agnostic: foot-ahead overstride detected for BOTH L->R and R->L travel.
    {
      const n = 30, idxs = [5, 12, 19, 26];
      const sR = osArrays(n, {              // rightward: dir=+1, foot ahead => ankleX > kneeX (30deg)
        leftAnkleX: 0.61547, rightAnkleX: 0.61547, leftAnkleY: 0.62, rightAnkleY: 0.62,
        leftKneeX: 0.50, leftKneeY: 0.42, rightKneeX: 0.50, rightKneeY: 0.42,
        leftKneeVis: 0.9, rightKneeVis: 0.9
      });
      sR.hipMidX = []; for (let i = 0; i < n; i++) sR.hipMidX.push(0.20 + 0.01 * i);
      const sL = osArrays(n, {              // leftward: dir=-1, foot ahead => ankleX < kneeX (mirror, 30deg)
        leftAnkleX: 0.38453, rightAnkleX: 0.38453, leftAnkleY: 0.62, rightAnkleY: 0.62,
        leftKneeX: 0.50, leftKneeY: 0.42, rightKneeX: 0.50, rightKneeY: 0.42,
        leftKneeVis: 0.9, rightKneeVis: 0.9
      });
      sL.hipMidX = []; for (let i = 0; i < n; i++) sL.hipMidX.push(0.80 - 0.01 * i);
      const rR = CadenceCore.computeOverstride(sR, { left: idxs, right: [] }, { w: 1000, h: 1000 });
      const rL = CadenceCore.computeOverstride(sL, { left: idxs, right: [] }, { w: 1000, h: 1000 });
      const ok = Math.abs(rR.shinAngleDeg - rL.shinAngleDeg) <= 0.2 &&
                 rR.classification === 'overstride' && rL.classification === 'overstride';
      console.log(`${ok ? 'PASS' : 'FAIL'} OSTest-g: R-travel=${rR.shinAngleDeg}(${rR.classification}) L-travel=${rL.shinAngleDeg}(${rL.classification})`);
    }

    // --- Touchdown window + per-stride median robustness (the real-jitter fix) ---
    // Realistic foot trajectory: stance plateau domed at midstance (foot lowest, shin~vertical=0),
    // swing lift, and the foot reaches FORWARD at touchdown -> rotates back through midstance.
    function genRun(n, fps, period, hipVel, reach) {
      const s = { leftAnkleY:[], rightAnkleY:[], hipMidY:[], leftVis:[], rightVis:[], tMs:[],
                  leftAnkleX:[], rightAnkleX:[], leftKneeX:[], leftKneeY:[], rightKneeX:[], rightKneeY:[],
                  hipMidX:[], leftKneeVis:[], rightKneeVis:[], frameW:1000, frameH:1000 };
      function ayOf(u){
        return (u < 0.4) ? (0.86 - 0.03 * Math.pow((u - 0.2)/0.2, 2))   // stance dome, max(lowest foot)@midstance
                         : (0.83 - 0.30 * Math.sin(Math.PI * (u - 0.4)/0.6)); // swing lift
      }
      function fxOf(u){ return (u < 0.4) ? (1 - 2*(u/0.4)) : (-1 + 2*((u-0.4)/0.6)); } // +1 td ->0 mid ->-1 toe-off ->+1
      for (let i=0;i<n;i++){
        const t=i/fps; s.tMs.push(t*1000);
        const hip = 0.20 + hipVel*t;                  // rightward -> dir=+1
        const uL = (((i/period)+0.0)%1+1)%1, uR = (((i/period)+0.5)%1+1)%1;
        s.leftAnkleY.push(ayOf(uL));  s.rightAnkleY.push(ayOf(uR));
        s.leftAnkleX.push(hip + reach*fxOf(uL)); s.rightAnkleX.push(hip + reach*fxOf(uR));
        s.leftKneeY.push(ayOf(uL) - 0.18);  s.rightKneeY.push(ayOf(uR) - 0.18);
        s.leftKneeX.push(hip);              s.rightKneeX.push(hip);
        s.hipMidX.push(hip);               s.hipMidY.push(0.50);
        s.leftVis.push(1); s.rightVis.push(1); s.leftKneeVis.push(0.9); s.rightKneeVis.push(0.9);
      }
      return s;
    }

    const sampleRun = genRun(300, 30, 20, 0.03, 0.045);
    const spmRun  = CadenceCore.computeSPM(sampleRun);
    const overRun = CadenceCore.computeOverstride(sampleRun, spmRun.peaks, { w:1000, h:1000 });

    // (TD-a) touchdown sits on the rising edge, not midstance: every td <= its peak, and the
    // touchdown-window angle is clearly larger than the raw peak-frame angle (peak ≈ midstance ~0).
    {
      const strides = overRun.debug.strides;
      const tdBeforePeak = strides.length > 0 && strides.every(s => s.touchdownFrame <= s.peakFrame);
      function rawPeakAngle(s, axA, ayA, kxA, kyA){
        const k=s.peakFrame;
        const dx=Math.max(0, (axA[k]-kxA[k])*1000);
        return Math.atan2(dx, Math.abs((ayA[k]-kyA[k])*1000))*180/Math.PI;
      }
      const peakAngles = strides.map(s => s.foot==='L'
        ? rawPeakAngle(s, sampleRun.leftAnkleX,  sampleRun.leftAnkleY,  sampleRun.leftKneeX,  sampleRun.leftKneeY)
        : rawPeakAngle(s, sampleRun.rightAnkleX, sampleRun.rightAnkleY, sampleRun.rightKneeX, sampleRun.rightKneeY));
      const sorted = peakAngles.slice().sort((a,b)=>a-b);
      const medPeak = sorted.length ? sorted[Math.floor(sorted.length/2)] : 0;
      const tdWins  = overRun.shinAngleDeg > medPeak + 3;
      const forward = overRun.shinAngleDeg >= 5 && overRun.shinAngleDeg <= 16;
      const ok = tdBeforePeak && tdWins && forward;
      console.log(`${ok ? 'PASS' : 'FAIL'} TD-a: final=${overRun.shinAngleDeg}° vs peakFrame≈${Math.round(medPeak*10)/10}° tdBeforePeak=${tdBeforePeak} strides=${overRun.strideCount}`);
    }

    // (TD-b) ±1-frame seek jitter: shifting all peak indices by +1/-1 keeps the final within ±2°.
    {
      function shiftP(peaks, d){
        return { left: peaks.left.map(i => Math.max(0, i+d)), right: peaks.right.map(i => Math.max(0, i+d)) };
      }
      const a0 = overRun.shinAngleDeg;
      const aP = CadenceCore.computeOverstride(sampleRun, shiftP(spmRun.peaks, +1), {w:1000,h:1000}).shinAngleDeg;
      const aM = CadenceCore.computeOverstride(sampleRun, shiftP(spmRun.peaks, -1), {w:1000,h:1000}).shinAngleDeg;
      const spread = Math.max(a0,aP,aM) - Math.min(a0,aP,aM);
      const ok = spread <= 2.0;
      console.log(`${ok ? 'PASS' : 'FAIL'} TD-b: jitter spread=${Math.round(spread*10)/10}° [${a0}, ${aP}, ${aM}]`);
    }

    // (TD-c) sign: foot ahead of knee -> positive; knee ahead of foot (propulsion) -> 0.
    {
      function osArrays(n, obj){ const o={}; for(const k in obj) o[k]=new Array(n).fill(obj[k]); return o; }
      const idxs=[5,12,19];
      const ahead  = osArrays(26, { leftAnkleX:0.62, leftAnkleY:0.62, leftKneeX:0.50, leftKneeY:0.42, leftKneeVis:0.9 });
      ahead.hipMidX=[];  for(let i=0;i<26;i++) ahead.hipMidX.push(0.3+0.01*i);
      const behind = osArrays(26, { leftAnkleX:0.40, leftAnkleY:0.62, leftKneeX:0.55, leftKneeY:0.42, leftKneeVis:0.9 });
      behind.hipMidX=[]; for(let i=0;i<26;i++) behind.hipMidX.push(0.3+0.01*i);
      const rA = CadenceCore.computeOverstride(ahead,  {left:idxs,right:[]}, {w:1000,h:1000});
      const rB = CadenceCore.computeOverstride(behind, {left:idxs,right:[]}, {w:1000,h:1000});
      const ok = rA.shinAngleDeg > 5 && rB.shinAngleDeg === 0;
      console.log(`${ok ? 'PASS' : 'FAIL'} TD-c: foot-ahead=${rA.shinAngleDeg}° knee-ahead=${rB.shinAngleDeg}°`);
    }

    // --- Overstride gating (demotion to experimental) ---
    function osFill(n, obj){ const o={}; for(const k in obj) o[k]=new Array(n).fill(obj[k]); return o; }

    // (GATE-discard) shin > 35deg is physically impossible at contact => discarded => 'unknown'.
    {
      // knee(0.5,0.4), ankle(0.8,0.62) -> dx=300,dy=200 -> atan2(300,200)=56.3deg (>35).
      const s = osFill(30, { leftAnkleX:0.80, leftAnkleY:0.62, leftKneeX:0.50, leftKneeY:0.42,
                             leftKneeVis:0.9, leftVis:0.9, hipMidX:0.5 });
      const r = CadenceCore.computeOverstride(s, { left:[5,12,19,26], right:[] }, { w:1000, h:1000 });
      const ok = r.classification === 'unknown' && r.shinAngleDeg === null;
      console.log(`${ok ? 'PASS' : 'FAIL'} GATE-discard: >35deg dropped -> class=${r.classification} shin=${r.shinAngleDeg}`);
    }

    // (GATE-lowconf) fewer than 6 valid strides -> confidence forced to 'low'.
    {
      // 20deg overstride, only 4 strides.
      const ax = (500 + 200*Math.tan(20*Math.PI/180))/1000;   // ~0.5728
      const s = osFill(30, { leftAnkleX:ax, leftAnkleY:0.60, leftKneeX:0.50, leftKneeY:0.40,
                             leftKneeVis:0.9, leftVis:0.9, hipMidX:0.5 });
      const r = CadenceCore.computeOverstride(s, { left:[5,12,19,26], right:[] }, { w:1000, h:1000 });
      const ok = r.confidence === 'low' && r.classification !== 'unknown' && r.strideCount < 6;
      console.log(`${ok ? 'PASS' : 'FAIL'} GATE-lowconf: strides=${r.strideCount} conf=${r.confidence} class=${r.classification}`);
    }

    // (GATE-nearleg) occluded far leg (low vis) excluded; the near (high-vis) leg is measured.
    {
      const n=30;
      const s = osFill(n, {
        leftAnkleX:0.62,  leftAnkleY:0.62,  leftKneeX:0.50,  leftKneeY:0.42,  leftKneeVis:0.15, leftVis:0.15,
        rightAnkleX:0.62, rightAnkleY:0.62, rightKneeX:0.50, rightKneeY:0.42, rightKneeVis:0.9, rightVis:0.9
      });
      s.hipMidX=[]; for(let i=0;i<n;i++) s.hipMidX.push(0.30+0.01*i);
      const r = CadenceCore.computeOverstride(s, { left:[5,12,19,26], right:[5,12,19,26] }, { w:1000, h:1000 });
      const ok = r.nearFoot === 'R' && r.classification !== 'unknown';
      console.log(`${ok ? 'PASS' : 'FAIL'} GATE-nearleg: nearFoot=${r.nearFoot} (far L vis 0.15 dropped) shin=${r.shinAngleDeg}`);
    }

    // --- Trunk forward lean (reliable metric) ---
    function trunkSample(n, dirSign, shoulderXoff, shoulderY) {
      const s = { shoulderMidX:[], shoulderMidY:[], hipMidX:[], hipMidY:[], shoulderVis:[], hipVis:[] };
      for (let i=0;i<n;i++){
        const hipx = dirSign >= 0 ? (0.30 + 0.005*i) : (0.80 - 0.005*i);  // +1 rightward, -1 leftward
        s.hipMidX.push(hipx); s.hipMidY.push(0.60);
        s.shoulderMidX.push(hipx + shoulderXoff);
        s.shoulderMidY.push(shoulderY);
        s.shoulderVis.push(0.9); s.hipVis.push(0.9);
      }
      return s;
    }

    // (TL-a) known 5deg forward lean: shoulder 17.5px ahead, 200px above -> atan2(17.5,200)=5deg.
    {
      const s = trunkSample(40, +1, 0.0175, 0.40);
      const r = CadenceCore.computeTrunkLean(s, { w:1000, h:1000 });
      const ok = r.leanDeg !== null && Math.abs(r.leanDeg - 5) <= 0.5 && r.classification === 'good' && r.confidence === 'high';
      console.log(`${ok ? 'PASS' : 'FAIL'} TL-a: leanDeg=${r.leanDeg}° class=${r.classification} conf=${r.confidence} frames=${r.frameCount}`);
    }

    // (TL-b) aspect-ratio correctness: same normalized coords, portrait vs landscape -> different angle.
    {
      const s = trunkSample(40, +1, 0.03, 0.40);
      const rP = CadenceCore.computeTrunkLean(s, { w:720,  h:1280 });
      const rL = CadenceCore.computeTrunkLean(s, { w:1280, h:720  });
      const expP = Math.atan2(0.03*720, 0.20*1280) * 180/Math.PI;   // ~4.8
      const expL = Math.atan2(0.03*1280, 0.20*720) * 180/Math.PI;   // ~14.9
      const ok = Math.abs(rP.leanDeg - Math.round(expP*10)/10) <= 0.2 &&
                 Math.abs(rL.leanDeg - Math.round(expL*10)/10) <= 0.2 &&
                 Math.abs(rP.leanDeg - rL.leanDeg) > 1;
      console.log(`${ok ? 'PASS' : 'FAIL'} TL-b: portrait=${rP.leanDeg}(exp≈${Math.round(expP*10)/10}) landscape=${rL.leanDeg}(exp≈${Math.round(expL*10)/10})`);
    }

    // (TL-c) sign: shoulders ahead in travel dir -> + (both L->R and R->L); shoulders behind -> -.
    {
      const fwdR = CadenceCore.computeTrunkLean(trunkSample(40, +1, +0.03, 0.40), { w:1000, h:1000 }); // rightward, ahead
      const fwdL = CadenceCore.computeTrunkLean(trunkSample(40, -1, -0.03, 0.40), { w:1000, h:1000 }); // leftward, ahead (mirror)
      const back = CadenceCore.computeTrunkLean(trunkSample(40, +1, -0.03, 0.40), { w:1000, h:1000 }); // rightward, behind
      const ok = fwdR.leanDeg > 0 && fwdL.leanDeg > 0 &&
                 Math.abs(fwdR.leanDeg - fwdL.leanDeg) <= 0.2 && back.leanDeg < 0;
      console.log(`${ok ? 'PASS' : 'FAIL'} TL-c: fwdR=${fwdR.leanDeg}° fwdL=${fwdL.leanDeg}° back=${back.leanDeg}°`);
    }

    // (TL-d) dead-zone: a near-vertical torso (-2.2deg) must classify 'good', not 'back'.
    {
      const n=40;
      const s={ shoulderMidX:[], shoulderMidY:[], hipMidX:[], hipMidY:[], shoulderVis:[], hipVis:[] };
      // tan(2.2deg)=0.03843 -> dx = -0.03843*200 = -7.686px -> (sx-hx) = -0.007686 at w=1000
      for(let i=0;i<n;i++){ const hipx=0.30+0.005*i;
        s.hipMidX.push(hipx); s.hipMidY.push(0.60);
        s.shoulderMidX.push(hipx - 0.007686); s.shoulderMidY.push(0.40);
        s.shoulderVis.push(0.9); s.hipVis.push(0.9); }
      const r = CadenceCore.computeTrunkLean(s, { w:1000, h:1000 });
      const ok = Math.abs(r.leanDeg - (-2.2)) <= 0.3 && r.classification === 'good';
      console.log(`${ok ? 'PASS' : 'FAIL'} TL-d deadzone: leanDeg=${r.leanDeg}° class=${r.classification}`);
    }

    // --- Vertical oscillation ---
    function genVO(n, fps, A, driftAmp, torsoNorm){
      const s={ hipMidY:[], hipMidX:[], shoulderMidY:[], shoulderMidX:[], hipVis:[], shoulderVis:[], tMs:[] };
      for(let i=0;i<n;i++){
        const t=i/fps; s.tMs.push(t*1000);
        const drift = driftAmp*(i/(n-1));
        const hy = 0.50 + A*Math.sin(2*Math.PI*1.5*t) + drift;
        s.hipMidY.push(hy); s.hipMidX.push(0.50);
        s.shoulderMidY.push(hy - torsoNorm); s.shoulderMidX.push(0.50);  // pure-vertical torso
        s.hipVis.push(0.9); s.shoulderVis.push(0.9);
      }
      return s;
    }

    // (VO-a) amplitude linearity: doubling hip bounce ~doubles the ratio.
    {
      const r1 = CadenceCore.computeVerticalOscillation(genVO(300,30,0.03,0,0.20), { w:1000, h:1000 });
      const r2 = CadenceCore.computeVerticalOscillation(genVO(300,30,0.06,0,0.20), { w:1000, h:1000 });
      const lin = Math.abs(r2.ratioPct - 2*r1.ratioPct) <= 0.05*2*r1.ratioPct;
      const ok = r1.ratioPct > 0 && lin && r1.classification !== 'unknown' && r2.classification !== 'unknown';
      console.log(`${ok ? 'PASS' : 'FAIL'} VO-a: r1=${r1.ratioPct}%(${r1.classification}) r2=${r2.ratioPct}% linear=${lin}`);
    }

    // (VO-b) drift invariance: detrend removes slow pan, ratio stays put.
    {
      const noDrift = CadenceCore.computeVerticalOscillation(genVO(300,30,0.03,0,   0.20), { w:1000, h:1000 });
      const drift   = CadenceCore.computeVerticalOscillation(genVO(300,30,0.03,0.10,0.20), { w:1000, h:1000 });
      const ok = Math.abs(noDrift.ratioPct - drift.ratioPct) <= 1.0;
      console.log(`${ok ? 'PASS' : 'FAIL'} VO-b: noDrift=${noDrift.ratioPct}% drift=${drift.ratioPct}% (Δ≤1%p)`);
    }

    // (VO-c) scale invariance: same normalized coords, different frameH -> same ratio.
    {
      const s = genVO(300,30,0.03,0,0.20);
      const rA = CadenceCore.computeVerticalOscillation(s, { w:720,  h:720  });
      const rB = CadenceCore.computeVerticalOscillation(s, { w:1440, h:1440 });
      const ok = Math.abs(rA.ratioPct - rB.ratioPct) <= 0.2;
      console.log(`${ok ? 'PASS' : 'FAIL'} VO-c: h720=${rA.ratioPct}% h1440=${rB.ratioPct}% (scale-free)`);
    }
  }

  return CadenceCore;
}));
