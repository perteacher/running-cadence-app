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
    let bestLag = minLag, bestAC = -Infinity;
    for (let lag = minLag; lag <= Math.min(maxLag, ac.length - 1); lag++) {
      if (ac[lag] > bestAC) { bestAC = ac[lag]; bestLag = lag; }
    }
    const spmAuto = bestLag > 0 ? (fps / bestLag) * 60 * 2 : 0;

    // 8. Choose method
    const enoughPeaks = totalSteps >= 6;
    const regularSpacing = combinedCV < 0.30;
    let method, spmRaw, confidence;

    if (enoughPeaks && regularSpacing) {
      method = 'peak';
      spmRaw = spmPeak;
      if (combinedCV < 0.15) confidence = 'high';
      else confidence = 'medium';
    } else {
      method = 'autocorr';
      spmRaw = spmAuto;
      const heavilyInterp = leftPrep.missingFrac > 0.4 && rightPrep.missingFrac > 0.4;
      confidence = heavilyInterp ? 'low' : 'medium';
    }

    // 9. Round
    const spm = Math.round(spmRaw);

    return {
      spm,
      confidence,
      method,
      stepCount: totalSteps,
      durationSec,
      fps,
      peaks: { left: leftPeaks, right: rightPeaks },
      debug: {
        fps, spmPeak, spmAuto, combinedCV, leftCV, rightCV,
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
      leftKneeVis, rightKneeVis
    } = samples;

    const strikes = [];

    function addSide(idxArr, ankleXArr, ankleYArr, kneeXArr, kneeYArr, kneeVisArr) {
      for (const i of idxArr) {
        const ax = ankleXArr && ankleXArr[i];
        const ay = ankleYArr && ankleYArr[i];
        const kx = kneeXArr  && kneeXArr[i];
        const ky = kneeYArr  && kneeYArr[i];
        const hx = hipMidX   && hipMidX[i];
        const vis= kneeVisArr && kneeVisArr[i];
        if (!isFinite(ax) || !isFinite(ay) || !isFinite(kx) ||
            !isFinite(ky) || !isFinite(hx) || !(vis >= 0.3)) continue;
        const axpx = ax * w, aypx = ay * h;
        const kxpx = kx * w, kypx = ky * h;
        const hxpx = hx * w;
        const shin   = Math.atan2(Math.abs(axpx - kxpx), Math.abs(aypx - kypx)) * 180 / Math.PI;
        const offset = axpx - hxpx;
        strikes.push({ shin, offset, vis });
      }
    }

    addSide(peaks.left,  leftAnkleX,  leftAnkleY,  leftKneeX,  leftKneeY,  leftKneeVis);
    addSide(peaks.right, rightAnkleX, rightAnkleY, rightKneeX, rightKneeY, rightKneeVis);

    if (strikes.length === 0) {
      return { shinAngleDeg: null, horizOffset: null, classification: 'unknown', confidence: 'low', strikeCount: 0 };
    }

    const shins   = strikes.map(s => s.shin);
    const offsets = strikes.map(s => s.offset);
    const visArr  = strikes.map(s => s.vis);

    const shinAngleDeg = Math.round(median(shins) * 10) / 10;
    const horizOffset  = Math.round(median(offsets));
    const strikeCount  = strikes.length;
    const avgVis       = visArr.reduce((a, b) => a + b, 0) / visArr.length;

    let classification;
    if (shinAngleDeg < 10)       classification = 'good';
    else if (shinAngleDeg < 15)  classification = 'borderline';
    else                         classification = 'overstride';

    let confidence;
    if (strikeCount >= 4 && avgVis >= 0.6)      confidence = 'high';
    else if (strikeCount >= 3 && avgVis >= 0.4) confidence = 'medium';
    else                                         confidence = 'low';

    return { shinAngleDeg, horizOffset, classification, confidence, strikeCount };
  }

  const CadenceCore = { computeSPM, computeOverstride, lowpass, detrend, findPeaks, autocorr };

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

    // (b) Known 45 deg: knee(0.5,0.4)->px(500,400), ankle(0.7,0.6)->px(700,600), dx=dy=200 -> 45deg
    {
      const idxs = [10, 20, 30, 40];
      const s = makeOSsamples(50, 0.7, 0.7, 0.5, 0.4, 0.5, 0.4, 0.5, 0.6, 0.6, 0.9, 0.9);
      const peaks = { left: idxs, right: [] };
      const r = CadenceCore.computeOverstride(s, peaks, { w: 1000, h: 1000 });
      const ok = r.shinAngleDeg !== null && Math.abs(r.shinAngleDeg - 45) <= 1 && r.confidence === 'high';
      console.log(`${ok ? 'PASS' : 'FAIL'} OSTest-b: 45deg shinAngleDeg=${r.shinAngleDeg} conf=${r.confidence}`);
    }

    // (c) Aspect-ratio correctness: same normalized coords, portrait vs landscape -> DIFFERENT pixel angles
    {
      // knee(0.5,0.4) ankle(0.6,0.6) with portrait w=720 h=1280:
      //   kxpx=360,kypx=512; axpx=432,aypx=768 -> dx=72,dy=256 -> atan2(72,256)
      // with landscape w=1280 h=720:
      //   kxpx=640,kypx=288; axpx=768,aypx=432 -> dx=128,dy=144 -> atan2(128,144)
      const expectedPortrait  = Math.atan2(Math.abs(0.6*720 - 0.5*720), Math.abs(0.6*1280 - 0.4*1280)) * 180 / Math.PI;
      const expectedLandscape = Math.atan2(Math.abs(0.6*1280 - 0.5*1280), Math.abs(0.6*720 - 0.4*720)) * 180 / Math.PI;
      const idxs = [5, 10, 15, 20];
      const s = makeOSsamples(30, 0.6, 0.6, 0.5, 0.4, 0.5, 0.4, 0.5, 0.6, 0.6, 0.9, 0.9);
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
  }

  return CadenceCore;
}));
