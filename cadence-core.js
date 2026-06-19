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
      debug: {
        fps, spmPeak, spmAuto, combinedCV, leftCV, rightCV,
        leftPeakCount: leftPeaks.length, rightPeakCount: rightPeaks.length,
        leftMissingFrac: leftPrep.missingFrac, rightMissingFrac: rightPrep.missingFrac,
        bestLag, bestAC
      }
    };
  }

  const CadenceCore = { computeSPM, lowpass, detrend, findPeaks, autocorr };

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
  }

  return CadenceCore;
}));
