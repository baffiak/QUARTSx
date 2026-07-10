//! Shared timestamped logging helper used across every pipeline stage.
//!
//! All output goes to stderr (never stdout). A single [`Stage`] value carries a
//! stage name and a start instant so every line is prefixed with a wall-clock
//! timestamp and the stage can report its own elapsed time when it finishes.
//! Heartbeats are rate-limited so long-running loops can call [`Stage::beat`]
//! per batch without ever emitting per-read output.
//!
//! [`Progress`] is the shared byte-progress-bar primitive that replaces the old
//! ~5 s `beat` heartbeats inside the long streaming loops (FILTERING §1 over input
//! bytes; COUNTING §3 over BGZF bytes of checkpoint 2). It is a lock-free counter
//! with a rate-limited, timestamped in-place redraw: workers call
//! [`Progress::set`] on the hot path (a relaxed atomic write plus a cheap time
//! check) and at most one thread ever repaints the line.

use std::io::Write;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::{Duration, Instant};

/// Wall-clock timestamp "2026-07-05 21:59:03" via libc (no chrono dependency).
pub fn ts() -> String {
    let t = unsafe { libc::time(std::ptr::null_mut()) };
    let mut tm: libc::tm = unsafe { std::mem::zeroed() };
    unsafe { libc::localtime_r(&t, &mut tm) };
    let mut buf = [0u8; 32];
    let n = unsafe {
        libc::strftime(
            buf.as_mut_ptr() as *mut libc::c_char,
            buf.len(),
            b"%Y-%m-%d %H:%M:%S\0".as_ptr() as *const libc::c_char,
            &tm,
        )
    };
    String::from_utf8_lossy(&buf[..n]).into_owned()
}

/// Human-friendly duration: "27.4s" / "12m03s" / "1h04m".
pub fn fmt_dur(d: Duration) -> String {
    let s = d.as_secs_f64();
    if s < 90.0 {
        format!("{s:.1}s")
    } else if s < 3600.0 {
        format!("{}m{:02}s", (s / 60.0) as u64, (s % 60.0) as u64)
    } else {
        format!("{}h{:02}m", (s / 3600.0) as u64, ((s % 3600.0) / 60.0) as u64)
    }
}

/// One pipeline stage. Emits a banner on [`Stage::begin`] and an elapsed line on
/// [`Stage::done`]; [`Stage::step`] and [`Stage::beat`] print progress in between.
pub struct Stage {
    name: &'static str,
    start: Instant,
    last_beat: Mutex<Instant>,
}

impl Stage {
    /// Print the stage banner and start the clock.
    pub fn begin(name: &'static str, msg: impl AsRef<str>) -> Stage {
        let now = Instant::now();
        eprintln!("[{}] {name}: {}", ts(), msg.as_ref());
        Stage { name, start: now, last_beat: Mutex::new(now) }
    }

    /// A discrete sub-step within the stage (always printed).
    pub fn step(&self, msg: impl AsRef<str>) {
        eprintln!("[{}] {}: {}", ts(), self.name, msg.as_ref());
    }

    /// Heartbeat: prints only if at least `every` has elapsed since the last
    /// emitted beat. Cheap; safe to call per batch (never per read).
    pub fn beat(&self, every: Duration, msg: impl FnOnce() -> String) {
        let mut lb = self.last_beat.lock().unwrap();
        if lb.elapsed() >= every {
            eprintln!("[{}] {}: {}", ts(), self.name, msg());
            *lb = Instant::now();
        }
    }

    /// Print the closing line with the stage's total elapsed time.
    pub fn done(&self, msg: impl AsRef<str>) {
        eprintln!("[{}] {}: done in {} ({})", ts(), self.name, fmt_dur(self.start.elapsed()), msg.as_ref());
    }

    /// Start a byte-denominated progress bar carrying this stage's name.
    /// `total` is the number of bytes that will be consumed (0 = unknown, in
    /// which case only a running count + rate + elapsed are shown, no bar/%/ETA).
    pub fn progress_bytes(&self, total: u64) -> Progress {
        Progress::new(self.name, total)
    }
}

/// Shared byte-progress-bar primitive: a lock-free counter with a rate-limited,
/// timestamped, in-place redraw to stderr. Replaces the old 5 s `beat`
/// heartbeats in the FILTERING (§1) and COUNTING (§3) streaming loops.
///
/// Concurrency model: [`set`](Progress::set) does one relaxed atomic write plus a
/// cheap elapsed-time comparison, so any number of rayon workers may report
/// progress on the hot path. When the redraw interval has elapsed exactly one
/// thread wins a `compare_exchange` on the last-draw timestamp and repaints the
/// line; all others return immediately. This keeps per-report overhead to a
/// couple of atomics and guarantees the redraw is never per-read output.
///
/// Rendering adapts to the sink: on a TTY the line is repainted in place with a
/// carriage return (trailing pad clears any shrinkage, e.g. a falling ETA); when
/// stderr is redirected to a log file the same content is emitted as discrete
/// timestamped lines on a slower cadence so the log is not flooded.
pub struct Progress {
    label: &'static str,
    /// Denominator; 0 means "unknown total" (no bar / percent / ETA).
    total: u64,
    current: AtomicU64,
    start: Instant,
    /// Millis-since-`start` of the last repaint; the redraw arbiter.
    last_draw_ms: AtomicU64,
    /// Minimum gap between repaints (small on a TTY for smoothness, large when
    /// writing to a log file to avoid thousands of lines).
    interval_ms: u64,
    /// Width in chars of the previously painted line, used to pad-clear on a TTY.
    last_len: AtomicUsize,
    tty: bool,
    done: AtomicBool,
}

impl Progress {
    /// Create a progress bar. Prefer [`Stage::progress_bytes`] /
    /// [`Stage::progress_reads`] so the bar inherits the stage name.
    pub fn new(label: &'static str, total: u64) -> Progress {
        let tty = stderr_is_tty();
        Progress {
            label,
            total,
            current: AtomicU64::new(0),
            start: Instant::now(),
            last_draw_ms: AtomicU64::new(0),
            // 150 ms → smooth on a terminal; 5 s → sparse in a redirected log,
            // matching the old heartbeat cadence the bar replaces.
            interval_ms: if tty { 150 } else { 5_000 },
            last_len: AtomicUsize::new(0),
            tty,
            done: AtomicBool::new(false),
        }
    }

    /// Set the absolute counter value and possibly repaint. Use when the consumer
    /// already tracks a cumulative position (e.g. bytes read from a BAM reader).
    pub fn set(&self, value: u64) {
        self.current.store(value, Ordering::Relaxed);
        self.maybe_draw(value);
    }

    /// Rate-limited repaint arbiter: at most one thread repaints per interval.
    fn maybe_draw(&self, cur: u64) {
        if self.done.load(Ordering::Relaxed) {
            return;
        }
        let now_ms = self.start.elapsed().as_millis() as u64;
        let last = self.last_draw_ms.load(Ordering::Relaxed);
        if now_ms.saturating_sub(last) < self.interval_ms {
            return;
        }
        // Claim the repaint slot; if another thread beat us, let it draw.
        if self
            .last_draw_ms
            .compare_exchange(last, now_ms, Ordering::Relaxed, Ordering::Relaxed)
            .is_err()
        {
            return;
        }
        self.draw(cur, now_ms);
    }

    /// Finish the bar: force a final repaint at the current value and terminate
    /// the in-place line (a newline on a TTY). Idempotent. Call before the
    /// stage's own [`Stage::done`] closing line so the bar does not clobber it.
    pub fn finish(&self) {
        if self.done.swap(true, Ordering::Relaxed) {
            return;
        }
        let cur = self.current.load(Ordering::Relaxed);
        let now_ms = self.start.elapsed().as_millis() as u64;
        self.draw(cur, now_ms);
        if self.tty {
            let mut err = std::io::stderr().lock();
            let _ = writeln!(err);
            let _ = err.flush();
        }
    }

    /// Build and emit one line. On a TTY: `\r` + content + pad-clear, no newline.
    /// Off a TTY: a plain timestamped line terminated with `\n`.
    fn draw(&self, cur: u64, elapsed_ms: u64) {
        let secs = elapsed_ms as f64 / 1000.0;
        let rate = if secs > 0.0 { cur as f64 / secs } else { 0.0 };
        let rate_str = format!("{}/s", fmt_bytes(rate as u64));
        let body = if self.total > 0 {
            let frac = (cur as f64 / self.total as f64).min(1.0);
            let pct = (frac * 100.0) as u64;
            let eta = if rate > 0.0 && cur < self.total {
                fmt_dur(Duration::from_secs_f64((self.total - cur) as f64 / rate))
            } else {
                "0s".to_string()
            };
            format!(
                "{} {:>3}%  {}/{}  {}  eta {}",
                bar(frac, 24), pct, fmt_bytes(cur), fmt_bytes(self.total), rate_str, eta
            )
        } else {
            // Unknown total: running count + rate + elapsed, no bar/percent/ETA.
            format!("{}  {}  ({})", fmt_bytes(cur), rate_str, fmt_dur(Duration::from_secs_f64(secs)))
        };
        let line = format!("[{}] {}: {}", ts(), self.label, body);
        let mut err = std::io::stderr().lock();
        if self.tty {
            let prev = self.last_len.swap(line.len(), Ordering::Relaxed);
            let pad = prev.saturating_sub(line.len());
            let _ = write!(err, "\r{}{}", line, " ".repeat(pad));
        } else {
            let _ = writeln!(err, "{line}");
        }
        let _ = err.flush();
    }
}

/// True when stderr (fd 2) is an interactive terminal, so in-place `\r` redraw is
/// safe; false when redirected to a file/pipe (use discrete lines instead).
fn stderr_is_tty() -> bool {
    unsafe { libc::isatty(2) == 1 }
}

/// Fixed-width `[===>   ]` bar for a fraction in `[0,1]`.
fn bar(frac: f64, width: usize) -> String {
    let filled = ((frac * width as f64).round() as usize).min(width);
    let mut s = String::with_capacity(width + 2);
    s.push('[');
    for i in 0..width {
        if i < filled {
            s.push('=');
        } else if i == filled {
            s.push('>');
        } else {
            s.push(' ');
        }
    }
    s.push(']');
    s
}

/// Binary byte formatting: "0B" / "512B" / "1.5K" / "3.2G".
fn fmt_bytes(n: u64) -> String {
    const U: [&str; 5] = ["B", "K", "M", "G", "T"];
    let mut v = n as f64;
    let mut i = 0;
    while v >= 1024.0 && i < U.len() - 1 {
        v /= 1024.0;
        i += 1;
    }
    if i == 0 {
        format!("{n}B")
    } else {
        format!("{v:.1}{}", U[i])
    }
}

/// Final unmistakable SUCCESS line for the whole run.
pub fn success(out_dir: &str, total: Duration) {
    eprintln!("[{}] SUCCESS: run complete in {} -> {out_dir}", ts(), fmt_dur(total));
}

/// Final unmistakable FAILURE line, optionally pointing at a full child log.
pub fn failure(stage: &str, reason: &str, logpath: Option<&std::path::Path>) {
    match logpath {
        Some(p) => eprintln!("[{}] FAILURE: {stage} failed — {reason}; full log: {}", ts(), p.display()),
        None => eprintln!("[{}] FAILURE: {stage} failed — {reason}", ts()),
    }
}
