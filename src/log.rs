//! Shared timestamped logging helper used across every pipeline stage.
//!
//! All output goes to stderr (never stdout). A single [`Stage`] value carries a
//! stage name and a start instant so every line is prefixed with a wall-clock
//! timestamp and the stage can report its own elapsed time when it finishes.
//! Heartbeats are rate-limited so long-running loops can call [`Stage::beat`]
//! per batch without ever emitting per-read output.

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

    /// Elapsed time since [`Stage::begin`].
    pub fn elapsed(&self) -> Duration {
        self.start.elapsed()
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
