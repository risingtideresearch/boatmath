//! Thread fan-out for the library's embarrassingly parallel inner loops.
//!
//! The crate stays dependency-free: fan-out is `std::thread::scope` with one
//! worker per core, each carrying its own clone of whatever scratch state the
//! work items need, pulling items off a shared counter so that fast and slow
//! cores (performance and efficiency cores on Apple silicon, say) each take
//! what they can. Results come back in input order, and every reduction the
//! callers perform over them runs in the original serial order — so a
//! parallel run is **bit-for-bit identical** to a serial one; only the wall
//! clock changes.
//!
//! # Thread budget
//!
//! How many workers a fan-out may use is a **per-thread** setting
//! ([`set_threads`] / [`with_threads`]), defaulting to the machine's
//! available parallelism. Workers spawned here run with a budget of one, so
//! nested library calls never fan out again; a caller that already spreads
//! independent jobs (sweep points, say) across its own threads should give
//! each of them a budget of `cores / jobs` — or one — so the machine is not
//! oversubscribed.
//!
//! The process-wide default can be pinned with the `BOATMATH_THREADS`
//! environment variable (falling back to `MICHELL_THREADS` for compatibility;
//! a positive integer, `1` disables threading), which is handy for
//! benchmarking and for sharing a machine.

use std::cell::Cell;
use std::num::NonZeroUsize;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

thread_local! {
    /// `None` = auto ([`available`]); `Some(n)` = at most `n` workers.
    static BUDGET: Cell<Option<NonZeroUsize>> = const { Cell::new(None) };
}

/// The default worker count: `BOATMATH_THREADS` or `MICHELL_THREADS` if set to a positive integer,
/// otherwise the machine's available parallelism (at least 1).
pub fn available() -> usize {
    static DEFAULT: OnceLock<usize> = OnceLock::new();
    *DEFAULT.get_or_init(|| {
        std::env::var("BOATMATH_THREADS")
            .or_else(|_| std::env::var("MICHELL_THREADS"))
            .ok()
            .and_then(|v| v.trim().parse::<usize>().ok())
            .filter(|&n| n > 0)
            .unwrap_or_else(|| {
                std::thread::available_parallelism()
                    .map(|n| n.get())
                    .unwrap_or(1)
            })
    })
}

/// The worker budget for fan-outs started from the current thread.
pub fn threads() -> usize {
    BUDGET.with(|b| b.get().map_or_else(available, |n| n.get()))
}

/// Set the worker budget for fan-outs started from the current thread.
/// `0` restores the default (available parallelism); `1` makes every
/// library call from this thread run serially.
pub fn set_threads(n: usize) {
    BUDGET.with(|b| b.set(NonZeroUsize::new(n)));
}

/// Run `f` with the current thread's budget temporarily set to `n`
/// (`0` = default), restoring the previous budget afterwards.
pub fn with_threads<R>(n: usize, f: impl FnOnce() -> R) -> R {
    let prev = BUDGET.with(|b| b.replace(NonZeroUsize::new(n)));
    struct Restore(Option<NonZeroUsize>);
    impl Drop for Restore {
        fn drop(&mut self) {
            BUDGET.with(|b| b.set(self.0));
        }
    }
    let _restore = Restore(prev);
    f()
}

/// Evaluate `f(state, i)` for `i in 0..n`, returning the results in index
/// order. Up to `threads()` workers pull indices from a shared counter, each
/// holding its own `state` (from `init`) — the per-worker scratch that makes
/// the items independent. Serial when the budget is one or `n` is small; in
/// either case the per-item arithmetic is exactly that of a plain loop, and
/// the caller's reduction over the returned `Vec` fixes the summation order.
pub fn map_indexed<S, R, I, F>(n: usize, init: I, f: F) -> Vec<R>
where
    S: Send,
    R: Send,
    I: Fn() -> S + Sync,
    F: Fn(&mut S, usize) -> R + Sync,
{
    // Below this many items per worker the spawn cost is not worth paying.
    const MIN_ITEMS_PER_WORKER: usize = 4;
    let workers = threads().min(n / MIN_ITEMS_PER_WORKER).max(1);
    if workers == 1 {
        let mut state = init();
        return (0..n).map(|i| f(&mut state, i)).collect();
    }
    let next = AtomicUsize::new(0);
    let (init, f, next) = (&init, &f, &next);
    let mut tagged: Vec<(usize, R)> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|_| {
                scope.spawn(move || {
                    // Nested fan-outs from a worker would oversubscribe.
                    set_threads(1);
                    let mut state = init();
                    let mut local: Vec<(usize, R)> = Vec::with_capacity(n / workers + 1);
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        if i >= n {
                            break;
                        }
                        local.push((i, f(&mut state, i)));
                    }
                    local
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().expect("worker panicked"))
            .collect()
    });
    tagged.sort_unstable_by_key(|(i, _)| *i);
    tagged.into_iter().map(|(_, r)| r).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn results_come_back_in_order_and_serial_matches_parallel() {
        let n = 1000;
        let work = |s: &mut usize, i: usize| {
            *s += 1;
            (i, *s)
        };
        let serial = with_threads(1, || map_indexed(n, || 0usize, work));
        let parallel = with_threads(8, || map_indexed(n, || 0usize, work));
        assert_eq!(serial.len(), n);
        assert_eq!(parallel.len(), n);
        assert!(serial.iter().enumerate().all(|(k, (i, _))| k == *i));
        assert!(parallel.iter().enumerate().all(|(k, (i, _))| k == *i));
        // Per-worker state restarts per worker; the serial state runs through.
        assert_eq!(serial.last().unwrap().1, n);
        assert!(
            parallel.iter().map(|(_, s)| *s).sum::<usize>() <= serial.iter().map(|(_, s)| *s).sum()
        );
    }

    #[test]
    fn budget_is_per_thread_and_restored() {
        assert_eq!(threads(), available());
        with_threads(3, || {
            assert_eq!(threads(), 3);
            // Workers see a budget of one.
            let seen = map_indexed(64, || (), |_, _| threads());
            assert!(seen.iter().all(|&t| t == 1));
        });
        assert_eq!(threads(), available());
        set_threads(2);
        assert_eq!(threads(), 2);
        set_threads(0);
        assert_eq!(threads(), available());
    }
}
