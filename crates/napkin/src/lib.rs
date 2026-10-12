//! # napkin
//!
//! A small language for engineering calculations where every number carries
//! **units** and **uncertainty**. Each uncertain value is a vector of samples
//! that share a "world" index, so arithmetic runs sample by sample and
//! correlations survive: `h - h` is exactly zero.
//!
//! This crate is the language: syntax, units, sampling, and the general
//! builtins (distributions, statistics, geometry, weight groups). A host adds
//! domain functions and opaque values through [`Interp::register`] and
//! [`value::Native`]; `napkin-boat` adds boatmath's hulls that way.

pub mod ast;
pub mod builtins;
pub mod dim;
pub mod interp;
pub mod lexer;
pub mod ops;
pub mod parser;
pub mod samples;
pub mod units;
pub mod value;

pub use dim::{Dim, Rational};
pub use interp::Interp;
pub use samples::Samples;
pub use units::UnitExpr;
pub use value::{Native, Quantity, Value};

use std::cell::RefCell;
use std::io::{BufRead, Write};
use std::rc::Rc;

/// Run `src` on `it`, returning what it printed instead of writing it out.
pub fn run_captured(it: &mut Interp, src: &str) -> (String, Result<(), String>) {
    let buf = Rc::new(RefCell::new(String::new()));
    let sink = buf.clone();
    let old = std::mem::replace(
        &mut it.output,
        Box::new(move |s| {
            let mut b = sink.borrow_mut();
            b.push_str(s);
            b.push('\n');
        }),
    );
    let r = it.run(src);
    it.output = old;
    let out = buf.borrow().clone();
    (out, r)
}

/// An interactive session on stdin; each complete input (brackets balanced)
/// runs as it is entered, keeping its bindings for the next.
pub fn repl(it: &mut Interp) {
    println!("napkin — units and uncertainty. Ctrl-D to exit.");
    let stdin = std::io::stdin();
    let mut buf = String::new();
    loop {
        print!("{}", if buf.is_empty() { "› " } else { "… " });
        std::io::stdout().flush().ok();
        let mut line = String::new();
        if stdin.lock().read_line(&mut line).unwrap_or(0) == 0 {
            println!();
            break;
        }
        buf.push_str(&line);
        let open = buf.chars().filter(|c| "([{".contains(*c)).count();
        let close = buf.chars().filter(|c| ")]}".contains(*c)).count();
        if open > close {
            continue;
        }
        if let Err(e) = it.run(&buf) {
            eprintln!("error: {}", e);
        }
        buf.clear();
    }
}
