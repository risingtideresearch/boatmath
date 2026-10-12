const USAGE: &str = "usage: napkin [--samples N] [--seed S] [file.nap]
  With no file, starts an interactive session.";

fn main() {
    let mut args = std::env::args().skip(1);
    let mut file = None;
    let mut n = 10_000;
    let mut seed = 1;
    while let Some(a) = args.next() {
        match a.as_str() {
            "--samples" | "-n" => {
                n = args
                    .next()
                    .and_then(|s| s.parse().ok())
                    .unwrap_or_else(|| die(USAGE))
            }
            "--seed" => {
                seed = args
                    .next()
                    .and_then(|s| s.parse().ok())
                    .unwrap_or_else(|| die(USAGE))
            }
            "-h" | "--help" => {
                println!("{}", USAGE);
                return;
            }
            _ => file = Some(a),
        }
    }
    let mut it = napkin_boat::interp(n, seed);
    match file {
        Some(path) => {
            let src =
                std::fs::read_to_string(&path).unwrap_or_else(|e| die(&format!("{}: {}", path, e)));
            if let Err(e) = it.run(&src) {
                die(&format!("{}: {}", path, e));
            }
        }
        None => napkin::repl(&mut it),
    }
}

fn die(msg: &str) -> ! {
    eprintln!("error: {}", msg);
    std::process::exit(1)
}
