//! Tokenizer. Tracks whether whitespace preceded each token, which the parser uses to
//! decide where a unit expression like `kg/m^3` ends (`5 kg/m^3` vs `x / m`).

#[derive(Clone, Debug, PartialEq)]
pub enum Tok {
    Num(f64),
    Ident(String),
    Str(String),
    Sym(&'static str),
    Newline,
    Eof,
}

#[derive(Clone, Debug)]
pub struct Token {
    pub tok: Tok,
    pub line: usize,
    pub space_before: bool,
}

const SYMS: &[&str] = &[
    "->", "=>", "==", "!=", "<=", ">=", "+-", "±", "·", "(", ")", "[", "]", "{", "}", ",", ":",
    ";", ".", "+", "-", "*", "/", "^", "<", ">", "=", "|", "%",
];

pub fn lex(src: &str) -> Result<Vec<Token>, String> {
    let chars: Vec<char> = src.chars().collect();
    let mut out = vec![];
    let mut i = 0;
    let mut line = 1;
    let mut space = true;
    // Newlines are insignificant inside () and [].
    let mut depth: Vec<char> = vec![];
    while i < chars.len() {
        let c = chars[i];
        if c == '#' {
            while i < chars.len() && chars[i] != '\n' {
                i += 1;
            }
            continue;
        }
        if c == '\n' {
            if !matches!(depth.last(), Some('(') | Some('[')) {
                out.push(Token {
                    tok: Tok::Newline,
                    line,
                    space_before: true,
                });
            }
            line += 1;
            i += 1;
            space = true;
            continue;
        }
        if c.is_whitespace() {
            i += 1;
            space = true;
            continue;
        }
        let start_line = line;
        let tok = if c.is_ascii_digit()
            || (c == '.' && i + 1 < chars.len() && chars[i + 1].is_ascii_digit())
        {
            let st = i;
            while i < chars.len()
                && (chars[i].is_ascii_digit() || chars[i] == '.' || chars[i] == '_')
            {
                i += 1;
            }
            if i < chars.len() && (chars[i] == 'e' || chars[i] == 'E') {
                let mut j = i + 1;
                if j < chars.len() && (chars[j] == '+' || chars[j] == '-') {
                    j += 1;
                }
                if j < chars.len() && chars[j].is_ascii_digit() {
                    i = j;
                    while i < chars.len() && chars[i].is_ascii_digit() {
                        i += 1;
                    }
                }
            }
            let s: String = chars[st..i].iter().filter(|&&c| c != '_').collect();
            Tok::Num(
                s.parse()
                    .map_err(|_| format!("line {}: bad number '{}'", line, s))?,
            )
        } else if c.is_alphabetic() || c == '_' || c == 'µ' || c == '°' {
            let st = i;
            while i < chars.len()
                && (chars[i].is_alphanumeric()
                    || chars[i] == '_'
                    || chars[i] == 'µ'
                    || chars[i] == '°')
            {
                i += 1;
            }
            Tok::Ident(chars[st..i].iter().collect())
        } else if c == '"' {
            i += 1;
            let mut s = String::new();
            while i < chars.len() && chars[i] != '"' {
                if chars[i] == '\n' {
                    line += 1;
                }
                if chars[i] == '\\' && i + 1 < chars.len() {
                    i += 1;
                    s.push(match chars[i] {
                        'n' => '\n',
                        't' => '\t',
                        c => c,
                    });
                } else {
                    s.push(chars[i]);
                }
                i += 1;
            }
            if i >= chars.len() {
                return Err(format!("line {}: unterminated string", start_line));
            }
            i += 1;
            Tok::Str(s)
        } else {
            let rest: String = chars[i..chars.len().min(i + 2)].iter().collect();
            let sym = SYMS
                .iter()
                .find(|s| rest.starts_with(**s))
                .ok_or_else(|| format!("line {}: unexpected character '{}'", line, c))?;
            i += sym.chars().count();
            match *sym {
                "(" | "[" | "{" => depth.push(sym.chars().next().unwrap()),
                ")" | "]" | "}" => {
                    depth.pop();
                }
                _ => {}
            }
            Tok::Sym(match *sym {
                "+-" => "±",
                s => s,
            })
        };
        out.push(Token {
            tok,
            line: start_line,
            space_before: space,
        });
        space = false;
    }
    out.push(Token {
        tok: Tok::Eof,
        line,
        space_before: true,
    });
    Ok(out)
}
