use crate::ast::{Expr, Stmt, StmtKind};
use crate::dim::Rational;
use crate::lexer::{Tok, Token};
use crate::units::UnitExpr;

const KEYWORDS: &[&str] = &[
    "let", "fn", "if", "then", "else", "and", "or", "not", "to", "true", "false",
];

pub struct Parser {
    toks: Vec<Token>,
    pos: usize,
}

type PResult<T> = Result<T, String>;

impl Parser {
    pub fn new(toks: Vec<Token>) -> Parser {
        Parser { toks, pos: 0 }
    }

    fn peek(&self) -> &Tok {
        &self.toks[self.pos].tok
    }

    fn peek_at(&self, k: usize) -> &Tok {
        &self.toks[(self.pos + k).min(self.toks.len() - 1)].tok
    }

    fn tok(&self) -> &Token {
        &self.toks[self.pos]
    }

    fn line(&self) -> usize {
        self.toks[self.pos].line
    }

    fn advance(&mut self) -> Tok {
        let t = self.toks[self.pos].tok.clone();
        if self.pos < self.toks.len() - 1 {
            self.pos += 1;
        }
        t
    }

    fn err<T>(&self, msg: impl Into<String>) -> PResult<T> {
        Err(format!("line {}: {}", self.line(), msg.into()))
    }

    fn is_sym(&self, s: &str) -> bool {
        matches!(self.peek(), Tok::Sym(x) if *x == s)
    }

    fn is_kw(&self, s: &str) -> bool {
        matches!(self.peek(), Tok::Ident(x) if x == s)
    }

    fn eat_sym(&mut self, s: &str) -> bool {
        if self.is_sym(s) {
            self.advance();
            true
        } else {
            false
        }
    }

    fn expect_sym(&mut self, s: &str) -> PResult<()> {
        if self.eat_sym(s) {
            Ok(())
        } else {
            self.err(format!("expected '{}', found {}", s, describe(self.peek())))
        }
    }

    fn skip_newlines(&mut self) {
        while matches!(self.peek(), Tok::Newline) || self.is_sym(";") {
            self.advance();
        }
    }

    fn ident(&mut self) -> PResult<String> {
        match self.peek().clone() {
            Tok::Ident(s) if !KEYWORDS.contains(&s.as_str()) => {
                self.advance();
                Ok(s)
            }
            t => self.err(format!("expected a name, found {}", describe(&t))),
        }
    }

    pub fn program(&mut self) -> PResult<Vec<Stmt>> {
        let mut stmts = vec![];
        self.skip_newlines();
        while !matches!(self.peek(), Tok::Eof) {
            stmts.push(self.stmt()?);
            if !matches!(self.peek(), Tok::Eof | Tok::Newline) && !self.is_sym(";") {
                return self.err(format!("unexpected {}", describe(self.peek())));
            }
            self.skip_newlines();
        }
        Ok(stmts)
    }

    fn stmt(&mut self) -> PResult<Stmt> {
        let line = self.line();
        let kind = if self.is_kw("let") {
            self.advance();
            let name = self.ident()?;
            self.expect_sym("=")?;
            self.skip_newlines();
            StmtKind::Let(name, self.expr()?)
        } else if self.is_kw("fn") && matches!(self.peek_at(1), Tok::Ident(_)) {
            self.advance();
            let name = self.ident()?;
            let params = self.params("(", ")")?;
            let body = if self.is_sym("{") {
                self.primary()?
            } else {
                self.expect_sym("=")?;
                self.skip_newlines();
                self.expr()?
            };
            StmtKind::Fn(name, params, body)
        } else {
            StmtKind::Expr(self.expr()?)
        };
        Ok(Stmt { kind, line })
    }

    fn params(&mut self, open: &str, close: &str) -> PResult<Vec<String>> {
        self.expect_sym(open)?;
        let mut ps = vec![];
        while !self.is_sym(close) {
            ps.push(self.ident()?);
            if !self.eat_sym(",") {
                break;
            }
        }
        self.expect_sym(close)?;
        Ok(ps)
    }

    pub fn expr(&mut self) -> PResult<Expr> {
        let mut e = self.or()?;
        while self.eat_sym("->") {
            let u = self.unit_expr(true)?;
            e = Expr::Convert(Box::new(e), u);
        }
        Ok(e)
    }

    fn or(&mut self) -> PResult<Expr> {
        let mut e = self.and()?;
        while self.is_kw("or") {
            self.advance();
            self.skip_newlines();
            e = Expr::Binary("or", Box::new(e), Box::new(self.and()?));
        }
        Ok(e)
    }

    fn and(&mut self) -> PResult<Expr> {
        let mut e = self.not()?;
        while self.is_kw("and") {
            self.advance();
            self.skip_newlines();
            e = Expr::Binary("and", Box::new(e), Box::new(self.not()?));
        }
        Ok(e)
    }

    fn not(&mut self) -> PResult<Expr> {
        if self.is_kw("not") {
            self.advance();
            return Ok(Expr::Unary("not", Box::new(self.not()?)));
        }
        self.cmp()
    }

    fn cmp(&mut self) -> PResult<Expr> {
        let e = self.pm()?;
        for op in ["==", "!=", "<=", ">=", "<", ">"] {
            if self.is_sym(op) {
                self.advance();
                self.skip_newlines();
                let r = self.pm()?;
                return Ok(Expr::Binary(op, Box::new(e), Box::new(r)));
            }
        }
        Ok(e)
    }

    fn pm(&mut self) -> PResult<Expr> {
        let mut l = self.add()?;
        let is_to = self.is_kw("to");
        if self.is_sym("±") || is_to {
            self.advance();
            self.skip_newlines();
            let r = self.add()?;
            // `12 ± 0.5 m` and `2 to 5 t`: a bare left literal shares the right literal's unit.
            if let (Expr::Num(_, lu @ None), Expr::Num(_, Some(ru))) = (&mut l, &r) {
                if ru.0 != vec![("%".to_string(), Rational::ONE)] {
                    *lu = Some(ru.clone());
                }
            }
            if is_to {
                return Ok(Expr::To(Box::new(l), Box::new(r)));
            }
            if let Expr::Num(p, Some(u)) = &r {
                if u.0 == vec![("%".to_string(), Rational::ONE)] {
                    return Ok(Expr::PlusMinusRel(Box::new(l), p / 100.0));
                }
            }
            return Ok(Expr::PlusMinus(Box::new(l), Box::new(r)));
        }
        Ok(l)
    }

    fn add(&mut self) -> PResult<Expr> {
        let mut e = self.mul()?;
        loop {
            let op = if self.is_sym("+") {
                "+"
            } else if self.is_sym("-") {
                "-"
            } else {
                break;
            };
            self.advance();
            self.skip_newlines();
            e = Expr::Binary(op, Box::new(e), Box::new(self.mul()?));
        }
        Ok(e)
    }

    fn mul(&mut self) -> PResult<Expr> {
        let mut e = self.unary()?;
        loop {
            let op = if self.is_sym("*") || self.is_sym("·") {
                "*"
            } else if self.is_sym("/") {
                "/"
            } else {
                break;
            };
            self.advance();
            self.skip_newlines();
            e = Expr::Binary(op, Box::new(e), Box::new(self.unary()?));
        }
        Ok(e)
    }

    fn unary(&mut self) -> PResult<Expr> {
        if self.eat_sym("-") {
            return Ok(Expr::Unary("-", Box::new(self.unary()?)));
        }
        if self.eat_sym("+") {
            return self.unary();
        }
        self.pow()
    }

    fn pow(&mut self) -> PResult<Expr> {
        let base = self.postfix()?;
        if self.eat_sym("^") {
            let exp = self.unary()?;
            return Ok(Expr::Binary("^", Box::new(base), Box::new(exp)));
        }
        Ok(base)
    }

    /// Is the next token the start of a unit (an identifier that isn't a keyword, or `%`)?
    fn at_unit(&self) -> bool {
        match self.peek() {
            Tok::Ident(s) => !KEYWORDS.contains(&s.as_str()),
            Tok::Sym("%") => true,
            _ => false,
        }
    }

    /// Parse a unit expression such as `kg/m^3`, `t·m` or `m^-1`. Unless `spaced` is set,
    /// the expression ends at the first whitespace, so `5 m / t` divides by a variable `t`.
    fn unit_expr(&mut self, spaced: bool) -> PResult<UnitExpr> {
        let mut u = self.unit_term()?;
        loop {
            let tight = !self.tok().space_before && !self.toks[self.pos + 1].space_before;
            let follows_unit = matches!(self.peek_at(1), Tok::Ident(_) | Tok::Sym("%"));
            if !(tight || spaced) || !follows_unit {
                break;
            }
            if self.is_sym("/") {
                self.advance();
                u = u.mul(&self.unit_term()?.inv());
            } else if self.is_sym("*") || self.is_sym("·") {
                self.advance();
                u = u.mul(&self.unit_term()?);
            } else {
                break;
            }
        }
        Ok(u)
    }

    fn unit_term(&mut self) -> PResult<UnitExpr> {
        let name = match self.advance() {
            Tok::Ident(s) => s,
            Tok::Sym("%") => "%".to_string(),
            t => return self.err(format!("expected a unit, found {}", describe(&t))),
        };
        if crate::units::lookup_unit(&name).is_none() {
            return self.err(format!(
                "unknown unit '{}' (to multiply by a variable, write '*' explicitly)",
                name
            ));
        }
        let mut exp = Rational::ONE;
        if self.is_sym("^") && !self.tok().space_before {
            self.advance();
            let neg = self.eat_sym("-");
            match self.advance() {
                Tok::Num(n) => {
                    exp = Rational::from_f64(if neg { -n } else { n })
                        .ok_or_else(|| format!("line {}: bad unit exponent", self.line()))?;
                }
                Tok::Sym("(") => {
                    // m^(1/2)
                    let a = match self.advance() {
                        Tok::Num(n) => n,
                        _ => return self.err("bad unit exponent"),
                    };
                    self.expect_sym("/")?;
                    let b = match self.advance() {
                        Tok::Num(n) => n,
                        _ => return self.err("bad unit exponent"),
                    };
                    self.expect_sym(")")?;
                    exp = Rational::new(a as i32, b as i32);
                    if neg {
                        exp = -exp;
                    }
                }
                _ => return self.err("bad unit exponent"),
            }
        }
        Ok(UnitExpr(vec![(name, exp)]))
    }

    fn postfix(&mut self) -> PResult<Expr> {
        let mut e = self.primary()?;
        loop {
            if self.is_sym("(") && !self.tok().space_before {
                self.advance();
                let mut args = vec![];
                while !self.is_sym(")") {
                    args.push(self.expr()?);
                    if !self.eat_sym(",") {
                        break;
                    }
                }
                self.expect_sym(")")?;
                e = Expr::Call(Box::new(e), args);
            } else if self.is_sym(".") {
                self.advance();
                let f = self.ident()?;
                e = Expr::Field(Box::new(e), f);
            } else if self.is_sym("[") && !self.tok().space_before {
                self.advance();
                let i = self.expr()?;
                self.expect_sym("]")?;
                e = Expr::Index(Box::new(e), Box::new(i));
            } else {
                break;
            }
        }
        // `(12 ± 0.5) m` and `[1, 2, 3] m`
        let prev = &self.toks[self.pos.saturating_sub(1)].tok;
        if matches!(prev, Tok::Sym(")") | Tok::Sym("]"))
            && self.at_unit()
            && !matches!(e, Expr::Call(..))
        {
            let u = self.unit_expr(false)?;
            e = Expr::ApplyUnit(Box::new(e), u);
        }
        Ok(e)
    }

    fn primary(&mut self) -> PResult<Expr> {
        let line = self.line();
        match self.peek().clone() {
            Tok::Num(n) => {
                self.advance();
                let unit = if self.at_unit() {
                    Some(self.unit_expr(false)?)
                } else {
                    None
                };
                Ok(Expr::Num(n, unit))
            }
            Tok::Str(s) => {
                self.advance();
                Ok(Expr::Str(s))
            }
            Tok::Ident(s) if s == "if" => {
                self.advance();
                let c = self.expr()?;
                self.skip_newlines();
                if !self.is_kw("then") {
                    return self.err("expected 'then'");
                }
                self.advance();
                self.skip_newlines();
                let a = self.expr()?;
                self.skip_newlines();
                if !self.is_kw("else") {
                    return self.err("expected 'else'");
                }
                self.advance();
                self.skip_newlines();
                let b = self.expr()?;
                Ok(Expr::If(Box::new(c), Box::new(a), Box::new(b)))
            }
            Tok::Ident(s) if s == "true" || s == "false" => {
                self.advance();
                Ok(Expr::Var(s))
            }
            Tok::Ident(s) if s == "fn" => {
                self.advance();
                let ps = self.params("(", ")")?;
                self.skip_newlines();
                if self.is_sym("=") || self.is_sym("=>") {
                    self.advance();
                }
                self.skip_newlines();
                Ok(Expr::Lambda(ps, Box::new(self.expr()?)))
            }
            Tok::Ident(_) => Ok(Expr::Var(self.ident()?)),
            Tok::Sym("|") => {
                let ps = self.params("|", "|")?;
                Ok(Expr::Lambda(ps, Box::new(self.expr()?)))
            }
            Tok::Sym("(") => {
                self.advance();
                let e = self.expr()?;
                self.expect_sym(")")?;
                Ok(e)
            }
            Tok::Sym("[") => {
                self.advance();
                let mut items = vec![];
                while !self.is_sym("]") {
                    items.push(self.expr()?);
                    if !self.eat_sym(",") {
                        break;
                    }
                }
                self.expect_sym("]")?;
                Ok(Expr::List(items))
            }
            Tok::Sym("{") => {
                self.advance();
                self.skip_newlines();
                let is_record = self.is_sym("}")
                    || (matches!(self.peek(), Tok::Ident(_) | Tok::Str(_))
                        && matches!(self.peek_at(1), Tok::Sym(":")));
                if is_record {
                    let mut fields = vec![];
                    while !self.is_sym("}") {
                        let name = match self.advance() {
                            Tok::Ident(s) | Tok::Str(s) => s,
                            t => {
                                return self
                                    .err(format!("expected field name, found {}", describe(&t)))
                            }
                        };
                        self.expect_sym(":")?;
                        self.skip_newlines();
                        fields.push((name, self.expr()?));
                        self.skip_newlines();
                        if !self.eat_sym(",") {
                            break;
                        }
                        self.skip_newlines();
                    }
                    self.skip_newlines();
                    self.expect_sym("}")?;
                    Ok(Expr::Record(fields))
                } else {
                    let mut stmts = vec![];
                    while !self.is_sym("}") {
                        stmts.push(self.stmt()?);
                        self.skip_newlines();
                    }
                    self.expect_sym("}")?;
                    Ok(Expr::Block(stmts))
                }
            }
            t => Err(format!("line {}: unexpected {}", line, describe(&t))),
        }
    }
}

fn describe(t: &Tok) -> String {
    match t {
        Tok::Num(n) => format!("number {}", n),
        Tok::Ident(s) => format!("'{}'", s),
        Tok::Str(s) => format!("string \"{}\"", s),
        Tok::Sym(s) => format!("'{}'", s),
        Tok::Newline => "end of line".into(),
        Tok::Eof => "end of input".into(),
    }
}
