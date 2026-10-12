use crate::units::UnitExpr;

#[derive(Clone, Debug)]
pub enum Expr {
    Num(f64, Option<UnitExpr>),
    Str(String),
    Var(String),
    List(Vec<Expr>),
    Record(Vec<(String, Expr)>),
    Block(Vec<Stmt>),
    Unary(&'static str, Box<Expr>),
    Binary(&'static str, Box<Expr>, Box<Expr>),
    /// `a ± b`: normal with mean a and standard deviation b.
    PlusMinus(Box<Expr>, Box<Expr>),
    /// `a ± p%`: normal with sd = p·|a|.
    PlusMinusRel(Box<Expr>, f64),
    /// `a to b`: lognormal (or normal if a ≤ 0) with a 90% interval of [a, b].
    To(Box<Expr>, Box<Expr>),
    ApplyUnit(Box<Expr>, UnitExpr),
    Convert(Box<Expr>, UnitExpr),
    Call(Box<Expr>, Vec<Expr>),
    Field(Box<Expr>, String),
    Index(Box<Expr>, Box<Expr>),
    If(Box<Expr>, Box<Expr>, Box<Expr>),
    Lambda(Vec<String>, Box<Expr>),
}

#[derive(Clone, Debug)]
pub enum StmtKind {
    Let(String, Expr),
    Fn(String, Vec<String>, Expr),
    Expr(Expr),
}

#[derive(Clone, Debug)]
pub struct Stmt {
    pub kind: StmtKind,
    pub line: usize,
}
