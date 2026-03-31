#include <bits/stdc++.h>
using namespace std;

// We implement a simple symbolic engine for expressions built from terms a*x^b*sin^c(x)*cos^d(x),
// combined by +,-,*,/, and parentheses. We output a single fraction p/q and its derivative h.

struct Term {
    long long a; int b,c,d;
};

struct Poly { // sum of terms
    vector<Term> t;
};

static void simplify_poly(Poly &p){
    // merge like terms (same b,c,d)
    sort(p.t.begin(), p.t.end(), [](const Term&x,const Term&y){
        if (x.b!=y.b) return x.b>y.b;
        if (x.c!=y.c) return x.c>y.c;
        if (x.d!=y.d) return x.d>y.d;
        return x.a>y.a;
    });
    vector<Term> out; out.reserve(p.t.size());
    for (auto &x: p.t){
        if (!out.empty() && out.back().b==x.b && out.back().c==x.c && out.back().d==x.d){
            out.back().a += x.a;
            if (out.back().a==0) out.pop_back();
        }else{
            if (x.a!=0) out.push_back(x);
        }
    }
    p.t.swap(out);
}

static Poly poly_add(const Poly&a,const Poly&b){ Poly r=a; r.t.insert(r.t.end(), b.t.begin(), b.t.end()); simplify_poly(r); return r; }
static Poly poly_sub(const Poly&a,const Poly&b){ Poly r=a; for(auto x:b.t){ Term y=x; y.a=-y.a; r.t.push_back(y);} simplify_poly(r); return r; }
static Poly poly_mul(const Poly&a,const Poly&b){
    Poly r; r.t.reserve((size_t)a.t.size()*b.t.size());
    for(auto &x:a.t) for(auto &y:b.t){
        Term z; z.a=x.a*y.a; z.b=x.b+y.b; z.c=x.c+y.c; z.d=x.d+y.d; if(z.a) r.t.push_back(z);
    }
    simplify_poly(r); return r;
}

static Poly poly_deriv(const Poly& a){
    Poly r;
    for (auto &x: a.t){
        // derivative of a*x^b*sin^c*cos^d = a*( b x^{b-1} sin^c cos^d + c x^b sin^{c-1} cos * + (-d) x^b sin * cos^{d-1} )
        if (x.b>0){ Term y{x.a*(long long)x.b, x.b-1, x.c, x.d}; if(y.a) r.t.push_back(y);}        
        if (x.c>0){ Term y{x.a*(long long)x.c, x.b, x.c-1, x.d+1}; if(y.a) r.t.push_back(y);}      
        if (x.d>0){ Term y{x.a*(long long)(-x.d), x.b, x.c+1, x.d-1}; if(y.a) r.t.push_back(y);}   
    }
    simplify_poly(r); return r;
}

struct Frac{ Poly p,q; }; // p/q

static Frac make_int(long long v){ Frac f; f.p.t.push_back({v,0,0,0}); f.q.t.push_back({1,0,0,0}); return f; }
static Frac make_term(Term t){ Frac f; f.p.t.push_back(t); f.q.t.push_back({1,0,0,0}); return f; }

static Frac frac_add(const Frac&a,const Frac&b){ return { poly_add(poly_mul(a.p,b.q), poly_mul(b.p,a.q)), poly_mul(a.q,b.q)}; }
static Frac frac_sub(const Frac&a,const Frac&b){ return { poly_sub(poly_mul(a.p,b.q), poly_mul(b.p,a.q)), poly_mul(a.q,b.q)}; }
static Frac frac_mul(const Frac&a,const Frac&b){ return { poly_mul(a.p,b.p), poly_mul(a.q,b.q)}; }
static Frac frac_div(const Frac&a,const Frac&b){ return { poly_mul(a.p,b.q), poly_mul(a.q,b.p)}; }

static Frac frac_deriv(const Frac& f){ // (p/q)' = (p'*q - q'*p)/(q*q)
    Poly p1 = poly_mul(poly_deriv(f.p), f.q);
    Poly p2 = poly_mul(poly_deriv(f.q), f.p);
    Poly num = poly_sub(p1, p2);
    Poly den = poly_mul(f.q, f.q);
    return {num, den};
}

// Pretty printing rules
static string term_to_str(const Term& t, bool first){
    // handle coefficient and signs
    long long a = t.a;
    string s;
    if (first){
        if (a<0) { s.push_back('-'); a = -a; }
    }else{
        if (a<0){ s += "-"; a=-a; }
        else s += "+";
    }
    bool isConst = (t.b==0 && t.c==0 && t.d==0);
    if (!(a==1 && !isConst)) s += to_string(a);
    if (!isConst){
        if (!(a==1 && !isConst)) ; // coefficient already printed, no '*'
        // x part
        if (t.b){ s += "x"; if (t.b!=1){ s += "^"+to_string(t.b);} }
        // sin part
        if (t.c){ s += "sinx"; if (t.c!=1){ s += "^"+to_string(t.c);} }
        // cos part
        if (t.d){ s += "cosx"; if (t.d!=1){ s += "^"+to_string(t.d);} }
    }
    if (s.empty()) s = "0"; // shouldn't happen
    return s;
}

static string poly_to_str(const Poly& p){
    if (p.t.empty()) return "0";
    // already simplified and sorted descending by (b,c,d)
    // ensure sorted order per spec
    vector<Term> v = p.t;
    sort(v.begin(), v.end(), [](const Term&x,const Term&y){
        if (x.b!=y.b) return x.b>y.b;
        if (x.c!=y.c) return x.c>y.c;
        if (x.d!=y.d) return x.d>y.d;
        return x.a>y.a;
    });
    string s;
    for (size_t i=0;i<v.size();++i){ s += term_to_str(v[i], i==0); }
    return s;
}

static string frac_to_str(const Frac& f){
    string p = poly_to_str(f.p);
    string q = poly_to_str(f.q);
    if (q=="1") return p;
    if (p=="0") return string("0");
    // s1/s2 format with parentheses unless q==1 handled
    return string("(")+p+")/"+q;
}

// Lexer / Parser
struct Lexer{
    const string s; int n; int i=0;
    Lexer(const string&str):s(str),n((int)str.size()){}
    void skip(){ while(i<n && isspace((unsigned char)s[i])) ++i; }
    bool match(char c){ skip(); if (i<n && s[i]==c){ ++i; return true;} return false; }
    bool peek(char c){ skip(); return i<n && s[i]==c; }
};

static Frac parse_expr(Lexer&L);

static Frac parse_factor(Lexer&L){
    L.skip();
    // number or term like x, sinx, cosx, with optional exponent ^k
    bool neg=false;
    if (L.match('(')){
        Frac inside = parse_expr(L);
        if (!L.match(')')) ;
        // exponent?
        if (L.match('^')){
            // integer exponent >=1 assumed small
            long long k=0; bool negpow=false; if (L.match('-')) negpow=true;
            while(L.i<L.n && isdigit((unsigned char)L.s[L.i])){ k = k*10 + (L.s[L.i++]-'0'); }
            // raise inside^k by repeated multiplication
            if (negpow){ // negative power translates to 1/(inside^{|k|})
                // compute inside^{|k|}
                Frac pow = make_int(1);
                long long e = k;
                for (long long t=0;t<e;++t) pow = frac_mul(pow, inside);
                // invert
                Frac one = make_int(1);
                inside = frac_div(one, pow);
            }else{
                Frac pow = make_int(1);
                for (long long t=0;t<k;++t) pow = frac_mul(pow, inside);
                inside = pow;
            }
        }
        return inside;
    }
    // unary + or -
    if (L.match('+')) return parse_factor(L);
    if (L.match('-')){
        Frac x = parse_factor(L);
        Frac z = make_int(0); return frac_sub(z, x);
    }
    // try sinx, cosx, x, integer
    if (L.i+3<=L.n && L.s.substr(L.i,3)=="sin"){
        L.i+=3; // expect x
        if (L.i<L.n && L.s[L.i]=='x') L.i++;
        Term t{1,0,1,0};
        // optional ^k
        if (L.match('^')){
            long long k=0; while(L.i<L.n && isdigit((unsigned char)L.s[L.i])){ k=k*10+(L.s[L.i++]-'0'); }
            t.c=(int)k;
        }
        return make_term(t);
    }
    if (L.i+3<=L.n && L.s.substr(L.i,3)=="cos"){
        L.i+=3; if (L.i<L.n && L.s[L.i]=='x') L.i++;
        Term t{1,0,0,1};
        if (L.match('^')){ long long k=0; while(L.i<L.n && isdigit((unsigned char)L.s[L.i])){ k=k*10+(L.s[L.i++]-'0'); } t.d=(int)k; }
        return make_term(t);
    }
    if (L.i<L.n && L.s[L.i]=='x'){
        L.i++;
        Term t{1,1,0,0};
        if (L.match('^')){ long long k=0; while(L.i<L.n && isdigit((unsigned char)L.s[L.i])){ k=k*10+(L.s[L.i++]-'0'); } t.b=(int)k; }
        return make_term(t);
    }
    // integer possibly with sign already handled above
    bool sign=false; long long val=0; bool has=false;
    if (L.i<L.n && (isdigit((unsigned char)L.s[L.i]))){ has=true; while(L.i<L.n && isdigit((unsigned char)L.s[L.i])){ val=val*10+(L.s[L.i++]-'0'); }}
    if (has){ return make_int(val); }
    // fallback zero
    return make_int(0);
}

static Frac parse_powered_atom(Lexer&L){
    // factor may already handle ^ for grouped, but for simple atoms like x, sinx, cosx handled too.
    return parse_factor(L);
}

static Frac parse_term(Lexer&L){
    // multiplication by juxtaposition or explicit *
    Frac cur = parse_powered_atom(L);
    for(;;){
        L.skip();
        if (L.i>=L.n) break;
        char c=L.s[L.i];
        if (c=='+'||c=='-'||c==')'||c=='/') break;
        if (c=='*'){ L.i++; }
        Frac nxt = parse_powered_atom(L);
        cur = frac_mul(cur, nxt);
    }
    return cur;
}

static Frac parse_expr(Lexer&L){
    Frac cur = parse_term(L);
    for(;;){
        L.skip(); if (L.i>=L.n) break;
        if (L.s[L.i]=='+'){ L.i++; Frac t=parse_term(L); cur=frac_add(cur,t); }
        else if (L.s[L.i]=='-'){ L.i++; Frac t=parse_term(L); cur=frac_sub(cur,t); }
        else if (L.s[L.i]=='/'){
            L.i++; Frac t=parse_term(L); cur=frac_div(cur,t);
        } else break;
    }
    return cur;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string s; if(!(cin>>s)) return 0;
    Lexer L(s);
    Frac f = parse_expr(L);
    string g = frac_to_str(f);
    Frac df = frac_deriv(f);
    string h = frac_to_str(df);
    cout<<g<<"\n"<<h<<"\n";
    return 0;
}

