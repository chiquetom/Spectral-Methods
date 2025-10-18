# ============================================================
# Tarefa complementar da Aula 4 — Itens (a) e (b) numéricos + comparação (c)
# ============================================================

using Printf

# EDO do item (a): u'(x) = (4x^2 - 2u)/x, com u(1) = 2, em x ∈ [1,3]
f(x,u) = (4*x^2 - 2*u)/x

# Passo de Runge–Kutta 4
function rk4_step(f, x, u, h)
    k1 = f(x, u)
    k2 = f(x + h/2, u + h*k1/2)
    k3 = f(x + h/2, u + h*k2/2)
    k4 = f(x + h,   u + h*k3)
    return u + (h/6)*(k1 + 2k2 + 2k3 + k4)
end

# Integra u'(x)=f(x,u) de x0→x1 (x1 ≥ x0) com passo fixo h
function integrate_to(f, x0, u0, x1; h=1e-4)
    @assert x1 ≥ x0 "integrador supõe x1 ≥ x0"
    x = x0
    u = u0
    while x < x1
        hstep = min(h, x1 - x)
        u = rk4_step(f, x, u, hstep)
        x += hstep
    end
    return u
end

# g(x) = u(x) - 4, computado numericamente integrando de 1 até x
g_of_x(x; h=1e-4) = integrate_to(f, 1.0, 2.0, x; h=h) - 4.0

# Bisseção para encontrar x* em [a,b] com precisão alvo em x
function bisect_root(a, b; h=1e-4, xtol=5e-13, maxit=200)
    ga = g_of_x(a; h=h)
    gb = g_of_x(b; h=h)
    @assert ga*gb ≤ 0 "g(a) e g(b) devem ter sinais opostos em [a,b]"
    for k in 1:maxit
        m = (a + b)/2
        gm = g_of_x(m; h=h)
        if (b - a) < xtol
            return m, k, gm
        end
        if ga*gm ≤ 0
            b = m; gb = gm
        else
            a = m; ga = gm
        end
    end
    m = (a + b)/2
    return m, maxit, g_of_x(m; h=h)
end

# ---------------- Execução (a) e (b) ----------------

# (a) Integração em malha (opcional, para inspeção)
xs = range(1.0, 3.0; length=2001)
us = similar(collect(xs))
u = 2.0
us[1] = u
for i in 2:length(xs)
    h = xs[i] - xs[i-1]
    u = rk4_step(f, xs[i-1], u, h)
    us[i] = u
end
@printf("Item (a): u(1)=%.6f, u(3)≈%.12f (numérico RK4)\n", us[1], us[end])

# (b) Encontrar x* com u(x*)=4 por bisseção (10+ dígitos)
xstar, it, gval = bisect_root(1.0, 3.0; h=1e-4, xtol=5e-13, maxit=200)
u_xstar = g_of_x(xstar; h=1e-5) + 4.0   # reavalia para conferir u(x*)

@printf("Item (b): x* ≈ %.12f (iters=%d, g=%.3e)\n", xstar, it, gval)
@printf("Verificação: u(x*) ≈ %.12f (alvo = 4.000000000000)\n", u_xstar)

# ---------------- Comparação (c) ----------------
# Apenas verificação/relato analítico — não usado no algoritmo
xk = sqrt(2 + sqrt(3))
err_abs = abs(xstar - xk)
err_rel = err_abs / abs(xk)
xk_alt = (sqrt(6) + sqrt(2))/2

@printf("Item (c): xk (analítico) = %.12f\n", xk)
@printf("Item (c): |x* - xk| = %.3e (erro relativo = %.3e)\n", err_abs, err_rel)
@printf("Item (c): xk_alt = %.12f, |xk - xk_alt| = %.3e\n", xk_alt, abs(xk - xk_alt))

# ---------------- Fim ----------------
