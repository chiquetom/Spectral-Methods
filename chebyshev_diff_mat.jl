using LinearAlgebra, Plots

# --------------------------------------------------------------------
# 1. Ferramentas necessárias (baseadas nos slides e exercícios)
# --------------------------------------------------------------------

"""
    chebyshev_diff_mat(N::Int)
Retorna a matriz de diferenciação espectral de Chebyshev D e os nodos x
de Gauss-Lobatto para um polinômio de grau N.
"""
function chebyshev_diff_mat(N::Int)
    x = [cos(j * pi / N) for j in 0:N] # Nodos de 1 a -1
    c = ones(N + 1); c[1] = 2; c[end] = 2
    
    D = zeros(N + 1, N + 1)
    for i in 0:N
        for j in 0:N
            if i == j
                if i == 0
                    D[i+1, j+1] = (2*N^2 + 1) / 6
                elseif i == N
                    D[i+1, j+1] = -(2*N^2 + 1) / 6
                else
                    D[i+1, j+1] = -x[j+1] / (2 * (1 - x[j+1]^2))
                end
            else
                D[i+1, j+1] = (c[i+1] / c[j+1]) * ((-1)^(i+j) / (x[i+1] - x[j+1]))
            end
        end
    end
    return D, x
end

# Funções dos exercícios anteriores para obter e avaliar os coeficientes
function montar_matrizes_chebyshev(n::Int)
    x = [-cos(k * pi / n) for k in 0:n]
    B = zeros(n + 1, n + 1)
    for i in 1:(n + 1)
        B[i, 1] = 1.0
        if n > 0; B[i, 2] = x[i]; end
        for j in 2:n
            B[i, j + 1] = 2 * x[i] * B[i, j] - B[i, j - 1]
        end
    end
    B_inv = inv(B) # Usar inv() é mais simples aqui
    return B, B_inv, x
end

function avaliar_serie_chebyshev(coeficientes::Vector, x_aval)
    n = length(coeficientes) - 1
    y = zeros(length(x_aval))
    for i in eachindex(x_aval)
        xi = x_aval[i]
        b_k2, b_k1 = 0.0, 0.0
        for k in n:-1:1
            b_k = coeficientes[k+1] + 2*xi*b_k1 - b_k2
            b_k2 = b_k1
            b_k1 = b_k
        end
        y[i] = coeficientes[1] + xi*b_k1 - b_k2
    end
    return y
end


# --------------------------------------------------------------------
# 2. Resolução da Equação Diferencial (Exercício III)
# --------------------------------------------------------------------

println("--- Iniciando Exercício III: Resolução Numérica da EDO ---")

# Grau da série de Chebyshev (conforme enunciado)
N = 8
println("Usando polinômio de grau N = $N\n")

# --- Passo 1: Obter matrizes de diferenciação e nodos ---
# Note que `cheb_diff_mat` gera os nodos de 1 a -1.
D, x_rev = chebyshev_diff_mat(N)
D2 = D * D
x = reverse(x_rev) # Invertemos para ter de -1 a 1, mais intuitivo

# --- Passo 2: Montar o operador linear L da EDO ---
# EDO: (x²/20)y'' + x*y' - y = 5x⁵ - 1
I_mat = I(N + 1)
L = (diagm(x.^2) / 20) * D2 + diagm(x) * D - I_mat

# --- Passo 3: Montar o vetor do lado direito F ---
F = 5 * x.^5 .- 1

# --- Passo 4: Impor as condições de contorno ---
# y(-1) = 2  (corresponde ao primeiro ponto, x[1])
# y(1) = 0   (corresponde ao último ponto, x[end])

# Modifica as linhas do sistema L*Y = F
# Linha para y(-1): a primeira equação se torna 1*Y[1] = 2
L[1, :] .= 0.0
L[1, 1] = 1.0
F[1] = 2.0

# Linha para y(1): a última equação se torna 1*Y[end] = 0
L[end, :] .= 0.0
L[end, end] = 1.0
F[end] = 0.0

# --- Passo 5: Resolver o sistema linear ---
# Método 1: O resultado são os valores da solução nos nodos
Y = L \ F

println("Solução nos nodos de Chebyshev (Y):")
display(Y)

# --- Passo 6: Análise e Comparação ---
# Método 2: Obter a representação espectral (coeficientes)
_, B_inv, _ = montar_matrizes_chebyshev(N)
coeficientes = B_inv * Y

println("\nCoeficientes da Série de Chebyshev para a solução Y(x):")
display(coeficientes)

# --- Passo 7: Visualizar a solução ---
# Criar um grid fino para plotar uma curva suave
x_fino = range(-1, 1, length=200)
# Avaliar o polinômio nos pontos finos usando os coeficientes
y_fino = avaliar_serie_chebyshev(coeficientes, x_fino)

# Plotar a solução
p = plot(x_fino, y_fino, label="Solução Interpolada (Método 2)", lw=2, legend=:topright)
scatter!(p, x, Y, label="Solução nos Nodos (Método 1)", markersize=4)
title!("Solução Numérica da EDO (N=$N)")
xlabel!("x")
ylabel!("y(x)")
display(p)