using LinearAlgebra

"""
    generalized_diff_mat(xs::Vector{Float64})

Calcula a matriz de diferenciação de Lagrange generalizada para um dado
vetor de pontos `xs`.

Esta é uma tradução para Julia do método compacto em MATLAB para
gerar matrizes de diferenciação.

# Argumentos
- `xs::Vector{Float64}`: Um vetor coluna com os pontos (nós) do grid.

# Retorna
- `D::Matrix{Float64}`: A matriz de diferenciação de dimensão (N+1) x (N+1).
"""
function generalized_diff_mat(xs::Vector{Float64})
    n = length(xs) - 1
    # Garante que xs seja um vetor coluna para as operações seguintes
    xs_col = reshape(xs, :, 1)

    # Replica o vetor de pontos para formar a matriz de diferenças
    X = repeat(xs_col, 1, n + 1)
    # Calcula a matriz de diferenças xi - xj, adicionando I para evitar zeros na diagonal
    dX = (X - X') + I

    # Calcula os pesos a_j = prod(xj - xk) de forma vetorial
    aj = vec(prod(dX', dims=1))

    # Calcula os elementos fora da diagonal D_ij = a_i / (a_j * (xi - xj))
    # A diagonal fica temporariamente com 1
    D = (aj ./ aj') ./ dX

    # Calcula os elementos da diagonal forçando a soma da linha a ser zero
    # D_jj = -sum(D_jk) para k != j
    D = D - Diagonal(vec(sum(D, dims=2)))
    
    return D
end

# --- CÓDIGO DO TESTE ---

# Definição das funções e suas derivadas exatas
f1(x) = exp(-x)
f1_prime(x) = -exp(-x)

f2(x) = sin(pi*x)
f2_prime(x) = pi*cos(pi*x)

f3(x) = 1e4 * cos(pi*x)
f3_prime(x) = -1e4 * pi * sin(pi*x)

f4(x) = 32x^6 - 48x^4 + 18x^2 - 1
f4_prime(x) = 192x^5 - 192x^3 + 36x

# Lista de funções para o teste
test_cases = [
    ("e⁻ˣ", f1, f1_prime),
    ("sin(πx)", f2, f2_prime),
    ("10⁴cos(πx)", f3, f3_prime),
    ("Polinômio T₆(x)", f4, f4_prime)
]

# Definição dos grids
grids = [
    ("Grid 1 (12 pontos)", -1.0 .+ (0:11) .* (2/11)),
    ("Grid 2 (121 pontos)", -1.0 .+ (0:120) .* (2/120))
]

println("--- Iniciando Teste de Exatidão da Matriz de Diferenciação Generalizada ---")

for (grid_name, xs) in grids
    println("\nAnalisando com $grid_name:")
    
    # Gera a matriz de diferenciação para o grid atual
    D = generalized_diff_mat(collect(xs))
    
    for (func_name, f, f_prime) in test_cases
        # Avalia a função e sua derivada analítica nos pontos do grid
        f_vals = f.(xs)
        f_prime_analytic_vals = f_prime.(xs)
        
        # Calcula a derivada numérica usando a matriz D
        f_prime_numeric_vals = D * f_vals
        
        # Calcula o erro máximo absoluto
        max_error = maximum(abs.(f_prime_numeric_vals - f_prime_analytic_vals))
        
        println("  Função: $(rpad(func_name, 18)) | Erro Máximo: $max_error")
    end
end