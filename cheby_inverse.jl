using LinearAlgebra

# --------------------------------------------------------------------
# 1. Função para montar as matrizes B e B⁻¹ (Slide 7)
# --------------------------------------------------------------------
"""
    montar_matrizes_chebyshev(n::Int)

Monta a matriz de base de Chebyshev `B` e sua inversa `B⁻¹` para um 
polinômio de grau `n`, usando os `n+1` nodos de Chebyshev-Lobatto.
"""
function montar_matrizes_chebyshev(n::Int)
    # Gerar os nodos de Chebyshev-Lobatto de -1 a 1.
    x = [-cos(k * pi / n) for k in 0:n]
    
    # Montar a Matriz de Base B, onde B[i, j] = T_{j-1}(x_i)
    B = zeros(n + 1, n + 1)
    for i in 1:(n + 1)
        B[i, 1] = 1.0
        if n > 0
            B[i, 2] = x[i]
        end
        for j in 2:n
            B[i, j + 1] = 2 * x[i] * B[i, j] - B[i, j - 1]
        end
    end

    # Montar a Matriz Inversa B⁻¹ analiticamente (via DCT-I)
    B_inv = zeros(n + 1, n + 1)
    w = ones(n + 1); w[1] = 0.5; w[end] = 0.5
    α = 2.0 * ones(n + 1); α[1] = 1.0; α[end] = 1.0
    for k in 1:(n + 1)
        for j in 1:(n + 1)
            B_inv[k, j] = (α[k] / n) * B[j, k] * w[j]
        end
    end
    
    return B, B_inv, x
end


# --------------------------------------------------------------------
# 2. Script de teste focado em Análise e Reconstrução
# --------------------------------------------------------------------

println("--- Exercício I: Teste de Análise e Reconstrução com Matrizes ---")

# Definir o grau do polinômio
n = 8

# Montar as matrizes e obter os nodos
B, B_inv, x_nodos = montar_matrizes_chebyshev(n)

# Verificar se B * B⁻¹ é a matriz identidade
identidade = I(n + 1)
erro_inversa = norm(B * B_inv - identidade)
println("Grau do Polinômio n = $n")
println("Verificação da Inversa: ||B * B⁻¹ - I|| = $erro_inversa")
@assert erro_inversa < 1e-14 "A matriz inversa não foi calculada corretamente!"
println("Verificação OK! As matrizes são inversas uma da outra.\n")

# Funções para testar
funcoes_teste = [
    (f = x -> exp(-x), nome = "e⁻ˣ"),
    (f = x -> sin(pi * x), nome = "sin(πx)"),
    (f = x -> 10^4 * cos(pi * x), nome = "10⁴cos(πx)"),
    (f = x -> 32x^6 - 48x^4 + 18x^2 - 1, nome = "T₆(x)")
]

for (f, nome) in funcoes_teste
    println("--- Testando a função: $nome ---")
    
    # 1. Obter os valores da função nos nodos
    f_valores_originais = f.(x_nodos)
    
    # 2. ANÁLISE: Transformar valores em coeficientes com B⁻¹
    coeficientes = B_inv * f_valores_originais
    
    # 3. SÍNTESE: Transformar coeficientes de volta em valores com B
    f_valores_reconstruidos = B * coeficientes
    
    # 4. Verificar o erro da reconstrução
    erro_reconstrucao = maximum(abs.(f_valores_originais - f_valores_reconstruidos))
    
    println("Primeiros 5 coeficientes calculados: ", round.(coeficientes[1:5], digits=6))
    println("Erro máximo da reconstrução: ", erro_reconstrucao)
    println()
end