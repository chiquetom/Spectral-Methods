using Plots
############################ 1. FUNÇÃO DE INTERPOLAÇÃO DE LAGRANGE
# (Inspirada no code do slide)
function lagrange_interp(xk, yk, xv)
"""
Calcula os valores interpolados `yv` para um conjunto de pontos `xv`,
usando a interpolação de Lagrange com base nos nós de interpolação (`xk`, `yk`).

A função implementa a fórmula P(x) = Σ [yk_j * L_j(x)], onde L_j(x) é o j-ésimo
polinômio da base de Lagrange.

# Argumentos
- `xk`: Vetor com as coordenadas x dos nós de interpolação.
- `yk`: Vetor com as coordenadas y dos nós de interpolação (valores de f(xk)).
- `xv`: Vetor de pontos onde a interpolação deve ser calculada.

# Retorna
- `yv`: Um vetor com os valores interpolados (`yv`) correspondentes a cada ponto em `xv`.
"""
    N = length(xk)
    M = length(xv)
    yv = zeros(M)

    for i in 1:M
        x_atual = xv[i]
        soma = 0.0
        for j in 1:N
            produtorio_base = 1.0
            for m in 1:N
                if j != m
                    produtorio_base *= (x_atual - xk[m]) / (xk[j] - xk[m])
                end
            end
            soma += yk[j] * produtorio_base
        end
        yv[i] = soma
    end
    return yv
end

############################ 2. PLOT E ERRO
function analisar_e_plotar(nodes, f, titulo_funcao, titulo_base)
"""
A função automatiza o processo de:
1. Calcular os valores da interpolação sobre um intervalo denso de pontos.
2. Calcular o erro máximo absoluto entre a função real e a curva interpolada.
3. Imprimir o resultado do erro no console.
4. Gerar e exibir um plot comparativo.

# Argumentos
- `nodes`: Vetor de nós (coordenadas x) para construir o polinômio interpolador.
- `f`: A função original a ser interpolada (ex: `x -> sin(pi * x)`).
- `titulo_funcao`: String contendo o nome da função para ser usado no título do gráfico.
- `titulo_base`: String descrevendo a base de nós utilizada (ex: "Base: 12 Nós").
"""
    println("\n--- ", titulo_base, " ---")

    # Gera os valores y para os nós (yk)
    y_nodes = f.(nodes)

    # Cria um conjunto de pontos para avaliar a interpolação (xv)
    pontos_teste = range(-1, 1, length=500)

    # Usa a nova função para obter os valores interpolados
    y_interpolado = lagrange_interp(nodes, y_nodes, pontos_teste)

    # Compara com os valores reais da função
    y_real = f.(pontos_teste)

    # Calcula o erro
    erro_maximo = maximum(abs.(y_real .- y_interpolado))
    println("Erro Máximo Absoluto: ", erro_maximo)
    titulo_plot = "$titulo_funcao\n$titulo_base\nErro Máximo: $(round(erro_maximo, sigdigits=4))"

    # Gera o gráfico
    p = plot(pontos_teste, y_real, label="Função Real", title=titulo_plot, lw=2)
    plot!(p, pontos_teste, y_interpolado, label="Interpolação", linestyle=:dash, lw=2)
    scatter!(p, nodes, y_nodes, label="Nós")
    display(p)
end

############################ 3. TESTES COM AS FUNÇÕES E AS SEQUENCIAS X_K E X_I
xk_nodes = [-1 + k * (2/11) for k in 0:11]
xi_nodes = [-1 + i * (2/119) for i in 0:119]

funcoes_para_testar = [
    (x -> exp(-x), "f(x) = exp(-x)"),
    (x -> sin(pi * x), "f(x) = sin(πx)"),
    (x -> 10 * cos(pi * x), "f(x) = 10*cos(πx)"),
    (x -> 32*x^6 - 48*x^4 + 18*x^2 - 1, "f(x) = 32x⁶ - 48x⁴ + 18x² - 1")
]

# ==========================================================
# 4. EXECUÇÃO PRINCIPAL
# ==========================================================
println("==========================================================")
println("Iniciando Análise Comparativa de Interpolação de Lagrange")
println("==========================================================")

for (f, titulo) in funcoes_para_testar
    println("\n\n**********************************************************")
    println("Analisando Função: ", titulo)
    println("**********************************************************")

    analisar_e_plotar(xk_nodes, f, titulo, "Base de Interpolação: 12 Nós (xk)")
    analisar_e_plotar(xi_nodes, f, titulo, "Base de Interpolação: 120 Nós (xi)")
end

println("\n------ Fim ------")