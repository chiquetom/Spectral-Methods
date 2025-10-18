using Plots
using Printf

############################ 1A. FUNÇÃO DE INTERPOLAÇÃO DE LAGRANGE
function lagrange_interp(xk, yk, xv)
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

############################ 1B. FUNÇÃO DE INTERPOLAÇÃO BARICÊNTRICA
function interpolacao_baricentrica(xk, yk, xj)
    n = length(xk) - 1
    wk = ones(Float64, n + 1)
    for k in 1:(n + 1)
        diferencas = [xk[k] - xk[i] for i in 1:(n+1) if i != k]
        wk[k] = 1.0 / prod(diferencas)
    end

    yj = zeros(Float64, length(xj))
    for i in 1:length(xj)
        x_eval = xj[i]
        indices_exatos = findall(x -> x == x_eval, xk)
        if !isempty(indices_exatos)
            yj[i] = yk[first(indices_exatos)]
            continue
        end
        termos = wk ./ (x_eval .- xk)
        numerador = sum(termos .* yk)
        denominador = sum(termos)
        yj[i] = numerador / denominador
    end
    return yj
end


############################ 2. PLOT E ERRO COMPARATIVO
function analisar_e_comparar(nodes, f, titulo_funcao, titulo_base)
    println("\n--- ", titulo_base, " ---")

    # Gera os valores y para os nós (yk)
    y_nodes = f.(nodes)

    # Cria um conjunto denso de pontos para avaliar a interpolação (xv)
    pontos_teste = range(-1, 1, length=500)
    y_real = f.(pontos_teste)

    # Calcula a interpolação com ambos os métodos
    y_lagrange = lagrange_interp(nodes, y_nodes, pontos_teste)
    y_baricentrico = interpolacao_baricentrica(nodes, y_nodes, pontos_teste)

    # Calcula os erros
    erro_lagrange = maximum(abs.(y_real .- y_lagrange))
    erro_baricentrico = maximum(abs.(y_real .- y_baricentrico))
    
    println("Erro Máximo Absoluto (Lagrange):     ", erro_lagrange)
    println("Erro Máximo Absoluto (Baricêntrico): ", erro_baricentrico)

    # Prepara o título do gráfico
    titulo_plot = "$titulo_funcao\n$titulo_base"
    
    # Gera o gráfico
    p = plot(pontos_teste, y_real, label="Função Real", title=titulo_plot, lw=3)
    plot!(p, pontos_teste, y_lagrange, label="Lagrange (Erro: $(@sprintf "%.2e" erro_lagrange))", linestyle=:dash, lw=2)
    plot!(p, pontos_teste, y_baricentrico, label="Baricêntrico (Erro: $(@sprintf "%.2e" erro_baricentrico))", linestyle=:dot, lw=2, color=:red)
    scatter!(p, nodes, y_nodes, label="Nós", markersize=3)
    display(p)
end

############################ 3. TESTES COM AS FUNÇÕES E AS SEQUENCIAS X_K E X_I
xk_nodes = [-1 + k * (2/11) for k in 0:11]
xi_nodes = [-1 + i * (2/119) for i in 0:119]

funcoes_para_testar = [
    (x -> exp(-x), "f(x) = exp(-x)"),
    (x -> sin(pi * x), "f(x) = sin(πx)"),
    (x -> 10^4 * cos(pi * x), "f(x) = 10⁴cos(πx)"),
    (x -> 32*x^6 - 48*x^4 + 18*x^2 - 1, "f(x) = T₆(x)")
]

# ==========================================================
# 4. EXECUÇÃO PRINCIPAL
# ==========================================================
println("==========================================================")
println("  Análise Comparativa: Lagrange vs. Baricêntrico")
println("==========================================================")

for (f, titulo) in funcoes_para_testar
    println("\n\n**********************************************************")
    println("Analisando Função: ", titulo)
    println("**********************************************************")

    analisar_e_comparar(xk_nodes, f, titulo, "Base de Interpolação: 12 Nós (xk)")
    analisar_e_comparar(xi_nodes, f, titulo, "Base de Interpolação: 120 Nós (xi)")
end

println("\n------ Fim da Análise ------")