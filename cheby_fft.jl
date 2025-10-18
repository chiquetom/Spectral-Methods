using FFTW, Plots

# --------------------------------------------------------------------
# 1. Função principal baseada no Slide 10
# --------------------------------------------------------------------
"""
    coeficientes_chebyshev_fft(f_valores)

Calcula os coeficientes da série de Chebyshev a partir dos valores da função `f_valores`
nos nodos de Chebyshev-Gauss-Lobatto (extremos inclusos), usando o método da FFT.

A implementação segue a lógica do código em MATLAB do slide 10.
"""
function coeficientes_chebyshev_fft(f_valores::Vector{Float64})
    N = length(f_valores)
    
    # Passo 1: Inverter o vetor para seguir a direção de θ (slide 10: "flipud(A)")
    # Os nodos de Chebyshev x = cos(θ) vão de 1 a -1 para θ de 0 a π.
    # O ifft espera a ordem natural de θ.
    A_rev = reverse(f_valores)
    
    # Passo 2: Construir o vetor estendido e simétrico para a FFT
    # Corresponde a: `[A(1:N,:); A(N-1:-1:2,:)]`
    v_ext = [A_rev; A_rev[end-1:-1:2]]
    
    # Passo 3: Aplicar a IFFT (Inverse Fast Fourier Transform)
    F = ifft(v_ext)
    
    # Passo 4: Extrair os coeficientes reais, ajustando a escala
    # Corresponde a: `B=([F(1,:); 2*F(2:(N-1),:); F(N,:)])`
    # Apenas a parte real nos interessa.
    c = real.([F[1]; 2 * F[2:N-1]; F[N]])
    
    return c
end


# --------------------------------------------------------------------
# 2. Função auxiliar para testar e visualizar os resultados
# --------------------------------------------------------------------
"""
    avaliar_serie_chebyshev(coeficientes, x)

Avalia a série de Chebyshev com dados `coeficientes` no ponto `x`
usando o algoritmo de Clenshaw para estabilidade e eficiência.
"""
function avaliar_serie_chebyshev(coeficientes::Vector{Float64}, x::Float64)
    n = length(coeficientes) - 1
    bk2, bk1 = 0.0, 0.0
    
    for k in n:-1:1
        bk = coeficientes[k+1] + 2 * x * bk1 - bk2
        bk2 = bk1
        bk1 = bk
    end
    # Caso especial para k=0
    b0 = coeficientes[1] + x * bk1 - bk2
    return b0
end


# --------------------------------------------------------------------
# 3. Script para testar as funções fornecidas
# --------------------------------------------------------------------

# Definição do grau do polinômio de aproximação
n = 16 # Grau do polinômio (resulta em n+1 pontos)

# Gerar os nodos de Chebyshev-Gauss-Lobatto (de -1 a 1)
# x_k = cos(kπ/n) para k=n..0, o que é igual a -cos(jπ/n) para j=0..n
x_nodos = [-cos(j * pi / n) for j in 0:n]

# Funções para testar (da imagem)
funcoes_teste = [
    (f = x -> exp(-x), nome = "e⁻ˣ"),
    (f = x -> sin(pi * x), nome = "sin(πx)"),
    (f = x -> 10^4 * cos(pi * x), nome = "10⁴cos(πx)"),
    (f = x -> 32x^6 - 48x^4 + 18x^2 - 1, nome = "T₆(x) = 32x⁶-48x⁴+18x²-1")
]

# Grid fino para plotar os resultados com alta resolução
x_plot = range(-1, 1, length=400)

println("--- Iniciando Testes de Aproximação com Polinômios de Chebyshev (Grau n=$n) ---")

# Iterar sobre cada função, calcular coeficientes, e plotar o resultado
for (f, nome) in funcoes_teste
    println("\n--- Testando a função: $nome ---")
    
    # a) Obter os valores da função nos nodos
    f_valores = f.(x_nodos)
    
    # b) Calcular os coeficientes de Chebyshev usando a FFT
    coeficientes = coeficientes_chebyshev_fft(f_valores)
    println("Primeiros 5 coeficientes: ", round.(coeficientes[1:5], digits=6))
    
    # c) Reconstruir a função no grid de plotagem usando os coeficientes
    y_aproximado = [avaliar_serie_chebyshev(coeficientes, xi) for xi in x_plot]
    y_exato = f.(x_plot)
    
    # d) Calcular e exibir o erro máximo da aproximação
    erro_max = maximum(abs.(y_exato - y_aproximado))
    println("Erro máximo da aproximação: ", erro_max)
    
    # e) Plotar os resultados
    p = plot(x_plot, y_exato, label="Função Original", lw=3, legend=:topleft)
    plot!(p, x_plot, y_aproximado, label="Aproximação de Chebyshev", ls=:dash, lw=2)
    scatter!(p, x_nodos, f_valores, label="Nodos (Pontos de Amostra)", markersize=4)
    title!("Aproximação para f(x) = $nome")
    xlabel!("x")
    ylabel!("f(x)")
    display(p)
end