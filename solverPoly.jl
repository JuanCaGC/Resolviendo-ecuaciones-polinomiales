using Oscar
using LinearAlgebra



function dimQI(v,lista)
    max = 0
    leadmons = Set([leading_monomial(f) for f in lista])
    for p in leadmons
        if v in vars(p)
            m = degree(p, v)
            if m> max
                max = m
            end
        end
    end
    return max
end

function kbasis(GB, VList)
    # devuelve una lista de monomios que forman una base del anillo cociente, donde GB es una base de Gröbner para un ideal de dimensión cero, y genera un mensaje de error si el ideal no es de dimensión cero.
    
    if iszero(dim(GB))  # verifica si el ideal es de dimensión cero
        Gb = collect(gens(GB))
        B = [(div(VList[1],VList[1])) ]  # inicializa la lista de monomios de la base con el monomio 1
        leadmons = Set([leading_monomial(f) for f in Gb])  # crea un conjunto de los monomios principales de cada polinomio en la base de Gröbner
        for v in VList
            m = dimQI(v, leadmons)  # determina el grado de la variable v en el anillo cociente
            C = copy(B)  # crea una copia de la lista de monomios actual como base para la siguiente iteración
            for t in C
                for l in 1:m-1
                    t = t*v  # multiplica el monomio actual por la variable v
                    if  !(any((divrem(t,u)[2] == 0) for u in leadmons))  # verifica si alguno de los monomios principales de la base de Gröbner divide al nuevo monomio
                        push!(B, t)  # si ningún monomio principal lo divide, agrega el nuevo monomio a la lista de monomios de la base
                    end
                end
            end
        end
        return B  # devuelve la lista de monomios de la base
    else
        error("el ideal no es de dimensión cero")
    end
end

function matrixM(variable,kbase,Gbase)
    kbx = kbase*variable
    reducida = reduce(collect(gens(kbx)),collect(Gbase))
    Mx = [coeff(kr,k) for k in collect(gens(kbase)) , kr in collect(reducida)] 
    M = [Int64.(numerator(i))//Int64.(denominator(i)) for i in Mx]
    return M
end

function buenaCombinacion(ListaM,E0,F0)
    n, m = size(ListaM[1])
    Ei = E0
    Fi = F0
    S0 = [diag(Fi*Mi*Ei) for Mi in ListaM]
    S = hvcat(length(ListaM), S0...)
    sigma = [ℯ^(1*im*(t-1)*π/n) for t in 1:n]
    Σ=diagm(sigma)
    alpha=pinv(S)*sigma
    M =sum([alpha[i]*ListaM[i] for i in 1:length(ListaM)])
    return M,Σ
end

function XYS(Σ,Z,Δ)
    n, m = size(Z)
    S = diagm(diag(Δ-Z*Σ))
    X = [i != j ? (-Δ[i,j]+ Z[i,j] * Σ[j,j])/(Σ[i,i]-Σ[j,j]) : 0 for i in 1:n , j in 1:n]
    Y = [i != j ? (Δ[i,j]- Z[i,j] * Σ[i,i])/(Σ[i,i]-Σ[j,j]) : -Z[i,i] for i in 1:n , j in 1:n]
    return X,Y,S
end

function iteracionNewton(E,F,Σ,M)
    n, m = size(M)
    I_n = Matrix{Float64}(I, n, n)
    Z = F*E-I_n
    Δ = F * M * E - Σ
    X,Y,S = XYS(Σ,Z,Δ)
    Em1 = E * (I_n + X)
    Fm1 = (I_n + Y) * F
    Σm1 = Σ + S
    return Em1,Fm1,Σm1
end

function newtonDiagSimul(M,E0,F0,Σ0)
    n, m = size(M)
    Ei = E0
    Fi = F0
    Σi = Σ0
    err = opnorm((Fi*M*Ei-Σi), Inf)
    errores = [err]
    numbern = 0
    while err > 10^(-15) 
        Ei,Fi,Σi = iteracionNewton(Ei,Fi,Σi,M)
        err = opnorm((Fi*M*Ei-Σi), Inf)
        numbern+=1
        append!(errores, err)
    end
    return Ei,Fi,numbern,Σi, errores
end


n = 3 # Numero de variables

# Inicializar el anillo
R, x = PolynomialRing(QQ, ["x$i" for i in 1:n])

# Polinomios del sistema
f1 = x[1]^2 - 2 * x[1] * x[3] + 5
f2 = x[1] * x[2]^2 + x[2] * x[3] + 1
f3 = 3 * x[2]^2 - 8 * x[1] * x[3]
Id = ideal(R,[f1,f2,f3])

G = groebner_basis(Id)

Gi = ideal(G) 

KB = kbasis(Gi,gens(R))

kb = ideal(R,KB)

# Lista de matrices de multiplicacion
ListaM = [matrixM(xi,kb,G) for xi in gens(R) ]

# Inicializaciones para la diagonalizacion silmultanea
L0,E0 = eigen(ListaM[1])
F0 = inv(E0)
M0,Σ0 = buenaCombinacion(ListaM,E0,F0)

# Metodo tipo Newton
E,F,niter,Sig,errores = newtonDiagSimul(M0,E0,F0,Σ0)

# Los puntos que solucionan el sistema 
# el punto ai = (puntosCompilados[i,j]) para 1<=j<=n
puntos = [diag(F*Mi*E) for Mi in ListaM]
puntosCompilados = hvcat(length(ListaM), puntos...)

println("Los puntos solucion del sistema son los siguientes: ")
[println("a" * string(i) * " =( " * join(puntosCompilados[i,:], " , ") * " ) ") for i in 1:size(puntosCompilados, 1)]

"Fin"