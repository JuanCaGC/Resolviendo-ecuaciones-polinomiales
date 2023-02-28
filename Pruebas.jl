using Singular

R, (x, y, z) = PolynomialRing(QQ, ["x","y","z"],ordering = :degrevlex)
f1 = x^2 - 2 * x * z + 5
f2 = x * y^2 + y * z + 1
f3 = 3 * y^2 - 8 * x * z
I = Ideal(R,[f1,f2,f3])

# Calcular la base de Groebner del ideal
@time G = std(I)


@time KB = kbase(G)

ge = collect(gens(G))
kb = collect(gens(KB))

function funcionPrueba(variable)
    Mx = []
    for k in kb
        arr=[]
        if !(variable*k in kb)
            r = variable*k
            monos = collect(monomials(r))
            coeffs = collect(coefficients(r))
            i = 1
            while !(all(b -> in(b, kb), monos)) 
                g = ge[i]
                
                println(r)
                
                q,r = divrem(r,g)
                i+=1
                if i > length(ge)
                    i = 1
                end
                monos = collect(monomials(r))
                coeffs = collect(coefficients(r))
            end
            
            println("Ultimo r")
            println(r)
            
            
            for ka2 in kb
                if ka2 in monos 
                    for j in 1:length(monos)
                        if ka2 == monos[j]
                            push!(arr,coeffs[j])
                        end
                    end
                else
                    push!(arr,0)
                end
            end
        else
            for ka in kb
                if variable*k == ka
                    push!(arr,1)
                else
                    push!(arr,0)
                end
            end
        end
        push!(Mx,arr)
        
         println(arr)
        
    end
    Mx = reduce(hcat,Mx)
    return Mx
end

Mx = funcionPrueba(x)

println(Mx)