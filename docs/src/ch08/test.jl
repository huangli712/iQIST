# To calculate the matrix representation of spin-orbit coupling in the
# complex orbital basis.
function calc_matrix(l::Int64)
    println("Construct complex orbital basis for 𝑙 = $l")
    COB = [] # To save the complex orbital basis
    # m = -l, -l+1, ..., l-1, l
    mlist = collect(-l:1:l)
    for m in mlist
    for s in ("up", "down")
        
            push!(COB, [m, s])
        end
    end
    #
    for i in eachindex(COB)
        m = COB[i][1]
        s = COB[i][2] == "up" ? "↑" : "↓"
        println("$i -> | $l, $m, $s ⟩")
    end

    println("Calculate matrix elements of spin-orbit coupling")
    for i in eachindex(COB)
        lᵢ = l
        mᵢ = COB[i][1]
        sᵢ = COB[i][2]
        for j in eachindex(COB)
            lⱼ = l
            mⱼ = COB[j][1]
            sⱼ = COB[j][2]

            # For l₊σ₋ term
            T₁ = (lⱼ - mⱼ) * (lⱼ + mⱼ + 1)
            m1ⱼ = mⱼ + 1
            #
            if sⱼ == "up"
                s1ⱼ = "down"
            else
                s1ⱼ = "null"
            end
            #
            if m1ⱼ == mᵢ && sᵢ == s1ⱼ
                println("SOC($i, $j) = sqrt(", T₁, ")")
            end

            # For l₋σ₊ term
            T₂ = (lⱼ + mⱼ)*(lⱼ - mⱼ + 1)
            m2ⱼ = mⱼ - 1
            #
            if sⱼ == "down"
                s2ⱼ = "up"
            else
                s2ⱼ = "null"
            end
            #
            if m2ⱼ == mᵢ && sᵢ == s2ⱼ
                println("SOC($i, $j) = sqrt(", T₂, ")")
            end
            
            # For lzσz term
            T₃ = mⱼ * (sⱼ == "up" ? 1 : -1)
            if mⱼ == mᵢ && sᵢ == sⱼ
                println("SOC($i, $j) = ", T₃)
            end
        end
    end
end

l = 1
calc_matrix(l)
l = 2
calc_matrix(l)
l = 3
calc_matrix(l)