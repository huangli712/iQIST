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

    for i in eachindex(COB)
        lᵢ = l
        mᵢ = COB[i][1]
        sᵢ = COB[i][2]
        for j in eachindex(COB)
            lⱼ = l
            mⱼ = COB[j][1]
            sⱼ = COB[j][2]

            T₁ = (lⱼ - mⱼ)*(lⱼ + mⱼ + 1)
            m1ⱼ = mⱼ + 1
            if sⱼ == "up"
                s1ⱼ = "down"
            else
                s1ⱼ = "null"
            end
            if m1ⱼ == mᵢ && sᵢ == s1ⱼ
                @show i, j, " sqrt ", T₁
            end

            T₂ = (lⱼ + mⱼ)*(lⱼ - mⱼ + 1)
            m2ⱼ = mⱼ - 1
            if sⱼ == "down"
                s2ⱼ = "up"
            else
                s2ⱼ = "null"
            end
            if m2ⱼ == mᵢ && sᵢ == s2ⱼ
                @show i, j, " sqrt ", T₂
            end
            
            T₃ = mⱼ * (sⱼ == "up" ? 1 : -1)
            if mⱼ == mᵢ && sᵢ == sⱼ
                @show i, j, T₃
            end
        end
    end
end

l = 2
calc_matrix(l)