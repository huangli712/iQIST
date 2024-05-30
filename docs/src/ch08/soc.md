# Spin-orbit coupling

**Spin-orbit interaction Hamiltonian**

The spin-orbit coupling (SOC) is implemented at atomic level,

```math
\begin{equation}
\hat{H}_{\text{SOC}}=
\lambda\sum_{i}\hat{l}_{i} \cdot \hat{s}_{i},
\end{equation}
```

where ``\hat{l}`` is orbital angular momentum, and ``\hat{s}`` is spin angular momentum. In second quantization form,

```math
\begin{equation}
\hat{H}_{\text{SOC}}=
\lambda\sum_{\alpha\sigma,\beta\sigma^{\prime}}
\left\langle
\alpha\sigma\left| \hat{l} \cdot \hat{s} \right|\beta\sigma^{\prime}
\right\rangle
\hat{f}_{\alpha\sigma}^{\dagger}\hat{f}_{\beta\sigma^{\prime}},
\end{equation}
```

where ``\alpha`` is orbital index and ``\sigma`` is spin index. We note that

```math
\begin{equation}
\hat{l} \cdot \hat{s} =
\frac{1}{2} \hat{l} \cdot \hat{\sigma}.
\end{equation}
```

where ``\hat{\sigma}`` is the Pauli operator:

```math
\begin{equation}
\hat{\sigma} = \hat{\sigma}_x \hat{x} +
               \hat{\sigma}_y \hat{y} +
               {\sigma}_z \hat{z}.
\end{equation}
```

```math
\begin{equation}
\hat{\sigma}_x = \left[
\begin{array}{cc}
0 & 1 \\
1 & 0 \\
\end{array}
\right],~
\hat{\sigma}_y = \left[
\begin{array}{cc}
0 & -i \\
i & 0 \\
\end{array}
\right],~
\hat{\sigma}_z = \left[
\begin{array}{cc}
1 & 0 \\
0 & -1 \\
\end{array}
\right].
\end{equation}
```

Now the question is how to write down the matrix elements for ``\hat{l} \cdot \hat{\sigma}``.

```math
\begin{equation}
\hat{l} \cdot \hat{\sigma} =
\hat{l}_{x}\hat{\sigma}_{x}+
\hat{l}_{y}\hat{\sigma}_{y}+
\hat{l}_{z}\hat{\sigma}_{z}
= \left[\begin{array}{cc}
\hat{l}_{z} & \hat{l}_{-}\\
\hat{l}_{+} & -\hat{l}_{z}
\end{array}\right],
\end{equation}
```

where ``\hat{l}_{\pm}=\hat{l}_{x}\pm\hat{l}_{y}``. It is easy to prove that the above equation can be rewritten as

```math
\begin{equation}
\hat{l} \cdot \hat{\sigma} =
\hat{l}_{+}\hat{\sigma}_{-}+
\hat{l}_{-}\hat{\sigma}_{+}+
\hat{l}_{z}\hat{\sigma}_{z},
\end{equation}
```

where

```math
\begin{equation}
\hat{\sigma}_+ = \left[
\begin{array}{cc}
0 & 1 \\
0 & 0 \\
\end{array}
\right],~
\hat{\sigma}_- = \left[
\begin{array}{cc}
0 & 0 \\
1 & 0 \\
\end{array}
\right].
\end{equation}
```

**Spin-orbit interaction in the complex orbital basis**

We just write down ``\hat{l}\cdot\hat{\sigma}`` in the complex spherical harmonics basis. The following relations are used to derive the matrix elements of the spin-orbit interaction:

```math
\begin{equation}
\hat{l}_{\pm}Y_{l}^{m}=\sqrt{(l\mp m)(l\pm m+1)}Y_{l}^{m\pm1}.
\end{equation}
```

```math
\begin{equation}
\hat{\sigma}_+ |\uparrow\rangle = 0,~
\hat{\sigma}_+ |\downarrow\rangle = |\uparrow\rangle.
\end{equation}
```

```math
\begin{equation}
\hat{\sigma}_- |\uparrow\rangle = |\downarrow\rangle,~
\hat{\sigma}_- |\downarrow\rangle = 0.
\end{equation}
```

```math
\begin{equation}
\hat{\sigma}_z |\uparrow\rangle = 1,~
\hat{\sigma}_z |\downarrow\rangle = -1.
\end{equation}
```

For ``p`` system,
```math
\begin{equation}
\hat{l}\cdot\hat{\sigma}=
\left[
\begin{array}{ccc|ccc}
-1 & 0 & 0 & 0 & \sqrt{2} & 0\\
0 & 0 & 0 & 0 & 0 & \sqrt{2}\\
0 & 0 & 1 & 0 & 0 & 0\\
\hline
0 & 0 & 0 & 1 & 0 & 0\\
\sqrt{2} & 0 & 0 & 0 & 0 & 0\\
0 & \sqrt{2} & 0 & 0 & 0 & -1\\
\end{array}
\right]
\end{equation}
```

For ``d`` system,
```math
\begin{equation}
\hat{l}\cdot\hat{\sigma}=
\left[
\begin{array}{ccccc|ccccc}
-2 & 0 & 0 & 0 & 0 & 0 & \sqrt{4} & 0 & 0 & 0\\
0 & -1 & 0 & 0 & 0 & 0 & 0 & \sqrt{6} & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{6} & 0\\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & \sqrt{4}\\
0 & 0 & 0 & 0 & 2 & 0 & 0 & 0 & 0 & 0\\
\hline
0 & 0 & 0 & 0 & 0 & 2 & 0 & 0 & 0 & 0\\
\sqrt{4} & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\
0 & \sqrt{6} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & \sqrt{6} & 0 & 0 & 0 & 0 & 0 & -1 & 0\\
0 & 0 & 0 & \sqrt{4} & 0 & 0 & 0 & 0 & 0 & -2\\
\end{array}
\right]
\end{equation}
```

For ``f`` system,
```math
\begin{equation}
\hat{l}\cdot\hat{\sigma}=
\left[
\begin{array}{ccccccc|ccccccc}
-3 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{6} & 0 & 0 & 0 & 0 & 0\\
0 & -2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{10} & 0 & 0 & 0 & 0\\
0 & 0 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{12} & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{12} & 0 & 0\\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{10} & 0\\
0 & 0 & 0 & 0 & 0 & 2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{6}\\
0 & 0 & 0 & 0 & 0 & 0 & 3 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
\hline
0 & 0 & 0 & 0 & 0 & 0 & 0 & 3 & 0 & 0 & 0 & 0 & 0 & 0\\
\sqrt{6} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 2 & 0 & 0 & 0 & 0 & 0\\
0 & \sqrt{10} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\
0 & 0 & \sqrt{12} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & \sqrt{12} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & -1 & 0 & 0\\
0 & 0 & 0 & 0 & \sqrt{10} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & -2 & 0\\
0 & 0 & 0 & 0 & 0 & \sqrt{6} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & -3\\
\end{array}
\right]
\end{equation}
```

The following Julia's script can be used to generate the matrix elements for the spin-orbit interaction.

```julia
# To calculate the matrix representation of spin-orbit coupling in the
# complex orbital basis.
function calc_matrix(l::Int64)
    println("Construct complex orbital basis for 𝑙 = $l")
    COB = [] # To save the complex orbital basis
    # m = -l, -l+1, ..., l-1, l
    mlist = collect(-l:1:l)
    for s in ("up", "down")
        for m in mlist
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
```

**Spin-orbit interaction in the ``\hat{j}^{2}-\hat{j}_{z}-\hat{l}^2-\hat{s}^2`` diagonal basis**

Since

```math
\begin{equation}
\hat{j} = \hat{l} + \hat{s},
\end{equation}
```

```math
\hat{j}^2 = \hat{l}^2 + \hat{s}^2 + 2\hat{l} \cdot \hat{s},
```

```math
\hat{l} \cdot \hat{\sigma} = 2 \hat{l} \cdot \hat{s}
                           = \hat{j}^2 - \hat{l}^2 - \hat{s}^2,
```

thus, in the ``\hat{j}^{2}-\hat{j}_{z}-\hat{l}^2-\hat{s}^2`` diagonal basis, the matrix of ``\hat{l} \cdot \hat{\sigma}`` is diagonal. The diagonal element reads ``j(j+1) - l(l+1) - s(s+1)``. Note that in this basis, ``j = l \pm \frac{1}{2}`` and ``s = \frac{1}{2}``.

For ``p`` system,
```math
\begin{equation}
\hat{l}\cdot\hat{\sigma}=
\left[
\begin{array}{ccc|ccc}
-2 & 0 & 0 & 0 & 0 & 0\\
0 & -2 & 0 & 0 & 0 & 0\\
0 & 0 & 1 & 0 & 0 & 0\\
\hline
0 & 0 & 0 & 1 & 0 & 0\\
0 & 0 & 0 & 0 & 1 & 0\\
0 & 0 & 0 & 0 & 0 & 1\\
\end{array}
\right]
\end{equation}
```

For ``d`` system,
```math
\begin{equation}
\hat{l}\cdot\hat{\sigma}=
\left[
\begin{array}{ccccc|ccccc}
-3 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & -3 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & -3 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & -3 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 2 & 0 & 0 & 0 & 0 & 0\\
\hline
0 & 0 & 0 & 0 & 0 & 2 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 2 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 2 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 2 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 2\\
\end{array}
\right]
\end{equation}
```

For ``f`` system,
```math
\begin{equation}
\hat{l}\cdot\hat{\sigma}=
\left[
\begin{array}{ccccccc|ccccccc}
-4 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & -4 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & -4 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & -4 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & -4 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & -4 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 3 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
\hline
0 & 0 & 0 & 0 & 0 & 0 & 0 & 3 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 3 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 3 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 3 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 3 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 3 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 3\\
\end{array}
\right]
\end{equation}
```

Next, it is easy to transform the above matrices from the ``\hat{j}^{2}-\hat{j}_{z}-\hat{l}^2-\hat{s}^2`` diagonal basis to complex orbital basis. See [Single particle basis](basis.md) for the transformation matrix.
