# iQIST (Interacting Quantum Impurity Solver Toolkit)

## JASMINE Test Cases 

### 1-band model

* 1a // Uc = 4.0, Uv = 4.0, Jz = 0.0, Js = 0.0, Jp = 0.0, NO GQNs 

* 1b // Uc = 4.0, Uv = 4.0, Jz = 0.0, Js = 0.0, Jp = 0.0, GQNs (N) 

* 1c // Uc = 4.0, Uv = 4.0, Jz = 0.0, Js = 0.0, Jp = 0.0, GQNs (N,Sz) 

* 1d // Uc = 4.0, Uv = 4.0, Jz = 0.0, Js = 0.0, Jp = 0.0, GQNs (N,Sz,PS) 

### 3-band model

* 3a // Uc = 4.0, Uv = 4.0, Jz = 0.0, Js = 0.0, Jp = 0.0, NO GQNs 

* 3b // Uc = 4.0, Uv = 2.0, Jz = 1.0, Js = 0.0, Jp = 0.0, NO GQNs 

* 3c // Uc = 4.0, Uv = 2.0, Jz = 1.0, Js = 1.0, Jp = 1.0, NO GQNs 

* 3d // Uc = 4.0, Uv = 2.0, Jz = 1.0, Js = 1.0, Jp = 1.0, diagonal CF, NO GQNs 

* 3e // Uc = 4.0, Uv = 2.0, Jz = 1.0, Js = 1.0, Jp = 1.0, non-diagonal CF, NO GQNs 

* 3f // Uc = 4.0, Uv = 2.0, Jz = 1.0, Js = 1.0, Jp = 1.0, mune=2.00, SOC, lambda = 0.50, NO GQNs 

* 3g // Uc = 4.0, Uv = 2.0, Jz = 1.0, Js = 1.0, Jp = 1.0, diagonal CF, SOC, lambda = 0.50, NO GQNs 

* 3h // Uc = 4.0, Uv = 2.0, Jz = 1.0, Js = 1.0, Jp = 1.0, diagonal CF, GQNs (N,Sz) 

* 3i // Uc = 4.0, Uv = 2.0, Jz = 1.0, Js = 1.0, Jp = 1.0, diagonal CF, GQNs (N,Sz,PS) 

* 3j // Uc = 4.0, Uv = 2.0, Jz = 1.0, Js = 1.0, Jp = 1.0, mune=2.00, SOC, lambda = 0.50, GQNs (N,Jz)
 
* 3k // Uc = 4.0, Uv = 2.0, Jz = 1.0, Js = 1.0, Jp = 1.0, diagonal CF, SOC, lambda = 0.50, GQNs (N) 

* 3l // Uc = 4.0, Uv = 2.0, Jz = 1.0, Js = 1.0, Jp = 1.0, diagonal CF, SOC, lambda = 0.50, GQNs (N), nmini=2, nmaxi=4

### 5-band model

* 5a // Uc = 4.0, Uv = 2.0, Jz = 1.0, Js = 1.0, Jp = 1.0, GQNs (N) 

* 5b // Ud = 4.0, Jh = 1.0, GQNs (N) 

* 5c // Ud = 4.0, Jh = 1.0, SOC, lambda = 0.50, GQNs (N)

* 5d // Ud = 4.0, Jh = 1.0, GQNs (N,Sz) 

* 5e // Uc = 4.0, Uv = 2.0, Jz = 1.0, Js = 1.0, Jp = 1.0, GQNs (N,Sz,PS)
 
* 5f // Ud = 4.0, Jh = 1.0, SOC, lambda = 0.50, GQNs (N,Jz) 

* 5g // Ud = 4.0, Jh = 1.0, diagonal CF, SOC, lambda = 0.50, GQNs (N) 

* 5h // Ud = 4.0, Jh = 1.0, CF, SOC, GQNs (N), make natural basis outside 

* 5i // Ud = 4.0, Jh = 1.0, non-diagonal CF, GQNs (N,Sz) 

* 5j // Ud = 4.0, Jh = 1.0, non-diagonal CF, GQNs (N,Sz), make natural basis outside 

* 5k // Ud = 4.0, Jh = 1.0, diagonal CF, SOC, lambda = 0.50, GQNs (N), nmini=3, nmaxi=7

### 7-band model

* 7a // Ud = 4.0, Jh = 1.0, GQNs (N,Sz), nmini=0, nmaxi=5

* 7b // Ud = 4.0, Jh = 1.0, SOC, lambda = 0.50, GQNs (N,Jz), nmini=0, nmaxi=5
