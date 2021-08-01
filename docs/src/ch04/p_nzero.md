### Parameter: nzero

**Definition**

Maximum allowed number of non-zero elements in $$F$$-matrix. In the iQIST software package, the terminology $$F$$-matrix denotes the $$d$$ operator (the destroy operator for impurity electron) matrix in the eigen-basis of the local Hamiltonian or in occupation number basis.

**Type**

Integer

**Default value**

128

**Component**

Only for the **BEGONIA**, **LAVENDER**, and **CAMELLIA** components.

**Behavior**

In the **BEGONIA**, **LAVENDER**, and **CAMELLIA** components, we used the sparse matrix technology to accelerate the matrix times vector operation (```spmv```) and matrix times matrix (```spmm```). In order to save the memory, the $$F$$-matrix, the local Hamiltonian matrix, and the product matrix of the above operations are stored in CSR format. Here, we used *nzero* parameter to describe the maximum allowed number of non-zero elements of these matrices, and allocate the corresponding memories.

In the ```spmm``` operations, sometimes the number of the non-zero elements in the resulted matrix exceeds the value of *nzero*, which will trigger a fatal error/exception. At that time, you have to increase *nzero* by twofold, and redo the calculation.

**Comment**

In the **BEGONIA** and **LAVENDER** components the $$F$$-matrix mean the $$d$$ operator matrix in the eigen-basis. However, in the **CAMELLIA** component, it means the $$d$$ operator matrix in the occupation number basis.
