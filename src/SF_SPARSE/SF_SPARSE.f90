MODULE SF_SPARSE
!:synopsis: Local real and complex sparse matrices with dynamic row or column storage
!
!Local sparse matrices with dynamic row or column storage, in double real
!(``D``) and double complex (``Z``) arithmetic. 
!
!The implementations are contained in explicit Fortran modules: ``SF_SPARSE_CSR_D``,
!``SF_SPARSE_CSR_Z``, ``SF_SPARSE_CSC_D``, and ``SF_SPARSE_CSC_Z``.
!
!Storage and insertion contract
!------------------------------
!
!These objects contain an array of allocatable sparse rows (CSR) or columns
!(CSC), rather than a single packed CSR/CSC buffer.
!
!.. warning::
!
!   **Insertion MUST preserve sorted indices.** For CSR, append new column
!   indices in strictly increasing order within each row. For CSC, append new
!   row indices in strictly increasing order within each column. Neither
!   ``insert`` nor ``fast_insert`` sorts indices or reorders values.
!
!``insert(value,i,j)`` adds to an existing coordinate, or appends a new one.
!An existing coordinate can be revisited without disturbing the ordering.
!
!``fast_insert(value,i,j)`` always appends: the caller must guarantee both
!uniqueness and increasing order. Rows/columns can be interleaved during
!construction; the ordering requirement applies independently within each.
!The same contract applies when directly filling the public storage components.
!Violating it can produce incorrect merge and multiplication results.
!
!``nnz()`` counts stored entries, including exact zeros. Algebra does not apply
!an implicit tolerance. Matrix multiplication omits exact-zero results to avoid performance losses.
!
!Dimensions and indices use default Fortran integers.
!
!Construction and access
!-----------------------
!
!The four public types retain their names:
!``sparse_dmatrix_csr``, ``sparse_zmatrix_csr``,
!``sparse_dmatrix_csc``, and ``sparse_zmatrix_csc``.
!
!.. code-block:: fortran
!
!   use SF_SPARSE
!   type(sparse_dmatrix_csr) :: A, B, C
!   real(8) :: dense(3,5), x(5)
!   real(8), allocatable :: y(:)
!
!   call A%init(3,5)
!   call A%insert(2d0,1,1)
!   call A%insert(3d0,1,4) ! increasing columns within row 1
!   call A%insert(1d0,2,2) ! row 2 has its own insertion order
!   x = 1d0
!   y = matmul(A,x)       ! length 3
!   call A%dump(dense)    ! replaces the destination contents
!   B = as_csr(dense)
!   C = transpose(A)     ! shape [5,3], CSR storage
!
!``init(n)`` constructs an empty square matrix; ``init(m,n)`` constructs an
!empty rectangular matrix. 
!  
!``free()`` is safe on fresh and already freed
!objects. 
!
!``shape(A)`` and ``A%shape()`` return the physical row/column counts.
!
!``get(i,j)`` returns zero for an absent coordinate.
!
!``as_csr(dense)`` and ``as_csc(dense)`` select storage explicitly and overload
!on real/complex input. Type-name constructors also accept dense input.
!
!``as_sparse`` and ``sparse`` remain available when importing a single
!orientation; they are intentionally absent from the combined facade because
!Fortran cannot choose an orientation from the result type alone.
!
!``as_matrix(A)`` / ``A%as_matrix()`` returns a dense array; ``dump`` overwrites
!an existing dense array of matching shape. ``copy`` and intrinsic matrix
!assignment perform deep copies. Assigning a scalar changes only stored values,
!not the sparsity pattern. The inherited ``Nelements`` field caches the stored-entry count and is maintained
!by the library. ``nnz()`` computes the count from the concrete storage. If filling
!storage arrays directly, keep their sizes and ``Nelements`` consistent.
!
!
!
!Common parent and polymorphism
!------------------------------
!
!All four types extend the concrete ``sparse_matrix`` parent from
!``SF_SPARSE_COMMON``; it is also exported by ``SF_SPARSE`` and ``SCIFOR``.
!Dimensions, ``Nelements``, and initialization status live in that parent.
!
!.. code-block:: fortran
!
!   use SF_SPARSE
!   class(sparse_matrix), allocatable :: A
!
!   allocate(sparse_zmatrix_csc :: A)
!   call A%init(3,5)        ! dispatches to CSC_Z initialization
!   print *, A%shape(), A%nnz(), A%bytes()
!   select type(A)
!   type is(sparse_zmatrix_csc)
!      call A%insert((2d0,1d0),1,2)
!   end select
!   call A%free()           ! frees the concrete CSC storage
!
!The common interface provides ``init``, ``free``, ``shape``, ``nnz``, ``bytes``,
!and (in MPI builds) ``bcast``. The concrete implementations override operations
!that depend on storage. Numerical insertion, extraction, and algebra remain
!on the concrete types; use ``select type`` when starting from a parent
!polymorphic declaration. Generic numerical operators do not dispatch from the
!parent's declared type.
!
!A plain ``type(sparse_matrix)`` object contains metadata only. Its initialization
!creates no numerical entries. An allocatable ``class(sparse_matrix)`` object can
!hold any concrete child. A polymorphic array has one dynamic type for the whole
!array; for heterogeneous collections, use a wrapper with an allocatable
!``class(sparse_matrix)`` component per element.
!
!For ``bcast`` through the parent interface, every rank must already have an
!object of the same concrete dynamic type. The method transfers its contents;
!it does not allocate or change the receiver's dynamic type.
!
!Algebra
!-------
!
!Supported operations for each type are addition, subtraction, scalar
!multiplication/division, transpose, Hermitian conjugate, matrix multiplication,
!both left and right vector multiplication, and Kronecker products.
!
!* ``matmul(A,B)`` and ``A .m. B`` multiply sparse matrices of the same type and
!  return the same storage orientation. The implementation intersects sorted
!  index streams against a transposed operand.
!* ``matmul(A,v)`` and ``A%dot(v)`` return a vector of length ``A%Nrow``.
!* ``matmul(v,A)`` returns a vector of length ``A%Ncol``. 
!  
!  .. warning::
!
!   **This is an ordinary product: it does not conjugate the complex vector.**
!
!* ``transpose(A)``, ``A%transpose()`` and ``A%t()`` transpose without conjugation.
!* ``hconjg(A)`` and ``A%dgr()`` return the Hermitian conjugate; for real matrices
!  these are the ordinary transpose. Complex types also provide ``A%conjg()``.
!* ``kron(A,B)`` and ``A .x. B`` use the convention
!  ``C((i-1)*B%Nrow+k,(j-1)*B%Ncol+l)=A(i,j)*B(k,l)``.
!* ``sp_eye_csr_d(n)``, ``sp_eye_csr_z(n)``, ``sp_eye_csc_d(n)``, and
!  ``sp_eye_csc_z(n)`` construct sparse identities directly.
!
!Mixed CSR/CSC and mixed real/complex sparse matrix products are not provided.
!Filtering, restricted Kronecker products, and parallel multiplication are
!outside this implementation.
!
!I/O and diagnostics
!-------------------
!
!``write(file=...)`` / ``read(file=...)`` serialize the dimensions followed by
!one size and its index/value entries per stored row or column. Files preserve
!the physical row/column counts but have orientation-specific payloads: read
!with the same orientation and numerical type. Exactly one of ``file`` or
!``unit`` must be supplied. Caller-owned units remain open. Reads require valid,
!strictly increasing indices.
!
!``show(fmt=...,file=...)`` prints the dense matrix; ``display()`` prints stored
!coordinates and values. ``spy(header)`` writes coordinate data and a Gnuplot
!command file; it does not execute an external plotting program.
!``bytes()`` returns the local stored payload in decimal kilobytes, excluding
!allocatable descriptors, allocation overhead, and padding. It is available
!with or without MPI.
!
!MPI helpers
!-----------
!
!With ``USE_MPI=ON``, each sparse matrix remains an ordinary local object.
!No distributed ownership of its rows or columns is implied.
!
!``A%bcast(comm,root)`` collectively replaces the object on receivers with the
!root's object, including dimensions and initialization status. The defaults
!are ``MPI_COMM_WORLD`` and root zero; a supplied root is a rank in the supplied
!communicator. An initialized empty matrix and an uninitialized object remain
!distinct. Only members of that communicator participate.
!
!``AllGather_MPI(comm,array)`` assembles rank-one or rank-two arrays of whole
!matrix objects. Array extents must agree on all participating ranks. Each
!entry must be initialized on at most one rank; an entry with no owner remains
!uninitialized. Inconsistent extents or multiple owners cause a collective
!``MPI_Abort``. After gathering, all ranks hold replicas: free replicas before
!using initialization status to express ownership in a subsequent gather.
!
!Validation
!----------
!
!Numerical regression tests have been carried out for all four types against
!small dense reference calculations, including rectangular and empty matrices,
!complex arithmetic, the polymorphic parent interface, and stored-entry counts.
!MPI tests cover split communicators, nonzero roots, whole-object array gathering,
!and rejection of duplicate owners or inconsistent array extents. Both MPI and
!serial library builds have been checked.
!
  USE SF_SPARSE_COMMON, only: sparse_matrix
  USE SF_SPARSE_CSR_D, only: &
    sparse_dmatrix_csr, &
    as_csr            , &
    sp_eye_csr_d      , &
    assignment(=)     , &
    operator(+)       , &
    operator(-)       , &
    operator(*)       , &
    operator(/)       , &
    operator(.m.)     , &
    operator(.x.)     , &
    kron              , &
    shape             , &
    transpose         , &
    hconjg            , &
    matmul            , &
    as_matrix
#ifdef _MPI
  USE SF_SPARSE_CSR_D, only: AllGather_MPI
#endif
  USE SF_SPARSE_CSR_Z, only: &
    sparse_zmatrix_csr, &
    as_csr            , &
    sp_eye_csr_z      , &
    assignment(=)     , &
    operator(+)       , &
    operator(-)       , &
    operator(*)       , &
    operator(/)       , &
    operator(.m.)     , &
    operator(.x.)     , &
    kron              , &
    shape             , &
    transpose         , &
    hconjg            , &
    matmul            , &
    as_matrix
#ifdef _MPI
  USE SF_SPARSE_CSR_Z, only: AllGather_MPI
#endif
  USE SF_SPARSE_CSC_D, only: &
    sparse_dmatrix_csc, &
    as_csc            , &
    sp_eye_csc_d      , &
    assignment(=)     , &
    operator(+)       , &
    operator(-)       , &
    operator(*)       , &
    operator(/)       , &
    operator(.m.)     , &
    operator(.x.)     , &
    kron              , &
    shape             , &
    transpose         , &
    hconjg            , &
    matmul            , &
    as_matrix
#ifdef _MPI
  USE SF_SPARSE_CSC_D, only: AllGather_MPI
#endif
  USE SF_SPARSE_CSC_Z, only: &
    sparse_zmatrix_csc, &
    as_csc            , &
    sp_eye_csc_z      , &
    assignment(=)     , &
    operator(+)       , &
    operator(-)       , &
    operator(*)       , &
    operator(/)       , &
    operator(.m.)     , &
    operator(.x.)     , &
    kron              , &
    shape             , &
    transpose         , &
    hconjg            , &
    matmul            , &
    as_matrix
#ifdef _MPI
  USE SF_SPARSE_CSC_Z, only: AllGather_MPI
#endif
  implicit none
END MODULE SF_SPARSE
