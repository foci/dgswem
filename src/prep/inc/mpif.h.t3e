
! "@(#) mpi_t3e/MPI.inc 10.1    07/29/96 15:45:41"
!
!       (C) COPYRIGHT CRAY RESEARCH, INC.
!       UNPUBLISHED PROPRIETARY INFORMATION.
!       ALL RIGHTS RESERVED.
!

! ***** This Fortran include file was installed for CRAY MPT T3E MPI *****

! ======================================================================*/
! (c) Copyright 1995, The University of Edinburgh.
! ======================================================================*/
!                                                                       */
!     File: MPI.inc                                                     */
!  Project: KTG-T3EMPI                                                  */
!   Author: K.Cameron & A.G. Smith                                      */
!  Created: 07/11/1994                                                  */
!  Descrip: The MPI Fortran interface header file                       */
!                                                                       */
! ======================================================================*/

! === Follows MPI (Message-passing Interface) Standard   Annex A. ===*/

! --- Return codes */
      integer*8    MPI_SUCCESS
      parameter (MPI_SUCCESS = 0)
      integer*8    MPI_ERR_BUFFER
      parameter (MPI_ERR_BUFFER = 1)
      integer*8    MPI_ERR_COUNT
      parameter (MPI_ERR_COUNT = 2)
      integer*8    MPI_ERR_TYPE
      parameter (MPI_ERR_TYPE = 3)
      integer*8    MPI_ERR_TAG
      parameter (MPI_ERR_TAG = 4)
      integer*8    MPI_ERR_COMM
      parameter (MPI_ERR_COMM = 5)
      integer*8    MPI_ERR_RANK
      parameter (MPI_ERR_RANK = 6)
      integer*8    MPI_ERR_REQUEST
      parameter (MPI_ERR_REQUEST = 7)
      integer*8    MPI_ERR_ROOT
      parameter (MPI_ERR_ROOT = 8)
      integer*8    MPI_ERR_GROUP
      parameter (MPI_ERR_GROUP = 9)
      integer*8    MPI_ERR_OP
      parameter (MPI_ERR_OP = 10)
      integer*8    MPI_ERR_TOPOLOGY
      parameter (MPI_ERR_TOPOLOGY = 11)
      integer*8    MPI_ERR_DIMS
      parameter (MPI_ERR_DIMS = 12)
      integer*8    MPI_ERR_ARG
      parameter (MPI_ERR_ARG = 13)
      integer*8    MPI_ERR_UNKNOWN
      parameter (MPI_ERR_UNKNOWN = 14)
      integer*8    MPI_ERR_TRUNCATE
      parameter (MPI_ERR_TRUNCATE = 15)
      integer*8    MPI_ERR_OTHER
      parameter (MPI_ERR_OTHER = 16)
      integer*8    MPI_ERR_INTERN
      parameter (MPI_ERR_INTERN = 17)
      integer*8    MPI_ERR_IN_STATUS
      parameter (MPI_ERR_IN_STATUS = 18)
      integer*8    MPI_ERR_PENDING
      parameter (MPI_ERR_PENDING = 19)
      integer*8    MPI_ERR_LASTCODE
      parameter (MPI_ERR_LASTCODE = 20)

! --- Assorted constants */
!     [MPI_BOTTOM must be recognised by its address:
!     it must reside in a common block, value is not used]
      integer*8 MPI_BOTTOM
      common /MPIPRIV/ MPI_BOTTOM
      integer*8    MPI_PROC_NULL
      parameter (MPI_PROC_NULL = -1)
      integer*8    MPI_ANY_SOURCE
      parameter (MPI_ANY_SOURCE = -101)
      integer*8    MPI_ANY_TAG
      parameter (MPI_ANY_TAG = -102)
      integer*8    MPI_UNDEFINED
      parameter (MPI_UNDEFINED = -1)
!     MPI_UB -- as elementary datatype below */
!     MPI_LB -- as elementary datatype below */
      integer*8    MPI_BSEND_OVERHEAD
      parameter (MPI_BSEND_OVERHEAD = 0)
      integer*8    MPI_KEYVAL_INVALID
      parameter (MPI_KEYVAL_INVALID = -1)

! --- Status size and reserved index values (Fortran) */
      integer*8    MPI_STATUS_SIZE
      parameter (MPI_STATUS_SIZE = 6)
      integer*8    MPI_SOURCE
      parameter (MPI_SOURCE = 1)
      integer*8    MPI_TAG
      parameter (MPI_TAG = 2)
      integer*8    MPI_ERROR
      parameter (MPI_ERROR = 3)

! --- Error-handling specifiers */
      integer*8    MPI_ERRORS_ARE_FATAL
      parameter (MPI_ERRORS_ARE_FATAL = 0)
      integer*8    MPI_ERRORS_RETURN
      parameter (MPI_ERRORS_RETURN = 1)

! --- Maximum sizes for strings */
      integer*8    MPI_MAX_PROCESSOR_NAME
      parameter (MPI_MAX_PROCESSOR_NAME = 256)
      integer*8    MPI_MAX_ERROR_STRING
      parameter (MPI_MAX_ERROR_STRING = 256)

! --- Elementary and optional datatypes (Fortran) */
      integer*8    MPI_integer
      parameter (MPI_integer = 14)
      integer*8    MPI_REAL
      parameter (MPI_REAL = 15)
      integer*8    MPI_DOUBLE_PRECISION
      parameter (MPI_DOUBLE_PRECISION = 16)
      integer*8    MPI_COMPLEX
      parameter (MPI_COMPLEX = 17)
      integer*8    MPI_DOUBLE_COMPLEX
      parameter (MPI_DOUBLE_COMPLEX = 18)
      integer*8    MPI_LOGICAL
      parameter (MPI_LOGICAL = 19)
      integer*8    MPI_CHARACTER
      parameter (MPI_CHARACTER = 20)
      integer*8    MPI_BYTE
      parameter (MPI_BYTE = 21)
      integer*8    MPI_PACKED
      parameter (MPI_PACKED = 22)
      integer*8    MPI_UB
      parameter (MPI_UB = 29)
      integer*8    MPI_LB
      parameter (MPI_LB = 30)

      integer*8    MPI_integer1
      parameter (MPI_integer1 = 23)
      integer*8    MPI_integer2
      parameter (MPI_integer2 = 24)
      integer*8    MPI_integer4
      parameter (MPI_integer4 = 25)
      integer*8    MPI_REAL4
      parameter (MPI_REAL4 = 27)
      integer*8    MPI_REAL8
      parameter (MPI_REAL8 = 28)

! --- Datatypes for reduction functions (Fortran) */

      integer*8    MPI_2REAL
      parameter (MPI_2REAL = 37)
      integer*8    MPI_2DOUBLE_PRECISION
      parameter (MPI_2DOUBLE_PRECISION = 38)
      integer*8    MPI_2integer
      parameter (MPI_2integer = 39)
!       MPI_2COMPLEX not defined [MPI Errata Oct 94]

! --- Reserved communicators */
      integer*8    MPI_COMM_WORLD
      parameter (MPI_COMM_WORLD = 0)
      integer*8    MPI_COMM_SELF
      parameter (MPI_COMM_SELF = 1)

! --- Results of communicator and group comparisons */
      integer*8    MPI_UNEQUAL
      parameter (MPI_UNEQUAL = -1)
      integer*8    MPI_SIMILAR
      parameter (MPI_SIMILAR = -2)
      integer*8    MPI_IDENT
      parameter (MPI_IDENT = -3)
      integer*8    MPI_CONGRUENT
      parameter (MPI_CONGRUENT = -4)

! --- Environmental inquiry keys */
      integer*8    MPI_TAG_UB
      parameter (MPI_TAG_UB = -50)
      integer*8    MPI_IO
      parameter (MPI_IO = -49)
      integer*8    MPI_HOST
      parameter (MPI_HOST = -48)
      integer*8    MPI_WTIME_IS_GLOBAL
      parameter (MPI_WTIME_IS_GLOBAL = -47)

! --- Collective operations */
      integer*8    MPI_MAX
      parameter (MPI_MAX = 0)
      integer*8    MPI_MIN
      parameter (MPI_MIN = 1)
      integer*8    MPI_SUM
      parameter (MPI_SUM = 2)
      integer*8    MPI_PROD
      parameter (MPI_PROD = 3)
      integer*8    MPI_MAXLOC
      parameter (MPI_MAXLOC = 4)
      integer*8    MPI_MINLOC
      parameter (MPI_MINLOC = 5)
      integer*8    MPI_BAND
      parameter (MPI_BAND = 6)
      integer*8    MPI_BOR
      parameter (MPI_BOR = 7)
      integer*8    MPI_BXOR
      parameter (MPI_BXOR = 8)
      integer*8    MPI_LAND
      parameter (MPI_LAND = 9)
      integer*8    MPI_LOR
      parameter (MPI_LOR = 10)
      integer*8    MPI_LXOR
      parameter (MPI_LXOR = 11)

! --- Null handles */
      integer*8    MPI_GROUP_NULL
      parameter (MPI_GROUP_NULL = -1)
      integer*8    MPI_COMM_NULL
      parameter (MPI_COMM_NULL = -1)
      integer*8    MPI_DATATYPE_NULL
      parameter (MPI_DATATYPE_NULL = -1)
      integer*8    MPI_REQUEST_NULL
      parameter (MPI_REQUEST_NULL = 0)
      integer*8    MPI_OP_NULL
      parameter (MPI_OP_NULL = -1)
      integer*8    MPI_ERRHANDLER_NULL
      parameter (MPI_ERRHANDLER_NULL = -1)

      integer*8    MPI_NULL_COPY_FN
      parameter (MPI_NULL_COPY_FN = -1)
      integer*8    MPI_DUP_FN
      parameter (MPI_DUP_FN = -2)
      integer*8    MPI_NULL_DELETE_FN
      parameter (MPI_NULL_DELETE_FN = -1)

! --- Empty group */
      integer*8    MPI_GROUP_EMPTY
      parameter (MPI_GROUP_EMPTY = 0)

! --- Topologies */
      integer*8    MPI_GRAPH
      parameter (MPI_GRAPH = 100)
      integer*8    MPI_CART
      parameter (MPI_CART = 101)

! --- Timing Functions */
! --- (should be DOUBLE PRECISION -- not on Cray T3E)
      external MPI_WTICK
      real*8   MPI_WTICK
      external MPI_WTIME
      real*8   MPI_WTIME

! ======================================================================*/
!  End of File: MPI.inc                                                 */
! ======================================================================*/
