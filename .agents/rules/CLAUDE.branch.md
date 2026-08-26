---
trigger: always_on
---

# What this branch is

DaGSWEM^2 (Discontinous Adaptive Galerkin Shallow Water Equation Model/Directed Acyclic Graph Shallow Water Execution Mode) is an experimental branch of the main DG-SWEM repository that implements an execution mode based on StarPU to allow for task-based parallelism. This branch is intended to explore the potential benefits of a task-based execution model for shallow water simulations, which may include the following: 

- Dynamic load balancing for dynamically heterogeneous problems

## v0.1.0 (Current)
This version will focus on implementing a basic task-based execution model for the shallow water equations, such that dynamic load balancing can be acheived on a multi-core CPU. 

### Orchestrator Architecture
This most primitive architecture will be almost 1-1 to the traditional flat MPI architecture of subdomains with ghost and resident data. A task will correspond roughly to doing some work on an MPI rank. 

#### Data
-  State variables: For each subdomain, each one will get a StarPU handle respresenting that subdomain's state variables, including ghost and resident elements/nodes. 
-  Recv vectors: For each subdomain's neighbors, a StarPU vector will be declared that is intended to be read by the subdomain and written to by the neighbor. 

#### Codelets
- `ADVANCE_STAGE_DISTRIBUTED`: Reads recv vectors, unpacks into state vectors, updates state vectors, and packs back into send vectors. 

#### Task Graph
The codelets will be submitted in a loop and dependencies will be inferred from the write/read patterns of the send vectors. 

#### What has been done so far
- To reuse the `DG.F` header, state variables that are not temporary to an advance stage task have been changed from ALLOCATABLE, to POINTER, so that each StarPU worker will have its own copy. The `ADVANCE_STATE_DISTRIBUTED` codelet contains a header subroutine that redirects the workers' pointers to the one handled by StarPU, so that each task effectively has its own private copies. Non state variables are left untouched. 
- `dagswem_cl.f90` has been filled out partially with a codelet stub and a `DG_STATE_ACTIVATE` subroutine that points the threadprivate declarations to the corresponding StarPU handles. 