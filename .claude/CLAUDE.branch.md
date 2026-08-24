# What this branch is

DaGSWEM^2 (Discontinous Adaptive Galerkin Shallow Water Equation Model/Directed Acyclic Graph Shallow Water Execution Mode) is an experimental branch of the main DG-SWEM repository that implements an execution mode based on StarPU to allow for task-based parallelism. This branch is intended to explore the potential benefits of a task-based execution model for shallow water simulations, which may include the following: 

- Dynamic load balancing for dynamically heterogeneous problems

## v0.1.0 (Current)
This version will focus on implementing a basic task-based execution model for the shallow water equations, such that dynamic load balancing can be acheived on a multi-core CPU. 

### Orchestrator Architecture
This most primitive architecture will be almost 1-1 to the traditional flat MPI architecture of subdomains with ghost and resident data. A task will correspond roughly to a MPI rank. 

#### Data
-  State vectors: StarPU vectors containing state data for each element/node, ghost and resident. Filtered using a single list-directed filter corresponding to the subdomain. Examples:
    - `ZE`: Water height
    - `QX`: Momentum in the x-direction
    - `QY`: Momentum in the y-direction
-  Send vectors: StarPU vectors containing data to be sent to neighboring subdomains. Filtered using a 2-tier list-directed filter corresponding to (subdomain, neighbor). This differs with MPI in that no recieving buffer will be declared, instead relies on StarPU to dynamically allocate the data. Examples: 
    - `ZE_send`: Water height to be sent to neighbors
    - `QX_send`: Momentum in the x-direction to be sent to neighbors
    - `QY_send`: Momentum in the y-direction to be sent to neighbors

#### Codelets
- `update`: Reads send vectors, unpacks into state vectors, updates state vectors, and packs back into send vectors. 

#### Task Graph
The codelets will be submitted in a loop and dependencies will be inferred from the write/read patterns of the send vectors. 

