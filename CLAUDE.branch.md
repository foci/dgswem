# What this branch is

DaGSWEM^2 (Discontinous Adaptive Galerkin Shallow Water Equation Model/Directed Acyclic Graph Shallow Water Execution Mode) is an experimental branch of the main DG-SWEM repository that implements an execution mode based on StarPU to allow for task-based parallelism. This branch is intended to explore the potential benefits of a task-based execution model for shallow water simulations, which may include the following: 

- Dynamic load balancing for dynamically heterogeneous problems
- Hardware-aware load balancing/assignment to heterogeneous architectures
- Improved compute/communication overlap
- Better cache reuse with weak loop tiling
- Speculative execution to relax stability constraints/enforcement