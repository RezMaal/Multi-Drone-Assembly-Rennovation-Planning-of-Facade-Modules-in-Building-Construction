# Multi-Drone Assembly Planning of Facade Rennovation and Module Installation in Building Construction

# Introduction

This repository provides a set of functions to optimize and visualize the assembly planning of prefabricated rectangular modules onto exterior building walls using robot drones. This planning framework can expedite not only the approval stage for such projects, but also the assembly stage of the project. The automatic planning framework can also serve as a basis to aid urban planners in renovation projects to calculate possible energy and CO2 reductions by adopting such solutions at scale. The methodology can also be extended to planning for general façade panel installation (i.e., architectural façade panels). Furthermore, the methodology can be applied to plan the installation sequence and work assignment of prefabricated elements to other robotic machinery, such as tower cranes, in building construction projects by adopting the same principles (but) in the horizontal floor planes instead of the vertical façade plane. Finally, the proposed multi-robot collaborative installation can be applied to coordinate the best renovation and rehabilitation plan to remedy façade damages, which impose problems in structural vulnerability, or architectural aesthetics. The latter is advantageous for efficient preservation and conservation of heritage buildings.

The methodology involves three main steps: (i) automated facade boundary detection; (ii) automated rectangular partitioning of facade regions; and (iii) automated collaborative multi-robot assembly planning of the prefabricated modules. The methodology was applied to the point cloud of an existing building acquired from OpenHeritage database. The building was chosen since it resembles common buildings that may require facade installation (e.g., architectural facade attachment/renovation in existing buildings, prefabricated facade installation in new projects, complete removal and energy retrofitting in aging buildings) with a complex boundary geometry. The following figure shows the stages of the processing (1) point cloud planar surface detection; (2) plane boundary detection and straightening/unfolding of the facade planes onto one plane to facilitate planning; (3) dividing the boundary region into rectangles to maximize coverage and minimize number of rectangles; and (4) planning the schedule optimized flight-path and sequencing of multiple robots/drones (in this example 5) to install all modules.

![Picture 1](https://github.com/user-attachments/assets/32950321-b8ba-475d-a419-8504424d6fb7)


# Methodology

In this study, the algorithm involved three stages. Stage 2: Rectangular Modularization and Stage 3: Multi-robot Assembly Planning were formulated as mixed integer linear programming (MILP). More specifically, Stage 2 was formulated as a lexicographic optimization to find the best rectangular sub-division that maximizes surface coverage and minimizes number of rectangular modules.

Stage 3 involved five variations of the multiple asymmetric traveling salesman problem (m-ATSP) as follows:
1. Variant 1: minimizing makespan with flexible sources and no connectivity constraint
2. Variant 2: minimizing makespan with flexible sources and with connectivity constraint
3. Variant 3: minimizing sum of absolute deviation with flexible sources and with connectivity constraint
4. Variant 4: minimizing makespan with connectivity constraint from a single source
5. Variant 5: minimizing sum of absolute deviation with connectivity constraint with a single source
   
Note: The MILP optimization was solved in MATLAB using the Gurobi solver. The MILP formulation of Stage 2 was generally easy to converge. The MILP formulations of Stage 3, particularly Variants 2–5 are NP-hard and require parameter tuning.


# Preliminary Results

The following image shows an example of the results of the optimal multi-robot collaborative assembly planning for Variant 3 above. The top image shows the Time vs. Module installation ID for a sample of five (5) robots/drones. The middle figure shows the sub-route paths taken by each drone/robot to complete the installation (based on the high-level planning)–rotated 90º to save space. The bottom image shows the segmented connectivity graph selected by each of the robots/drones.

![Picture 2](https://github.com/user-attachments/assets/a0b894a4-ef4f-409a-b08a-57ff6b627856)


# Explanation on Code and Provided Files

As part of this Repository, three folders are provided:

1- Input, which is a point cloud of major planar surfaces of the exterior walls of the building, decomposed into 19 files of less than 25MB for upload in GitHub.

2- Dependent Functions, which includes the new functions called in the main script.

3- Demonstration, which includes the main script along with the code for the general plotting.

To successfully run the code, all items must be accessible to Matlab at the time of execution. 

Disclaimer: The backbone of the code was generated with the help of ChatGPT 5.0 (Thinking) and debugged/thoroughly reviewed and tested by the author. While it works, there is no guarantee of efficiency or fast convergence of the algorithms (particularly the MILP solutions in Stage 3).

# Citation
The study is published in the proceedings of the Smart and Sustainable Built Environment (SASBE), where the algorithms and background are explained in more detail. You may cite the study using the following information:

R. Maalek (2025). Multi-Drone Collaborative Planning for Prefabricated Façade Installation on Heritage Buildings. Lecture Notes in Civil Engineering: Proceedings of the Smart and Sustainable Built Environment (SASBE), Lille, France, November 2025.

# Acknowledgements
The author wishes to acknowledge the generous endowment provided by GOLDBECK GmbH to the Karlsruhe Institute of Technology (KIT) for the establishment of the Professorship in Digital Engineering and Construction at the Institute of Technology and Management in Construction (TMB).
