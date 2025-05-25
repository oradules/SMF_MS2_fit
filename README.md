# Code used to infer finite state Markov process models for transcription using MS2 and SMF data

## Content

- directory generate_models contains the code to generate objective functions needed for inference

    - To generate objective functions needed for 5 states models run in Matlab inverse_problem_symbolic5onemoviekini_mergep4p5ponkinishort.m
    - The output files resulting from this should be moved to fit_models directory (this is already done)
    - To generatate objective functions needed for 4 states models launch in Matlab inverse_problem_symbolic4onemovieponkini.m
    - The output files should be moved to fit_models directory (this is already done)

- directory fit_models contains the code for parameter model inference
    - To compute MS2 survival functions launch in Matlab Survival_function_short_long.m
    - To fit models launch in Matlab the following scripts
          fit_5expMODEL67_merge_onemoviekini_mergep4p5ponkinishort.m
          fit_5expMODEL45merge_onemoviekini_mergep4p5ponkinishort.m
          fit_5expMODEL23merge_onemoviekini_mergep4p5ponkinishort.m
          fit_5expMODEL17merge_onemoviekini_mergep4p5ponkinishort.m
          fit_5expMODEL14merge_onemoviekini_mergep4p5ponkinishort.m
          fit_5expMODEL1merge_onemoviekini_mergep4p5ponkinishort.m
          fit_4expMODEL91516_merge_onemovieponkini.m
          fit_4expMODEL911_merge_onemovieponkini.m
          fit_4expMODEL12_merge_onemovieponkini.m
  - To read results and generate figures launch in Matlab read_distribution_parsmerge_onemoviekini_mergep4p5ponkinishort.m
 The results are written in the directory results_mergekinimergep4p5ponkiniallshort45 (already done)
## List of models (the models with * were used in the paper, with a different label)
Model1.  5 states, branched, no k8

Model17*(named Model2 in the paper). 5 states, branched

Model14. 5 states, cyclic, all transitions

Model 2*(named Model4 in the paper). 5 states, cyclic, k9 from state 4 to state 2

Model 3*(named Model5 in the paper). 5 states, cyclic, k9 from state 4 to state 1

Model 13*(named Model3 in the paper). 5 states, cyclic, k9 from state 4 to state 3

Model 4*(named Model7 in the paper). 5 states, cyclic, k9 from state 4 to state 2, no k7

Model 5*(named Model8 in the paper). 5 states, cyclic, k9 from state 4 to state 1, no k7

Model 18*(named Model6 in the paper). 5 states, cyclic, k9 from state 4 to state 3, no k7

Model 6*(named Model10 in the paper). 5 states, cyclic, k9 from state 4 to state 2, no k7, no k8

Model 7*(named Model11 in the paper). 5 states, cyclic, k9 from state 4 to state 1, no k7, no k8

Model 8*(named Model9 in the paper). 5 states, cyclic, k9 from state 4 to state 3, no k7, no k8

Model 9. 4 states, cyclic, k9 from state 4 to state 2, no k7, no k6

Model 10. 4 states, cyclic, k9 from state 4 to state 1, no k7, no k6

Model 11*(named Model1 in the paper). 4 states, chain

Model 12. 4 states, cyclic, k9 from state 4 to state 1, k11 from state 4 to state 2, k12 from state 4 to state 1

Model 15. 4 states, cyclic, k9 from state 4 to state 2

Model 16. 4 states, cyclic, k9 from state 4 to state 1

## How to cite this code
Vera Saninova, Flavia Mazzarda, Lasha Dalakishvili, David Depierre, Christina J. I. Moene, Rachel Topno, Jean-Marc Escudier, Marie-Cécile Robert, Olivier Cuvier, Arnaud Krebs, Nacho Molina, Ovidiu Radulescu, Edouard Bertrand. Single-Molecule DNA Footprinting and Transcription Imaging Reveal the Molecular Mechanisms of Promoter Dynamics.    





