# RUM-DFT Model

This repository presents the RUM-DFT model, which integrates Random Utility Maximization (RUM) and Decision Field Theory (DFT) to capture the dynamic and cognitive nature of decision-making processes. By incorporating breadth-first information search, attribute re-attendance, and evolving utilities over the deliberation time, the RUM-DFT model addresses limitations of DFT (e.g., identification issues, or lack of RUM compatibility).

## Contents

1. DGP.R:
Generates simulated data following the RUM-DFT model. This script allows you to create synthetic datasets for testing and validating the model’s estimation capabilities.

2. RUM_DFT_ISP.R:
Loads simulated data, estimates the RUM-DFT-ISP model (with observed information search sequences), and returns estimates and model performance.

3. functions.R:
Contains the core functions representing the derivations the RUM-DFT model. Each function is referenced with the corresponding equations number in the manuscript.
