<h1 align="center">
  <b>Supervised brain node and network construction under <br /> voxel-level functional imaging</b><br>
</h1>

This repo is the official R implementation for the "Supervised brain node and network construction under voxel-level functional imaging" submission. 

**Abstract**: 
> Recent advancements in understanding the brain's functional organization related to behavior have been pivotal, particularly in the development of predictive models based on brain connectivity. Traditional methods in this domain often involve a two-step process by first constructing a connectivity matrix from predefined brain regions, and then linking these connections to behaviors or clinical outcomes. However, these approaches with unsupervised node partitions predict outcomes inefficiently with independently established connectivity.  In this paper, we introduce the Supervised Brain Parcellation (SBP), a brain node parcellation scheme informed by the downstream predictive task. With voxel-level functional time courses generated under resting-state or cognitive tasks as input, our approach clusters voxels into nodes in a manner that maximizes the correlation between inter-node connections and the behavioral outcome, while also accommodating intra-node homogeneity.  We rigorously evaluate the SBP approach using resting-state and task-based fMRI data from both the Adolescent Brain Cognitive Development (ABCD) study and the Human Connectome Project (HCP). Our analyses show that SBP significantly improves out-of-sample connectome-based predictive performance compared to conventional step-wise methods under various brain atlases. This advancement holds promise for enhancing our understanding of brain functional architectures with behavior and establishing more informative network neuromarkers for clinical applications.


### Usage
There are two scripts included, with th required libraries listed on the top of each file. 

"SBP.R" contains the main algorithm described in the paper. 

"Simulation.R" can reproduce Figure 3. By running this script, generated figures will be automaticated saved in the same folder. 

