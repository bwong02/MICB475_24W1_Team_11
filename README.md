# MICB475_24W1_Team_11

Team 11: Brandon Wong, Gabi Reznik, McKenna Postles, Miles Silva, Nicolas Taraviras

**Hi Team 11! Here are the zoom link and passcode we will use for our next meeting (October 10th). Thank you and see you soon!**

Here is the link you can use: https://ubc.zoom.us/j/67359231965?pwd=1QOWWuFTTjXQcla2T32GQ5MoyAUHhG.1

Meeting ID: 673 5923 1965

Passcode: 574502


# Meeting Minutes and Agenda

## October 3, 2024

### Agenda

- Go through proposed papers and databases to decide on one to move forward with
- Brainstorm research question depending on final database chosen

### Meeting Minutes

2 Paper Choices: Gastrointestinal Cancer or Breast Cancer
- Breast Cancer does not have much metadata on the cancer patients
- **CHOSE GASTROINTESTINAL CANCER DATA**

Research Questions Brainstorm
- 2 question routes:
  - Taxonomic route: what species or communities exist or is associated with something
  - Functional route: predict what metabolic pathways exist in those communities
- How does H. pylori affect the microbiome community?
- Look at H. pylori both taxonomically and functionally?
- Is H. pylori associated with different histopathology stages?
- BEFORE MAKING A RESEARCH QUESTION, MAKE A CORRELATION PLOT OR CREATE A NEW METADATA CATEGORY (go with group and H. pylori positive/negative like HC+, HC-, CG+, CG-, etc.)

Next Steps
1. Create a column with the Group + H. pylori status - Nicholas
2. Make heat map comparing the different categories (UJEMI+ anemia group did a heat map) - Miles
4. Outline proposal before next meeting (aims and research question) - McKenna, Brandon, Gabi
5. Find background papers that discuss H. pylori, cancer, and microbiome
6. Eventually check with Hans if trimming or denoising is appropriate 

Other Notes:
 - If we find an non-annotated ASV, we can BLAST it
 - Finding an ASV does not mean it is biologically significant in the community




## October 10, 2024

### Agenda
**1. Discuss research question and specific project aims.**
- Potential questions:
  - Is there a difference in microbiome diversity and function comparing  +/- H. pylori status within and between disease stages of gastric cancer?
- Our thoughts/questions 
  - Is this okay for us to move forward with H. pylori analysis, even though they also examined
      
**2. Discuss project rational (discuss heat map, etc...)**
- Rationale: created a heat map looking at (+/-) HP samples in each stage:
- We noticed little difference between collection sites, which led us to look into differences between groups and not distinguish 
      between collection sites. 
- We noticed a difference in relative percentages between groups
- We discussed looking into differences between subtypes within each stage but decided against it due to sample constraints. 
- The heatmap is located in the R folder already. The R project and script file are also there. 

**3. Discuss project aims** 
Aim 1: Is there an H. Pylori status difference between different sites within a disease stage (supported by heat map)?
Aim 2: Does H. pylori status and disease stage affect beta diversity and microbiome taxonomy?
Aim 3: What functional differences exist between the communities at each disease stage and between disease stages? 

**4. Go through our qiime2 workflow** 
- Trimming, rarefaction, diversity plots, etc... 

### Meeting Minutes
**Hans Points from Assignment 5**
- To make things reproducible, always include all input files, output files, and scripts in the middle
  - No lines of code should be needed to change
- Avoid absolute paths and include relative paths from the folder you're sharing with the person
- Always comment preceding the code --> describe what code does
  - Very helpful!

 **Deadline for Proposal is October 20**

**For talking point 1:**
 - Paper looked within cohorts and only looked at the abundance of H. pylori, grouping it into low (below 0.01) and high abundance
   - 0.01 is already pretty undetectable. We can change this or create our own categories to i.e. Low can be 0.01 - 1% and high above 1%
 - Paper already did basic analysis of H. pylori and did basic microbiome analysis, BUT the paper did not go past the PCoA analysis
   
 - How about we make **new metadata**. New category we look at: Fusobacteria
   - Explore how the presence or absence of this microbe changes the community composition
   - Potentially, redo what they did and look at H. pylori and Fusobacteria (**OR only Fusobacteria**)
     - Assign undetectable, low, and high for abundance
       - **Undetectable** should be half of the error rate of the machine they use (find error rate in the seq they use)
       - I.e. Undetectable H. pylori against low fusobacteria
     - Though we may not have enough sampels for some of the categories, it would be interesting to look at these new categories
   - Should we also look at Functional too? How would Hans approach it?
     - We accept negative data!
     - We can look at both taxonomic and functional then change this 1 month later (proposal does not need to match final outcome)
   - Should we look at the low abundant species as there was unweighted significance
     
  - **FINAL DECISION (include it in the proposal)**:
    - Categorize the abundance of H. pylori to undetectable, low, and high (new categorization)
    - Categorize the abundance of Fusobacteria to undetectable, low, and high (novel analysis)
    - Look at the taxonomy across the stages and progression to gastric cancer (GC)
    - Also do a functional analysis of the bacteria at the stages
  
  - Action items:
    - Find another paper that categorizes undetectable, low, and high
    - Have aims, objectives, and approach set before the next meeting (points 2-5)


## October 17, 2024

### Agenda

Go through **Project Proposal**: https://docs.google.com/document/d/1onvvpNUt2O39NBdv-u6NiFVkD_dlgdf6OKItRVRdAC0/edit?usp=sharing
- Confirm the abdunance levels and convention with Hans
- Confirm Aims with Hans
  - Do we combine Aims 2 and 3?
  - What happens if Aim 2 goes wrong? Can it and how do you justify this? 
- Confirm Approaches with Hans
- Confirm the Data Set Overview (which ones to do)

### Meeting Minutes

Notes
- Keeping the nucleotide reads longer will give better resolution
- **Abundance Levels Question**
  - 0.01 is a commonly accepted number for relative abundance and the notation we used is good with the %
  - We don't need a specific paper that used our abundance ranges
  - Say we defined 1 as low abundance then justify why we did this
  - For the study, look at how many species are present above 1% abundance and how many are below
- **Data Set Overview Question**
  - Use the original data set overview in the canvas
- **Research Objective Questions**
  - It's not upregulation, it's an increase in representation/presence
  - Pie Crest: Takes input of sequences --> sees what genes are present --> sees what pathways COULD be expressed or used (functional POTENTIAL and PREDICTION)
    - With 16S, we are looking at prediction of potential
    - Future directions can say we do transcriptomics
- **Aims Questions**
  - Paper Overall: You should be redundant
  - Most comments are on the google doc
  - For Aim 3, if we care more about disease stage, do Indicator Species Analysis (ISA) on progression
    - NOT FOR PROPOSAL BUT GOOD TO KNOW
    - Read a paper by morton et al. 2019 - reference frames
    - Knight Lab
      - Interesting way to look at changes in abundance
      - Log 10 ratio of two taxa with 1 being constant and then the other changing across samples???
  - What happens if Aim 2 goes wrong?
    - Hans says you can have a paper of negative results --> if this happens, focus on indicator taxa part and see if there is any correlation with function etc.
    - Just say there's nothing interested we found about H. pylori and Fusobacterium
- **Proposed Timeline**
  - Include learning QIIME (everything since September)

## October 21, 2024

### Agenda
- Go over proposal

### Meeting Minutes
Notes
- We have the chance to resubmit the proposal
  - Automatic 5% more if we resubmit (probably)
- **Question**: Confused about Functional.
  - Our thoughts: Do whole thing --> reclassification of H. pylori & Fusobactera --> 2 new columns (bacteria + disease stage) --> run PiCRUSt2 based on these groups
  - Hans:
    - PiCRUSt2 output: Gives pathway representation (how much is present and how much is not) --> raw count of pathway present in one sample
    - Calculate the log fold change of representation of pathways based on the healthy control
    - DESEQ2 does normalization
    - In the 3 categories in each disease stage, meaning 15 groups, we can make a stacked bar plot --> each panel has 1 pathway (do top 10 overrepresented and top 10 underrepresented, each one is going to have 5 stacked bar plots)
- **Question**: Do we run Indicator Analysis on each disease stage then find correlation with abundance OR do we incorporate abundance into the analysis itself
  - Hans
    - We can't include abundance in the analysis
    - Should simplify analysis
    - We were right in our proposal
      
- **GOAL for this Week: Complete Aim 1**

## October 31, 2024

### Agenda


### Meeting Minutes


## November 7, 2024

### Agenda


### Meeting Minutes


## November 14, 2024

### Agenda


### Meeting Minutes
