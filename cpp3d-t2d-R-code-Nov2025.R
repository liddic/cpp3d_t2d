#########################
#
# cpp3d: Compound processing potential, at the scale of individual compounds, with 3-d chemical/bioenergetic mapping visualisation (O:C, H:C, N:C)
# Craig Liddicoat | Flinders University, South Australia 

# cpp3d-t2d: Compound processing potential: translates functional potential profile 
# relative abundance data to potential metabolism at the scale of individual compounds, 
# with 3-d chemical/bioenergetic mapping-visualisation using O:C, H:C, N:C ratios.
# This study compares case study metagenome datasets of soil microbiomes from ecosystem 
# restoration and disturbed vs natural soils; and gut microbiomes in type 2 diabetes (T2D)
# versus normal health.
# Code by Craig Liddicoat, Flinders University, South Australia
#
#########################

# record library and version info
.libPaths() # "/Library/Frameworks/R.framework/Versions/4.2/Resources/library"

R.Version()
# "R version 4.2.2 (2022-10-31)"
citation()
# R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL
# https://www.R-project.org/.

library(readxl); packageVersion("readxl") # '1.4.1'
library(plyr); packageVersion("plyr") # '1.8.8'
library(dplyr); packageVersion("dplyr") # '1.0.10'
library(vegan);packageVersion("vegan") # '2.6.4'
library(phyloseq); packageVersion("phyloseq") # '1.42.0'
library(ggplot2); packageVersion("ggplot2") # '3.4.0'
library(grid); packageVersion("grid") #  '4.2.2'
library(reshape2); packageVersion("reshape2") # '1.4.4'
library(tidyr); packageVersion("tidyr") # '1.2.1'
library(corrr); packageVersion("corrr") # '0.4.4'
library(ggforce); packageVersion("ggforce") # '0.4.1'
library(ggrepel); packageVersion("ggrepel") # '0.9.2'
library(stringdist); packageVersion("stringdist") # ‘0.9.10’
library(stringr); packageVersion("stringr") # ‘1.5.0’
library(doParallel); packageVersion("doParallel") # '1.0.17'
library(RColorBrewer); packageVersion("RColorBrewer") # '1.1.3'
library(ggpp); packageVersion("ggpp") # ‘0.5.0’ # https://cran.r-project.org/web/packages/ggpp/vignettes/grammar-extensions.html
library(MASS)                     ;packageVersion("MASS") # ‘7.3.58.1’
library(ggsignif); packageVersion("ggsignif") # '0.6.4'
library(moments)                  ;packageVersion("moments") # ‘0.14.1’
library(grDevices); packageVersion("grDevices") #  '4.2.2'
library(ggbiplot); packageVersion("ggbiplot") #  ‘0.55’
library(viridis); packageVersion("viridis") #  ‘0.6.2’
library(FSA); packageVersion("FSA") # '0.9.3'
library(rcompanion); packageVersion("rcompanion") # '2.4.18'
library(fields); packageVersion("fields") # ‘14.1’
library(car); packageVersion("car") # ‘3.1.1’
library(multcompView); packageVersion("multcompView") # ‘0.1.8’
library(gtools); packageVersion("gtools") # ‘3.9.4’
library(igraph); packageVersion("igraph") #  '1.4.2'
library(pheatmap); packageVersion("pheatmap") # '1.0.12'
library(colorspace); packageVersion("colorspace") # ‘2.1.0’

#########################
## save.image("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/cpp3d-indiv-resto-vs-t2d-WORKSPACE-v6.RData")
## load("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/cpp3d-indiv-resto-vs-t2d-WORKSPACE-v6.RData")
#########################

workdir <- "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R"
setwd(workdir)
getwd()

par.default <- par()


#### ModelSEED lookup tables
#-------------------------

seed_db_dir <- "/Users/lidd0026/WORKSPACE/DATA/ModelSEEDDB/select_files"

## Unique_ModelSEED_Reaction_ECs

rxn_ECs.lut <- read_excel(path = paste0(seed_db_dir,"/","Unique_ModelSEED_Reaction_ECs.xlsx") , range = "A1:C30834")
rxn_ECs.lut <- as.data.frame(rxn_ECs.lut)
head(rxn_ECs.lut)
names(rxn_ECs.lut) # "ModelSEED ID" "External ID"  "Source"
names(rxn_ECs.lut) <- c("rxn_id", "EC_id",  "Source")

dim(rxn_ECs.lut) # 30833     3
length(unique(rxn_ECs.lut$rxn_id)) # 25771
length(unique(rxn_ECs.lut$EC_id)) # 7354


## Unique_ModelSEED_Reaction_Pathways

# rxn_pathways.lut <- read_excel(path = paste0(seed_db_dir,"/","Unique_ModelSEED_Reaction_Pathways.xlsx") , range = "A1:C121445")
# rxn_pathways.lut <- as.data.frame(rxn_pathways.lut)
rxn_pathways.lut <- read.csv(file = paste0(seed_db_dir,"/","Unique_ModelSEED_Reaction_Pathways.txt"), header = TRUE, sep = "\t")
class(rxn_pathways.lut) # "data.frame"

head(rxn_pathways.lut)
names(rxn_pathways.lut) # "ModelSEED.ID" "External.ID"  "Source"
names(rxn_pathways.lut) <- c("rxn_id", "External_rxn_name",  "Source")

dim(rxn_pathways.lut) # 121444      3
length(unique(rxn_pathways.lut$rxn_id)) # 21734
length(unique(rxn_pathways.lut$External_rxn_name)) # 3819

saveRDS(object = rxn_pathways.lut, file = "/Users/lidd0026/WORKSPACE/DATA/ModelSEEDDB/select_files/rxn_pathways.lut.RDS")



## Model SEED Subsystems

subsys.lut <- read_excel(path = paste0(seed_db_dir,"/","ModelSEED_Subsystems.xlsx") , range = "A1:E9821")
subsys.lut <- as.data.frame(subsys.lut)
dim(subsys.lut) # 9820    5

head(subsys.lut)
tail(subsys.lut)
names(subsys.lut) # "Class"     "Sub-class" "Name"      "Role"      "Reaction" 
names(subsys.lut) <- c("Class",    "Subclass", "Name",  "Role",     "Reaction")

subsys.lut[sample(1:dim(subsys.lut)[1], 10,replace = FALSE), ]
#                                     Class                                    Subclass                                                Name                                                                           Role Reaction
# 2984          Clustering-based subsystems                                           -                              CBSS-196620.1.peg.2477                                      Maltose O-acetyltransferase (EC 2.3.1.79) rxn01133
# 6834 Fatty Acids, Lipids, and Isoprenoids                                 Fatty acids                       Fatty_Acid_Biosynthesis_FASII          (3R)-hydroxymyristoyl-[acyl carrier protein] dehydratase (EC 4.2.1.-) rxn05427
# 7057 Fatty Acids, Lipids, and Isoprenoids                                 Fatty acids                   Unsaturated_Fatty_Acid_Metabolism                  3-oxoacyl-[acyl-carrier-protein] synthase, KASI (EC 2.3.1.41) rxn05350
# 9384                      Stress Response                              Osmotic stress Choline_and_Betaine_Uptake_and_Betaine_Biosynthesis                   Glycine betaine ABC transport system, permease protein OpuAB rxn05181
# 5904              Experimental Subsystems                                           -                              Transporters_In_Models                 Oligopeptide transport ATP-binding protein oppF (TC 3.A.1.5.1) rxn05539
# 5267              Experimental Subsystems                                           -              Sugar_catabolome_in_Shewanella_species                                              Alpha-galactosidase (EC 3.2.1.22) rxn00818
# 2626                Cell Wall and Capsule   Capsular and extracellular polysacchrides                                 Alginate_metabolism                                      Acetoin (diacetyl) reductase (EC 1.1.1.5) rxn01685
# 8538                   Protein Metabolism                        Protein biosynthesis                       Translation_factors_bacterial                                       Methionine aminopeptidase (EC 3.4.11.18) rxn12635
# 1853                        Carbohydrates                                Fermentation                 Acetyl-CoA_fermentation_to_Butyrate                     Acyl-CoA dehydrogenase, short-chain specific (EC 1.3.99.2) rxn10012
# 776           Amino Acids and Derivatives Lysine, threonine, methionine, and cysteine                     Lysine_Biosynthesis_DAP_Pathway 2,3,4,5-tetrahydropyridine-2,6-dicarboxylate N-acetyltransferase (EC 2.3.1.89) rxn03030
## Correspondences:

# Superfocus (imported 'tab' object): Subsystem Level 1 = subsys.lut$Class
#                                     Subsystem Level 2 = subsys.lut$Subclass
#    (skip the subsys.lut$Name field)
#                                     Function = subsys.lut$Role  # replace " " with "_" to match Superfocus output

dim(subsys.lut) # 9820    5
length(unique(subsys.lut$Class )) # 29
length(unique(subsys.lut$Subclass )) # 143
length(unique(subsys.lut$Name )) # 666
length(unique(subsys.lut$Role )) # 1885
length(unique(subsys.lut$Reaction )) # 1986

saveRDS(object = subsys.lut, file = "/Users/lidd0026/WORKSPACE/DATA/ModelSEEDDB/select_files/subsys.lut.RDS")


## Reactions

rxns.lut <- read_excel(path = paste0(seed_db_dir,"/","reactions.xlsx") , range = "A1:V43775")
rxns.lut <- as.data.frame(rxns.lut)
dim(rxns.lut) # 43774    22
names(rxns.lut)
# [1] "id"                "abbreviation"      "name"              "code"              "stoichiometry"     "is_transport"      "equation"          "definition"       
# [9] "reversibility"     "direction"         "abstract_reaction" "pathways"          "aliases"           "ec_numbers"        "deltag"            "deltagerr"        
# [17] "compound_ids"      "status"            "is_obsolete"       "linked_reaction"   "notes"             "source" 

head(rxns.lut)

length(unique(rxns.lut$id)) # 43774
length(unique(rxns.lut$name)) # 27668
length(unique(rxns.lut$code)) # 36212
length(unique(rxns.lut$equation)) # 37310
length(unique(rxns.lut$aliases)) # 32718    # SEARCH aliases ???????
length(unique(rxns.lut$ec_numbers)) # 7607

saveRDS(object = rxns.lut, file = "/Users/lidd0026/WORKSPACE/DATA/ModelSEEDDB/select_files/rxns.lut.RDS")


## Compounds

compounds.lut <- read_excel(path = paste0(seed_db_dir,"/","compounds.xlsx") , range = "A1:T33993")
warnings()
compounds.lut <- as.data.frame(compounds.lut)
dim(compounds.lut) # 33992    20
names(compounds.lut)
# [1] "id"                "abbreviation"      "name"              "formula"           "mass"              "source"            "inchikey"          "charge"           
# [9] "is_core"           "is_obsolete"       "linked_compound"   "is_cofactor"       "deltag"            "deltagerr"         "pka"               "pkb"              
# [17] "abstract_compound" "comprised_of"      "aliases"           "smiles"

head(compounds.lut)

length(unique(compounds.lut$id)) # 33992
length(unique(compounds.lut$formula)) # 16763

saveRDS(object = compounds.lut, file = "/Users/lidd0026/WORKSPACE/DATA/ModelSEEDDB/select_files/compounds.lut.RDS")


#-------------------------

#### Atomic ratio lookup function? O, C, H, N & map all ModelSEED compounds
#    Include multi-elements & element counts per 'beyond vK' paper
#    Include all critical nutrients for microbes
#-------------------------

# from earlier
## Compounds

compounds.lut <- read_excel(path = paste0(seed_db_dir,"/","compounds.xlsx") , range = "A1:T33993") # excludes column 'U'
# warnings()
compounds.lut <- as.data.frame(compounds.lut)
dim(compounds.lut) # 33992    20
names(compounds.lut)
# [1] "id"                "abbreviation"      "name"              "formula"           "mass"              "source"            "inchikey"          "charge"           
# [9] "is_core"           "is_obsolete"       "linked_compound"   "is_cofactor"       "deltag"            "deltagerr"         "pka"               "pkb"              
# [17] "abstract_compound" "comprised_of"      "aliases"           "smiles"

head(compounds.lut)
head(compounds.lut[ ,1:7])
#         id abbreviation  name       formula mass           source                    inchikey
# 1 cpd00001          h2o   H2O           H2O   18 Primary Database XLYOFNOQVPJJNP-UHFFFAOYSA-N
# 2 cpd00002          atp   ATP C10H13N5O13P3  504 Primary Database ZKHQWZAMYRWXGA-KQYNXXCUSA-K
# 3 cpd00003          nad   NAD C21H26N7O14P2  662 Primary Database BAWFJGJZGIEFAR-NNYOXOHSSA-M
# 4 cpd00004         nadh  NADH C21H27N7O14P2  663 Primary Database BOPGDPNILDQYTO-NNYOXOHSSA-L
# 5 cpd00005        nadph NADPH C21H26N7O17P3  742 Primary Database ACFIXJIJDZMPPO-NNYOXOHSSA-J
# 6 cpd00006         nadp  NADP C21H25N7O17P3  741 Primary Database XJLXINKUBYWONI-NNYOXOHSSA-K

length(unique(compounds.lut$id)) # 33992
length(unique(compounds.lut$formula)) # 16763
length(unique(compounds.lut$name)) # 33844


## get all combinations of 'O?' 'C?' 'H?' - to isolate oxygen, carbon, hydrogen containing compounds

# https://pmc.ncbi.nlm.nih.gov/articles/PMC4100946/
#   Elemental Economy: microbial strategies for optimizing growth in the face of nutrient limitation
# 
# Merchant, S.S. & Helmann, J.D. Chapter 2 - Elemental Economy: Microbial Strategies for Optimizing Growth in the Face of Nutrient Limitation. In: Poole, R.K. (ed). Advances in Microbial Physiology, vol. 60. Academic Press, 2012, pp 91-210.
# 
# Table 1 - include 'Required for all cells' and 'Required for most cells'
# 
# Required for All Cells
# C	basis of all organic molecules
# H	H2O, organic molecules
# N	organic molecules, esp. proteins and nucleic acids
# O	H2O, organic molecules
# # Already done
# 
# P	nucleic acids, NTPs, metabolites, phospholipids
# S	proteins, glutathione and LMW thiols, biotin, lipoic acid, thiamin
# Mg	major cation; cofactor for phosphotransferase reactions
# Zn	enzyme cofactor, protein folding
# 
# Required for Most Cells
# K	major cation, common in cells
# Ca	major cation, required by many eukaryotes
# Mn	enzyme cofactor, ribonucleotide reductase, SOD, PS II
# Fe	heme, iron-sulfur cluster, non-heme enzymes
# Co	enzyme cofactor, B12-dependent enzymes
# Cu	enzyme cofactor, electron carrier, respiration, SOD
# Mo	FeMoCo cofactor (nitrogenase), Mo cofactor enzymes



temp <- compounds.lut

sel <- grep(pattern = "Na", x = temp$formula)

sel <- grep(pattern = "P", x = temp$formula)
sel <- grep(pattern = "Si", x = temp$formula)
sel <- grep(pattern = "Mg", x = temp$formula)
sel <- grep(pattern = "Zn", x = temp$formula)
sel <- grep(pattern = "K", x = temp$formula)
sel <- grep(pattern = "Ca", x = temp$formula)
sel <- grep(pattern = "Mn", x = temp$formula)
sel <- grep(pattern = "Fe", x = temp$formula)
sel <- grep(pattern = "Co", x = temp$formula)
sel <- grep(pattern = "Cu", x = temp$formula)
sel <- grep(pattern = "Mo", x = temp$formula)


temp$formula[sel[1:100]]


#.          1.      2.      3.       4.       5.       6.        7.        8.      9.     10.          11             12.        13.   14.       15.               16
form <- c("NaCl", "H2O", "CaOH", "CaOCl2Rh","HgCl","C5H8O7PR", "MoOH2", "CoCH3", "OsR", "ScO2", "C10H12N4O12P2", "C15H27N5O5", "HgR", "Zn",
          
          "C25H26N9NaO8S2", "C12H7Cl2NNaO2", "C21H26N7O14P2", "H9AlFeMgO15Si4", "C55H75N4O6Zn", "C16H16K3O6PS", "C23H30CaN3Na3O11",
          
          "C4H6MnN2R2S4", "C33H30FeN4O4", "C55H83CoN14O9R2", "C7H4ClCuNO3S", "C20H22MoN10O15P2S2", "")

atom <- "O" # avoid: Os, Og
atom <- "C" # avoid: Cs, Ca, Ce, Cr, Co, Cm, Cu, Cd, Cn, Cf, Cl, 
atom <- "H" # avoid: Hf, Hs, Hg, Ho, He
atom <- "N" # avoid: Na, Nd, Ne, Np, Ni, Nh, Nb, No

atom <- "P" # avoid: Pd, Pt, Pb, Po
atom <- "S" # avoid: Si, Sc, Se, Sr, Sn, Sb
atom <- "Mg" #
atom <- "Zn" # 

atom <- "K" # avoid Kr 
atom <- "Ca" # 
atom <- "Mn" # 
atom <- "Fe" # 
atom <- "Co" # 
atom <- "Cu" # 
atom <- "Mo" # 



atomic_no <- function(form, atom) {
  ##form=form
  ##atom="O"
  ##atom="N"
  coefs <- list()
  for (c in 1:length(form)) {
    #c<-1
    pos <- str_locate(string = form[c], pattern = atom)[1]
    if (is.na(pos)) { # not present
      ##print("zero")
      coef<-0
    } else if (pos==nchar(form[c])) { # at the end, i.e. coef = 1.    # A 
      ##print("A")
      coef<-1
    } else {
      # check for off-targets ... otherwise decide coef = 1, 2, etc?       # B
      ##print("B")
      checkstring <- substring(text = form[c], first = pos, last = pos+1)
      
      if (atom == "O" & checkstring %in% c("Os", "Og")) {   # C
        ##print("C")
        coef<-0
      } else if (atom == "C" & checkstring %in% c("Cs", "Ca", "Ce", "Cr", "Co", "Cm", "Cu", "Cd", "Cn", "Cf", "Cl")) { # D
        ##print("D")
        coef<-0
      } else if (atom == "H" & checkstring %in% c("Hf", "Hs", "Hg", "Ho", "He")) { # E 
        ##print("E")
        coef<-0
        
      } else if (atom == "N" & checkstring %in% c("Na", "Nd", "Ne", "Np", "Ni", "Nh", "Nb", "No")) { # F 
        ##print("F")
        coef<-0
        
      } else if (atom == "P" & checkstring %in% c("Pd", "Pt", "Pb", "Po")) { # G 
        ##print("G")
        coef<-0
        
      } else if (atom == "S" & checkstring %in% c("Si", "Sc", "Se", "Sr", "Sn", "Sb")) { # H 
        ##print("H")
        coef<-0
        
      } else if (atom == "Mg" & checkstring %in% c("")) { # I 
        ##print("I")
        coef<-0
        
      } else if (atom == "Zn" & checkstring %in% c("")) { # J 
        ##print("J")
        coef<-0
        
      } else if (atom == "K" & checkstring %in% c("Kr")) { # K 
        ##print("K")
        coef<-0
        
      } else if (atom == "Ca" & checkstring %in% c("")) { # L 
        ##print("L")
        coef<-0
        
      }else if (atom == "Mn" & checkstring %in% c("")) { # M 
        ##print("M")
        coef<-0
        
      }else if (atom == "Fe" & checkstring %in% c("")) { # N 
        ##print("N")
        coef<-0
        
      }else if (atom == "Co" & checkstring %in% c("")) { # O 
        ##print("O")
        coef<-0
        
      }else if (atom == "Cu" & checkstring %in% c("")) { # P 
        ##print("P")
        coef<-0
        
      }else if (atom == "Mo" & checkstring %in% c("")) { # Q 
        ##print("Q")
        coef<-0
        
      } else { # R
        ##print("R")
        coef<-1
        coef.try <- coef
        pos.end <- pos+1
        
        while (is.numeric(coef.try) & !is.na(coef.try) & pos.end <= nchar(form[c])+1 ) {
          # substring allows last index to be too long & need to allow 'pos.end <= nchar(form[c])+1' for atomic numbers at end of formula
          coef <- coef.try
          coef.try <- as.numeric( substring(text = form[c], first = pos+1, last = pos.end) )
          pos.end <- pos.end+1
          ##print("G")
        } # END while loop
      } # END else
      
    } # END else
    
    coefs[[c]] <- coef
    
  } # END c loop
  return( unlist(coefs) )
}

form
# [1] "NaCl"               "H2O"                "CaOH"               "CaOCl2Rh"           "HgCl"               "C5H8O7PR"           "MoOH2"             
# [8] "CoCH3"              "OsR"                "ScO2"               "C10H12N4O12P2"      "C15H27N5O5"         "HgR"                "Zn"                
# [15] "C25H26N9NaO8S2"     "C12H7Cl2NNaO2"      "C21H26N7O14P2"      "H9AlFeMgO15Si4"     "C55H75N4O6Zn"       "C16H16K3O6PS"       "C23H30CaN3Na3O11"  
# [22] "C4H6MnN2R2S4"       "C33H30FeN4O4"       "C55H83CoN14O9R2"    "C7H4ClCuNO3S"       "C20H22MoN10O15P2S2" "" 

atomic_no(form, "O")
# [1]  0  1  1  1  0  7  1  0  0  2 12  5  0  0  8  2 14 15  6  6 11  0  4  9  3 15  0

atomic_no(form="NaCl",atom="C")
# [1] 0
atomic_no(form=form,atom="C")
# [1]  0  0  0  0  0  5  0  0  0  0 10 15  0  0 25 12 21  0 55 16 23  4 33 55  7 20  0
atomic_no(form=form,atom="H")
# 0  2  1  0  0  8  2  3  0  0 12 27  0  0 26  7 26  9 75 16 30  6 30 83  4 22  0
atomic_no(form=form,atom="N")
# 0 0 0 0 0 0 0 0 0 0 4 5 0 0 9 1
atomic_no(form=form,atom="P")
# [1] 0 0 0 0 0 1 0 0 0 0 2 0 0 0 0 0 2 0 0 1 0 0 0 0 0 2 0
atomic_no(form=form,atom="S")
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 1 0 4 0 0 1 2 0
atomic_no(form=form,atom="Mg")
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0
atomic_no(form=form,atom="Zn")
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0
atomic_no(form=form,atom="K")
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0
atomic_no(form=form,atom="Ca")
# [1] 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0
atomic_no(form=form,atom="Mn")
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0
atomic_no(form=form,atom="Fe")
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0
atomic_no(form=form,atom="Co")
# [1] 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0
atomic_no(form=form,atom="Cu")
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0
atomic_no(form=form,atom="Mo")
# [1] 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0


df.comp <- data.frame(id=compounds.lut$id, 
                      abbrev=compounds.lut$abbreviation,
                      name=compounds.lut$name,
                      form=compounds.lut$formula,
                      OC_ratio=NA,
                      HC_ratio=NA,
                      NC_ratio=NA,
                      
                      PC_ratio=NA,
                      
                      # new
                      NP_ratio = NA,
                      O_count = NA,
                      N_count = NA,
                      P_count = NA,
                      S_count = NA,
                      mass = NA,
                      
                      SC_ratio=NA,
                      MgC_ratio=NA,
                      ZnC_ratio=NA,
                      
                      KC_ratio=NA,
                      CaC_ratio=NA,
                      MnC_ratio=NA,
                      FeC_ratio=NA,
                      CoC_ratio=NA,
                      CuC_ratio=NA,
                      MoC_ratio=NA
                      
)

# uses function:

atomic_no
# function(form, atom) {
#   ##form=form
#   ##atom="O"
#   ##atom="N"
#   coefs <- list()
#   for (c in 1:length(form)) {
#     #c<-1
#     pos <- str_locate(string = form[c], pattern = atom)[1]
#     if (is.na(pos)) { # not present
#       ##print("zero")
#       coef<-0
#     } else if (pos==nchar(form[c])) { # at the end, i.e. coef = 1.    # A 
#       ##print("A")
#       coef<-1
#     } else {
#       # check for off-targets ... otherwise decide coef = 1, 2, etc?       # B
#       ##print("B")
#       checkstring <- substring(text = form[c], first = pos, last = pos+1)
#       
#       if (atom == "O" & checkstring %in% c("Os", "Og")) {   # C
#         ##print("C")
#         coef<-0
#       } else if (atom == "C" & checkstring %in% c("Cs", "Ca", "Ce", "Cr", "Co", "Cm", "Cu", "Cd", "Cn", "Cf", "Cl")) { # D
#         ##print("D")
#         coef<-0
#       } else if (atom == "H" & checkstring %in% c("Hf", "Hs", "Hg", "Ho", "He")) { # E 
#         ##print("E")
#         coef<-0
#         
#       } else if (atom == "N" & checkstring %in% c("Na", "Nd", "Ne", "Np", "Ni", "Nh", "Nb", "No")) { # F 
#         ##print("F")
#         coef<-0
#         
#       } else if (atom == "P" & checkstring %in% c("Pd", "Pt", "Pb", "Po")) { # G 
#         ##print("G")
#         coef<-0
#         
#       } else if (atom == "S" & checkstring %in% c("Si", "Sc", "Se", "Sr", "Sn", "Sb")) { # H 
#         ##print("H")
#         coef<-0
#         
#       } else if (atom == "Mg" & checkstring %in% c("")) { # I 
#         ##print("I")
#         coef<-0
#         
#       } else if (atom == "Zn" & checkstring %in% c("")) { # J 
#         ##print("J")
#         coef<-0
#         
#       } else if (atom == "K" & checkstring %in% c("Kr")) { # K 
#         ##print("K")
#         coef<-0
#         
#       } else if (atom == "Ca" & checkstring %in% c("")) { # L 
#         ##print("L")
#         coef<-0
#         
#       }else if (atom == "Mn" & checkstring %in% c("")) { # M 
#         ##print("M")
#         coef<-0
#         
#       }else if (atom == "Fe" & checkstring %in% c("")) { # N 
#         ##print("N")
#         coef<-0
#         
#       }else if (atom == "Co" & checkstring %in% c("")) { # O 
#         ##print("O")
#         coef<-0
#         
#       }else if (atom == "Cu" & checkstring %in% c("")) { # P 
#         ##print("P")
#         coef<-0
#         
#       }else if (atom == "Mo" & checkstring %in% c("")) { # Q 
#         ##print("Q")
#         coef<-0
#         
#       } else { # R
#         ##print("R")
#         coef<-1
#         coef.try <- coef
#         pos.end <- pos+1
#         
#         while (is.numeric(coef.try) & !is.na(coef.try) & pos.end <= nchar(form[c])+1 ) {
#           # substring allows last index to be too long & need to allow 'pos.end <= nchar(form[c])+1' for atomic numbers at end of formula
#           coef <- coef.try
#           coef.try <- as.numeric( substring(text = form[c], first = pos+1, last = pos.end) )
#           pos.end <- pos.end+1
#           ##print("G")
#         } # END while loop
#       } # END else
#       
#     } # END else
#     
#     coefs[[c]] <- coef
#     
#   } # END c loop
#   return( unlist(coefs) )
# }
# <bytecode: 0x7f778d8cb2e8>



for (i in 1:length(df.comp$id)) {
  #i<-1
  formx.char <- df.comp$form[i]
  
  # O:C ratio (replace Inf with NA when C not present)
  df.comp$OC_ratio[i] <- atomic_no(form = formx.char, atom = "O")/atomic_no(form = formx.char, atom = "C")
  df.comp$OC_ratio[i][is.infinite(df.comp$OC_ratio[i])] <- NA
  
  # H:C ratio (replace Inf with NA when C not present)
  df.comp$HC_ratio[i] <- atomic_no(form = formx.char, atom = "H")/atomic_no(form = formx.char, atom = "C")
  df.comp$HC_ratio[i][is.infinite(df.comp$HC_ratio[i])] <- NA
  
  # N:C ratio (replace Inf with NA when C not present)
  df.comp$NC_ratio[i] <- atomic_no(form = formx.char, atom = "N")/atomic_no(form = formx.char, atom = "C")
  df.comp$NC_ratio[i][is.infinite(df.comp$NC_ratio[i])] <- NA
  
  # PC_ratio
  df.comp$PC_ratio[i] <- atomic_no(form = formx.char, atom = "P")/atomic_no(form = formx.char, atom = "C")
  df.comp$PC_ratio[i][is.infinite(df.comp$PC_ratio[i])] <- NA
  
  # new
  # NP_ratio
  df.comp$NP_ratio[i] <- atomic_no(form = formx.char, atom = "N")/atomic_no(form = formx.char, atom = "P")
  df.comp$NP_ratio[i][is.infinite(df.comp$NP_ratio[i])] <- NA
  
  # O_count
  df.comp$O_count[i] <- atomic_no(form = formx.char, atom = "O")
  
  # N_count
  df.comp$N_count[i] <- atomic_no(form = formx.char, atom = "N")
  
  # P_count
  df.comp$P_count[i] <- atomic_no(form = formx.char, atom = "P")
  
  # S_count
  df.comp$S_count[i] <- atomic_no(form = formx.char, atom = "S")
  
  # mass
  df.comp$mass[i] <- compounds.lut$mass[i]
  
  # SC_ratio
  df.comp$SC_ratio[i] <- atomic_no(form = formx.char, atom = "S")/atomic_no(form = formx.char, atom = "C")
  df.comp$SC_ratio[i][is.infinite(df.comp$SC_ratio[i])] <- NA
  
  # MgC_ratio
  df.comp$MgC_ratio[i] <- atomic_no(form = formx.char, atom = "Mg")/atomic_no(form = formx.char, atom = "C")
  df.comp$MgC_ratio[i][is.infinite(df.comp$MgC_ratio[i])] <- NA
  
  # ZnC_ratio
  df.comp$ZnC_ratio[i] <- atomic_no(form = formx.char, atom = "Zn")/atomic_no(form = formx.char, atom = "C")
  df.comp$ZnC_ratio[i][is.infinite(df.comp$ZnC_ratio[i])] <- NA
  
  # KC_ratio
  df.comp$KC_ratio[i] <- atomic_no(form = formx.char, atom = "K")/atomic_no(form = formx.char, atom = "C")
  df.comp$KC_ratio[i][is.infinite(df.comp$KC_ratio[i])] <- NA
  
  # CaC_ratio
  df.comp$CaC_ratio[i] <- atomic_no(form = formx.char, atom = "Ca")/atomic_no(form = formx.char, atom = "C")
  df.comp$CaC_ratio[i][is.infinite(df.comp$CaC_ratio[i])] <- NA
  
  # MnC_ratio
  df.comp$MnC_ratio[i] <- atomic_no(form = formx.char, atom = "Mn")/atomic_no(form = formx.char, atom = "C")
  df.comp$MnC_ratio[i][is.infinite(df.comp$MnC_ratio[i])] <- NA
  
  # FeC_ratio
  df.comp$FeC_ratio[i] <- atomic_no(form = formx.char, atom = "Fe")/atomic_no(form = formx.char, atom = "C")
  df.comp$FeC_ratio[i][is.infinite(df.comp$FeC_ratio[i])] <- NA
  
  # CoC_ratio
  df.comp$CoC_ratio[i] <- atomic_no(form = formx.char, atom = "Co")/atomic_no(form = formx.char, atom = "C")
  df.comp$CoC_ratio[i][is.infinite(df.comp$CoC_ratio[i])] <- NA
  
  # CuC_ratio
  df.comp$CuC_ratio[i] <- atomic_no(form = formx.char, atom = "Cu")/atomic_no(form = formx.char, atom = "C")
  df.comp$CuC_ratio[i][is.infinite(df.comp$CuC_ratio[i])] <- NA
  
  # MoC_ratio
  df.comp$MoC_ratio[i] <- atomic_no(form = formx.char, atom = "Mo")/atomic_no(form = formx.char, atom = "C")
  df.comp$MoC_ratio[i][is.infinite(df.comp$MoC_ratio[i])] <- NA
  
  print(paste0("completed ",i))
  
}



saveRDS(object = df.comp, file = "df.comp__multi-element-v2.RDS")



## Exploring ... Add additional Compound classes from Rivas-Ubach 2018

df.comp2 <- df.comp
head(df.comp2)
# id abbrev  name          form  OC_ratio HC_ratio  NC_ratio  PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 1 cpd00001    h2o   H2O           H2O        NA       NA       NaN       NaN      NaN       1       0       0       0   18      NaN       NaN       NaN      NaN       NaN
# 2 cpd00002    atp   ATP C10H13N5O13P3 1.3000000 1.300000 0.5000000 0.3000000 1.666667      13       5       3       0  504        0         0         0        0         0
# 3 cpd00003    nad   NAD C21H26N7O14P2 0.6666667 1.238095 0.3333333 0.0952381 3.500000      14       7       2       0  662        0         0         0        0         0
# 4 cpd00004   nadh  NADH C21H27N7O14P2 0.6666667 1.285714 0.3333333 0.0952381 3.500000      14       7       2       0  663        0         0         0        0         0
# 5 cpd00005  nadph NADPH C21H26N7O17P3 0.8095238 1.238095 0.3333333 0.1428571 2.333333      17       7       3       0  742        0         0         0        0         0
# 6 cpd00006   nadp  NADP C21H25N7O17P3 0.8095238 1.190476 0.3333333 0.1428571 2.333333      17       7       3       0  741        0         0         0        0         0
# MnC_ratio FeC_ratio CoC_ratio CuC_ratio MoC_ratio
# 1       NaN       NaN       NaN       NaN       NaN
# 2         0         0         0         0         0
# 3         0         0         0         0         0
# 4         0         0         0         0         0
# 5         0         0         0         0         0
# 6         0         0         0         0         0


names(df.comp2)
# [1] "id"        "abbrev"    "name"      "form"      "OC_ratio"  "HC_ratio"  "NC_ratio"  "PC_ratio"  "NP_ratio"  "O_count"   "N_count"   "P_count"   "S_count"   "mass"     
# [15] "SC_ratio"  "MgC_ratio" "ZnC_ratio" "KC_ratio"  "CaC_ratio" "MnC_ratio" "FeC_ratio" "CoC_ratio" "CuC_ratio" "MoC_ratio"

# all start as unspecified
df.comp2$class <- "Unspecified"

# define classes using "OC_ratio"  "HC_ratio"  "NC_ratio"  "PC_ratio"  "NP_ratio"  "O_count"   "N_count"   "P_count"   "S_count"   "mass"
# based on Rivas-Ubach et al 2018 beyond vK paper

# Lipid: qty 1727
sel <- which(df.comp2$OC_ratio <= 0.6 &
               df.comp2$HC_ratio >= 1.32 & 
               df.comp2$NC_ratio <= 0.126 & 
               df.comp2$PC_ratio < 0.35 & 
               df.comp2$NP_ratio <= 5 )  # df.comp2$O_count   df.comp2$N_count   df.comp2$P_count   df.comp2$S_count   df.comp2$mass )
df.comp2$class[sel] <- "Lipid"

# Protein (constraints 1): qty 4597
sel <- which(df.comp2$OC_ratio > 0.12 & df.comp2$OC_ratio <= 0.6 &
               df.comp2$HC_ratio > 0.9 & df.comp2$HC_ratio < 2.5 &
               df.comp2$NC_ratio >= 0.126 & df.comp2$NC_ratio <= 0.7 &
               df.comp2$PC_ratio < 0.17 &
               
               df.comp2$N_count >= 1 )  # df.comp2$O_count   df.comp2$N_count   df.comp2$P_count   df.comp2$S_count   df.comp2$mass )
df.comp2$class[sel] <- "Protein"

# Protein (constraints 2): qty 1199
sel <- which(df.comp2$OC_ratio > 0.6 & df.comp2$OC_ratio <= 1 &
               df.comp2$HC_ratio > 1.2 & df.comp2$HC_ratio < 2.5 &
               df.comp2$NC_ratio > 0.2 & df.comp2$NC_ratio <= 0.7 &
               df.comp2$PC_ratio < 0.17 &
               
               df.comp2$N_count >= 1 )
df.comp2$class[sel] <- "Protein"

# A-Sugar: qty 317
sel <- which(df.comp2$OC_ratio >= 0.61 &
               df.comp2$HC_ratio >= 1.45 & 
               df.comp2$NC_ratio > 0.07 & df.comp2$NC_ratio <= 0.2 & 
               df.comp2$PC_ratio < 0.3 &
               df.comp2$NP_ratio <= 2 &
               df.comp2$O_count >= 3 &
               df.comp2$N_count >= 1
               )  # df.comp2$O_count   df.comp2$N_count   df.comp2$P_count   df.comp2$S_count   df.comp2$mass )
df.comp2$class[sel] <- "Amino sugar"

# Carbohydrate: qty 1214
sel <- which(df.comp2$OC_ratio >= 0.8 &
               df.comp2$HC_ratio >= 1.65 & df.comp2$HC_ratio < 2.7 & 
               
               df.comp2$N_count == 0
)  # df.comp2$O_count   df.comp2$N_count   df.comp2$P_count   df.comp2$S_count   df.comp2$mass )
df.comp2$class[sel] <- "Carbohydrate"

# Nucleotide: qty 169 
sel <- which(df.comp2$OC_ratio >= 0.5 & df.comp2$OC_ratio < 1.7 &
               df.comp2$HC_ratio > 1 & df.comp2$HC_ratio < 1.8 &
               df.comp2$NC_ratio >= 0.2 & df.comp2$NC_ratio <= 0.5 &
               df.comp2$PC_ratio >= 0.1 & df.comp2$PC_ratio < 0.35 &
               df.comp2$NP_ratio > 0.6 & df.comp2$NP_ratio <= 5 &

               df.comp2$N_count >= 2 &
               df.comp2$P_count >= 1 &
               df.comp2$S_count == 0 &
               df.comp2$mass > 305 & df.comp2$mass < 523 
               )  # df.comp2$O_count   df.comp2$N_count   df.comp2$P_count   df.comp2$S_count   df.comp2$mass )
df.comp2$class[sel] <- "Nucleotide"

# Phytochemical (oxy-aromatic): qty 114
sel <- which(df.comp2$OC_ratio <= 1.15 &
               df.comp2$HC_ratio < 1.32 & 
               df.comp2$NC_ratio < 0.126 & 
               df.comp2$PC_ratio <= 0.2 & 
               df.comp2$NP_ratio <= 3 )  # df.comp2$O_count   df.comp2$N_count   df.comp2$P_count   df.comp2$S_count   df.comp2$mass )
df.comp2$class[sel] <- "Phytochemical"

#df.comp2[sel, ]

dim(df.comp2) # 33992    25
length(which(df.comp2$class=="Unspecified")) # 24710

table(df.comp2$class)
# Amino sugar  Carbohydrate         Lipid    Nucleotide Phytochemical       Protein   Unspecified 
#         317          1214          1727           169           114          5741         24710 

saveRDS(object = df.comp2, file = "df.comp__multi-element-v2b.RDS")


shapes.class <- c("Amino sugar" = 2 ,
                  "Carbohydrate" = 0 ,
                  "Lipid" = 3,
                  "Nucleotide" = 5,
                  "Phytochemical" = 8  ,
                  "Protein" = 4,
                  "Unspecified" = 1 )


#-------------------------

#### Sun & Badgley - post-mining - RE-ANALYSIS of SUPER-FOCUS data from Bioenergetic mapping 
#### /CPP paper - https://doi.org/10.1016/j.scitotenv.2024.173543
#### now at individual compound level !

#### Sun & Badgley - post-mining - build reaction search in parallel - get_reactions & compounds
#-------------------------

phy <- readRDS("phy-phyloseq-object-sunbad-resto.RDS")

## convert each row in functional tax_table to "mean van Krevelen distance to health-associated compounds"

df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 30125     4


get_rxns_and_compounds_indiv <- function( df.tax, subsys.lut, rxns.lut, rxn_pathways.lut ) {
  
  rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
  rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
  
  sub1 <- df.tax$subsys_L1[i]
  sub2 <- df.tax$subsys_L2[i]
  sub3 <- df.tax$subsys_L3[i]
  
  fxn.temp <- df.tax$fxn[i]
  fxn.superfocus.rowlabel <- row.names(df.tax)[i]
  
  # store results corresponding to each Superfocus row
  fxn.list <- list()
  fxn.list[[ fxn.superfocus.rowlabel  ]] <- list()
  
  # check for multiple functions/reactions?
  flag1 <- grepl(pattern = "_/_|/", x = fxn.temp)
  flag2 <- grepl(pattern = "_@_", x = fxn.temp)
  if (!any(flag1,flag2)==TRUE) {
    # no multiples
    fxns <- fxn.temp
  } else if (flag1==TRUE) {
    fxns <- unlist( strsplit(fxn.temp, split = "_/_") )  ###### WHAT ABOUT SPLIT FOR "/" WITHOUT UNDERSCORES ??
  } else {
    fxns <- unlist( strsplit(fxn.temp, split = "_@_") )
  }
  # remove underscores
  ( fxns <- gsub(pattern = "_", replacement = " ", x = fxns) )
  
  # process each fxn & store attributes
  #df.fxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=1:length(fxns),`f__in`=fxns, matching_method=NA, min_adist_modelSEED=NA, min_amatch_modelSEED=NA, rxns=NA, tot_mean_OC_x=NA, tot_mean_HC_y=NA , tot_mean_NC_z=NA )
  
  df.fxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=1:length(fxns),`f__in`=fxns, matching_method=NA, rxns=NA) #, compound_name=NA, compound_form=NA, compound_molar_prop=NA, OC_x=NA, HC_y=NA , NC_z=NA )
  
  # # do round brackets interfere with search? - YES
  # lookfor <- "option 4 (this one)"
  # lookuplist <- c("option 1", "option 2", "option 3 (this one)", "option 4 (this one)")
  # grep(pattern = lookfor, x = lookuplist)
  
  # Identify '/' separators with no '_'  ??
  
  for (f in 1:length(fxns)) {  # this accounts for multiple functions/reactions reported in Superfocus outputs
    #f<-1
    #f<-2
    f.in <- fxns[f]
    
    # these concatenated expressions will be used to look for exact match using hierarchy in ModelSEED Subsystem table
    full_hier_target <- paste0(sub1,"__",sub2,"__",sub3,"__",f.in)
    full_hier_list <- paste0(subsys.lut$Class,"__",subsys.lut$Subclass,"__",gsub("_"," ",subsys.lut$Name),"__",subsys.lut$Role)
    
    ## data cleaning
    
    # trim off '_#' and '_##' tags
    trim_nchar <- str_locate(string = f.in, pattern = " # | ## ")[1]
    if (!is.na(trim_nchar) & length(trim_nchar)==1) {
      f.in <- substring(text = f.in , first = 1, last = trim_nchar-1)
    }
    
    # Eliminate unwanted parsing of regular expressions: '[', ']','***', '(', ')'
    f.in <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\} ", replacement ="." , x = f.in) # used later
    
    #rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
    #rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
    
    full_hier_target <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_target)
    full_hier_list <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_list)
    
    sel.rx <- grep(pattern = full_hier_target, x = full_hier_list)
    
    ## ALTERNATIVE #1 == FULL HIERACHICAL MATCH
    if (length(sel.rx)>=1) {
      df.fxns$matching_method[f] <- "Exact hierachy match"
      df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
      
    } else if (str_detect(string = fxns[f], pattern = " \\(EC ")) {  ## ALTERNATIVE #2 == MATCHING ECs
      # search by EC id if present
      
      f.in <- fxns[f] # this goes back to string with brackets for EC
      ## LOOK FOR MULTIPLE ECs ??????????
      # 22889
      # 22894
      # 31768
      
      how_many_ECs <- str_count(string = f.in, pattern = "\\(EC.*?\\)")
      
      ECs <- as.character( str_extract_all(string = f.in, pattern = "\\(EC.*?\\)", simplify = TRUE) )
      #class(ECs)
      ECs <- gsub(pattern = "\\(EC |\\)", replacement = "", x = ECs)
      ECs.collapse <- paste0(ECs, collapse = "|")
      
      sel.rx <- which(rxns.lut$ec_numbers == ECs.collapse)
      
      if (length(how_many_ECs)==0 | length(ECs)==0) {
        # there was a glitch, database typo, or some error in identifying the EC number
        df.fxns$matching_method[f] <- "No match found"
        df.fxns$rxns[f] <- NA
        
      } else if (length(sel.rx)>=1) {
        # combined EC hits identified
        df.fxns$matching_method[f] <- "EC number"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if (length(which(rxns.lut$ec_numbers %in% ECs)) >=1) {
        # treat EC hits individually
        sel.rx <- which(rxns.lut$ec_numbers %in% ECs) # look 1st where ECs are exact matches for EC numbers in Reactions lookup table
        
        df.fxns$matching_method[f] <- "EC number"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if (length(grep(pattern = ECs, x = rxns.lut$ec_numbers)) >=1) {
        # this allows EC to be part of a combination of EC numbers that are listed in Reactions lookup table
        sel.rx <- grep(pattern = ECs, x = rxns.lut$ec_numbers)
        
        df.fxns$matching_method[f] <- "EC number"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else {
        # it had an EC number but couldn't find a match in the EC numbers listed in Reaction lookup table
        df.fxns$matching_method[f] <- "No match found"
        df.fxns$rxns[f] <- NA
        
      }
      # END EC matching
      
      
    } else {  ## ALTERNATIVE 3 == FXN NAME MATCHING
      ## otherwise attempt to match function name - a) first look for exact matches   ########## then b) closest match above a threshold
      # 1. 'reactions' table by name: rxns.lut$name
      # 2. 'reactions' table by aliases: rxns.lut$aliases
      # 3. 'Model SEED Subsystems' table by Role: subsys.lut$Role
      # 4. 'Unique_ModelSEED_Reaction_Pathways' table by External ID: rxn_pathways.lut$External_rxn_name
      
      if ( length( grep(pattern = f.in, x = rxns.lut$name) )>=1 ) {
        # 1a - exact match - rxns.lut$name
        sel.rx <- grep(pattern = f.in, x = rxns.lut$name)
        #rxns.lut$name[sel.rx]
        df.fxns$matching_method[f] <- "Matched Reactions name"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if ( length( grep(pattern = f.in, x = rxns.lut$aliases) )>=1 ) {
        # 2a - exact match - rxns.lut$aliases
        sel.rx <- grep(pattern = f.in, x = rxns.lut$aliases)
        #rxns.lut$aliases[sel.rx]
        #rxns.lut$name[sel.rx]
        
        df.fxns$matching_method[f] <- "Matched Reactions aliases"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if ( length( grep(pattern = f.in, x = subsys.lut$Role) )>=1 ) {
        # 3a - exact match - subsys.lut$Role
        sel.rx <- grep(pattern = f.in, x = subsys.lut$Role)
        #subsys.lut$Role[sel.rx]
        #subsys.lut$Reaction[sel.rx]
        
        df.fxns$matching_method[f] <- "Matched Subsytem role"
        df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
        
      } else if ( length( grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name) )>=1 ) {
        # 4a - exact match - rxn_pathways.lut$External_rxn_name
        sel.rx <- grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name)
        
        df.fxns$matching_method[f] <- "Matched ModelSEED Reaction pathways"
        df.fxns$rxns[f] <- paste0( unique(rxn_pathways.lut$rxn_id[sel.rx]), collapse = ";")
        
        
      } else {
        df.fxns$matching_method[f] <- "No match found"
        df.fxns$rxns[f] <- NA
        
      }
      
      
    } # END function - reaction search
    
    #fxn.list[[ fxn.superfocus.rowlabel  ]][[ f ]][[ "fxns" ]] <- df.fxns
    
    print(paste0("completed fxn ", f))
    
    
    ## now investigate these reactions ...
    # Reactions lookup table: 
    # - "equation": Definition of reaction expressed using compound IDs and after protonation
    # Compounds lookup table:
    # - "formula": Standard chemical format (using Hill system) in protonated form to match reported charge
    #df.fxns
    
    
    #if (df.fxns$matching_method == "No match found") {
    if (df.fxns$rxns[f] == "" | is.na(df.fxns$rxns[f])) {
      
      df.Rxns <- NA
      df.Compounds <- NA
      
    } else { # reaction(s) were identified
      
      # consider reactions for this f.in only (possibly > 1 f.in per Superfocus row)
      f.in.rxns <- unique(unlist(str_split(string = df.fxns$rxns[f], pattern = ";")))
      
      df.Rxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel, f=f, f__in=f.in,rxn_id= f.in.rxns,
                            rxn_name=NA, rxn_eqn=NA, rxn_defn=NA,compds=NA,compd_coef=NA, chem_formx=NA ) #, OC_ratios=NA, HC_ratios=NA, NC_ratios=NA, coefwtmean_OC_x=NA, coefwtmean_HC_y=NA, coefwtmean_NC_z=NA)
      
      #df.Compounds <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=f, f__in=f.in,rxn_id= f.in.rxns, compound_name=NA, compound_form=NA, compound_molar_prop=NA, OC_x=NA, HC_y=NA , NC_z=NA )
      
      for (r in 1:dim(df.Rxns)[1]) {
        #r<-1
        #this_rxn <- "rxn00004"
        this_rxn <- df.Rxns$rxn_id[r]
        sel <- which(rxns.lut$id == this_rxn)
        ( df.Rxns$rxn_name[r] <- rxns.lut$name[sel] )
        ( df.Rxns$rxn_eqn[r] <- rxns.lut$equation[sel] )
        ( df.Rxns$rxn_defn[r] <- rxns.lut$definition[sel] )
        
        # extract compound info
        
        #df.Rxns$rxn_eqn[r]
        #[1] "(1) cpd00010[0] + (1) cpd29672[0] <=> (1) cpd00045[0] + (1) cpd11493[0]"
        #[1] "(45) cpd00144[0] + (45) cpd00175[0] <=> (45) cpd00014[0] + (45) cpd00091[0] + (1) cpd15634[0]"
        
        ( compds.idx <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd")[[1]][,"start"] )
        # 5 23 43 61
        # 6 25 46 65 83
        
        ( compds <- as.character( str_extract_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd.....", simplify = TRUE) ) )
        # "cpd00010" "cpd29672" "cpd00045" "cpd11493"
        
        if (length(compds)>=1) {
          
          df.Rxns$compds[r] <- paste0(compds, collapse = ";")
          
          ## get compound coefficients?
          start_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\(")[[1]][,"start"]
          end_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\)")[[1]][,"start"]
          ( compd.coeff <- as.numeric( substring(text = df.Rxns$rxn_eqn[r], first = start_brackets+1, last = end_brackets-1)) )
          
          df.Rxns$compd_coef[r] <- paste0(compd.coeff, collapse = ";")
          
          # get formulas of compounds
          
          formx <-filter(compounds.lut, id %in% compds )
          row.names(formx) <- formx$id
          ( formx.char <- formx[compds, ]$formula )
          # "C21H32N7O16P3S" "HOR"            "C10H11N5O10P2"  "C11H22N2O7PRS" 
          # "C15H19N2O18P2"      "C17H25N3O17P2"      "C9H12N2O12P2"       "C9H11N2O9P"         "C630H945N45O630P45"
          # "C7H7O7" "H2O"    "C7H5O6"
          df.Rxns$chem_formx[r] <- paste0(formx.char, collapse = ";")
          
          ( compd.names <- formx[compds, ]$name )
          # "2-methyl-trans-aconitate" "cis-2-Methylaconitate"
          
          
          temp.df.Compounds <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=f, f__in=f.in,rxn_id= f.in.rxns[r], 
                                     cpd_id=compds, cpd_name=compd.names, cpd_form=formx.char, cpd_molar_prop=compd.coeff #, 
                                     #OC_x=OC_ratio, HC_y=HC_ratio , NC_z=NC_ratio 
                                     )
          
        } else {
          # No specified reaction equation or chemical formula info
          df.Rxns$compds[r] <- NA
          df.Rxns$compd_coef[r] <- NA
          df.Rxns$chem_formx[r] <- NA
          
          temp.df.Compounds <- NA
         
        }
        
        if (r==1) { df.Compounds <- temp.df.Compounds }
        
        if (r>1 & is.data.frame(df.Compounds) & is.data.frame(temp.df.Compounds)) { df.Compounds <- rbind(df.Compounds, temp.df.Compounds) }
        
        # clean up - if there are additional reactions?
        temp.df.Compounds <- NA
        
      } # END loop for r - rxn_id's per f/f.in
      
    } # END else loop when reactions identified
    
    # store results corresponding to each sub-reaction of each Superfocus row
    fxn.list[[ fxn.superfocus.rowlabel  ]][[ "fxns" ]] <- df.fxns
    
    if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]] <- list() } # set this only once
    fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]][[ f ]] <- df.Rxns
    
    if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]] <- list() } # set this only once
    fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]][[ f ]] <- df.Compounds
    
  
  } # END loop - f in 1:length(fxns)) - to account for multiple functions/reactions reported in each row of Superfocus outputs
  
  
  #return(fxn.list)
  
  saveRDS(object = fxn.list, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/fxn-list-",fxn.superfocus.rowlabel,".rds") ) # use readRDS()
  
  #print(paste0("COMPLETED ROW ",i," OF SUPERFOCUS FUNCTIONAL TAXA  # # # # # # # # # # # # # # # # # # # # #"))
  
} # END function to be run in parallel for each superfocus row


# # # # # # # # # # # # # # # # # #


time.start <- Sys.time()
cl<-makeCluster( detectCores()-1 )
registerDoParallel(cl)

#foreach(i=1:100 , .packages=c('stringr', 'dplyr')) %dopar%
foreach(i=1:dim(df.tax)[1] , .packages=c('stringr', 'dplyr')) %dopar%  #
  get_rxns_and_compounds_indiv( df.tax=df.tax, subsys.lut=subsys.lut, rxns.lut=rxns.lut, rxn_pathways.lut=rxn_pathways.lut )

stopCluster(cl)
time.finish <- Sys.time()





## assemble results

# results saved here
#saveRDS(object = fxn.list, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/AMI-R-working-files-indiv/fxn-list-",fxn.superfocus.rowlabel,".rds") ) # use readRDS()


modelSEED_rxn_result_dir <- "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv"


dim(df.tax)
# 30125     4




# read first output
i<-1
#temp <- readRDS(paste0(modelSEED_rxn_result_dir,"/fxn-list-fxn_",i,".rds"))
fxn.superfocus.rowlabel <- row.names(df.tax)[i]
temp <- readRDS(paste0(modelSEED_rxn_result_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))

print( length(temp) )
print( names(temp) )
# "fxn_1"
print( paste0( class(temp[[1]])," ",names(temp)," of length ", length(temp[[1]])," named ",paste0(names(temp[[1]]), collapse = " & ") ))
# "list fxn_1 of length 3 named fxns & rxns & compounds"
class ( temp[[1]][["fxns"]] ) # data.frame
class ( temp[[1]][["rxns"]] ) # list
class ( temp[[1]][["compounds"]] ) # list

length( temp[[1]][["compounds"]] ) # 1
names( temp[[1]][["compounds"]] ) # NULL
class( temp[[1]][["compounds"]][[1]] ) # "logical"
is.na( temp[[1]][["compounds"]][[1]] )

i<-2
fxn.superfocus.rowlabel <- row.names(df.tax)[i]
temp <- readRDS(paste0(modelSEED_rxn_result_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))

length(temp) # 1
names(temp) # "fxn_2"
print( paste0( class(temp[[1]])," ",names(temp)," of length ", length(temp[[1]])," named ",paste0(names(temp[[1]]), collapse = " & ") ))
# "list fxn_2 of length 3 named fxns & rxns & compounds"
class ( temp[[1]][["fxns"]] ) # data.frame
class ( temp[[1]][["rxns"]] ) # list
class ( temp[[1]][["compounds"]] ) # list

length( temp[[1]][["compounds"]] ) # 1
names( temp[[1]][["compounds"]] ) # NULL
class( temp[[1]][["compounds"]][[1]] ) # "data.frame"

df.out <- temp[[1]][["compounds"]][[1]]

names(df.out) #
# [1] "superfocus_fxn" "f"              "f__in"          "rxn_id"         "cpd_id"         "cpd_name"       "cpd_form"       "cpd_molar_prop" 




# total results written to disk - 30125
dim(df.tax) # 30125     4
num_results_files <- dim(df.tax)[1]
#num_results_files <- 100


# assemble all compound data outputs
# start with blank row

df.out <- data.frame(superfocus_fxn=NA, f=NA, f__in=NA, rxn_id=NA, cpd_id=NA, cpd_name=NA, cpd_form=NA, cpd_molar_prop=NA #, 
                     #OC_x=NA, HC_y=NA, NC_z=NA
                     )

for (i in 1:num_results_files) {

  fxn.superfocus.rowlabel <- row.names(df.tax)[i]
  temp <- readRDS(paste0(modelSEED_rxn_result_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))
  
  f_no <- length( temp[[1]][["compounds"]] )
  
  for (f in 1:f_no) {
    #f<-2
    # only add non-NA results
    if (is.data.frame( temp[[1]][["compounds"]][[f]] )) {
      
      df.temp <- temp[[1]][["compounds"]][[f]]
      ok <- complete.cases(df.temp)
      df.temp <- df.temp[ which(ok==TRUE), ] # updated version will include some compounds with vK coordinates that are NA. vK coordinates are considered later
      df.out <- rbind(df.out,df.temp)
    }
  }
 
  print(paste0("added df ",i," of ",num_results_files ))
   
}



str(df.out)
# 'data.frame':	1154548 obs. of  8 variables:

saveRDS(object = df.out, file = "df.out--get_rxns_and_compounds_indiv--sunbad-resto.RDS")
df.out <- readRDS(file = "df.out--get_rxns_and_compounds_indiv--sunbad-resto.RDS")

# remove NA first row
head(df.out)
#   superfocus_fxn  f                                                                         f__in   rxn_id   cpd_id
# 1            <NA> NA                                                                          <NA>     <NA>     <NA>
# 2           fxn_2  1                                                   2-methylaconitate isomerase rxn25278 cpd25681
# 3           fxn_2  1                                                   2-methylaconitate isomerase rxn25278 cpd02597
# 11          fxn_3  1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd24620
# 31          fxn_3  1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd25681
# 14          fxn_4  1                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79) rxn03060 cpd01501
#                                        cpd_name cpd_form cpd_molar_prop      OC_x      HC_y NC_z
# 1                                          <NA>     <NA>             NA        NA        NA   NA
# 2                      2-methyl-trans-aconitate   C7H5O6              1 0.8571429 0.7142857    0
# 3                         cis-2-Methylaconitate   C7H5O6              1 0.8571429 0.7142857    0
# 11 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7              1 1.0000000 1.0000000    0
# 31                     2-methyl-trans-aconitate   C7H5O6              1 0.8571429 0.7142857    0
# 14                              2-Methylcitrate   C7H7O7              1 1.0000000 1.0000000    0

df.out <- df.out[-1, ]


# check for different cpd_molar_prop ??
hist(df.out$cpd_molar_prop)

dim(df.out) # 1154547       8


# normalise molar_prop to cpd_relabun so total of 1 per superfocus function !!

df.out$cpd_molar_prop_norm <- NA

length(unique(df.out$superfocus_fxn)) # 16518

phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 30125 taxa and 15 samples ]
# sample_data() Sample Data:       [ 15 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 30125 taxa by 4 taxonomic ranks ]

100*(length(unique(df.out$superfocus_fxn)) / ntaxa(phy)) # 54.83154 % of functions represented
100*(16518/30125) # 54.83154

fxns_found <- unique(df.out$superfocus_fxn)

for (k in 1:length(fxns_found)) {
  #k<-1
  this_fxn <- fxns_found[k]
  sel <- which(df.out$superfocus_fxn == this_fxn)
  
  sum_molar_prop <- sum( df.out$cpd_molar_prop[sel], na.rm = TRUE)
  # calculate 
  
  df.out$cpd_molar_prop_norm[sel] <- df.out$cpd_molar_prop[sel]/sum_molar_prop
  
  print(paste0("completed ",k))
  
}

sum(df.out$cpd_molar_prop_norm) # 16518


sample_sums(phy)
# 20C 30B 30A UMA 10B  5A 10C 20A 20B  5B UMB 10A  5C UMC 30C 
# 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 

dim(df.out) # 1154547       9




getwd() # "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R"

saveRDS(object = df.out, file = "df.out--tidy-compounds_indiv--cpp3d-sunbad-resto.RDS")

df.out <-readRDS( "df.out--tidy-compounds_indiv--cpp3d-sunbad-resto.RDS" )

#-------------------------

#### Sun & Badgley - post-mining
#    Collate CPP (compound rel abun %)
#    Test trending with restoration
#-------------------------

this_study <- "-sunbad-resto-"
header <- "cpp3d-cpdtrend"
phy <- readRDS("phy-phyloseq-object-sunbad-resto.RDS")
df.out <-readRDS( "df.out--tidy-compounds_indiv--cpp3d-sunbad-resto.RDS" )
dim(df.out) # 1154547       9

df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 30125     4


phy@sam_data
# Sample Data:        [15 samples by 8 sample variables]:
#   sample_name mgrast_id metagenome_id metagenome_name investigation_type       seq_meth file_name age
# 20C         20C mgl422999  mgm4679658.3       20C.fastq         metagenome pyrosequencing 20C.fastq  22
# 30B         30B mgl423005  mgm4679659.3       30B.fastq         metagenome pyrosequencing 30B.fastq  31
# 30A         30A mgl423002  mgm4679660.3       30A.fastq         metagenome pyrosequencing 30A.fastq  31
# UMA         UMA mgl423011  mgm4679661.3       UMA.fastq         metagenome pyrosequencing UMA.fastq  UM
# 10B         10B mgl422987  mgm4679662.3       10B.fastq         metagenome pyrosequencing 10B.fastq  12
# 5A           5A mgl422975  mgm4679663.3        5A.fastq         metagenome pyrosequencing  5A.fastq   6
# 10C         10C mgl422990  mgm4679664.3       10C.fastq         metagenome pyrosequencing 10C.fastq  12
# 20A         20A mgl422993  mgm4679665.3       20A.fastq         metagenome pyrosequencing 20A.fastq  22
# 20B         20B mgl422996  mgm4679666.3       20B.fastq         metagenome pyrosequencing 20B.fastq  22
# 5B           5B mgl422978  mgm4679667.3        5B.fastq         metagenome pyrosequencing  5B.fastq   6
# UMB         UMB mgl423014  mgm4679668.3       UMB.fastq         metagenome pyrosequencing UMB.fastq  UM
# 10A         10A mgl422984  mgm4679669.3       10A.fastq         metagenome pyrosequencing 10A.fastq  12
# 5C           5C mgl422981  mgm4679670.3        5C.fastq         metagenome pyrosequencing  5C.fastq   6
# UMC         UMC mgl423017  mgm4679671.3       UMC.fastq         metagenome pyrosequencing UMC.fastq  UM
# 30C         30C mgl423008  mgm4679672.3       30C.fastq         metagenome pyrosequencing 30C.fastq  31


# NOTE ADDITIONAL AVERAGE PLOT (n = 3) LEVEL DATA, available from Avera et al 2015
# New Forests (2015) 46:683–702
# DOI 10.1007/s11056-015-9502-8


sample_names(phy)
identical( sample_names(phy), colnames( as.matrix( phy@otu_table)) ) # TRUE

df.OTU <- as.data.frame( phy@otu_table ) # this is Superfocus functional relative abundance data represented in phyloseq OTU abundance table
dim(df.OTU) # 30125    15
df.OTU[1:5, 1:8]
# 20C          30B          30A          UMA          10B           5A          10C          20A
# fxn_1 5.233125e-04 5.972645e-05 1.551835e-04 1.195234e-04 2.501741e-04 5.509386e-04 3.930307e-04 4.019769e-04
# fxn_2 2.012740e-05 0.000000e+00 0.000000e+00 0.000000e+00 2.204177e-05 1.756616e-05 0.000000e+00 3.122151e-05
# fxn_3 1.363632e-04 2.526888e-05 4.877197e-05 2.988084e-05 1.146172e-04 1.250163e-04 9.707384e-05 1.951344e-05
# fxn_4 3.975162e-05 1.148586e-05 6.650723e-05 8.964252e-05 2.204177e-05 7.984617e-05 8.286791e-05 5.854033e-05
# fxn_5 3.371340e-03 2.216770e-03 2.194738e-03 1.942255e-03 2.182136e-03 2.626939e-03 2.379493e-03 2.566018e-03

mean( df.OTU[ , "20C"] ) # 0.003319502
sum( df.OTU[ , "20C"] ) # 100

sample_sums(phy) # all values of 100


1154547*15 # 17318205 = up to 17,318,205 rows





# loop through each sample

# add grouping variables

# for each function, assign relative abundance across selected compounds



get_cpd_relabun_per_sample <- function(phy_in, dat.cpd) {
  #i<-1
  #phy_in = phy
  #dat.cpd = df.out
  
  this_samp <- sample_names(phy_in)[i]
  df.OTU <- as.data.frame( phy_in@otu_table[ ,this_samp] )
  
  dat.cpd$sample <- this_samp
  
  dat.cpd$cpd_rel_abun_norm <- NA
  
  fxns_all <- row.names(df.OTU)
  
  for (k in 1:length(fxns_all)) {
    #k<-1
    this_fxn <- fxns_all[k]
    sel <- which(dat.cpd$superfocus_fxn == this_fxn)
    
    if (length(sel)>=1) {
      dat.cpd$cpd_rel_abun_norm[sel] <- df.OTU[this_fxn, ]*dat.cpd$cpd_molar_prop_norm[sel]
      
    }
  } # END rel abun values for all relevant functions added
  
  saveRDS(object = dat.cpd, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
  
} # END


time.start <- Sys.time()
cl<-makeCluster( detectCores()-1 )
registerDoParallel(cl)

foreach(i=1: length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
  get_cpd_relabun_per_sample( phy_in = phy, dat.cpd = df.out)

stopCluster(cl)
time.finish <- Sys.time()

# output 1
i<-1
this_samp <- sample_names(phy)[i]
#saveRDS(object = dat.cpd, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
dat <- readRDS ( file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
head(dat)

for ( i in 2:length(sample_names(phy)) ) {
  #i<-1
  this_samp <- sample_names(phy)[i]
  #saveRDS(object = dat.cpd, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
  temp <- readRDS ( file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
  
  dat <- rbind(dat, temp)
  
  print(paste0("completed ",i))
}


saveRDS(object = dat, file = "dat.cpd-long-all-samps-cpp3d-sunbad-resto.rds" )
dat <- readRDS("dat.cpd-long-all-samps-cpp3d-sunbad-resto.rds")

rm(temp)

str(dat)
# 'data.frame':	17318205 obs. of  11 variables:
# $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_3" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylaconitate isomerase" "2-methylaconitate isomerase" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" ...
# $ rxn_id             : chr  "rxn25278" "rxn25278" "rxn25279" "rxn25279" ...
# $ cpd_id             : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ cpd_name           : chr  "2-methyl-trans-aconitate" "cis-2-Methylaconitate" "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" ...
# $ cpd_form           : chr  "C7H5O6" "C7H5O6" "C7H7O7" "H2O" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.333 0.333 0.333 ...
# $ sample             : chr  "20C" "20C" "20C" "20C" ...
# $ cpd_rel_abun_norm  : num  1.01e-05 1.01e-05 4.55e-05 4.55e-05 4.55e-05 ...

#dat$combos <- paste0(dat$OC_x,"__",dat$HC_y,"__",dat$NC_z)

sum(dat$cpd_rel_abun_norm) # 1005.557
sum(dat$cpd_rel_abun_norm)/nsamples(phy) # 67.03712 = average 67% functional relative abundance per sample

length(which(is.na(dat$cpd_rel_abun_norm))==TRUE) # 0
length(which( dat$cpd_rel_abun_norm > 0) == TRUE) # 11960630
length(which( dat$cpd_rel_abun_norm == 0) == TRUE) # 5357575

# so this step does collect some zero relative abundances from the 'otu-table' 

names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"               "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm" "sample"              "cpd_rel_abun_norm" 

## create dat.cpd.distil
## later collapse to unique combos of "OC_x__HC_y__NC_z", add rel_abun, collate unique cpd_form, collate rxn_id, collate cpd_id 

length(unique(dat$cpd_id)) # 8370



## Collate compounds within each sample 


unique_cpd <- unique(dat$cpd_id)
samp_names <- sample_names(phy)




collate_compounds <- function(dat.cpd, unique_cpd, samp) {
  #i<-1
  #samp = samp_names[i]
  #dat.cpd = dat[which(dat$sample == samp_names[i]), ]
  
  this_samp <- samp
  
  cpd_data <- data.frame(cpd_id = unique_cpd, sample=this_samp, #OC_x=NA, HC_y=NA, NC_z=NA, 
                         cpd_rel_abun=NA)
  
  for (c in 1:length(unique_cpd)) {
    #c<-1
    this_cpd <- unique_cpd[c]
    sel.cpd <- which(dat.cpd$cpd_id == this_cpd)
    
    # cpd_data$OC_x[c] <- unique(dat.cpd$OC_x[sel.cpd])
    # cpd_data$HC_y[c] <- unique(dat.cpd$HC_y[sel.cpd])
    # cpd_data$NC_z[c] <- unique(dat.cpd$NC_z[sel.cpd])
    
    # # now collate rel abun for this sample only
    # sel <- which(dat.cpd$sample == this_samp)
    # dat.cpd <- dat.cpd[sel, ]
    
    # # again select in this sample-specific dataset
    # sel.cpd <- which(dat.cpd$cpd_id == this_cpd)
    
    if (length(sel.cpd) >=1) {
      cpd_data$cpd_rel_abun[c] <- sum(dat.cpd$cpd_rel_abun_norm[sel.cpd])
    }
    
  } # END all compounds
  
  saveRDS(object = cpd_data, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/cpd_data.collate-",this_samp,".rds") ) # use readRDS()
  
} # END


time.start <- Sys.time()
cl<-makeCluster( detectCores()-1 )
registerDoParallel(cl)

foreach(i=1: length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
  collate_compounds(dat.cpd = dat[which(dat$sample == samp_names[i]), ], unique_cpd = unique_cpd, samp = samp_names[i])

stopCluster(cl)
time.finish <- Sys.time()




# output 1
i<-1
this_samp <- sample_names(phy)[i]
dat.cpd.collate <- readRDS ( file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/cpd_data.collate-",this_samp,".rds") ) # use readRDS()
head(dat.cpd.collate)

for ( i in 2:length(sample_names(phy)) ) {
  #i<-1
  this_samp <- sample_names(phy)[i]
  temp <- readRDS ( file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/Sunbad-R-working-files-indiv/cpd_data.collate-",this_samp,".rds") ) # use readRDS()
  
  dat.cpd.collate <- rbind(dat.cpd.collate, temp)
  
  print(paste0("completed ",i))
}


str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  3 variables:
# $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "20C" "20C" "20C" "20C" ...
# $ cpd_rel_abun: num  0.000435 0.022738 0.000456 5.135411 0.015705 ...

sum(dat.cpd.collate$cpd_rel_abun) # 1005.557.   929.1978
sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample)) # 67.03712.   61.94652

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-sunbad-resto.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-sunbad-resto.rds")

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000013 0.000165 0.008009 0.001455 8.107427
# PREVIOUSLY, BASED ON C,O,H,N
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000020 0.000210 0.008008 0.002060 4.046109 

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf -4.9017 -3.7832    -Inf -2.8372  0.9089
# PREVIOUSLY, BASED ON C,O,H,N
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf  -4.693  -3.678    -Inf  -2.686   0.607


# log10 abun
dat.cpd.collate$log10_abun <- dat.cpd.collate$cpd_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(dat.cpd.collate$log10_abun == 0) # qty 8897 7945
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(dat.cpd.collate$log10_abun[ -subsel.zero ])
  dat.cpd.collate$log10_abun[ subsel.zero ] <- zero_replace
}
dat.cpd.collate$log10_abun <- log10(dat.cpd.collate$log10_abun)

hist(dat.cpd.collate$log10_abun); summary( dat.cpd.collate$log10_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -8.9875 -4.9017 -3.7832 -4.1549 -2.8372  0.9089

# PREVIOUSLY, BASED ON C,O,H,N
# Indiv
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -7.950  -4.693  -3.678  -3.943  -2.686   0.607
# Fxn-Mean
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -5.3909 -5.3909 -5.3909 -4.8501 -5.3909  0.3853




# make group variable from sample name

dat.cpd.collate$group <- NA


for (i in 1:length(sample_names(phy))) {
  #i<-1
  this_samp <- sample_names(phy)[i]
  sel <- which(dat.cpd.collate$sample == this_samp)
  dat.cpd.collate$group[sel] <- phy@sam_data$age[i]
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group) # "22" "31" "UM" "12" "6"
dat.cpd.collate$group <- factor(dat.cpd.collate$group, levels = c("6", "12", "22", "31", "UM"), ordered = TRUE)

dat.cpd.collate$group_label <- factor(dat.cpd.collate$group, 
                                      levels = c("6","12", "22", "31", "UM"),
                                      labels = c("6 yr","12 yr", "22 yr", "31 yr", "Unmined"),ordered = TRUE)



levels(dat.cpd.collate$group) # "6"  "12" "22" "31" "UM"

dat.cpd.collate$ord_group <- NA
sel <- which(dat.cpd.collate$group == "6") # qty Compound: 25110   23208  Indiv 31500 ; Fxn-MEan 12240
dat.cpd.collate$ord_group[sel] <- 1
sel <- which(dat.cpd.collate$group == "12") # qty 25110
dat.cpd.collate$ord_group[sel] <- 2
sel <- which(dat.cpd.collate$group == "22") # qty 25110
dat.cpd.collate$ord_group[sel] <- 3
sel <- which(dat.cpd.collate$group == "31") # qty 25110
dat.cpd.collate$ord_group[sel] <- 4
sel <- which(dat.cpd.collate$group == "UM") # qty 25110
dat.cpd.collate$ord_group[sel] <- 5


saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")



str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:
# $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "20C" "20C" "20C" "20C" ...
# $ cpd_rel_abun: num  0.000435 0.022738 0.000456 5.135411 0.015705 ...
# $ log10_abun  : num  -3.362 -1.643 -3.341 0.711 -1.804 ...
# $ group       : Ord.factor w/ 5 levels "6"<"12"<"22"<..: 3 3 3 3 3 3 3 3 3 3 ...
# $ group_label : Ord.factor w/ 5 levels "6 yr"<"12 yr"<..: 3 3 3 3 3 3 3 3 3 3 ...
# $ ord_group   : num  3 3 3 3 3 3 3 3 3 3 ...


length( unique(dat.cpd.collate$cpd_id) ) # 8370
8370*15 # 125550


dat.test <- data.frame(cpd = unique(dat.cpd.collate$cpd_id), data_for_this_cpd=NA , p_val = NA, kendall_tau = NA, trend_with_age = NA )


for (i in 1:dim(dat.test)[1]) {
  #i<-1
  this_cpd <- dat.test$cpd[i]
  sel <- which(dat.cpd.collate$cpd_id == this_cpd)
  # Kendall Tau correlation
  
  x = dat.cpd.collate$ord_group[sel]
  y = dat.cpd.collate$log10_abun[sel]
  df = as.data.frame(cbind(x,y))
  sel.ok = which(complete.cases(df)==TRUE)
  df = df[sel.ok, ]
  if (sd(df$y)==0 & length(unique(df$y))==1 & unique(df$y)[1]==min(df$y)) {  
    # disqualified: if they were all zeros before zero-replacement
    dat.test$data_for_this_cpd[i] <- "disqualified"
    
  } else if ( length(which(df$y > min(df$y))) < 0.25*length(sel) ) { # i.e. 3.75 ; 3 or less
    # low data: if non-zero cases do not comprise at least a replicate / a quarter of available scenario data
    dat.test$data_for_this_cpd[i] <- "low data"
    
  } else {
    ktcor<- cor.test(x = df$x, y = df$y, method = "kendall")
    dat.test$p_val[i] <- ktcor$p.value
    dat.test$kendall_tau[i] <- ktcor$statistic
    
  }
  
  
  if (!(is.na(dat.test$p_val[i])|is.na(dat.test$kendall_tau[i]))) {
    if (dat.test$p_val[i] <= 0.05 & dat.test$kendall_tau[i] > 0) { dat.test$trend_with_age[i] <- "Increasing" }
    if (dat.test$p_val[i] <= 0.05 & dat.test$kendall_tau[i] < 0) { dat.test$trend_with_age[i] <- "Decreasing" }
  }
  print(paste0("Completed ",i))
}

sel.disq <- which(dat.test$data_for_this_cpd == "disqualified") # empty
sel.low <- which(dat.test$data_for_this_cpd == "low data") # 381

length(unique(dat.cpd.collate$cpd_id)) # 8370
length(unique(dat.cpd.collate$cpd_id)) - (length(sel.disq) + length(sel.low) ) # 7989

sel.nonNA <- which(!is.na(dat.test$p_val)) # 7989 applicable tests

sel.sig <- which(dat.test$p_val <= 0.05) # 2958

# only keep applicable tests; the remainder are likely due to zero replacement:
dat.test <- dat.test[ sel.nonNA, ]


summary(dat.test$p_val)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000133 0.0113463 0.1879858 0.2947939 0.5434239 1.0000000

hist(dat.test$p_val)


# ## False Discovery Rate correction (Benjamini & Hochberg 1995)
# http://www.statisticshowto.com/benjamini-hochberg-procedure/

# what is m?

m <- dim(dat.test)[1]  # 7989 applicable compounds tested
alpha <- 0.05

p_values <- dat.test[ order(dat.test$p_val, decreasing = FALSE) ,  "p_val" ]

summary(p_values)


plot(x=1:m, y=p_values, xlab="k", ylab="P(k)")

# criteria values for Benjamini & Hochberg test
test_values <- rep(NA, times=length(p_values))

for (i in 1:m) { test_values[i] <- (i/m)*alpha }
test_values

test_results <- rep(NA, times=length(p_values))

## calculate test results
for (i in 1:m) {
  if ( p_values[i] < test_values[i] ) { test_results[i] <- "yes" }
}
test_results
# get index of largest ranked p-value with "yes" result (i.e. p-value is smaller than test criteria)

idx <- rev(which(!is.na(test_results)))[1] # 2122.  compounds: 1953 ; cf  313

dat.test$sigBH <- NA
if (length(1:idx) >0) {
  dat.test[ order(dat.test$p_val,decreasing = FALSE)[1:idx] , "sigBH" ] <- "sig"
}

length(1:idx) # 2122.  compounds 1953 ;  heatmap cell 313


p_values[idx] # 0.01311466;  0.01219892 ;  0.01984315


plot(x=1:m, y=p_values, xlab="k (index of ranked P-values)", ylab="P-value(k)" ) #, xlim=c(0,100), ylim=c(0,0.05))
# title("(a) 16S", adj=0)
abline(a=0, b=(alpha/m), col="red" )
text(x = 4000, y = 0.08, adj = 0,labels = "slope = alpha/m", col = "red")

points(x=c(1:m)[1:idx], y=p_values[1:idx], col="purple" )
text(x = 1000, y = 0.1, adj=0.5, labels = "P(k) < \n(k/m)*alpha", col = "purple")


dev.print(tiff, filename = paste0(workdir,"/plots/","3d-indiv-compound-Benjamini-Hochberg-significant-p-values-",this_study,header,".tiff"),
          width = 14, height = 14, units = "cm", res=600, compression = "lzw",type="cairo" )


# extract sig results

#sel.sig <- which(dat.test$sigBH == "sig") # 60
sel.sig <- which(dat.test$p_val <= 0.05) # 2958.  compounds: 2742  ; heatmap cells 381

dat.test.sig <- dat.test[sel.sig, ]

dat.test.sig$minuslog10_p_val <- -log10(dat.test.sig$p_val)

plot(x = dat.test.sig$kendall_tau , y =dat.test.sig$minuslog10_p_val , xlab="Kendall Tau", ylab="-log10(P-value)")

sel.sigBH <- which(dat.test.sig$sigBH == "sig")
points(x=dat.test.sig$kendall_tau[sel.sigBH], y=dat.test.sig$minuslog10_p_val[sel.sigBH], col="purple" )

dev.print(tiff, filename = paste0(workdir,"/plots/","3d-indiv-compound-VolcanoPlot-Benjamini-Hochberg-significant-p-values-",this_study,header,".tiff"),
          width = 12, height = 14, units = "cm", res=600, compression = "lzw",type="cairo" )


# join compound info??

length(which(dat.test.sig$sigBH=="sig")) # 2122


dat.test.sig$cpd_names <- NA
dat.test.sig$cpd_forms <- NA

dat.test.sig$OC_x <- NA
dat.test.sig$HC_y <- NA
dat.test.sig$NC_z <- NA

dat.test.sig$mass <- NA
dat.test.sig$class <- NA

for (i in 1:dim(dat.test.sig)[1]) {
  #i<-1
  this_cpd <- dat.test.sig$cpd[i]
  
  sel.cpd <- which(df.comp2$id == this_cpd)
  
  dat.test.sig$cpd_names[i] <- df.comp2$name[sel.cpd]
  dat.test.sig$cpd_forms[i] <- df.comp2$form[sel.cpd]
  
  dat.test.sig$OC_x[i] <- df.comp2$OC_ratio[sel.cpd]
  dat.test.sig$HC_y[i] <- df.comp2$HC_ratio[sel.cpd]
  dat.test.sig$NC_z[i] <- df.comp2$NC_ratio[sel.cpd]
  
  dat.test.sig$mass[i] <- df.comp2$mass[sel.cpd]
  dat.test.sig$class[i] <- df.comp2$class[sel.cpd]
  
  print(paste0("completed ",i))
}

#write.table(x = dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(x = dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto-v2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto.rds")
saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto-v2.rds")
# updated save later



## plot as Increasing or Decreasing?? in vK space


names(dat.test.sig)
# [1] "cpd"               "data_for_this_cpd" "p_val"             "kendall_tau"       "trend_with_age"    "sigBH"             "minuslog10_p_val" 
# [8] "cpd_names"         "cpd_forms"         "OC_x"              "HC_y"              "NC_z"


p <- ggplot(data = filter(dat.test.sig, sigBH == "sig")) +
  coord_equal()+
  ggtitle("Microbiota compound processing potential\n- Post-mining restoration case study")+
  xlim(0,3.4)+ ylim(0,4.1)+
  geom_point(aes(x = OC_x, y = HC_y, color = trend_with_age), size = 1, alpha = 0.3 ) + # 
  xlab("O:C molar ratio")+ ylab("H:C molar ratio")+
  guides(color = guide_legend(title = "Trend with age\nin functional\ncapacity (%)\nallocated to\ncompounds"))+
  
  geom_mark_rect(data= vkgrouprect, aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ 
  
  annotate(geom="text", x= 0+0.01, y= 2.3+0.02, label = "Lipid", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0.2+0.01, y= 2.2+0.02, label = "Protein", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  #annotate(geom="text", x= 0.52+0.01, y= 2.2+0.02, label = "Amino\nsugar", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0.52+0.01, y= 2.2+0.02, label = "Amino\nsugar", hjust=0, vjust=0, size = 2 , col="#737373", lineheight = 0.8) + # top-left
  annotate(geom="text", x= 0.7+0.01, y= 2.4+0.02, label = "Carbohydrate", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0+0.01, y= 0.5-0.02, label = "Condensed\naromatics", hjust=0, vjust=1, size = 2 , col="#737373", lineheight = 0.8) + # bottom-left
  annotate(geom="text", x= 0.25+0.01, y= 0.75-0.02, label = "Lignin", hjust=0, vjust=1, size = 2 , col="#737373" ) + # bottom-left
  annotate(geom="text", x= 0.67+0.01, y= 0.53-0.02, label = "Tannin", hjust=0, vjust=1, size = 2 , col="#737373") + # bottom-left
  
  theme_bw()+
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.text = element_text(size = rel(0.9)) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","Compounds-indiv-vKSpace-Trend-with-Age-",this_study,header,".tiff"), width = 14, height = 14, units = "cm", res=350, compression="lzw",type="cairo")
# Removed 147 rows containing missing values or values outside the scale range (`geom_point()`). 


hist(dat.test.sig$NC_z)


min(dat.test.sig$NC_z[ dat.test.sig$NC_z > 0 ], na.rm = TRUE) # 0.01449275


dim(dat.test.sig) # 2958   15
sel.ok <- which(!is.na(dat.test.sig$NC_z) ) # qty 2757

dat.test.sig$z_layer <- NA

subsel <- which(dat.test.sig$NC_z[sel.ok] == 0) # 1184
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C = 0"

subsel <- which(dat.test.sig$NC_z[sel.ok] > 0 & dat.test.sig$NC_z[sel.ok] <= 0.2 ) # 845
# Rivas-Ubach et al 2018 Table 1 highlights N/C breaks at 0.126, 0.2, 0.5, 0.7
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C >0 to 0.2"

subsel <- which(dat.test.sig$NC_z[sel.ok] > 0.2) # 728
max(dat.test.sig$NC_z[sel.ok]) # 2
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C >0.2 to 2"

unique(dat.test.sig$z_layer[sel.ok]) # "N:C = 0"       "N:C >0 to 0.2" "N:C >0.2 to 2"

dat.test.sig$z_layer <- factor(dat.test.sig$z_layer, levels = c("N:C = 0",
                                                                "N:C >0 to 0.2",
                                                                "N:C >0.2 to 2"), ordered = TRUE)


write.table(x = dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto-v2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto.rds")
saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto-v2.rds")
#dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto-v2.rds")

dim(dat.test.sig) #  2958   15
head(dat.test.sig)
dim(dat.test.sig[ which(dat.test.sig$sigBH == "sig"), ]) # 2122   15
sel <- which(dat.test.sig$sigBH == "sig" & !is.na(dat.test.sig$OC_x) ) # 1975
sel <- which(dat.test.sig$sigBH == "sig" & dat.test.sig$trend_with_age == "Decreasing") # 1365
sel <- which(dat.test.sig$sigBH == "sig" & dat.test.sig$trend_with_age == "Decreasing" & !is.na(dat.test.sig$OC_x) ) # 1255
sel <- which(dat.test.sig$sigBH == "sig" & dat.test.sig$trend_with_age == "Increasing") # 757
sel <- which(dat.test.sig$sigBH == "sig" & dat.test.sig$trend_with_age == "Increasing" & !is.na(dat.test.sig$OC_x) ) # 720

# build vk layer for each z-layer facet AND age group
for (i in 1:length(levels(dat.test.sig$z_layer))) {
  #i<-1
  temp <- vkgrouprect
  temp$z_layer <- levels(dat.test.sig$z_layer)[i]
  if (i==1) { keep <- temp }
  if (i>1) { keep <- rbind(keep, temp)}
  print(paste0("completed ",i))
}
vkgrouprect.facets <- keep
rm(keep)

str( vkgrouprect.facets )
vkgrouprect.facets$z_layer <- factor( vkgrouprect.facets$z_layer )




p <- ggplot(data = filter(dat.test.sig[sel.ok, ], sigBH == "sig")) +
  coord_equal()+
  ggtitle("Compound processing potential of microbiota - Post-mining restoration case study")+
  #xlim(0,3.4)+ ylim(0,4.1)+
  xlim(0,2.6)+ ylim(0,3.1)+
  geom_point(aes(x = OC_x, y = HC_y, color = trend_with_age), size = 1, alpha = 0.3 ) + # 
  xlab("O:C ratio")+ ylab("H:C ratio")+
  guides(color = guide_legend(title = "Trend with age in functional capacity\n(%) allocated to compounds"))+
  
  facet_wrap(facets = vars(z_layer))+
  
  #geom_mark_rect(data= vkgrouprect, aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ 
  
  geom_mark_rect(data= filter(vkgrouprect.facets, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets, z_layer == "N:C >0.2 to 2" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  
  annotate(geom="text", x= 0+0.01, y= 2.3+0.02, label = "Lipid", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0.2+0.01, y= 2.2+0.02, label = "Protein", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  #annotate(geom="text", x= 0.52+0.01, y= 2.2+0.02, label = "Amino\nsugar", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0.52+0.01, y= 2.2+0.02, label = "Amino\nsugar", hjust=0, vjust=0, size = 2 , col="#737373", lineheight = 0.8) + # top-left
  annotate(geom="text", x= 0.7+0.01, y= 2.4+0.02, label = "Carbohydrate", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0+0.01, y= 0.5-0.02, label = "Condensed\naromatics", hjust=0, vjust=1, size = 2 , col="#737373", lineheight = 0.8) + # bottom-left
  annotate(geom="text", x= 0.25+0.01, y= 0.75-0.02, label = "Lignin", hjust=0, vjust=1, size = 2 , col="#737373" ) + # bottom-left
  annotate(geom="text", x= 0.67+0.01, y= 0.53-0.02, label = "Tannin", hjust=0, vjust=1, size = 2 , col="#737373") + # bottom-left
  
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.text = element_text(size = rel(0.9)) ,
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    title = element_text(size = rel(0.9), lineheight = 1),
    strip.background = element_rect(fill = "transparent")
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-Trend-with-Age-",this_study,header,"-v2.tiff"), width = 20, height = 12, units = "cm", res=350, compression="lzw",type="cairo")


# Exploring possible groupings based on literature ...

subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0.2 to 2" & dat.test.sig$sigBH[sel.ok] == "sig" &
                  #dat.test.sig$OC_x[sel.ok] > 0.52 & dat.test.sig$OC_x[sel.ok] <= 0.7 &
                  dat.test.sig$OC_x[sel.ok] > 0.52 & dat.test.sig$OC_x[sel.ok] <= 1.05 &
                  dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.2 )
dat.test.sig[sel.ok[subsel] , ] # Proteins!
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Nucleotide     Protein Unspecified 
# 9         122           8
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" & dat.test.sig$sigBH[sel.ok] == "sig" &
                  dat.test.sig$OC_x[sel.ok] > 0.25 & dat.test.sig$OC_x[sel.ok] <= 1.05 &
                  #dat.test.sig$OC_x[sel.ok] > 0.52 & dat.test.sig$OC_x[sel.ok] <= 1.05 &
                  dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.2 )
#dat.test.sig[sel.ok[subsel] , ] # 
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Amino sugar       Lipid     Protein Unspecified 
# 19         116          71          94
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" & dat.test.sig$sigBH[sel.ok] == "sig" &
                  dat.test.sig$OC_x[sel.ok] > 0.35 & dat.test.sig$OC_x[sel.ok] <= 0.61 &
                  #dat.test.sig$OC_x[sel.ok] > 0.61 & dat.test.sig$OC_x[sel.ok] <= 1.05 & # 0.52  1.05
                  dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.2 )
#dat.test.sig[sel.ok[subsel] , ] # 
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Lipid     Protein Unspecified 
# 46          11          14

## Check dominant compound classes in each vK zone?
## "N:C = 0"
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C = 0" & dat.test.sig$sigBH[sel.ok] == "sig" & dat.test.sig$OC_x[sel.ok] >= 0 & dat.test.sig$OC_x[sel.ok] <= 0.6 & dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.3 )
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Lipid Unspecified - ok Lipids!
# 58         251 

subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C = 0" & dat.test.sig$sigBH[sel.ok] == "sig" & dat.test.sig$OC_x[sel.ok] >= 0.7 & dat.test.sig$OC_x[sel.ok] <= 1.25 & dat.test.sig$HC_y[sel.ok] > 1.5 & dat.test.sig$OC_x[sel.ok] <= 2.4 )
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Carbohydrate  Unspecified - ok Carbohydrates
# 92           19 

# "N:C >0 to 0.2"
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" & dat.test.sig$sigBH[sel.ok] == "sig" & dat.test.sig$OC_x[sel.ok] >= 0 & dat.test.sig$OC_x[sel.ok] <= 0.5 & dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.3 )
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Lipid     Protein Unspecified 
# 143          55          66 
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" & dat.test.sig$sigBH[sel.ok] == "sig" & dat.test.sig$OC_x[sel.ok] >= 0 & dat.test.sig$OC_x[sel.ok] <= 0.6 & dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.3 )
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Lipid     Protein Unspecified - Lipids with Protein overlap
# 147          85          74
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" & dat.test.sig$sigBH[sel.ok] == "sig" & dat.test.sig$OC_x[sel.ok] >= 0.12 & dat.test.sig$OC_x[sel.ok] <= 0.6 & dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.3 )
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Lipid     Protein Unspecified 
# 147          85          58 
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" & dat.test.sig$sigBH[sel.ok] == "sig" & dat.test.sig$OC_x[sel.ok] >= 0.12 & dat.test.sig$OC_x[sel.ok] <= 1.05 & dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.3 )
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Amino sugar       Lipid     Protein Unspecified 
#          19         147          85         121

# "N:C >0.2 to 2"
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0.2 to 2" & dat.test.sig$sigBH[sel.ok] == "sig" & dat.test.sig$OC_x[sel.ok] >= 0.3 & dat.test.sig$OC_x[sel.ok] <= 1.05 & dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.3 )
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Nucleotide     Protein Unspecified 
# 9         191          17 
subsel <- which(dat.test.sig$z_layer[sel.ok] == "N:C >0.2 to 2" & dat.test.sig$sigBH[sel.ok] == "sig" & dat.test.sig$OC_x[sel.ok] >= 0.12 & dat.test.sig$OC_x[sel.ok] <= 1.05 & dat.test.sig$HC_y[sel.ok] > 1.32 & dat.test.sig$OC_x[sel.ok] <= 2.3 )
table( dat.test.sig[sel.ok[subsel] ,"class" ] )
# Nucleotide     Protein Unspecified 
# 9         220          20



## Use adjusted compound classes (adapted from Wu 2018, D'Andrilli, Rivas-Ubach 2018, and Minor et al 2015)

vkgrouprect.facets2 <- read.table(file = "cpp3d-compound-classes.tsv", header = TRUE, sep = "\t" )
vkgrouprect.facets2.labels <- read.table(file = "cpp3d-compound-classes-labels.tsv", header = TRUE, sep = "\t" )

vkgrouprect.facets2.labels$label <- gsub(pattern = "Other p", replacement = "P", x = vkgrouprect.facets2.labels$label)




p <- ggplot(data = filter(dat.test.sig[sel.ok, ], sigBH == "sig")) +
  coord_equal()+
  
  #ggtitle("Compound processing potential of microbiota - Post-mining restoration case study")+
  
  ##xlim(0,3.4)+ ylim(0,4.1)+
  xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_age), size = 1, alpha = 0.3 ) + # 
  
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_age, shape = class), size = 1, alpha = 0.3 ) + # 
  #scale_shape_manual(values = shapes.class, name = "Compound type")+
  
  xlab("O:C ratio")+ ylab("H:C ratio")+
  guides(color = guide_legend(title = "Trend with age in functional capacity\n(%) allocated to compounds"))+
  
  facet_wrap(facets = vars(z_layer))+
  
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, group = label), color="#969696", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373", "grey"
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, group = label), color="#969696", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0.2 to 2" ), aes(x = OC_x, y = HC_y, group = label), color="#969696", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_with_age), size = 1, alpha = 0.3 ) + # 
  
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T,  size = 2.75 , col="#252525" , lineheight = 0.8)+ # "#737373"
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2.75 , col="#252525" , lineheight = 0.8)+
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0.2 to 2" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2.75 , col="#252525" , lineheight = 0.8)+
  
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.text = element_text(size = rel(0.9)) ,
    
    legend.box.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt"),
    title = element_text(size = rel(0.9), lineheight = 1),
    strip.background = element_rect(fill = "transparent")
  )

p
#dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-Trend-with-Age-",this_study,header,"-v2c.tiff"), width = 20, height = 12, units = "cm", res=350, compression="lzw",type="cairo")
#dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-Trend-with-Age-Sunbad-resto-NoGGtitle-v2c.tiff"), width = 20, height = 12, units = "cm", res=600, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-Trend-with-Age-Sunbad-resto-NoGGtitle-v4.tiff"), width = 20, height = 12, units = "cm", res=600, compression="lzw",type="cairo")




### also visualise all compounds for a single sample

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:

# isolate a single sample
df1 <- dat.cpd.collate
unique(df1$sample)
#  "20C" "30B" "30A" "UMA" "10B" "5A"  "10C" "20A" "20B" "5B"  "UMB" "10A" "5C"  "UMC" "30C"
sel <- which(df1$sample == "UMA")
df1 <- df1[sel, ]
names(df1)
# [1] "cpd_id"       "sample"       "cpd_rel_abun" "log10_abun"   "group"        "group_label"  "ord_group"    "cpd_names"    "cpd_forms"   
# [10] "OC_x"         "HC_y"         "NC_z" 

df1$cpd_names <- NA
df1$cpd_forms <- NA
df1$OC_x <- NA
df1$HC_y <- NA
df1$NC_z <- NA

for (i in 1:dim(df1)[1]) {
  #i<-1
  this_cpd <- df1$cpd_id[i]
  
  sel.cpd <- which(df.comp2$id == this_cpd)
  
  df1$cpd_names[i] <- df.comp2$name[sel.cpd]
  df1$cpd_forms[i] <- df.comp2$form[sel.cpd]
  
  df1$OC_x[i] <- df.comp2$OC_ratio[sel.cpd]
  df1$HC_y[i] <- df.comp2$HC_ratio[sel.cpd]
  df1$NC_z[i] <- df.comp2$NC_ratio[sel.cpd]
  
  print(paste0("completed ",i))
}

hist(df1$NC_z)


min(df1$NC_z[ df1$NC_z > 0 ], na.rm = TRUE) # 0.004201681

dim(df1) # 8370 12
dim(df1[df1$cpd_rel_abun > 0, ]) # 7826   13

sel.ok <- which(!is.na(df1$NC_z) ) # qty 7736

df1$z_layer <- NA

subsel <- which(df1$NC_z[sel.ok] == 0) # 3661
df1$z_layer[sel.ok [subsel]] <- "N:C = 0"

subsel <- which(df1$NC_z[sel.ok] > 0 & df1$NC_z[sel.ok] <= 0.2 ) # 2290
# Rivas-Ubach et al 2018 Table 1 highlights N/C breaks at 0.126, 0.2, 0.5, 0.7
df1$z_layer[sel.ok [subsel]] <- "N:C >0 to 0.2"

subsel <- which(df1$NC_z[sel.ok] > 0.2) # 1785
max(df1$NC_z[sel.ok]) # 3
#df1$z_layer[sel.ok [subsel]] <- "N:C >0.2 to 2"
df1$z_layer[sel.ok [subsel]] <- "N:C >0.2 to 3"

unique(df1$z_layer[sel.ok]) # "N:C = 0"       "N:C >0 to 0.2" "N:C >0.2 to 3"

df1$z_layer <- factor(df1$z_layer, levels = c("N:C = 0",
                                                                "N:C >0 to 0.2",
                                                                "N:C >0.2 to 3"), ordered = TRUE)


# compound class guides
vkgrouprect.facets3 <- vkgrouprect.facets2
vkgrouprect.facets3.labels <- vkgrouprect.facets2.labels

sel <- which(vkgrouprect.facets3$z_layer == "N:C >0.2 to 2")
vkgrouprect.facets3$z_layer[sel] <- "N:C >0.2 to 3"

sel <- which(vkgrouprect.facets3.labels$z_layer == "N:C >0.2 to 2")
vkgrouprect.facets3.labels$z_layer[sel] <- "N:C >0.2 to 3"

vkgrouprect.facets3.labels$label <- gsub(pattern = "Other p", replacement = "P", x = vkgrouprect.facets3.labels$label)

unique(vkgrouprect.facets3.labels$z_layer)

p <- ggplot(data = df1[sel.ok, ]) +
  coord_equal()+
  
  #ggtitle("Compound processing potential of microbiota - Post-mining restoration case study")+
  
  #xlim(0,2.6)+ ylim(0,3.1)+
  
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_age), size = 1, alpha = 0.3 ) + # 
  
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_age, shape = class), size = 1, alpha = 0.3 ) + # 
  #scale_shape_manual(values = shapes.class, name = "Compound type")+
  
  xlab("O:C ratio")+ ylab("H:C ratio")+
  
  #guides(color = guide_legend(title = "Log10 functional capacity (%)\nallocated to compounds"))+
  
  scale_color_continuous(type = "viridis", direction = -1) + 
  #guides(color = guide_colorbar(direction = "horizontal", barheight = 0.65, barwidth = 3.5, vjust = 0, hjust = 0.5, title = "Log10 functional capacity (%)\nallocated to compounds"))+
  guides(color = guide_colorbar(direction = "horizontal", barheight = 0.65, barwidth = 3.5, vjust = 0, hjust = 0.5, title = "Log10 functional rel abun (%)\nallocated to compounds"))+   
  
  facet_wrap(facets = vars(z_layer))+
  
  geom_mark_rect(data= filter(vkgrouprect.facets3, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, group = label), color= "#969696", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373", "grey"
  geom_mark_rect(data= filter(vkgrouprect.facets3, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, group = label), color="#969696", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets3, z_layer == "N:C >0.2 to 3" ), aes(x = OC_x, y = HC_y, group = label), color="#969696", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  
  geom_point(aes(x = OC_x, y = HC_y, color = log10_abun), size = 1, alpha = 0.3 ) + # 
  
  # geom_text(data = filter(vkgrouprect.facets3.labels, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T,  size = 2.75 , col="#252525" , lineheight = 0.8)+ # "#737373"
  # geom_text(data = filter(vkgrouprect.facets3.labels, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2.75 , col="#252525" , lineheight = 0.8)+
  # geom_text(data = filter(vkgrouprect.facets3.labels, z_layer == "N:C >0.2 to 3" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2.75 , col="#252525" , lineheight = 0.8)+ # "#252525"
  # 
  theme_bw()+
  theme(
    #legend.position = "right",
    
    legend.position = c(0.5, 0.05), # as fraction of plot dimensions
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(0.7)),
    
    #legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #legend.background = element_rect(fill = "transparent"),
    
    #legend.title = element_text(size = rel(1)),
    legend.text = element_text(size = rel(0.6)) ,
    
    legend.box.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt"),
    #title = element_text(size = rel(0.9), lineheight = 1),
    strip.background = element_rect(fill = "transparent")
  )

p
#dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-Trend-with-Age-",this_study,header,"-v2c.tiff"), width = 20, height = 12, units = "cm", res=350, compression="lzw",type="cairo")

#grid.text(label = "(a)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
grid.text(label = "(a)", x = unit(0.02, "npc") , y = unit(0.85,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-1Sample-UMA-Sunbad-resto-NoGGtitle-v2c.tiff"), width = 20, height = 12, units = "cm", res=450, compression="lzw",type="cairo")

grid.text(label = "(b)", x = unit(0.02, "npc") , y = unit(0.9,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-1Sample-UMA-Sunbad-resto-NoGGtitle-v2c-Fullrange.tiff"), width = 20, height = 12, units = "cm", res=450, compression="lzw",type="cairo")


#-------------------------

#### Sun & Badgley - post-mining - CPP as phyloseq object
#    PCoA using CPP vs Functions
#-------------------------

dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:

data_in <- dat.cpd.collate

length( unique(data_in$cpd_id) ) # 8370
8370*15 # 125550


### get data into phyloseq object ...

head(data_in)
#     cpd_id sample cpd_rel_abun log10_abun group group_label ord_group
# 1 cpd25681    20C 0.0004345842 -3.3619261    22       22 yr         3
# 2 cpd02597    20C 0.0227382574 -1.6432428    22       22 yr         3
# 3 cpd24620    20C 0.0004564208 -3.3406346    22       22 yr         3
# 4 cpd00001    20C 5.1354106008  0.7105752    22       22 yr         3
# 5 cpd01501    20C 0.0157048650 -1.8039658    22       22 yr         3
# 6 cpd00851    20C 0.0121620567 -1.9149930    22       22 yr         3

df.wide <- dcast(data_in, formula = sample + group ~ cpd_id , value.var = "cpd_rel_abun" )

df.wide[1:5, 1:10]
#   sample group cpd00001 cpd00002 cpd00003 cpd00004  cpd00005  cpd00006  cpd00007 cpd00008
# 1    10A    12 5.126634 2.287831 1.058127 1.031673 0.7565038 0.7592450 0.3678576 1.393837
# 2    10B    12 5.109487 2.333880 1.058379 1.032674 0.7348372 0.7370765 0.3472805 1.392368
# 3    10C    12 5.103210 2.319224 1.070352 1.044885 0.7498325 0.7520871 0.3398280 1.399435
# 4    20A    22 5.138510 2.323322 1.069859 1.043509 0.7560155 0.7585376 0.3598548 1.399336
# 5    20B    22 5.173074 2.355075 1.064497 1.037668 0.7483448 0.7511583 0.3623488 1.419401

unique(paste0(df.wide$sample,"--",df.wide$group))
# [1] "10A--12" "10B--12" "10C--12" "20A--22" "20B--22" "20C--22" "30A--31" "30B--31" "30C--31" "5A--6"  
# [11] "5B--6"   "5C--6"   "UMA--UM" "UMB--UM" "UMC--UM"

# save group variable
samp <- df.wide[ ,1:2]
row.names(samp) <- samp$sample

# transpose
df.wide <- t(df.wide[ ,-2]) # minus 'group' column

head(df.wide)
# [,1]        [,2]        [,3]        [,4]        [,5]        [,6]        [,7]       
# sample   "10A"       "10B"       "10C"       "20A"       "20B"       "20C"       "30A"      
# cpd00001 "5.126634"  "5.109487"  "5.103210"  "5.138510"  "5.173074"  "5.135411"  "5.036718" 
# cpd00002 "2.287831"  "2.333880"  "2.319224"  "2.323322"  "2.355075"  "2.318166"  "2.309780" 
# cpd00003 "1.058127"  "1.058379"  "1.070352"  "1.069859"  "1.064497"  "1.060200"  "1.067475" 
# cpd00004 "1.031673"  "1.032674"  "1.044885"  "1.043509"  "1.037668"  "1.033088"  "1.041951" 
# cpd00005 "0.7565038" "0.7348372" "0.7498325" "0.7560155" "0.7483448" "0.7564190" "0.7424559"
# [,8]        [,9]        [,10]       [,11]       [,12]       [,13]       [,14]      
# sample   "30B"       "30C"       "5A"        "5B"        "5C"        "UMA"       "UMB"      
# cpd00001 "4.977622"  "5.054107"  "5.136126"  "5.128531"  "5.155902"  "4.988787"  "5.045806" 
# cpd00002 "2.272223"  "2.332804"  "2.335298"  "2.345242"  "2.334887"  "2.298596"  "2.329824" 
# cpd00003 "1.093609"  "1.077584"  "1.054323"  "1.056899"  "1.075898"  "1.098387"  "1.096881" 
# cpd00004 "1.068711"  "1.051745"  "1.026151"  "1.028856"  "1.047329"  "1.072772"  "1.071776" 
# cpd00005 "0.7597370" "0.7483150" "0.7482304" "0.7517447" "0.7519912" "0.7931363" "0.7843985"
# [,15]      
# sample   "UMC"      
# cpd00001 "5.064984" 
# cpd00002 "2.339306" 
# cpd00003 "1.089745" 
# cpd00004 "1.065067" 
# cpd00005 "0.7746962"

samp_names <- df.wide[1, ]
tax_names <- row.names(df.wide[-1, ])
head(tax_names) # "cpd00001" "cpd00002" "cpd00003" "cpd00004" "cpd00005" "cpd00006"
otu.df <- df.wide[-1, ] # remove sample labels in 1st row
# this is necessary to create numeric matrix

colnames(otu.df) <- samp_names

# convert OTU table to matrix
class(otu.df) # "matrix" "array" 
#otu.df <- as.matrix(otu.df)

# convert to numeric matrix
# https://stackoverflow.com/questions/20791877/convert-character-matrix-into-numeric-matrix
otu.df <- apply(otu.df, 2, as.numeric)

rownames(otu.df) # NULL
dim(otu.df) #  8370   15
rownames(otu.df) <- tax_names

## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
OTU <- otu_table(otu.df, taxa_are_rows = TRUE)

# Error in validObject(.Object) : invalid class “otu_table” object: 
#   Non-numeric matrix provided as OTU table.
# Abundance is expected to be numeric.



# # convert Taxonomy table to matrix  

tax <- data.frame(cpd_id = tax_names)
row.names(tax) <- tax_names

tax <- as.matrix(tax)

identical( row.names(otu.df), row.names(tax) ) # TRUE


## Create 'taxonomyTable'
#  tax_table - Works on any character matrix.
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
TAX <- tax_table(tax)


## Create a phyloseq object, merging OTU & TAX tables
phy.cpp = phyloseq(OTU, TAX)
phy.cpp
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8370 taxa and 15 samples ]
# tax_table()   Taxonomy Table:    [ 8370 taxa by 1 taxonomic ranks ]


sample_names(phy.cpp)
# "10A" "10B" "10C" "20A" "20B" "20C" "30A" "30B" "30C" "5A"  "5B"  "5C"  "UMA" "UMB" "UMC"

identical(sample_names(phy.cpp), samp$sample) # TRUE


# row.names need to match sample_names() from phyloseq object
row.names(samp) <- samp$sample




### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

SAMP <- sample_data(samp)


### Combine SAMPDATA into phyloseq object
phy.cpp <- merge_phyloseq(phy.cpp, SAMP)
phy.cpp
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8370 taxa and 15 samples ]
# sample_data() Sample Data:       [ 15 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 8370 taxa by 1 taxonomic ranks ]

phy.cpp@sam_data
# Sample Data:        [15 samples by 2 sample variables]:
#   sample group
# 10A    10A    12
# 10B    10B    12
# 10C    10C    12
# 20A    20A    22
# 20B    20B    22
# 20C    20C    22
# 30A    30A    31
# 30B    30B    31
# 30C    30C    31
# 5A      5A     6
# 5B      5B     6
# 5C      5C     6
# UMA    UMA    UM
# UMB    UMB    UM
# UMC    UMC    UM


saveRDS(object = phy.cpp, file = "phy.cpp-sunbad-resto.RDS")


phy_in <- phy.cpp

phy_in
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8370 taxa and 15 samples ]
# sample_data() Sample Data:       [ 15 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 8370 taxa by 1 taxonomic ranks ]

min(taxa_sums(phy_in)) # 8.136846e-08
sum(sample_sums(phy_in)) # 1005.557
sample_sums(phy_in)
# 10A      10B      10C      20A      20B      20C      30A      30B      30C       5A       5B 
# 66.37627 67.28058 67.43371 67.14898 67.44508 67.38783 67.00958 66.68865 67.16063 67.19161 67.09898 
# 5C      UMA      UMB      UMC 
# 67.64787 66.59001 66.44038 66.65667

summary( sample_sums(phy_in) )
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 66.38   66.67   67.15   67.04   67.33   67.65 

sd( sample_sums(phy_in) )
# 0.3958827

max(taxa_sums(phy_in)) # 120.7394



# don't rarefy - already in form of relative abundance %


table(phy_in@sam_data$group)
# 6 12 22 31 UM 
# 3  3  3  3  3



## ordination plot
## PCoA + Bray-Curtis

set.seed(123)
#ord <- ordinate(r1.ps, "NMDS", "bray")
ord <- ordinate(phy_in, "PCoA", "bray")


ord

unique(phy_in@sam_data$group)
# [1] 12 22 31 6  UM
# Levels: 6 < 12 < 22 < 31 < UM


p <- plot_ordination(phy_in, ord, type="samples", color="group")
p

p$labels$x # "Axis.1   [61.5%]"
x_lab <- "PCo1 (61.5%)"

p$labels$y # "Axis.2   [17.9%]"
y_lab <- "PCo2 (17.9%)"

#temp <- r1.ps
p_df <- p$data


cols.group <- c("6" = "#9e0142",
                "12" = "#d53e4f",
                "22" = "#fdae61",
                "31" = "#abdda4",
                "UM" = "#3288bd")



p <- #plot_ordination(temp, ord, type="samples", color="group") +
  #ggplot(data = p_df, aes(x = NMDS1, y = NMDS2, color = group))+
  ggplot(data = p_df, aes(x = Axis.1, y = Axis.2, color = group))+
  theme_bw()+
  geom_point()+
  xlab(x_lab) + ylab(y_lab)+
  scale_color_manual(values = cols.group, name = "Reveg\nage (yr)") +
  #annotate(geom="text", x= -1.2, y= 1.3, label = paste0("Stress = ",round(ord$stress,digits=4)),size = 3, hjust=0, vjust=1) +
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = "Restoration: CPP\n\n", size = 3 )+
  theme(
    panel.grid.major = element_blank(),
    
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.85)),
    legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = -8, unit = "pt"), # randomForest::margin() interferes !!
    
    panel.grid.minor = element_blank())
p

grid.text(label = "(b)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )

dev.print(tiff, file = paste0(workdir,"/plots/","PCoA-Cpp3d-sunbad-resto-v2.tiff"), width = 9, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## PERMANOVA

# Calculate bray curtis distance matrix
set.seed(123)
bray <- phyloseq::distance(phy_in, method = "bray")
sampledf <- data.frame(sample_data(phy_in))
str(sampledf)

names(phy_in@sam_data)
# "sample" "group" 

# Adonis test
set.seed(123)
adonis2(bray ~ group , data = sampledf)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ group, data = sampledf)
#          Df  SumOfSqs      R2      F Pr(>F)    
# group     4 0.0044641 0.75881 7.8651  0.001 ***
# Residual 10 0.0014190 0.24119                  
# Total    14 0.0058831 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



beta <- betadisper(bray, sampledf$group)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#           Df     Sum Sq    Mean Sq      F N.Perm Pr(>F)
# Groups     4 4.4289e-05 1.1072e-05 0.3636    999  0.842
# Residuals 10 3.0451e-04 3.0452e-05




## Functions

phy_in <- readRDS("phy-phyloseq-object-sunbad-resto.RDS")

phy_in
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 30125 taxa and 15 samples ]
# sample_data() Sample Data:       [ 15 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 30125 taxa by 4 taxonomic ranks ]

min(taxa_sums(phy_in)) # 2.521458e-06
sum(sample_sums(phy_in)) # 1500
sample_sums(phy_in)
# 20C 30B 30A UMA 10B  5A 10C 20A 20B  5B UMB 10A  5C UMC 30C 
# 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100

summary( sample_sums(phy_in) )

max(taxa_sums(phy_in)) # 9.941053




# do not rarefy as already normalized using rel abun (%)


table(phy_in@sam_data$age)
# 12 22 31  6 UM 
# 3  3  3  3  3 

phy_in@sam_data$age <- factor(phy_in@sam_data$age,
                              levels = c("6", "12", "22", "31", "UM"),
                              ordered = TRUE)



## ordination plot
## PCoA + Bray-Curtis

set.seed(123)
ord <- ordinate(phy_in, "PCoA", "bray")

ord

unique(phy_in@sam_data$age)
# [1] 12 22 31 6  UM
# Levels: 6 < 12 < 22 < 31 < UM


p <- plot_ordination(phy_in, ord, type="samples", color="age")
p

p$labels$x # "Axis.1   [49.6%]"
x_lab <- "PCo1 (49.6%)"

p$labels$y # "Axis.2   [18.8%]"
y_lab <- "PCo2 (18.8%)"

#temp <- r1.ps
p_df <- p$data


cols.group <- c("6" = "#9e0142",
                "12" = "#d53e4f",
                "22" = "#fdae61",
                "31" = "#abdda4",
                "UM" = "#3288bd")



p <- #plot_ordination(temp, ord, type="samples", color="group") +
  #ggplot(data = p_df, aes(x = NMDS1, y = NMDS2, color = group))+
  ggplot(data = p_df, aes(x = Axis.1, y = Axis.2, color = age))+ # NOTE change from group to age
  theme_bw()+
  geom_point()+
  xlab(x_lab) + ylab(y_lab)+
  scale_color_manual(values = cols.group, name = "Reveg\nage (yr)") +
  #annotate(geom="text", x= -1.2, y= 1.3, label = paste0("Stress = ",round(ord$stress,digits=4)),size = 3, hjust=0, vjust=1) +
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = "Restoration: Fxns", size = 3 )+
  theme(
    panel.grid.major = element_blank(),
    
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.85)),
    legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = -8, unit = "pt"), # randomForest::margin() interferes !!
    
    panel.grid.minor = element_blank())
p

grid.text(label = "(a)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )

dev.print(tiff, file = paste0(workdir,"/plots/","PCoA-Functions-sunbad-resto-v2.tiff"), width = 9, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## PERMANOVA

# Calculate bray curtis distance matrix
set.seed(123)
bray <- phyloseq::distance(phy_in, method = "bray")
sampledf <- data.frame(sample_data(phy_in))
str(sampledf)

names(phy_in@sam_data)
# [1] "sample_name"        "mgrast_id"          "metagenome_id"      "metagenome_name"    "investigation_type" "seq_meth"          
# [7] "file_name"          "age" 

# Adonis test
set.seed(123)
adonis2(bray ~ age , data = sampledf)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ age, data = sampledf)
#          Df SumOfSqs      R2      F Pr(>F)   
# age       4 0.033804 0.67039 5.0846  0.002 **
# Residual 10 0.016621 0.32961                 
# Total    14 0.050425 1.00000                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



beta <- betadisper(bray, sampledf$age)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#           Df     Sum Sq    Mean Sq      F N.Perm Pr(>F)
# Groups     4 0.00032954 8.2385e-05 0.3447    999  0.838
# Residuals 10 0.00238990 2.3899e-04 



#-------------------------


##########################
##########################
##########################
##########################


#### T2D Swedish (SWE) cohort - T2D Met- and Normal only - with human genome host-removal

#### Forslund T2D-SWE - identify sequence datasets
#-------------------------

## Prepare T2D - Swedish cohort dataset

# List of samples and diagnosis (Forslund et al SItableS1)
# "Samples-and-diagnosis-SItableS1.xlsx"
subjects <- read_excel(path= "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/data/Forslund-2015-Type2-diabetes/Samples-and-diagnosis-SItableS1.xlsx",
                       sheet=1, range="A1:C785", col_names = TRUE)
subjects <- as.data.frame(subjects)
str(subjects)
sel <- which(subjects$`Country subset` == "SWE")
subjects <- subjects[sel, ]
str(subjects)
# 'data.frame':	145 obs. of  3 variables:
# $ Sample        : chr  "NG-5425_137" "NG-5425_249" "NG-5425_342" "NG-5425_485" ...
# $ Country subset: chr  "SWE" "SWE" "SWE" "SWE" ...
# $ Status        : chr  "ND CTRL" "ND CTRL" "ND CTRL" "ND CTRL" ...
length(unique( subjects$Sample )) # 145
table(subjects$Status) # , useNA = TRUE
# ND CTRL T2D metformin- T2D metformin+ 
#   92             33             20 
92 + 33 + 20 # 145

run_details <- read_excel(path= "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/data/Forslund-2015-Type2-diabetes/Samples-and-diagnosis-SItableS1.xlsx",
                          sheet="Sequencing statistics SWE sampl", range="A2:K149", col_names = TRUE)
run_details <- as.data.frame(run_details)
str(run_details)
# 'data.frame':	147 obs. of  11 variables:
# $ Sample           : chr  "Sum" "Average" "NG-5425_137" "NG-5425_249" ...
# $ Raw bases        : num  4.53e+11 3.13e+09 7.98e+08 9.67e+08 1.06e+09 ...
# $ Raw reads        : num  4.49e+09 3.10e+07 8.31e+06 1.01e+07 1.11e+07 ...
# $ Trimmed bases    : num  2.78e+11 1.92e+09 3.87e+08 4.74e+08 5.31e+08 ...
# $ Trimmed reads    : num  3.31e+09 2.28e+07 5.04e+06 6.14e+06 6.78e+06 ...
# $ Trimmed inserts  : num  1.89e+09 1.30e+07 3.03e+06 3.68e+06 4.06e+06 ...
# $ Unassembled bases: num  2.78e+11 1.92e+09 3.87e+08 4.74e+08 5.31e+08 ...
# $ Unassembled reads: num  3.31e+09 2.28e+07 5.04e+06 6.14e+06 6.78e+06 ...
# $ Inserts          : num  1.89e+09 1.30e+07 3.03e+06 3.68e+06 4.06e+06 ...
# $ Max read length  : chr  "-" "97.855172413793099" "94" "92" ...
# $ Avg read length  : chr  "-" "83.303448275862067" "77" "77" ...

head(run_details)
# Sample    Raw bases  Raw reads Trimmed bases Trimmed reads Trimmed inserts Unassembled bases
# 1         Sum 453465327668 4494103168  277759781846    3312028283      1889248942      277759781846
# 2     Average   3127347087   30993815    1915584702      22841574        13029303        1915584702
# 3 NG-5425_137    797799936    8310416     387280402       5036087         3027036         387280402
# 4 NG-5425_249    966913152   10072012     474310006       6141684         3679946         474310006
# 5 NG-5425_286   1061488704   11057174     530656068       6778413         4062863         530656068
# 6 NG-5425_297    761805504    7935474     379222624       4895014         2914433         379222624
# Unassembled reads    Inserts    Max read length    Avg read length
# 1        3312028283 1889248942                  -                  -
#   2          22841574   13029303 97.855172413793099 83.303448275862067
# 3           5036087    3027036                 94                 77
# 4           6141684    3679946                 92                 77
# 5           6778413    4062863                 96                 78
# 6           4895014    2914433                 93                 77

# delete first 2 rows:
run_details <- run_details[-c(1,2), ]
head(run_details)
#        Sample  Raw bases Raw reads Trimmed bases Trimmed reads Trimmed inserts Unassembled bases Unassembled reads Inserts Max read length Avg read length
# 3 NG-5425_137  797799936   8310416     387280402       5036087         3027036         387280402           5036087 3027036              94              77
# 4 NG-5425_249  966913152  10072012     474310006       6141684         3679946         474310006           6141684 3679946              92              77
# 5 NG-5425_286 1061488704  11057174     530656068       6778413         4062863         530656068           6778413 4062863              96              78
# 6 NG-5425_297  761805504   7935474     379222624       4895014         2914433         379222624           4895014 2914433              93              77
# 7 NG-5425_342  789010752   8218862     378848377       4968844         2988398         378848377           4968844 2988398              92              76
# 8 NG-5425_470  820824384   8550254     402254631       5219247         3151422         402254631           5219247 3151422              93              77

length(unique(run_details$`Raw bases`)) # 145


# SRA Runs Metadata: 
# "Metadata-Forslund-ERP002469-SWE-samples--SraRunTable"
meta <- read_excel(path= "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/data/Forslund-2015-Type2-diabetes/Metadata-Forslund-ERP002469-SWE-samples--SraRunTable.xlsx",
                   sheet=1, range="A1:AI146", col_names = TRUE)
meta <- as.data.frame(meta)
str(meta)
# 'data.frame':	145 obs. of  35 variables:
names(meta)
# [1] "Run"                    "Assay Type"             "AvgSpotLen"             "Bases"                  "BioProject"             "BioSample"             
# [7] "Bytes"                  "Center Name"            "Consent"                "DATASTORE filetype"     "DATASTORE provider"     "DATASTORE region"      
# [13] "ENA-FIRST-PUBLIC (run)" "ENA-FIRST-PUBLIC"       "ENA-LAST-UPDATE (run)"  "ENA-LAST-UPDATE"        "Experiment"             "External_Id"           
# [19] "INSDC_center_alias"     "INSDC_center_name"      "INSDC_first_public"     "INSDC_last_update"      "INSDC_status"           "Instrument"            
# [25] "Library Name"           "LibraryLayout"          "LibrarySelection"       "LibrarySource"          "Organism"               "Platform"              
# [31] "ReleaseDate"            "Sample Name"            "Sample_Name"            "SRA Study"              "Submitter_Id" 
length(unique(meta$Bases)) # 145
length(meta$Bases %in% run_details$`Raw bases`) # 145 
# i.e. meta$Run - has meta$Bases 
# >> match to run_details$`Raw bases` and run_details$Sample
# link to subjects$Sample
dim(subjects) # 145   3
dim(meta)     # 145  35
dim(run_details) # 145  11

subjects$Bases <- NA
subjects$Run <- NA

for (i in 1:length(subjects$Sample)) {
  #i<-1
  this_samp <- subjects$Sample[i]
  sel1 <- which(run_details$Sample == this_samp)
  subjects$Bases[i] <- run_details$`Raw bases`[sel1]
  
  sel2 <- which(meta$Bases == run_details$`Raw bases`[sel1])
  
  subjects$Run[i] <- meta$Run[sel2]
  
  print(paste0("completed ",i))
}

length(meta$Run) # 145

head(subjects)
#          Sample Country subset  Status      Bases       Run
# 640 NG-5425_137            SWE ND CTRL  797799936 ERR260268
# 641 NG-5425_249            SWE ND CTRL  966913152 ERR260269
# 642 NG-5425_342            SWE ND CTRL  789010752 ERR260272
# 643 NG-5425_485            SWE ND CTRL 1090824000 ERR260275
# 644 NG-5636_112            SWE ND CTRL 1954234860 ERR260249
# 645 NG-5636_118            SWE ND CTRL 3307295096 ERR260250
tail(subjects)
#         Sample Country subset         Status      Bases       Run
# 779 NG-5636_530            SWE T2D metformin+ 3860516334 ERR260192
# 780 NG-5636_536            SWE T2D metformin+ 4292398798 ERR260194
# 781 NG-5636_543            SWE T2D metformin+ 4253796800 ERR260196
# 782 NG-5636_563            SWE T2D metformin+ 3464254550 ERR260202
# 783 NG-5636_591            SWE T2D metformin+ 3833173412 ERR260208
# 784 NG-5636_604            SWE T2D metformin+ 3923443172 ERR260214
names(subjects)
# "Sample"         "Country subset" "Status"         "Bases"          "Run"  
getwd() # "/Users/lidd0026/WORKSPACE/PROJ/PCaN-NZ/nz-city-resto/modelling/R"
write.csv(x = subjects, file = "Forslund-T2D-SWE-subjects-status-sra-run-list.txt", quote = FALSE, row.names = FALSE)

str(subjects)
# 'data.frame':	145 obs. of  5 variables:
# $ Sample        : chr  "NG-5425_137" "NG-5425_249" "NG-5425_342" "NG-5425_485" ...
# $ Country subset: chr  "SWE" "SWE" "SWE" "SWE" ...
# $ Status        : chr  "ND CTRL" "ND CTRL" "ND CTRL" "ND CTRL" ...
# $ Bases         : num  7.98e+08 9.67e+08 7.89e+08 1.09e+09 1.95e+09 ...
# $ Run           : chr  "ERR260268" "ERR260269" "ERR260272" "ERR260275" ...

table(subjects$Status)
# ND CTRL T2D metformin- T2D metformin+ 
#      92             33             20 
length(unique(subjects$Run)) # 145
sampid <- sort( subjects$Run )


# Sequences for these runs were downloaded from SRA, quality trimmed, sequences mapping to the human genome (GRCh38.p14 or hg38) were removed,
# and then R1 sequences were run through SUPER-FOCUS functional profiling software



# but additional data from Karlsson et al 2013 includes impaired glucose tolerance (IGT; n = 49) or normal glucose tolerance (NGT; n = 43)

# this lookup has IGT / NGT classification
t2d.class <- read_excel(path= "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/data/Forslund-2015-Type2-diabetes/Karlsson et al 2013--41586_2013_BFnature12198_MOESM507_ESM.xlsx",
                        sheet="Supplementary Table 3", range="A2:AD147", col_names = TRUE)
t2d.class <- as.data.frame(t2d.class)
str(t2d.class)
# 'data.frame':	145 obs. of  30 variables:
# $ Sample ID                                                                             : num  51 53 54 58 59 60 77 80 88 92 ...
# $ Age (years)                                                                           : num  69.1 70.3 69.9 70.2 69.4 ...
# $ Classification                                                                        : chr  "IGT" "NGT" "IGT" "NGT" ...
# etc ...

# this lookup has no of reads / bases
t2d.nbases <- read_excel(path= "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/data/Forslund-2015-Type2-diabetes/Karlsson et al 2013--41586_2013_BFnature12198_MOESM507_ESM.xlsx",
                         sheet="Supplementary Table 4", range="A2:C147", col_names = TRUE)
t2d.nbases <- as.data.frame(t2d.nbases)
str(t2d.nbases)
# 'data.frame':	145 obs. of  3 variables:
# $ Sample ID                 : num  51 53 54 58 59 60 77 80 88 92 ...
# $ Number of reads           : num  10653221 13971868 13561331 13213300 12460015 ...
# $ Total number of bases (bp): num  2.15e+09 2.82e+09 2.74e+09 2.67e+09 2.52e+09 ...

identical(t2d.class$`Sample ID`,t2d.nbases$`Sample ID`) # TRUE

length(unique(meta$Bases)) # 145

df.samp <- subjects
unique_samps <- unique(df.samp$Run) # qty 145

temp <- df.samp

df.samp$group_label <- NA
df.samp$age <- NA

for (i in 1:dim(df.samp)[1]) {
  #i<-1
  this_run <- unique_samps[i]
  sel.bases.row <- which(t2d.nbases$`Total number of bases (bp)` == df.samp$Bases[i])
  df.samp$group_label[i] <- paste0(df.samp$Status[i],"__", t2d.class$Classification[sel.bases.row])
  df.samp$age[i] <- t2d.class$`Age (years)`[sel.bases.row]
  print(paste0("completed ", i))
}

unique(df.samp$group_label) # "ND CTRL__NGT"        "ND CTRL__IGT"        "T2D metformin-__T2D" "T2D metformin+__T2D"

table(subjects$Status)
# ND CTRL T2D metformin- T2D metformin+ 
#   92             33             20 
table(df.samp$group_label)
# ND CTRL__IGT        ND CTRL__NGT T2D metformin-__T2D T2D metformin+__T2D 
# 49                  43                  33                  20 

# IGT = impaired glucose tolerance, NGT = normal glucose tolerance

# keep only NGT and untreated T2D (T2D metformin-)
sel <- which(df.samp$group_label %in% c("ND CTRL__NGT", "T2D metformin-__T2D") ) # qty 76
df.samp <- df.samp[sel, ]
df.samp$group_label <- factor(df.samp$group_label, 
                              levels = c("T2D metformin-__T2D" ,"ND CTRL__NGT"),
                              labels = c("T2D met-","Normal"),
                              ordered = TRUE)

saveRDS(object = df.samp, file = "df.samp.with-t2dclass-age-Forslund-SWE-T2D-qty76-v6.RDS")

keep_samps <- df.samp$Run # qty 76

#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - compare no reads w/ vs w/out host-removal
#-------------------------

r1.clean <- read.table(file = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe/4_compare_read_counts/clean_r1_sequence_counts_swe.txt", header = FALSE, sep = " ")
head(r1.clean)
names(r1.clean) <- c("sample","n_reads_cleaned")
r1.nonhost <- read.table(file = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe/4_compare_read_counts/nonhost_r1_sequence_counts_swe.txt", header = FALSE, sep = " ")
head(r1.nonhost)
names(r1.nonhost) <- c("sample","n_reads_nonhost")
identical(r1.clean$sample, r1.nonhost$sample) # TRUE
r1.swe <- r1.clean
r1.swe$n_reads_nonhost <- r1.nonhost$n_reads_nonhost
dim(r1.swe) # 145   3

r1.swe$percent_host_removed <- 100*(r1.swe$n_reads_cleaned - r1.swe$n_reads_nonhost)/r1.swe$n_reads_cleaned
r1.swe$n_reads_host_removed <- r1.swe$n_reads_cleaned - r1.swe$n_reads_nonhost

# what were subset of samples finally used?
keep_samps

sel <- which(r1.swe$sample %in% keep_samps) # qty 76

r1.swe <- r1.swe[sel, ]
dim(r1.swe) # 76 5


hist(r1.swe$percent_host_removed);summary(r1.swe$percent_host_removed)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001006 0.013148 0.031051 0.081442 0.077532 0.840383

median(r1.swe$percent_host_removed) # 0.03105059
# For the Swedish T2D cohort study samples contained a median (and interquartile range, IQR): 0.03 % (0.01 - 0.08 %) of reads that were classified as human sequences.

hist(r1.swe$n_reads_cleaned);summary(r1.swe$n_reads_cleaned)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 1225835  5579366  7822743  9080544 12870602 22474687

# For the Swedish T2D cohort study data we used median 7,822,743 (IQR: 5,579,366 - 12,870,602) cleaned sequences

hist(r1.swe$n_reads_host_removed);summary(r1.swe$n_reads_host_removed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  173    1203    2152    6970    6067  106552 

# median 2,152 (IQR 1203 - 6,067)

hist(r1.swe$n_reads_nonhost);summary(r1.swe$n_reads_nonhost)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 1223102  5572690  7820878  9073574 12868662 22466068 

median(r1.swe$n_reads_nonhost) # 7,820,878

# No-host sequences: median 7,820,878 IQR 5,572,690 - 12,868,662

#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - read in superfocus - fxn potential outputs
#-------------------------

# SUPER-FOCUS results copied here ...

superfocus_out_dir <- "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe/3_fxn_superfocus_copy_hostremoval"
list.dirs(superfocus_out_dir)
head( list.dirs(superfocus_out_dir) )

# don't keep 1st directory
( results_dirs <- list.dirs(superfocus_out_dir)[-c(1)] )
length(results_dirs) # 145

sel <- grep(pattern = paste0("*",keep_samps,collapse = "|"), x = results_dirs) # qty 76
results_dirs <- results_dirs[sel]
length(results_dirs) # 76

head(results_dirs)
# [1] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe/3_fxn_superfocus_copy_hostremoval/superfocus_out_ERR260139"
# [2] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe/3_fxn_superfocus_copy_hostremoval/superfocus_out_ERR260140"
# [3] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe/3_fxn_superfocus_copy_hostremoval/superfocus_out_ERR260144"
# [4] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe/3_fxn_superfocus_copy_hostremoval/superfocus_out_ERR260147"
# [5] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe/3_fxn_superfocus_copy_hostremoval/superfocus_out_ERR260151"
# [6] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe/3_fxn_superfocus_copy_hostremoval/superfocus_out_ERR260152"

names(results_dirs) <- gsub(pattern = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe/3_fxn_superfocus_copy_hostremoval/superfocus_out_", replacement = "", x = results_dirs)
head(results_dirs)
# ERR260139 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe/3_fxn_superfocus_copy_hostremoval/superfocus_out_ERR260139" 
# ERR260140 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe/3_fxn_superfocus_copy_hostremoval/superfocus_out_ERR260140" 
# ERR260144 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe/3_fxn_superfocus_copy_hostremoval/superfocus_out_ERR260144" 
# ERR260147 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe/3_fxn_superfocus_copy_hostremoval/superfocus_out_ERR260147" 
# ERR260151 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe/3_fxn_superfocus_copy_hostremoval/superfocus_out_ERR260151" 
# ERR260152 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe/3_fxn_superfocus_copy_hostremoval/superfocus_out_ERR260152"


sampid <- sort(keep_samps)

# check identical order
identical(sampid, names(results_dirs)) # TRUE
length(results_dirs) # 76


# In this data one Run corresponds to a single Sample_ID !!!

# collate results into a long-format table

sfx.long <- data.frame(sampleID=NA, subsys_L1=NA, subsys_L2=NA, subsys_L3=NA,fxn=NA,percent_abun=NA)

for (i in 1:length(sampid)) {
  #i<-1
  this_samp <- sampid[i]
  sel.folder <- grep(pattern = this_samp, x = results_dirs)
  this_folder <- results_dirs[sel.folder]
  
  #tab1 <- read_excel(path = paste0(this_folder,"/output_all_levels_and_function.xlsx"), skip = 4, col_names = TRUE)
  
  tab <- read.csv(file = paste0(this_folder,"/output_all_levels_and_function.xls"), sep = "\t", skip = 4 )
  # names(tab)
  # [1] "Subsystem.Level.1"
  # [2] "Subsystem.Level.2"
  # [3] "Subsystem.Level.3"
  # [4] "Function"
  # [5] "X.scratch.user.lidd0026.nz.pilot.nz_2_fastp_qc.10923_R1.good.fastq"
  # [6] "X.scratch.user.lidd0026.nz.pilot.nz_2_fastp_qc.10923_R1.good.fastq.."
  
  
  # [1] "Subsystem.Level.1"
  # [2] "Subsystem.Level.2"
  # [3] "Subsystem.Level.3"
  # [4] "Function"
  # [5] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq"
  # [6] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq"
  # [7] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq"
  # [8] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq"
  # [9] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq.." # this is %
  # [10] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq.." # this is %
  # [11] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq.." # this is %
  # [12] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq.." # this is %
  
  
  tab$sampid <- this_samp
  names(tab)
  
  #tab <- tab[,c(7,1,2,3,4,6)]
  
  # last column is sampid
  # take average of percentages
  
  sel.col.percent <- grep(pattern = "_non_host.1.fastq..$", x = names(tab))
  if (length(sel.col.percent)>1) {
    tab$percent_abun <- apply(X = tab[ ,sel.col.percent], MARGIN = 1, FUN = mean )
  } else {
    tab$percent_abun <- tab[ ,sel.col.percent]
  }
  
  # sum(tab$percent_abun) # 100
  # mean(tab$percent_abun) # 0.004338583
  
  names(sfx.long) # "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"    "percent_abun"
  # names(tab)
  # [1] "Subsystem.Level.1"
  # [2] "Subsystem.Level.2"
  # [3] "Subsystem.Level.3"
  # [4] "Function"
  # ...
  # [13] "sampid"
  # [14] "percent_abun"
  
  tab <- tab[ ,c("sampid","Subsystem.Level.1","Subsystem.Level.2","Subsystem.Level.3","Function","percent_abun")]
  names(tab) <- names(sfx.long)
  
  sfx.long <- rbind(sfx.long,tab)
  
  print(paste0("completed ",i," - sample ID: ",sampid[i]))
}

head(sfx.long)
# remove empty 1st row
sfx.long <- sfx.long[-1, ]
dim(sfx.long) # 595371      6
head(sfx.long)
# sampleID                   subsys_L1 subsys_L2           subsys_L3
# 2 ERR260139 Amino Acids and Derivatives         - Amino acid racemase
# 3 ERR260139 Amino Acids and Derivatives         - Amino acid racemase
# 4 ERR260139 Amino Acids and Derivatives         - Amino acid racemase
# 5 ERR260139 Amino Acids and Derivatives         - Amino acid racemase
# 6 ERR260139 Amino Acids and Derivatives         - Amino acid racemase
# 7 ERR260139 Amino Acids and Derivatives         - Amino acid racemase
# fxn percent_abun
# 2                                                                            Alanine_racemase_(EC_5.1.1.1) 9.245298e-03
# 3                                                                         Aspartate_racemase_(EC_5.1.1.13) 1.260722e-03
# 4                                                                   Diaminopimelate_epimerase_(EC_5.1.1.7) 9.455418e-03
# 5                                                                          Glutamate_racemase_(EC_5.1.1.3) 8.909105e-03
# 6 UDP-N-acetylmuramoyl-tripeptide--D-alanyl-D-alanine_ligase_(EC_6.3.2.10)_/_Alanine_racemase_(EC_5.1.1.1) 4.584445e-05
# 7                                                                                     hypothetical_protein 1.260722e-04

sfx.long$full_fxn_tax <- paste0(sfx.long$subsys_L1,"___", sfx.long$subsys_L2,"___", sfx.long$subsys_L3,"___", sfx.long$fxn)

## translate from long to wide format
names(sfx.long)
# "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"          "percent_abun" "full_fxn_tax"

sfx.wide <- dcast(sfx.long, formula = full_fxn_tax ~ sampleID, value.var = "percent_abun")
dim(sfx.wide) # 17912    77

sel.na <- which(is.na(sfx.wide),arr.ind = TRUE)
sfx.wide[sel.na] <- 0

# function taxonomy
full_fxn_names <- sfx.wide$full_fxn_tax

length(full_fxn_names) # 17912
length(unique(full_fxn_names)) # 17912

names(full_fxn_names) <- paste0("fxn_",c(1:length(full_fxn_names)))
head(full_fxn_names)
# fxn_1 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylaconitate_cis-trans_isomerase" 
# fxn_2 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylcitrate_dehydratase_(2-methyl-trans-aconitate_forming)_(EC_4.2.1.117)" 
# fxn_3 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylcitrate_dehydratase_FeS_dependent_(EC_4.2.1.79)" 
# fxn_4 
# "Amino Acids and Derivatives___-___Amino acid racemase___4-hydroxyproline_epimerase_(EC_5.1.1.8)" 
# fxn_5 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)" 
# fxn_6 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic"


tax.fxn <- separate(sfx.wide, full_fxn_tax, c("subsys_L1", "subsys_L2", "subsys_L3", "fxn"), sep= "___", remove=TRUE)
# remove sample ids
tax.fxn <- tax.fxn[ ,-which(names(tax.fxn) %in% sampid)]

row.names(tax.fxn) <- names(full_fxn_names)

head(sfx.wide)

names(sfx.wide)
# [1] "full_fxn_tax" "ERR260139"    "ERR260140"    "ERR260144"    "ERR260147"    "ERR260151"    "ERR260152"    "ERR260153"    "ERR260159"   
# [10] "ERR260161"    "ERR260162"    "ERR260163"    "ERR260165"    "ERR260166"    "ERR260167"    "ERR260169"    "ERR260170"    "ERR260171"   
# [19] "ERR260173"    "ERR260174"    "ERR260175"    "ERR260179"    "ERR260180"    "ERR260181"    "ERR260185"    "ERR260186"    "ERR260188"   
# [28] "ERR260189"    "ERR260190"    "ERR260193"    "ERR260198"    "ERR260199"    "ERR260201"    "ERR260203"    "ERR260204"    "ERR260205"   
# [37] "ERR260206"    "ERR260207"    "ERR260209"    "ERR260210"    "ERR260215"    "ERR260216"    "ERR260217"    "ERR260218"    "ERR260221"   
# [46] "ERR260223"    "ERR260224"    "ERR260225"    "ERR260226"    "ERR260227"    "ERR260230"    "ERR260231"    "ERR260234"    "ERR260241"   
# [55] "ERR260242"    "ERR260243"    "ERR260244"    "ERR260246"    "ERR260250"    "ERR260251"    "ERR260252"    "ERR260253"    "ERR260255"   
# [64] "ERR260256"    "ERR260258"    "ERR260259"    "ERR260260"    "ERR260263"    "ERR260264"    "ERR260265"    "ERR260266"    "ERR260267"   
# [73] "ERR260268"    "ERR260271"    "ERR260273"    "ERR260276"    "ERR275252"

#names(sfx.wide) <- gsub(pattern = "-", replacement = "_", x = names(sfx.wide))

identical(as.character(full_fxn_names), sfx.wide$full_fxn_tax) # TRUE

row.names(sfx.wide) <- names(full_fxn_names)
sfx.wide <- sfx.wide[ ,-1]

names(sfx.wide)

head(sampid)
# "ERR260139" "ERR260140" "ERR260144" "ERR260147" "ERR260151" "ERR260152"

length(sampid) # 76

names(sampid) # NULL - in this case there is NOT an alternative sample name being used

# check alignment of sample IDs and sample names
identical(names(sfx.wide) , as.character(sampid)) # TRUE

#NOT RUN THIS TIME
#names(sfx.wide) <- names(sampid)


names(tax.fxn) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn"
dim(tax.fxn) # 17912     4

length(unique(tax.fxn$subsys_L1)) # 35
length(unique(tax.fxn$subsys_L2)) # 185
length(unique(tax.fxn$subsys_L3)) # 1077
length(unique(tax.fxn$fxn)) # 9560


#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - get into Phyloseq object
#-------------------------

# sfx.wide - is equiv to OTU table

# tax.fxn - is equiv to TAX table

# meta - is equiv to sample table

## Create 'taxonomyTable'
#  tax_table - Works on any character matrix. 
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
tax.m <- as.matrix( tax.fxn )
dim(tax.m) # 17912     4

TAX <- tax_table( tax.m )


## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
otu.m <- as.matrix( sfx.wide )
dim(otu.m)
# 17912    76

OTU <- otu_table(otu.m, taxa_are_rows = TRUE)


## Create a phyloseq object, merging OTU & TAX tables
phy = phyloseq(OTU, TAX)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 17912 taxa and 76 samples ]
# tax_table()   Taxonomy Table:    [ 17912 taxa by 4 taxonomic ranks ]

sample_names(phy)
# [1] "ERR260139" "ERR260140" "ERR260144" "ERR260147" "ERR260151" "ERR260152" "ERR260153" "ERR260159" "ERR260161" "ERR260162" "ERR260163"
# [12] "ERR260165" "ERR260166" "ERR260167" "ERR260169" "ERR260170" "ERR260171" "ERR260173" "ERR260174" "ERR260175" "ERR260179" "ERR260180"
# [23] "ERR260181" "ERR260185" "ERR260186" "ERR260188" "ERR260189" "ERR260190" "ERR260193" "ERR260198" "ERR260199" "ERR260201" "ERR260203"
# [34] "ERR260204" "ERR260205" "ERR260206" "ERR260207" "ERR260209" "ERR260210" "ERR260215" "ERR260216" "ERR260217" "ERR260218" "ERR260221"
# [45] "ERR260223" "ERR260224" "ERR260225" "ERR260226" "ERR260227" "ERR260230" "ERR260231" "ERR260234" "ERR260241" "ERR260242" "ERR260243"
# [56] "ERR260244" "ERR260246" "ERR260250" "ERR260251" "ERR260252" "ERR260253" "ERR260255" "ERR260256" "ERR260258" "ERR260259" "ERR260260"
# [67] "ERR260263" "ERR260264" "ERR260265" "ERR260266" "ERR260267" "ERR260268" "ERR260271" "ERR260273" "ERR260276" "ERR275252"

### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

head(df.samp)
#          Sample Country subset  Status      Bases       Run group_label      age
# 640 NG-5425_137            SWE ND CTRL  797799936 ERR260268      Normal 70.04658
# 645 NG-5636_118            SWE ND CTRL 3307295096 ERR260250      Normal 70.01644
# 646 NG-5636_121            SWE ND CTRL 2503143600 ERR260251      Normal 69.27397
# 647 NG-5636_126            SWE ND CTRL 2981401426 ERR260252      Normal 69.48493
# 648 NG-5636_127            SWE ND CTRL 2480018034 ERR260253      Normal 70.00822
# 650 NG-5636_141            SWE ND CTRL 1999859186 ERR260255      Normal 69.42466

dim(df.samp) # 76 7

head(row.names(df.samp)) # "640" "645" "646" "647" "648" "650"
row.names(df.samp) <- df.samp$Run

identical(row.names(df.samp), sample_names(phy)) # FALSE
length(which(row.names(df.samp) %in% sample_names(phy))==TRUE) # 76

# rearrange order
df.samp2 <- df.samp[ sample_names(phy), ]
identical(row.names(df.samp2), sample_names(phy)) # TRUE

SAMP <- sample_data(df.samp2)


### Combine SAMPDATA into phyloseq object
phy <- merge_phyloseq(phy, SAMP)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 17912 taxa and 76 samples ]
# sample_data() Sample Data:       [ 76 samples by 5 sample variables ]
# tax_table()   Taxonomy Table:    [ 17912 taxa by 4 taxonomic ranks ]

head(taxa_names(phy))
# "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"

head(phy@tax_table)
# Taxonomy Table:     [6 taxa by 4 taxonomic ranks]:
#   subsys_L1                     subsys_L2 subsys_L3             fxn                                                                            
# fxn_1 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "2-methylaconitate_cis-trans_isomerase"                                        
# fxn_2 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "2-methylcitrate_dehydratase_(2-methyl-trans-aconitate_forming)_(EC_4.2.1.117)"
# fxn_3 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "2-methylcitrate_dehydratase_FeS_dependent_(EC_4.2.1.79)"                      
# fxn_4 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "4-hydroxyproline_epimerase_(EC_5.1.1.8)"                                      
# fxn_5 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)"                                                
# fxn_6 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic"  


getwd()  # "/Users/lidd0026/WORKSPACE/PROJ/PCaN-NZ/nz-city-resto/modelling/R"

saveRDS(object = phy, file = "phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval.RDS")

#phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval.RDS")

# get stats?

head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 76
NonZeroFxns

mean(NonZeroFxns) # 7833.829
sd(NonZeroFxns) # 1699.389


#-------------------------

#### Forslund T2D-SWE - w/ Host removal - COPY of R code to run CPP steps on HPC
#    1) build reaction search - get reactions and compounds
#    2) get cpd rel abun per sample
#    3) collate compounds for each sample
#-------------------------

# # # # # # # # # # # # #
# #
# # R script for cpp3d
# # - build reaction search in parallel - get_reactions & compounds
# # - get cpd rel abun per sample
# # - collate_compounds
# #
# # Forslund T2D - SWE - Host-removal
# # Craig Liddicoat - Flinders University
# # Running on Pawsey Setonix
# # # # # # # # # # # # #
# 
# # Add a new path
# .libPaths(c("/software/projects/pawsey1216/cliddicoat/setonix/2024.05/r/4.4.1",
#             "/software/projects/pawsey1216/cliddicoat/setonix/2024.05/r/4.3", .libPaths()))
# 
# R.Version()
# 
# # load packages
# #library(readxl); packageVersion("readxl")
# library(parallel); packageVersion("parallel")
# library(doParallel); packageVersion("doParallel")
# library(dplyr); packageVersion("dplyr")
# library(stringr); packageVersion("stringr")
# library(phyloseq); packageVersion("phyloseq") # '1.44.0'
# 
# message("\n# establish folders and input files")
# 
# message("\nworkdir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/5_rcode_cpp'")
# workdir <- "/scratch/pawsey1216/cliddicoat/ft2d_swe/5_rcode_cpp"
# message("\nsetwd(workdir)")
# setwd(workdir)
# message("\ntemp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/5_rcode_cpp/working'")
# temp_dir <- "/scratch/pawsey1216/cliddicoat/ft2d_swe/5_rcode_cpp/working"
# 
# message("\nthis_study <- '-forslund-t2d-swe-hostremoval-pawsey'")
# this_study <- "-forslund-t2d-swe-hostremoval-pawsey"
# message("\nphy <- readRDS('phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval.RDS')")
# phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval.RDS")
# 
# 
# subsys.lut <- readRDS("subsys.lut.RDS")
# rxns.lut <- readRDS("rxns.lut.RDS")
# rxn_pathways.lut <- readRDS("rxn_pathways.lut.RDS")
# compounds.lut <- readRDS("compounds.lut.RDS")
# 
# 
# 
# message("\n### 1) build reaction search in parallel - get_reactions & compounds")
# message("\n# # # # # # # # # #")
# message("\ndf.tax <- as.data.frame(phy@tax_table)")
# df.tax <- as.data.frame(phy@tax_table)
# message("\nhead(row.names(df.tax))")
# head(row.names(df.tax))
# message("\ndim(df.tax)")
# dim(df.tax)
# 
# 
# get_rxns_and_compounds_indiv <- function( df.tax, subsys.lut, rxns.lut, rxn_pathways.lut ) {
#   
#   rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
#   rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
#   
#   sub1 <- df.tax$subsys_L1[i]
#   sub2 <- df.tax$subsys_L2[i]
#   sub3 <- df.tax$subsys_L3[i]
#   
#   fxn.temp <- df.tax$fxn[i]
#   fxn.superfocus.rowlabel <- row.names(df.tax)[i]
#   
#   # store results corresponding to each Superfocus row
#   fxn.list <- list()
#   fxn.list[[ fxn.superfocus.rowlabel  ]] <- list()
#   
#   # check for multiple functions/reactions?
#   flag1 <- grepl(pattern = "_/_|/", x = fxn.temp)
#   flag2 <- grepl(pattern = "_@_", x = fxn.temp)
#   if (!any(flag1,flag2)==TRUE) {
#     # no multiples
#     fxns <- fxn.temp
#   } else if (flag1==TRUE) {
#     fxns <- unlist( strsplit(fxn.temp, split = "_/_") )  ###### WHAT ABOUT SPLIT FOR "/" WITHOUT UNDERSCORES ??
#   } else {
#     fxns <- unlist( strsplit(fxn.temp, split = "_@_") )
#   }
#   # remove underscores
#   ( fxns <- gsub(pattern = "_", replacement = " ", x = fxns) )
#   
#   # process each fxn & store attributes
#   df.fxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=1:length(fxns),`f__in`=fxns, matching_method=NA, rxns=NA)
#   
#   # Identify '/' separators with no '_'  ??
#   
#   for (f in 1:length(fxns)) {  # this accounts for multiple functions/reactions reported in Superfocus outputs
#     #f<-1
#     #f<-2
#     f.in <- fxns[f]
#     
#     # these concatenated expressions will be used to look for exact match using hierarchy in ModelSEED Subsystem table
#     full_hier_target <- paste0(sub1,"__",sub2,"__",sub3,"__",f.in)
#     full_hier_list <- paste0(subsys.lut$Class,"__",subsys.lut$Subclass,"__",gsub("_"," ",subsys.lut$Name),"__",subsys.lut$Role)
#     
#     ## data cleaning
#     
#     # trim off '_#' and '_##' tags
#     trim_nchar <- str_locate(string = f.in, pattern = " # | ## ")[1]
#     if (!is.na(trim_nchar) & length(trim_nchar)==1) {
#       f.in <- substring(text = f.in , first = 1, last = trim_nchar-1)
#     }
#     
#     # Eliminate unwanted parsing of regular expressions: '[', ']','***', '(', ')'
#     f.in <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\} ", replacement ="." , x = f.in) # used later
#     
#     #rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
#     #rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
#     
#     full_hier_target <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_target)
#     full_hier_list <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_list)
#     
#     sel.rx <- grep(pattern = full_hier_target, x = full_hier_list)
#     
#     ## ALTERNATIVE #1 == FULL HIERACHICAL MATCH
#     if (length(sel.rx)>=1) {
#       df.fxns$matching_method[f] <- "Exact hierachy match"
#       df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
#       
#     } else if (str_detect(string = fxns[f], pattern = " \\(EC ")) {  ## ALTERNATIVE #2 == MATCHING ECs
#       # search by EC id if present
#       
#       f.in <- fxns[f] # this goes back to string with brackets for EC
#       ## LOOK FOR MULTIPLE ECs ?
#       
#       how_many_ECs <- str_count(string = f.in, pattern = "\\(EC.*?\\)")
#       
#       ECs <- as.character( str_extract_all(string = f.in, pattern = "\\(EC.*?\\)", simplify = TRUE) )
#       #class(ECs)
#       ECs <- gsub(pattern = "\\(EC |\\)", replacement = "", x = ECs)
#       ECs.collapse <- paste0(ECs, collapse = "|")
#       
#       sel.rx <- which(rxns.lut$ec_numbers == ECs.collapse)
#       
#       if (length(how_many_ECs)==0 | length(ECs)==0) {
#         # there was a glitch, database typo, or some error in identifying the EC number
#         df.fxns$matching_method[f] <- "No match found"
#         df.fxns$rxns[f] <- NA
#         
#       } else if (length(sel.rx)>=1) {
#         # combined EC hits identified
#         df.fxns$matching_method[f] <- "EC number"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if (length(which(rxns.lut$ec_numbers %in% ECs)) >=1) {
#         # treat EC hits individually
#         sel.rx <- which(rxns.lut$ec_numbers %in% ECs) # look 1st where ECs are exact matches for EC numbers in Reactions lookup table
#         
#         df.fxns$matching_method[f] <- "EC number"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if (length(grep(pattern = ECs, x = rxns.lut$ec_numbers)) >=1) {
#         # this allows EC to be part of a combination of EC numbers that are listed in Reactions lookup table
#         sel.rx <- grep(pattern = ECs, x = rxns.lut$ec_numbers)
#         
#         df.fxns$matching_method[f] <- "EC number"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else {
#         # it had an EC number but couldn't find a match in the EC numbers listed in Reaction lookup table
#         df.fxns$matching_method[f] <- "No match found"
#         df.fxns$rxns[f] <- NA
#         
#       }
#       # END EC matching
#       
#       
#     } else {  ## ALTERNATIVE 3 == FXN NAME MATCHING
#       ## otherwise attempt to match function name - a) first look for exact matches   ########## then b) closest match above a threshold
#       # 1. 'reactions' table by name: rxns.lut$name
#       # 2. 'reactions' table by aliases: rxns.lut$aliases
#       # 3. 'Model SEED Subsystems' table by Role: subsys.lut$Role
#       # 4. 'Unique_ModelSEED_Reaction_Pathways' table by External ID: rxn_pathways.lut$External_rxn_name
#       
#       if ( length( grep(pattern = f.in, x = rxns.lut$name) )>=1 ) {
#         # 1a - exact match - rxns.lut$name
#         sel.rx <- grep(pattern = f.in, x = rxns.lut$name)
#         #rxns.lut$name[sel.rx]
#         df.fxns$matching_method[f] <- "Matched Reactions name"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if ( length( grep(pattern = f.in, x = rxns.lut$aliases) )>=1 ) {
#         # 2a - exact match - rxns.lut$aliases
#         sel.rx <- grep(pattern = f.in, x = rxns.lut$aliases)
#         #rxns.lut$aliases[sel.rx]
#         #rxns.lut$name[sel.rx]
#         
#         df.fxns$matching_method[f] <- "Matched Reactions aliases"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if ( length( grep(pattern = f.in, x = subsys.lut$Role) )>=1 ) {
#         # 3a - exact match - subsys.lut$Role
#         sel.rx <- grep(pattern = f.in, x = subsys.lut$Role)
#         #subsys.lut$Role[sel.rx]
#         #subsys.lut$Reaction[sel.rx]
#         
#         df.fxns$matching_method[f] <- "Matched Subsytem role"
#         df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
#         
#       } else if ( length( grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name) )>=1 ) {
#         # 4a - exact match - rxn_pathways.lut$External_rxn_name
#         sel.rx <- grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name)
#         
#         df.fxns$matching_method[f] <- "Matched ModelSEED Reaction pathways"
#         df.fxns$rxns[f] <- paste0( unique(rxn_pathways.lut$rxn_id[sel.rx]), collapse = ";")
#         
#         
#       } else {
#         df.fxns$matching_method[f] <- "No match found"
#         df.fxns$rxns[f] <- NA
#         
#       }
#       
#       ## DON'T RUN PARTIAL MATCHING AT THIS STAGE
#       
#       
#     } # END function - reaction search
#     
#     #fxn.list[[ fxn.superfocus.rowlabel  ]][[ f ]][[ "fxns" ]] <- df.fxns
#     
#     #print(paste0("completed fxn ", f))
#     
#     
#     ## now investigate these reactions ...
#     # Reactions lookup table: 
#     # - "equation": Definition of reaction expressed using compound IDs and after protonation
#     # Compounds lookup table:
#     # - "formula": Standard chemical format (using Hill system) in protonated form to match reported charge
#     #df.fxns
#     
#     
#     #if (df.fxns$matching_method == "No match found") {
#     if (df.fxns$rxns[f] == "" | is.na(df.fxns$rxns[f])) {
#       
#       df.Rxns <- NA
#       df.Compounds <- NA
#       
#     } else { # reaction(s) were identified
#       
#       # consider reactions for this f.in only (possibly > 1 f.in per Superfocus row)
#       f.in.rxns <- unique(unlist(str_split(string = df.fxns$rxns[f], pattern = ";")))
#       
#       df.Rxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel, f=f, f__in=f.in,rxn_id= f.in.rxns,
#                             rxn_name=NA, rxn_eqn=NA, rxn_defn=NA,compds=NA,compd_coef=NA, chem_formx=NA )
#       
#       for (r in 1:dim(df.Rxns)[1]) {
#         #r<-1
#         #this_rxn <- "rxn00004"
#         this_rxn <- df.Rxns$rxn_id[r]
#         sel <- which(rxns.lut$id == this_rxn)
#         ( df.Rxns$rxn_name[r] <- rxns.lut$name[sel] )
#         ( df.Rxns$rxn_eqn[r] <- rxns.lut$equation[sel] )
#         ( df.Rxns$rxn_defn[r] <- rxns.lut$definition[sel] )
#         
#         # extract compound info
#         
#         #df.Rxns$rxn_eqn[r]
#         #[1] "(1) cpd00010[0] + (1) cpd29672[0] <=> (1) cpd00045[0] + (1) cpd11493[0]"
#         #[1] "(45) cpd00144[0] + (45) cpd00175[0] <=> (45) cpd00014[0] + (45) cpd00091[0] + (1) cpd15634[0]"
#         
#         ( compds.idx <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd")[[1]][,"start"] )
#         # 5 23 43 61
#         # 6 25 46 65 83
#         
#         ( compds <- as.character( str_extract_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd.....", simplify = TRUE) ) )
#         # "cpd00010" "cpd29672" "cpd00045" "cpd11493"
#         
#         if (length(compds)>=1) {
#           
#           df.Rxns$compds[r] <- paste0(compds, collapse = ";")
#           
#           ## get compound coefficients?
#           start_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\(")[[1]][,"start"]
#           end_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\)")[[1]][,"start"]
#           ( compd.coeff <- as.numeric( substring(text = df.Rxns$rxn_eqn[r], first = start_brackets+1, last = end_brackets-1)) )
#           
#           df.Rxns$compd_coef[r] <- paste0(compd.coeff, collapse = ";")
#           
#           # get formulas of compounds
#           
#           formx <-filter(compounds.lut, id %in% compds )
#           row.names(formx) <- formx$id
#           ( formx.char <- formx[compds, ]$formula )
#           # "C21H32N7O16P3S" "HOR"            "C10H11N5O10P2"  "C11H22N2O7PRS" 
#           # "C15H19N2O18P2"      "C17H25N3O17P2"      "C9H12N2O12P2"       "C9H11N2O9P"         "C630H945N45O630P45"
#           # "C7H7O7" "H2O"    "C7H5O6"
#           df.Rxns$chem_formx[r] <- paste0(formx.char, collapse = ";")
#           
#           ( compd.names <- formx[compds, ]$name )
#           # "2-methyl-trans-aconitate" "cis-2-Methylaconitate"
#           
#           temp.df.Compounds <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=f, f__in=f.in,rxn_id= f.in.rxns[r], 
#                                           cpd_id=compds, cpd_name=compd.names, cpd_form=formx.char, cpd_molar_prop=compd.coeff #, 
#                                           #OC_x=OC_ratio, HC_y=HC_ratio , NC_z=NC_ratio 
#           )
#           
#         } else {
#           # No specified reaction equation or chemical formula info
#           df.Rxns$compds[r] <- NA
#           df.Rxns$compd_coef[r] <- NA
#           df.Rxns$chem_formx[r] <- NA
#           
#           temp.df.Compounds <- NA
#           
#         }
#         
#         if (r==1) { df.Compounds <- temp.df.Compounds }
#         
#         if (r>1 & is.data.frame(df.Compounds) & is.data.frame(temp.df.Compounds)) { df.Compounds <- rbind(df.Compounds, temp.df.Compounds) }
#         
#         # clean up - if there are additional reactions?
#         temp.df.Compounds <- NA
#         
#       } # END loop for r - rxn_id's per f/f.in
#       
#     } # END else loop when reactions identified
#     
#     # store results corresponding to each sub-reaction of each Superfocus row
#     fxn.list[[ fxn.superfocus.rowlabel  ]][[ "fxns" ]] <- df.fxns
#     
#     if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]] <- list() } # set this only once
#     fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]][[ f ]] <- df.Rxns
#     
#     if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]] <- list() } # set this only once
#     fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]][[ f ]] <- df.Compounds
#     
#     
#   } # END loop - f in 1:length(fxns)) - to account for multiple functions/reactions reported in each row of Superfocus outputs
#   
#   saveRDS(object = fxn.list, file = paste0(temp_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds") )
#   
# } # END function to be run in parallel for each superfocus row
# 
# 
# # # # # # # # # # # # # # # # # # #
# 
# no_forks <- 8
# 
# # this makes clusters on Unix-like system (may need to use other alternative for Windows)
# cl<-makeForkCluster(nnodes = no_forks)      # no of nodes will depend on your HPC facility
# registerDoParallel(cl)
# 
# foreach(i=1:dim(df.tax)[1] , .packages=c('stringr', 'dplyr')) %dopar%  #
#   get_rxns_and_compounds_indiv( df.tax=df.tax, subsys.lut=subsys.lut, rxns.lut=rxns.lut, rxn_pathways.lut=rxn_pathways.lut )
# 
# stopCluster(cl)
# 
# 
# message("\n## assemble results")
# 
# message("\n(num_results_files <- dim(df.tax)[1])")
# (num_results_files <- dim(df.tax)[1])
# 
# # assemble all compound data outputs
# # start with blank row
# 
# df.out <- data.frame(superfocus_fxn=NA, f=NA, f__in=NA, rxn_id=NA, cpd_id=NA, cpd_name=NA, cpd_form=NA, cpd_molar_prop=NA )
# 
# for (i in 1:num_results_files) {
#   fxn.superfocus.rowlabel <- row.names(df.tax)[i]
#   temp <- readRDS(paste0(temp_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))
#   
#   f_no <- length( temp[[1]][["compounds"]] )
#   
#   for (f in 1:f_no) {
#     #f<-2
#     # only add non-NA results
#     if (is.data.frame( temp[[1]][["compounds"]][[f]] )) {
#       
#       df.temp <- temp[[1]][["compounds"]][[f]]
#       ok <- complete.cases(df.temp)
#       df.temp <- df.temp[ which(ok==TRUE), ] # updated version will include some compounds with vK coordinates that are NA. vK coordinates are considered later
#       df.out <- rbind(df.out,df.temp)
#     }
#   }
#   print(paste0("added df ",i," of ",num_results_files ))
#   
# }
# 
# 
# message("\nstr(df.out)")
# str(df.out)
# 
# 
# saveRDS(object = df.out, file = paste0("df.out--get_rxns_and_compounds_indiv-",this_study,".RDS"))
# 
# # remove NA first row
# message("\nhead(df.out)")
# head(df.out)
# 
# df.out <- df.out[-1, ]
# 
# message("\ndim(df.out)")
# dim(df.out)
# 
# 
# message("\n## normalise molar_prop to cpd_relabun so total of 1 per superfocus function")
# 
# df.out$cpd_molar_prop_norm <- NA
# 
# message("\nlength(unique(df.out$superfocus_fxn))")
# length(unique(df.out$superfocus_fxn))
# 
# message("\nphy")
# phy
# 
# message("\n% of functions represented - with compound information")
# 100*(length(unique(df.out$superfocus_fxn)) / ntaxa(phy))
# 
# 
# fxns_found <- unique(df.out$superfocus_fxn)
# 
# for (k in 1:length(fxns_found)) {
#   #k<-1
#   this_fxn <- fxns_found[k]
#   sel <- which(df.out$superfocus_fxn == this_fxn)
#   
#   sum_molar_prop <- sum( df.out$cpd_molar_prop[sel], na.rm = TRUE)
#   # calculate 
#   
#   df.out$cpd_molar_prop_norm[sel] <- df.out$cpd_molar_prop[sel]/sum_molar_prop
#   
#   print(paste0("completed ",k))
#   
# }
# 
# message("\nsum(df.out$cpd_molar_prop_norm)")
# sum(df.out$cpd_molar_prop_norm)
# 
# message("\nsample_sums(phy)")
# sample_sums(phy)
# 
# message("\ngetwd()")
# getwd()
# 
# saveRDS(object = df.out, file = paste0("df.out--tidy-compounds_indiv-cpp3d-",this_study,".RDS"))
# 
# 
# 
# message("\n### 2) get cpd rel abun per sample")
# message("\n# # # # # # # # # #")
# 
# 
# df.OTU <- as.data.frame( phy@otu_table ) # this is Superfocus functional relative abundance data represented in phyloseq OTU abundance table
# message("\ndim(df.OTU)")
# dim(df.OTU)
# 
# 
# get_cpd_relabun_per_sample <- function(phy_in, dat.cpd) {
#   
#   this_samp <- sample_names(phy_in)[i]
#   df.OTU <- as.data.frame( phy_in@otu_table[ ,this_samp] )
#   
#   dat.cpd$sample <- this_samp
#   
#   dat.cpd$cpd_rel_abun_norm <- NA
#   
#   fxns_all <- row.names(df.OTU)
#   
#   for (k in 1:length(fxns_all)) {
#     #k<-1
#     this_fxn <- fxns_all[k]
#     sel <- which(dat.cpd$superfocus_fxn == this_fxn)
#     
#     if (length(sel)>=1) {
#       dat.cpd$cpd_rel_abun_norm[sel] <- df.OTU[this_fxn, ]*dat.cpd$cpd_molar_prop_norm[sel]
#       
#     }
#   } # END rel abun values for all relevant functions added
#   
#   saveRDS(object = dat.cpd, file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
#   
# } # END
# 
# 
# no_forks <- 8
# 
# # this makes clusters on Unix-like system
# cl<-makeForkCluster(nnodes = no_forks)      # no of nodes will depend on your HPC facility
# registerDoParallel(cl)
# 
# foreach(i=1: length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
#   get_cpd_relabun_per_sample( phy_in = phy, dat.cpd = df.out)
# 
# stopCluster(cl)
# 
# 
# message("\n## assemble results")
# 
# # output 1
# i<-1
# this_samp <- sample_names(phy)[i]
# dat <- readRDS( file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
# head(dat)
# 
# for ( i in 2:length(sample_names(phy)) ) {
#   this_samp <- sample_names(phy)[i]
#   temp <- readRDS( file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
#   dat <- rbind(dat, temp)
#   print(paste0("completed ",i))
# }
# 
# 
# saveRDS(object = dat, file = paste0("dat.cpd-long-all-samps-cpp3d-",this_study,".rds") )
# 
# rm(temp)
# 
# message("\nstr(dat)")
# str(dat)
# 
# message("\nsum(dat$cpd_rel_abun_norm)")
# sum(dat$cpd_rel_abun_norm)
# 
# message("\naverage functional relative abundance per sample")
# message("\nsum(dat$cpd_rel_abun_norm)/nsamples(phy)")
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# 
# message("\nnames(dat)")
# names(dat)
# 
# message("\nlength(unique(dat$cpd_id))")
# length(unique(dat$cpd_id))
# 
# 
# 
# 
# message("\n### 3) collate_compounds within each sample")
# message("\n# # # # # # # # # #")
# 
# 
# unique_cpd <- unique(dat$cpd_id)
# samp_names <- sample_names(phy)
# 
# 
# collate_compounds <- function(dat.cpd, unique_cpd, samp) {
#   #i<-1
#   #samp = samp_names[i]
#   #dat.cpd = dat[which(dat$sample == samp_names[i]), ]
#   
#   this_samp <- samp
#   
#   cpd_data <- data.frame(cpd_id = unique_cpd, sample=this_samp, cpd_rel_abun=NA)
#   
#   for (c in 1:length(unique_cpd)) {
#     #c<-1
#     this_cpd <- unique_cpd[c]
#     sel.cpd <- which(dat.cpd$cpd_id == this_cpd)
#     
#     if (length(sel.cpd) >=1) {
#       cpd_data$cpd_rel_abun[c] <- sum(dat.cpd$cpd_rel_abun_norm[sel.cpd])
#     }
#     
#   } # END all compounds
#   
#   saveRDS(object = cpd_data, file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
#   
# } # END
# 
# 
# 
# no_forks <- 4
# 
# # this makes clusters on Unix-like system
# cl<-makeForkCluster(nnodes = no_forks)   # no of nodes will depend on your HPC facility
# registerDoParallel(cl)
# 
# foreach(i=1:length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
#   collate_compounds(dat.cpd = dat[which(dat$sample == samp_names[i]), ], unique_cpd = unique_cpd, samp = samp_names[i])
# 
# stopCluster(cl)
# 
# 
# message("\n## assemble results")
# 
# # output 1
# i<-1
# this_samp <- sample_names(phy)[i]
# dat.cpd.collate <- readRDS( file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
# head(dat.cpd.collate)
# 
# for ( i in 2:length(sample_names(phy)) ) {
#   this_samp <- sample_names(phy)[i]
#   temp <- readRDS( file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
#   
#   dat.cpd.collate <- rbind(dat.cpd.collate, temp)
#   
#   print(paste0("completed ",i))
# }
# 
# 
# message("\nstr(dat.cpd.collate)")
# str(dat.cpd.collate)
# 
# message("\nsum(dat.cpd.collate$cpd_rel_abun)")
# sum(dat.cpd.collate$cpd_rel_abun)
# 
# message("\nsum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))")
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# 
# saveRDS(object = dat.cpd.collate, file = paste0("dat.cpd.collate-all-samps-cpp3d-",this_study,".rds" ))
# 
# # END


#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - COPY of OUTOUTS from R code after running CPP steps on HPC
#-------------------------

# $platform
# [1] "x86_64-pc-linux-gnu"
# 
# $arch
# [1] "x86_64"
# 
# $os
# [1] "linux-gnu"
# 
# $system
# [1] "x86_64, linux-gnu"
# 
# $status
# [1] ""
# 
# $major
# [1] "4"
# 
# $minor
# [1] "4.1"
# 
# $year
# [1] "2024"
# 
# $month
# [1] "06"
# 
# $day
# [1] "14"
# 
# $`svn rev`
# [1] "86737"
# 
# $language
# [1] "R"
# 
# $version.string
# [1] "R version 4.4.1 (2024-06-14)"
# 
# $nickname
# [1] "Race for Your Life"
# 
# [1] ‘4.4.1’
# Loading required package: foreach
# Loading required package: iterators
# [1] ‘1.0.17’
# 
# Attaching package: ‘dplyr’
# 
# The following objects are masked from ‘package:stats’:
#   
#   filter, lag
# 
# The following objects are masked from ‘package:base’:
#   
#   intersect, setdiff, setequal, union
# 
# [1] ‘1.1.4’
# [1] ‘1.5.2’
# [1] ‘1.46.0’
# 
# # establish folders and input files
# 
# workdir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/5_rcode_cpp'
# 
# setwd(workdir)
# 
# temp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/5_rcode_cpp/working'
# 
# this_study <- '-forslund-t2d-swe-hostremoval-pawsey'
# 
# phy <- readRDS('phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval.RDS')
# 
# ### 1) build reaction search in parallel - get_reactions & compounds
# 
# # # # # # # # # # #
# 
# df.tax <- as.data.frame(phy@tax_table)
# 
# head(row.names(df.tax))
# [1] "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"
# 
# dim(df.tax)
# [1] 17912     4
# 
# ...
# 
# ## assemble results
# 
# (num_results_files <- dim(df.tax)[1])
# [1] 17912
# 
# ...
# 
# 
# str(df.out)
# 'data.frame':	507392 obs. of  8 variables:
#   $ superfocus_fxn: chr  NA "fxn_2" "fxn_2" "fxn_2" ...
# $ f             : int  NA 1 1 1 1 1 1 1 1 1 ...
# $ f__in         : chr  NA "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" ...
# $ rxn_id        : chr  NA "rxn25279" "rxn25279" "rxn25279" ...
# $ cpd_id        : chr  NA "cpd24620" "cpd00001" "cpd25681" ...
# $ cpd_name      : chr  NA "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" "2-methyl-trans-aconitate" ...
# $ cpd_form      : chr  NA "C7H7O7" "H2O" "C7H5O6" ...
# $ cpd_molar_prop: num  NA 1 1 1 1 1 1 1 1 1 ...
# 
# head(df.out)
# superfocus_fxn  f
# 1           <NA> NA
# 2          fxn_2  1
# 3          fxn_2  1
# 4          fxn_2  1
# 5          fxn_3  1
# 6          fxn_3  1
# f__in
# 1                                                                          <NA>
#   2 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# 3 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# 4 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# 5                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 6                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# rxn_id   cpd_id                                     cpd_name cpd_form
# 1     <NA>     <NA>                                         <NA>     <NA>
#   2 rxn25279 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 3 rxn25279 cpd00001                                          H2O      H2O
# 4 rxn25279 cpd25681                     2-methyl-trans-aconitate   C7H5O6
# 5 rxn03060 cpd01501                              2-Methylcitrate   C7H7O7
# 6 rxn03060 cpd00001                                          H2O      H2O
# cpd_molar_prop
# 1             NA
# 2              1
# 3              1
# 4              1
# 5              1
# 6              1
# 
# dim(df.out)
# [1] 507391      8
# 
# ## normalise molar_prop to cpd_relabun so total of 1 per superfocus function
# 
# length(unique(df.out$superfocus_fxn))
# [1] 10009
# 
# phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 17912 taxa and 76 samples ]
# sample_data() Sample Data:       [ 76 samples by 5 sample variables ]
# tax_table()   Taxonomy Table:    [ 17912 taxa by 4 taxonomic ranks ]
# 
# % of functions represented - with compound information
# [1] 55.87874
# 
# ...
# 
# 
# 
# sum(df.out$cpd_molar_prop_norm)
# [1] 10009
# 
# sample_sums(phy)
# ERR260139 ERR260140 ERR260144 ERR260147 ERR260151 ERR260152 ERR260153 ERR260159 
# 100       100       100       100       100       100       100       100 
# ERR260161 ERR260162 ERR260163 ERR260165 ERR260166 ERR260167 ERR260169 ERR260170 
# 100       100       100       100       100       100       100       100 
# ERR260171 ERR260173 ERR260174 ERR260175 ERR260179 ERR260180 ERR260181 ERR260185 
# 100       100       100       100       100       100       100       100 
# ERR260186 ERR260188 ERR260189 ERR260190 ERR260193 ERR260198 ERR260199 ERR260201 
# 100       100       100       100       100       100       100       100 
# ERR260203 ERR260204 ERR260205 ERR260206 ERR260207 ERR260209 ERR260210 ERR260215 
# 100       100       100       100       100       100       100       100 
# ERR260216 ERR260217 ERR260218 ERR260221 ERR260223 ERR260224 ERR260225 ERR260226 
# 100       100       100       100       100       100       100       100 
# ERR260227 ERR260230 ERR260231 ERR260234 ERR260241 ERR260242 ERR260243 ERR260244 
# 100       100       100       100       100       100       100       100 
# ERR260246 ERR260250 ERR260251 ERR260252 ERR260253 ERR260255 ERR260256 ERR260258 
# 100       100       100       100       100       100       100       100 
# ERR260259 ERR260260 ERR260263 ERR260264 ERR260265 ERR260266 ERR260267 ERR260268 
# 100       100       100       100       100       100       100       100 
# ERR260271 ERR260273 ERR260276 ERR275252 
# 100       100       100       100 
# 
# getwd()
# [1] "/scratch/pawsey1216/cliddicoat/ft2d_swe/5_rcode_cpp"
# 
# ### 2) get cpd rel abun per sample
# 
# # # # # # # # # # #
# 
# dim(df.OTU)
# [1] 17912    76
# 
# ...
# 
# 
# 
# 
# 
# ## assemble results
# superfocus_fxn f
# 2          fxn_2 1
# 3          fxn_2 1
# 4          fxn_2 1
# 5          fxn_3 1
# 6          fxn_3 1
# 7          fxn_3 1
# f__in
# 2 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# 3 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# 4 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# 5                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 6                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 7                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# rxn_id   cpd_id                                     cpd_name cpd_form
# 2 rxn25279 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 3 rxn25279 cpd00001                                          H2O      H2O
# 4 rxn25279 cpd25681                     2-methyl-trans-aconitate   C7H5O6
# 5 rxn03060 cpd01501                              2-Methylcitrate   C7H7O7
# 6 rxn03060 cpd00001                                          H2O      H2O
# 7 rxn03060 cpd02597                        cis-2-Methylaconitate   C7H5O6
# cpd_molar_prop cpd_molar_prop_norm    sample cpd_rel_abun_norm
# 2              1          0.33333333 ERR260139                 0
# 3              1          0.33333333 ERR260139                 0
# 4              1          0.33333333 ERR260139                 0
# 5              1          0.05555556 ERR260139                 0
# 6              1          0.05555556 ERR260139                 0
# 7              1          0.05555556 ERR260139                 0
# 
# ...
# 
# 
# 
# str(dat)
# 'data.frame':	38561716 obs. of  11 variables:
#   $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_2" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" ...
# $ rxn_id             : chr  "rxn25279" "rxn25279" "rxn25279" "rxn03060" ...
# $ cpd_id             : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ cpd_name           : chr  "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" "2-methyl-trans-aconitate" "2-Methylcitrate" ...
# $ cpd_form           : chr  "C7H7O7" "H2O" "C7H5O6" "C7H7O7" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.3333 0.3333 0.3333 0.0556 0.0556 ...
# $ sample             : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun_norm  : num  0 0 0 0 0 0 0 0 0 0 ...
# 
# sum(dat$cpd_rel_abun_norm)
# [1] 5359.577
# 
# average functional relative abundance per sample
# 
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# [1] 70.52075
# 
# names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"              
# [4] "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm"
# [10] "sample"              "cpd_rel_abun_norm"  
# 
# length(unique(dat$cpd_id))
# [1] 7024
# 
# ### 3) collate_compounds within each sample
# 
# # # # # # # # # # #
# 
# ...
# 
# 
# 
# ## assemble results
# cpd_id    sample cpd_rel_abun
# 1 cpd24620 ERR260139 0.0000000000
# 2 cpd00001 ERR260139 4.9823891271
# 3 cpd25681 ERR260139 0.0000000000
# 4 cpd01501 ERR260139 0.0000000000
# 5 cpd02597 ERR260139 0.0001882679
# 6 cpd00851 ERR260139 0.0012152908
# 
# ...
# 
# 
# 
# str(dat.cpd.collate)
# 'data.frame':	533824 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 4.982389 0 0 0.000188 ...
# 
# sum(dat.cpd.collate$cpd_rel_abun)
# [1] 5359.577
# 
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# [1] 70.52075


#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - continue CPP analysis
#-------------------------

# copy output file from HPC
dat.cpd.collate <- readRDS("/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe/run_cpp_r_on_hpc_gcb_r1_swe/dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-pawsey.rds")

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000000 0.000158 0.010040 0.001541 7.135613 

length(unique(dat.cpd.collate$cpd_id)) # 7024
length(unique(dat.cpd.collate$sample)) # 76
str(dat.cpd.collate)
# 'data.frame':	533824 obs. of  3 variables:
# $ cpd_id      : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 4.982389 0 0 0.000188 ...
7024*76 # 533824

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf -7.4308 -3.8023    -Inf -2.8121  0.8534 

# log10 abun
dat.cpd.collate$log10_abun <- dat.cpd.collate$cpd_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(dat.cpd.collate$log10_abun == 0) #
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(dat.cpd.collate$log10_abun[ -subsel.zero ])
  dat.cpd.collate$log10_abun[ subsel.zero ] <- zero_replace
}
dat.cpd.collate$log10_abun <- log10(dat.cpd.collate$log10_abun)

hist(dat.cpd.collate$log10_abun); summary( dat.cpd.collate$log10_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -8.1872 -7.4308 -3.8023 -4.6013 -2.8121  0.8534 

# make group variable from sample name

dat.cpd.collate$group_label <- NA

# from above
df.samp <- as(phy@sam_data, "data.frame")

identical( phy@sam_data$Run , df.samp$Run ) # TRUE
identical( sample_names(phy), df.samp$Run ) # TRUE
unique(df.samp$group_label)
# [1] T2D met- Normal  
# Levels: T2D met- < Normal

#for (i in 1:length(sample_names(phy))) {
for (i in 1:length( df.samp$Run )) {
  #i<-1
  #this_samp <- sample_names(phy)[i]
  this_samp <- df.samp$Run[i]
  sel <- which(dat.cpd.collate$sample == this_samp)
  #dat.cpd.collate$group[sel] <- phy@sam_data$age[i]
  dat.cpd.collate$group_label[sel] <- as.character( df.samp$group_label[i] )
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group_label) # "T2D met-" "Normal" 
dat.cpd.collate$group_label <- factor(dat.cpd.collate$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)
head(dat.cpd.collate)

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-v6.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-v6.rds")

str(dat.cpd.collate)
# 'data.frame':	533824 obs. of  5 variables:
# $ cpd_id      : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 4.982389 0 0 0.000188 ...
# $ log10_abun  : num  -8.187 0.697 -8.187 -8.187 -3.725 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...


#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - CPP as phyloseq object
#    PCoA using CPP vs Functions
#-------------------------

data_in <- dat.cpd.collate
str(data_in)
# 'data.frame':	533824 obs. of  5 variables:

length( unique(data_in$cpd_id) ) # 7024
length( unique(data_in$cpd_id[data_in$cpd_rel_abun > 0]) ) # 7024
length( unique(data_in$sample) ) # 76


### get data into phyloseq object ...

head(data_in)
#     cpd_id    sample cpd_rel_abun log10_abun group_label
# 1 cpd24620 ERR260139 0.0000000000 -8.1872270    T2D met-
# 2 cpd00001 ERR260139 4.9823891271  0.6974376    T2D met-
# 3 cpd25681 ERR260139 0.0000000000 -8.1872270    T2D met-
# 4 cpd01501 ERR260139 0.0000000000 -8.1872270    T2D met-
# 5 cpd02597 ERR260139 0.0001882679 -3.7252238    T2D met-
# 6 cpd00851 ERR260139 0.0012152908 -2.9153198    T2D met-

df.wide <- dcast(data_in, formula = sample + group_label ~ cpd_id , value.var = "cpd_rel_abun" )

df.wide[1:5, 1:10]
#      sample group_label cpd00001 cpd00002  cpd00003  cpd00004  cpd00005  cpd00006   cpd00007 cpd00008
# 1 ERR260139    T2D met- 4.982389 3.228089 0.5465259 0.4748061 0.6432422 0.6426122 0.08357366 1.851411
# 2 ERR260140    T2D met- 5.851583 2.461925 0.5040640 0.4577167 0.5016121 0.5017730 0.06366115 1.415858
# 3 ERR260144    T2D met- 5.069301 3.202583 0.5370572 0.4688221 0.5935540 0.5943987 0.07326647 1.787144
# 4 ERR260147      Normal 4.786099 2.078081 0.6370613 0.6022029 0.4735150 0.4778028 0.14751635 1.296487
# 5 ERR260151    T2D met- 5.013970 2.988423 0.4906686 0.4340136 0.5758568 0.5763056 0.06677339 1.617121

sort(unique(paste0(df.wide$sample,"--",df.wide$group_label)))
# [1] "ERR260139--T2D met-" "ERR260140--T2D met-" "ERR260144--T2D met-" "ERR260147--Normal"   "ERR260151--T2D met-" "ERR260152--T2D met-" "ERR260153--Normal"  
# [8] "ERR260159--T2D met-" "ERR260161--T2D met-" "ERR260162--T2D met-" "ERR260163--Normal"   "ERR260165--T2D met-" "ERR260166--T2D met-" "ERR260167--T2D met-"
# [15] "ERR260169--T2D met-" "ERR260170--Normal"   "ERR260171--Normal"   "ERR260173--T2D met-" "ERR260174--T2D met-" "ERR260175--Normal"   "ERR260179--T2D met-"
# [22] "ERR260180--Normal"   "ERR260181--T2D met-" "ERR260185--T2D met-" "ERR260186--T2D met-" "ERR260188--T2D met-" "ERR260189--T2D met-" "ERR260190--T2D met-"
# [29] "ERR260193--Normal"   "ERR260198--T2D met-" "ERR260199--T2D met-" "ERR260201--T2D met-" "ERR260203--T2D met-" "ERR260204--Normal"   "ERR260205--Normal"  
# [36] "ERR260206--T2D met-" "ERR260207--T2D met-" "ERR260209--Normal"   "ERR260210--T2D met-" "ERR260215--Normal"   "ERR260216--Normal"   "ERR260217--Normal"  
# [43] "ERR260218--Normal"   "ERR260221--Normal"   "ERR260223--Normal"   "ERR260224--Normal"   "ERR260225--Normal"   "ERR260226--Normal"   "ERR260227--Normal"  
# [50] "ERR260230--Normal"   "ERR260231--Normal"   "ERR260234--Normal"   "ERR260241--T2D met-" "ERR260242--Normal"   "ERR260243--Normal"   "ERR260244--Normal"  
# [57] "ERR260246--Normal"   "ERR260250--Normal"   "ERR260251--Normal"   "ERR260252--Normal"   "ERR260253--Normal"   "ERR260255--Normal"   "ERR260256--Normal"  
# [64] "ERR260258--Normal"   "ERR260259--Normal"   "ERR260260--Normal"   "ERR260263--Normal"   "ERR260264--Normal"   "ERR260265--Normal"   "ERR260266--Normal"  
# [71] "ERR260267--Normal"   "ERR260268--Normal"   "ERR260271--T2D met-" "ERR260273--T2D met-" "ERR260276--T2D met-" "ERR275252--T2D met-"

# save group variable
samp <- df.wide[ ,1:2]
row.names(samp) <- samp$sample

# transpose
df.wide <- t(df.wide[ ,-2]) # minus 'group' column

head(df.wide)
# [,1]        [,2]        [,3]        [,4]        [,5]        [,6]        [,7]        [,8]        [,9]        [,10]       [,11]       [,12]       [,13]      
# sample   "ERR260139" "ERR260140" "ERR260144" "ERR260147" "ERR260151" "ERR260152" "ERR260153" "ERR260159" "ERR260161" "ERR260162" "ERR260163" "ERR260165" "ERR260166"
# cpd00001 "4.982389"  "5.851583"  "5.069301"  "4.786099"  "5.013970"  "4.988688"  "5.323399"  "5.879231"  "5.145938"  "5.110849"  "5.144012"  "5.435803"  "5.926055" 
# cpd00002 "3.228089"  "2.461925"  "3.202583"  "2.078081"  "2.988423"  "3.106154"  "2.809472"  "2.548482"  "3.056357"  "3.027847"  "3.044203"  "2.890642"  "2.633544" 
# cpd00003 "0.5465259" "0.5040640" "0.5370572" "0.6370613" "0.4906686" "0.5455284" "0.5234037" "0.4817143" "0.5773655" "0.5263570" "0.5622112" "0.5548784" "0.5020394"
# cpd00004 "0.4748061" "0.4577167" "0.4688221" "0.6022029" "0.4340136" "0.4862296" "0.4703171" "0.4308653" "0.5199245" "0.4694344" "0.5064305" "0.4990025" "0.4596358"
# cpd00005 "0.6432422" "0.5016121" "0.5935540" "0.4735150" "0.5758568" "0.6145899" "0.5330685" "0.5004035" "0.6432077" "0.5756384" "0.5971796" "0.6007156" "0.5327750"
# etc ...

samp_names <- df.wide[1, ]
tax_names <- row.names(df.wide[-1, ])
head(tax_names) # "cpd00001" "cpd00002" "cpd00003" "cpd00004" "cpd00005" "cpd00006"
otu.df <- df.wide[-1, ] # remove sample labels in 1st row
# this is necessary to create numeric matrix

colnames(otu.df) <- samp_names

# convert OTU table to matrix
class(otu.df) # "matrix" "array" 
#otu.df <- as.matrix(otu.df)

# convert to numeric matrix
# https://stackoverflow.com/questions/20791877/convert-character-matrix-into-numeric-matrix
otu.df <- apply(otu.df, 2, as.numeric)

rownames(otu.df) # NULL
dim(otu.df) #  7024   76
rownames(otu.df) <- tax_names

## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
OTU <- otu_table(otu.df, taxa_are_rows = TRUE)

# # convert Taxonomy table to matrix  
tax <- data.frame(cpd_id = tax_names)
row.names(tax) <- tax_names

tax <- as.matrix(tax)

identical( row.names(otu.df), row.names(tax) ) # TRUE


## Create 'taxonomyTable'
#  tax_table - Works on any character matrix.
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
TAX <- tax_table(tax)


## Create a phyloseq object, merging OTU & TAX tables
phy.cpp = phyloseq(OTU, TAX)
phy.cpp
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7024 taxa and 76 samples ]
# tax_table()   Taxonomy Table:    [ 7024 taxa by 1 taxonomic ranks ]


sample_names(phy.cpp)
# [1] "ERR260139" "ERR260140" "ERR260144" "ERR260147" ... etc.

identical(sample_names(phy.cpp), samp$sample) # TRUE

# row.names need to match sample_names() from phyloseq object
#row.names(samp) <- samp$sample
identical(row.names(samp), samp$sample) # TRUE


### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object
SAMP <- sample_data(samp)


### Combine SAMPDATA into phyloseq object
phy.cpp <- merge_phyloseq(phy.cpp, SAMP)
phy.cpp
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7024 taxa and 76 samples ]
# sample_data() Sample Data:       [ 76 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 7024 taxa by 1 taxonomic ranks ]

min(taxa_sums(phy.cpp)) # 5.040179e-07

phy.cpp@sam_data
# Sample Data:        [76 samples by 2 sample variables]:
# sample group_label
# ERR260139 ERR260139    T2D met-
# ERR260140 ERR260140    T2D met-
# ERR260144 ERR260144    T2D met-
# ERR260147 ERR260147      Normal
# ERR260151 ERR260151    T2D met-
# ERR260152 ERR260152    T2D met-
# ERR260153 ERR260153      Normal
# etc. ...


saveRDS(object = phy.cpp, file = "phy.cpp-cleaned-Forslund-SWE-T2D-Hostremoval-v6.RDS")

phy_in <- phy.cpp

sum(sample_sums(phy_in)) # 5359.577
sort(sample_sums(phy_in))
# ERR260147 ERR260265 ERR260234 ERR260227 ERR260221 ERR260167 ERR260140 ERR260216 ERR260230 ERR260203 ERR260159 ERR260218 ERR260250 ERR260170 ERR260204 ERR260246 
# 59.85823  63.87817  64.50890  65.04760  65.42071  65.58304  66.13222  66.21866  66.36343  67.35876  67.44448  67.73393  67.87818  68.20336  68.89130  68.89405 
# ERR260255 ERR260166 ERR260181 ERR260264 ERR260276 ERR260251 ERR260256 ERR260153 ERR260242 ERR260231 ERR260171 ERR260151 ERR260258 ERR260268 ERR260252 ERR260224 
# 68.89614  68.90639  68.91962  68.92277  69.08120  69.14342  69.15113  69.25359  69.65289  69.69328  69.79413  69.94410  70.07215  70.10489  70.18464  70.25143 
# ERR260169 ERR260190 ERR260225 ERR260217 ERR260193 ERR260223 ERR260259 ERR260273 ERR260174 ERR260199 ERR260165 ERR260260 ERR260241 ERR260263 ERR260201 ERR260163 
# 70.35783  70.56365  70.58628  70.92997  70.96271  71.08933  71.10379  71.30520  71.38072  71.43197  71.43581  71.46987  71.62479  71.74410  71.80095  71.92468 
# ERR275252 ERR260198 ERR260186 ERR260206 ERR260210 ERR260189 ERR260162 ERR260244 ERR260253 ERR260188 ERR260209 ERR260152 ERR260266 ERR260267 ERR260179 ERR260205 
# 72.03489  72.13153  72.14033  72.18534  72.20040  72.21808  72.41453  72.46792  72.53755  72.54218  72.73864  72.82765  72.99417  73.03861  73.07521  73.15618 
# ERR260243 ERR260185 ERR260215 ERR260226 ERR260161 ERR260173 ERR260207 ERR260175 ERR260144 ERR260139 ERR260271 ERR260180 
# 73.15886  73.20355  73.32795  73.37034  73.73351  73.84853  73.96212  74.04600  74.31994  74.73856  74.76905  75.29727 

summary( sample_sums(phy_in) )
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 59.86   68.92   71.10   70.52   72.54   75.30 

mean( sample_sums(phy_in) )
# 70.52075
sd( sample_sums(phy_in) )
# 2.881688

max(taxa_sums(phy_in)) # 499.8058

fxns <- as.data.frame( phy_in@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 76
NonZeroFxns

mean(NonZeroFxns) # 5287.329
sd(NonZeroFxns) # 584.9439


# don't rarefy - already in form of relative abundance %

table(phy_in@sam_data$group_label)
# T2D met-   Normal 
# 33       43 



## ordination plot
## PCoA + Bray-Curtis

set.seed(123)
ord <- ordinate(phy_in, "PCoA", "bray")

ord

unique(phy_in@sam_data$group_label)
# [1] T2D met- Normal  
# Levels: T2D met- < Normal

p <- plot_ordination(phy_in, ord, type="samples", color="group_label")
p

p$labels$x # "Axis.1   [40.6%]"
x_lab <- "PCo1 (40.6%)"

p$labels$y # "Axis.2   [22.9%]"
y_lab <- "PCo2 (22.9%)"

40.6 + 22.9 # 63.5

#temp <- r1.ps
p_df <- p$data

cols.group <- c("T2D met-" = "#af8dc3",
                "Normal" = "#7fbf7b")

p <- #plot_ordination(temp, ord, type="samples", color="group") +
  #ggplot(data = p_df, aes(x = NMDS1, y = NMDS2, color = group))+
  ggplot(data = p_df, aes(x = Axis.1, y = Axis.2, color = group_label))+
  theme_bw()+
  geom_point()+
  
  stat_ellipse(linetype = "dashed")+
  
  xlab(x_lab) + ylab(y_lab)+
  scale_color_manual(values = cols.group, name = "Diagnosis") +
  #annotate(geom="text", x= -1.2, y= 1.3, label = paste0("Stress = ",round(ord$stress,digits=4)),size = 3, hjust=0, vjust=1) +
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = "T2D: CPP", size = 3 )+
  theme(
    panel.grid.major = element_blank(),
    
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.85)),
    legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = -8, unit = "pt"), # randomForest::margin() interferes !!
    
    panel.grid.minor = element_blank())
p

grid.text(label = "(d)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )

dev.print(tiff, file = paste0(workdir,"/plots/","PCoA-Cpp3d-Forslund-SWE-T2D-Hostremoval-v6.tiff"), width = 9.5, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## PERMANOVA

# Calculate bray curtis distance matrix
set.seed(123)
bray <- phyloseq::distance(phy_in, method = "bray")
sampledf <- data.frame(sample_data(phy_in))
str(sampledf)

names(phy_in@sam_data)
# "sample" "group_label" 

# Adonis test
set.seed(123)
adonis2(bray ~ group_label , data = sampledf)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ group_label, data = sampledf)
#             Df SumOfSqs      R2      F Pr(>F)
# group_label  1  0.00814 0.01959 1.4784  0.175
# Residual    74  0.40763 0.98041              
# Total       75  0.41578 1.00000 

beta <- betadisper(bray, sampledf$group_label)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#           Df   Sum Sq    Mean Sq     F N.Perm Pr(>F)
# Groups     1 0.000001 0.00000066 5e-04    999  0.981
# Residuals 74 0.088995 0.00120264 



## Functions

phy_in <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval.RDS")

phy_in
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 17912 taxa and 76 samples ]
# sample_data() Sample Data:       [ 76 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 17912 taxa by 4 taxonomic ranks ]

sample_sums(phy_in)
# all 100

summary( sample_sums(phy_in) )

max(taxa_sums(phy_in)) # 233.1107


# do not rarefy as already normalized using rel abun (%)

table(phy_in@sam_data$group_label)
# T2D met-   Normal 
# 33       43 



## ordination plot
## PCoA + Bray-Curtis

set.seed(123)
ord <- ordinate(phy_in, "PCoA", "bray")

ord

p <- plot_ordination(phy_in, ord, type="samples", color="group_label")
p

p$labels$x # "Axis.1   [28.2%]"
x_lab <- "PCo1 (28.2%)"

p$labels$y # "Axis.2   [21.4%]"
y_lab <- "PCo2 (21.4%)"

p_df <- p$data

cols.group <- c("T2D met-" = "#af8dc3",
                "Normal" = "#7fbf7b")

p <- #plot_ordination(temp, ord, type="samples", color="group") +
  ggplot(data = p_df, aes(x = Axis.1, y = Axis.2, color = group_label))+
  theme_bw()+
  geom_point()+
  
  stat_ellipse(linetype = "dashed")+
  
  xlab(x_lab) + ylab(y_lab)+
  scale_color_manual(values = cols.group, name = "Diagnosis") +
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = "T2D: Fxns", size = 3 )+
  theme(
    panel.grid.major = element_blank(),
    
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.85)),
    legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = -8, unit = "pt"), # randomForest::margin() interferes !!
    
    panel.grid.minor = element_blank())
p

grid.text(label = "(c)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","PCoA-Functions-Forslund-SWE-T2D-Hostremoval-v6.tiff"), width = 9.5, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## PERMANOVA

# Calculate bray curtis distance matrix
set.seed(123)
bray <- phyloseq::distance(phy_in, method = "bray")
sampledf <- data.frame(sample_data(phy_in))
str(sampledf)
names(phy_in@sam_data)
#  "Sample"         "Country subset" "Status"         "Bases"          "Run"            "group_label"    "age"  

# Adonis test
set.seed(123)
adonis2(bray ~ group_label , data = sampledf)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ group_label, data = sampledf)
#             Df SumOfSqs      R2      F Pr(>F)
# group_label  1   0.0530 0.02051 1.5493  0.136
# Residual    74   2.5313 0.97949              
# Total       75   2.5843 1.00000 


beta <- betadisper(bray, sampledf$group_label)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.00087 0.0008726 0.1575    999  0.709
# Residuals 74 0.41000 0.0055405 

#-------------------------

#### Wilcoxon-Mann-Whitney Tests - Forslund-SWE-T2D - w/ Host-removal - CPP differences for ALL compounds
#    Test sig diff Normal vs T2D - consider all p <= 0.05
#-------------------------

dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-v6.rds")
str(dat.cpd.collate)
# 'data.frame':	533824 obs. of  5 variables:
# $ cpd_id      : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 4.982389 0 0 0.000188 ...
# $ log10_abun  : num  -8.187 0.697 -8.187 -8.187 -3.725 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

unique(dat.cpd.collate$group_label)
# [1] T2D met- Normal  
# Levels: T2D met- < Normal

data_in <- dat.cpd.collate

dat.test <- data.frame(cpd = unique(data_in$cpd_id), data_for_this_cpd=NA , 
                       median_t2d = NA, median_normal = NA,  
                       alt = NA,
                       p_val = NA, W_statistic = NA, hl_effect_wilcox = NA,
                       trend_with_disease = NA
)

for (i in 1:dim(dat.test)[1]) {
  #i<-1
  this_cpd <- dat.test$cpd[i]
  sel <- which(data_in$cpd_id == this_cpd)
  # prepare data
  df = data.frame(group = as.character( data_in$group[sel]), 
                  value = data_in$log10_abun[sel], 
                  value_perc = data_in$cpd_rel_abun[sel] )
  
  x <- df$value[df$group == "T2D met-"]
  y <- df$value[df$group == "Normal"]
  
  if ( length(which(x == min(x) )) > 0.5*length(x) & length(which(y == min(y) )) > 0.5*length(y) ) {
    # low data: if at least one dataset does not have at least 50% non-zero cases
    dat.test$data_for_this_cpd[i] <- "low data"
    
  } else {
    
    # # test for homogeneity of variances
    # var.test(x,y) # e.g., j<-100: F = 2.3477, num df = 32, denom df = 42, p-value = 0.009872
    
    dat.test$median_t2d[i] <- median( x ) # log10 abun  
    dat.test$median_normal[i] <- median( y ) # log10 abun
    
    alt <- NA
    if (median( x ) > median( y ) ) {
      alt <- "greater"
    } else {
      alt <- "less"
    }
    dat.test$alt[i] <- alt
    
    # Wilcoxon-Mann-Whitney Test
    wmw.test <- wilcox.test(x, y, alternative = alt, paired = FALSE, conf.int = TRUE) # based on log10 abun
    
    dat.test$p_val[i] <- wmw.test$p.value
    dat.test$W_statistic[i] <- wmw.test$statistic
    
    # https://search.r-project.org/CRAN/refmans/DescTools/html/HodgesLehmann.html
    # https://aakinshin.net/posts/r-hodges-lehmann-problems
    
    dat.test$hl_effect_wilcox[i] <- wmw.test$estimate
    #dat.test$hl_effect_hlfxn[i] <- hl(x, y)
    
    if (!(is.na(dat.test$p_val[i])|is.na(dat.test$W_statistic[i]))) {
      if (dat.test$p_val[i] <= 0.05 & alt == "greater") { dat.test$trend_with_disease[i] <- "Increasing" }
      else if (dat.test$p_val[i] <= 0.05 & alt == "less") { dat.test$trend_with_disease[i] <- "Decreasing" }
      else { dat.test$trend_with_disease[i] <- "No trend" }
    }
  }
  print(paste0("Completed ",i))
}


sel.low <- which(dat.test$data_for_this_cpd == "low data") # 1474

length(unique(dat.cpd.collate$cpd_id)) # 7024
length(unique(dat.cpd.collate$cpd_id)) - length(sel.low) # 5550

sel.nonNA <- which(!is.na(dat.test$p_val)) # 5550 applicable tests

dat.test %>% arrange(., p_val) %>% slice(1:50)


write.csv(x = dat.test, file = "dat.test-cpp3d-indiv-Compounds-Forslund-SWE-T2D-Hostremoval-All-Compounds-SELECT-T2D-VS-NORM--Wilcox-v6.csv")

saveRDS(dat.test, file = "dat.test-cpp3d-indiv-Compounds-Forslund-SWE-T2D-Hostremoval-All-Compounds-SELECT-T2D-VS-NORM--Wilcox-v6.RDS")

# extract sig results ... (no p-adjustment)

sel.sig <- which(dat.test$p_val <= 0.05) # 1269

dat.test.sig <- dat.test[sel.sig, ]

dat.test.sig$minuslog10_p_val <- -log10(dat.test.sig$p_val)

plot(x = dat.test.sig$W_statistic , y =dat.test.sig$minuslog10_p_val , xlab="W statistic", ylab="-log10(P-value)")

dev.print(tiff, filename = paste0(workdir,"/plots/","cpp3d-indiv-compound-VolcanoPlot-P-values--All-Compounds-SELECT-T2D-SWE-vs-NORM--Wilcox-v6.tiff"),
          width = 12, height = 14, units = "cm", res=600, compression = "lzw",type="cairo" )


# join compound info??

dat.test.sig$cpd_names <- NA
dat.test.sig$cpd_forms <- NA

dat.test.sig$OC_x <- NA
dat.test.sig$HC_y <- NA
dat.test.sig$NC_z <- NA

for (i in 1:dim(dat.test.sig)[1]) {
  #i<-1
  this_cpd <- dat.test.sig$cpd[i]

  sel.cpd <- which(df.comp2$id == this_cpd)

  dat.test.sig$cpd_names[i] <- df.comp2$name[sel.cpd]
  dat.test.sig$cpd_forms[i] <- df.comp2$form[sel.cpd]

  dat.test.sig$OC_x[i] <- df.comp2$OC_ratio[sel.cpd]
  dat.test.sig$HC_y[i] <- df.comp2$HC_ratio[sel.cpd]
  dat.test.sig$NC_z[i] <- df.comp2$NC_ratio[sel.cpd]

  print(paste0("completed ",i))
}

write.table(x = dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D-Hostremoval--All-Compounds-SELECT-T2D-vs-NORM--Wilcox-v6.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D-Hostremoval--All-Compounds-SELECT-T2D-vs-NORM--Wilcox-v6.rds")

hist(dat.test.sig$NC_z); summary(dat.test.sig$NC_z)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.0000  0.0000  0.1033  0.1538  1.0000      88 

dim(dat.test.sig) # 1269   17
sel.ok <- which(!is.na(dat.test.sig$NC_z) ) # qty 1181

dat.test.sig$z_layer <- NA

subsel <- which(dat.test.sig$NC_z[sel.ok] == 0) # 639
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C = 0"

subsel <- which(dat.test.sig$NC_z[sel.ok] > 0 & dat.test.sig$NC_z[sel.ok] <= 0.2 ) # 311
# Rivas-Ubach et al 2018 Table 1 highlights N/C breaks at 0.126, 0.2, 0.5, 0.7
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C >0 to 0.2"

subsel <- which(dat.test.sig$NC_z[sel.ok] > 0.2) # 231
max(dat.test.sig$NC_z[sel.ok]) # 1
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C >0.2 to 1"

unique(dat.test.sig$z_layer[sel.ok]) # "N:C = 0"       "N:C >0 to 0.2" "N:C >0.2 to 1"

dat.test.sig$z_layer <- factor(dat.test.sig$z_layer, levels = c("N:C = 0",
                                                                "N:C >0 to 0.2",
                                                                "N:C >0.2 to 1"), ordered = TRUE)

saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D-Hostremoval--All-Compounds-SELECT-T2D-vs-NORM--Wilcox-v6.rds")

dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D-Hostremoval--All-Compounds-SELECT-T2D-vs-NORM--Wilcox-v6.rds")


dim(dat.test.sig) #  1269   18
head(dat.test.sig)
dim(dat.test.sig[ which(dat.test.sig$p_val < 0.05), ]) # 1269   18
sel <- which(!is.na(dat.test.sig$OC_x) ) # 1181
sel <- which(dat.test.sig$trend_with_disease == "Decreasing") # 911
sel <- which(dat.test.sig$trend_with_disease == "Decreasing" & !is.na(dat.test.sig$OC_x) ) # 860 - displayed in O:C, H:C, N:C plot
sel <- which(dat.test.sig$trend_with_disease == "Increasing") # 358
sel <- which(dat.test.sig$trend_with_disease == "Increasing" & !is.na(dat.test.sig$OC_x) ) # 321 - displayed in O:C, H:C, N:C plot



## plot as Increasing or Decreasing?? in vK space
## Use adjusted compound classes (adapted from Wu 2018, D'Andrilli, Rivas-Ubach 2018, and Minor et al 2015)

# # Zones and labels as above
# vkgrouprect.facets2 <- read.table(file = "cpp3d-compound-classes.tsv", header = TRUE, sep = "\t" )
# vkgrouprect.facets2.labels <- read.table(file = "cpp3d-compound-classes-labels.tsv", header = TRUE, sep = "\t" )

vkgrouprect.facets2.labels


p <- ggplot(data = dat.test.sig[sel.ok, ]) +
  coord_equal()+
  #ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study")+
  #xlim(0,3.4)+ ylim(0,4.1)+
  xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + #
  xlab("O:C ratio")+ ylab("H:C ratio")+
  
  guides(color = guide_legend(title = "Trend with disease in functional capacity\n(%) allocated to compounds"))+
  
  facet_wrap(facets = vars(z_layer))+
  
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, group = label), color="#969696", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373", "grey"
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, group = label), color="#969696", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  #geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0.2 to 2" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0.2 to 1" ), aes(x = OC_x, y = HC_y, group = label), color="#969696", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 )+ # 
  
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T,  size = 2.75 , col="#252525" , lineheight = 0.8)+ # "#737373"
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2.75 , col="#252525" , lineheight = 0.8)+
  #geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0.2 to 2" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2 , col="#737373" , lineheight = 0.8)+
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0.2 to 1" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2.75 , col="#252525" , lineheight = 0.8)+
  
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.text = element_text(size = rel(0.9)) ,
    title = element_text(size = rel(0.9), lineheight = 1),
    strip.background = element_rect(fill = "transparent")
  )

p

dev.print(tiff, file = paste0(workdir,"/plots/","cpp3d-Compounds-indiv-vKSpace-Trend-with-Disease--All-Compounds-SELECT-T2D-SWE-Hostremoval-vs-NORM--Wilcox-v6-axis-lim.tiff"), width = 20, height = 12, units = "cm", res=500, compression="lzw",type="cairo")

#-------------------------

#### Wilcoxon-Mann-Whitney Tests - Forslund-SWE-T2D - w/ Host-removal 
#    Test sig diff Normal vs T2D
#    FOCUS ONLY ON B-H ADJUSTED SIG DIFF COMPOUNDS SELECTED FROM ECOSYSTEM RESTORATION
#-------------------------

#saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto.rds")
dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto.rds")
max(dat.test.sig$p_val) # 0.04552315

str(dat.test.sig)
# 'data.frame':	2958 obs. of  14 variables:
# $ cpd              : chr  "cpd25681" "cpd24620" "cpd00001" "cpd00851" ...
# $ data_for_this_cpd: chr  NA NA NA NA ...
# $ p_val            : num  0.000575 0.015075 0.015075 0.006249 0.006249 ...
# $ kendall_tau      : num  -3.44 -2.43 -2.43 -2.73 -2.73 ...
# $ trend_with_age   : chr  "Decreasing" "Decreasing" "Decreasing" "Decreasing" ...
# $ sigBH            : chr  "sig" NA NA "sig" ...
# $ minuslog10_p_val : num  3.24 1.82 1.82 2.2 2.2 ...
# $ cpd_names        : chr  "2-methyl-trans-aconitate" "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" "trans-4-Hydroxy-L-proline" ...
# $ cpd_forms        : chr  "C7H5O6" "C7H7O7" "H2O" "C5H9NO3" ...
# $ OC_x             : num  0.857 1 NA 0.6 0.6 ...
# $ HC_y             : num  0.714 1 NA 1.8 1.8 ...
# $ NC_z             : num  0 0 NaN 0.2 0.2 ...
# $ z_layer          : Ord.factor w/ 3 levels "N:C = 0"<"N:C >0 to 0.2"<..: 1 1 NA 2 2 3 3 3 3 NA ...
# $ mass             : num  188 206 18 130 130 89 175 147 147 1 ...

sel.sigBH <- which( dat.test.sig$sigBH == "sig") # 2122

resto.sig.cpds <- unique(dat.test.sig$cpd[sel.sigBH]) # qty 2122

dat.test.sig.resto <- dat.test.sig

saveRDS(resto.sig.cpds, file = "resto.sig.cpds-sigBH-from-sunbad-resto-v2.rds")


# compare to data from T2D-SWE
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-v6.rds")

str(dat.cpd.collate)
# 'data.frame':	533824 obs. of  5 variables:
# $ cpd_id      : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 4.982389 0 0 0.000188 ...
# $ log10_abun  : num  -8.187 0.697 -8.187 -8.187 -3.725 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

sel <- which(dat.cpd.collate$cpd_id %in% resto.sig.cpds)
dat.cpd.collate <- dat.cpd.collate[sel, ]

length( unique(dat.cpd.collate$cpd_id) ) # 1873
table(dat.cpd.collate$group)/1873
# T2D met-   Normal 
# 33       43 

length( unique(dat.cpd.collate$sample[ dat.cpd.collate$group == "T2D met-" ]) ) # 33
length( unique(dat.cpd.collate$sample[ dat.cpd.collate$group == "Normal" ]) ) # 43

data_in <- dat.cpd.collate

dat.test <- data.frame(cpd = unique(data_in$cpd_id), data_for_this_cpd=NA , 
                       median_t2d = NA, median_normal = NA,  
                       alt = NA,
                       p_val = NA, W_statistic = NA, hl_effect_wilcox = NA, #hl_effect_hlfxn = NA,
                       trend_with_disease = NA, 
                       
                       incResto_incT2D = NA, decResto_decT2D = NA, # same direction
                       incResto_decT2D = NA, decResto_incT2D = NA, # opposing direction
                       incResto_notT2D = NA, decResto_notT2D = NA,  # no trend in T2D
                       notResto_incT2D = NA, notResto_decT2D = NA, notResto_notT2D = NA, # not significant in resto
                       anyResto_incT2D = NA, anyResto_decT2D = NA, anyResto_notT2D = NA # any/all ecosystem scenarios
)

# # define fxn: Hodges-Lehmann estimator of effect size - https://aakinshin.net/posts/r-hodges-lehmann-problems/
# # this is same result as wmw.test$estimate
# hl <- function(x, y = NULL) {
#   if (is.null(y)) {
#     walsh <- outer(x, x, "+") / 2
#     median(walsh[lower.tri(walsh, diag = TRUE)])
#   } else {
#     median(outer(x, y, "-"))
#   }
# }



for (i in 1:dim(dat.test)[1]) {
  #i<-100 i<-1
  this_cpd <- dat.test$cpd[i]
  sel <- which(data_in$cpd_id == this_cpd)
  # prepare data
  df = data.frame(group = as.character( data_in$group[sel]), 
                  value = data_in$log10_abun[sel], 
                  value_perc = data_in$cpd_rel_abun[sel] )
  
  x <- df$value[df$group == "T2D met-"]
  y <- df$value[df$group == "Normal"]
  
  if ( length(which(x == min(x) )) > 0.5*length(x) & length(which(y == min(y) )) > 0.5*length(y) ) {
    # low data: if at least one dataset does not have at least 50% non-zero cases
    dat.test$data_for_this_cpd[i] <- "low data"
    
  } else {
    
    # # test for homogeneity of variances
    # var.test(x,y) # e.g., j<-100: F = 2.3477, num df = 32, denom df = 42, p-value = 0.009872
    
    dat.test$median_t2d[i] <- median( x ) # log10 abun  
    dat.test$median_normal[i] <- median( y ) # log10 abun
    
    alt <- NA
    if (median( x ) > median( y ) ) {
      alt <- "greater"
    } else {
      alt <- "less"
    }
    dat.test$alt[i] <- alt
    
    # Wilcoxon-Mann-Whitney Test
    wmw.test <- wilcox.test(x, y, alternative = alt, paired = FALSE, conf.int = TRUE) # based on log10 abun
    
    dat.test$p_val[i] <- wmw.test$p.value
    dat.test$W_statistic[i] <- wmw.test$statistic
    
    # https://search.r-project.org/CRAN/refmans/DescTools/html/HodgesLehmann.html
    # https://aakinshin.net/posts/r-hodges-lehmann-problems
    
    dat.test$hl_effect_wilcox[i] <- wmw.test$estimate
    #dat.test$hl_effect_hlfxn[i] <- hl(x, y)
    
    if (!(is.na(dat.test$p_val[i])|is.na(dat.test$W_statistic[i]))) {
      if (dat.test$p_val[i] <= 0.05 & alt == "greater") { dat.test$trend_with_disease[i] <- "Increasing" }
      else if (dat.test$p_val[i] <= 0.05 & alt == "less") { dat.test$trend_with_disease[i] <- "Decreasing" }
      else { dat.test$trend_with_disease[i] <- "No trend" }
    }
    
    # detect if same or opposing pattern to restoration?
    
    # trend with restoration?
    sel.resto <- which( dat.resto$cpd == this_cpd ) # which( dat.resto$cpd == "blahblah" )
    #dat.test.sig.resto[sel.resto, ]
    
    if (length(sel.resto) == 0) {
      
      # not significant in restoration
      if ( dat.test$trend_with_disease[i] == "Increasing") { dat.test$notResto_incT2D[i] <- 1 }
      if ( dat.test$trend_with_disease[i] == "Decreasing") { dat.test$notResto_decT2D[i] <- 1 }
      if ( dat.test$trend_with_disease[i] == "No trend") { dat.test$notResto_notT2D[i] <- 1 }
      
    } else {
      
      # same direction
      if ( dat.resto$trend_with_age[sel.resto] == "Increasing" & dat.test$trend_with_disease[i] == "Increasing") { dat.test$incResto_incT2D[i] <- 1 }
      if ( dat.resto$trend_with_age[sel.resto] == "Decreasing" & dat.test$trend_with_disease[i] == "Decreasing") { dat.test$decResto_decT2D[i] <- 1 }
      
      # opposing direction
      if ( dat.resto$trend_with_age[sel.resto] == "Increasing" & dat.test$trend_with_disease[i] == "Decreasing") { dat.test$incResto_decT2D[i] <- 1 }
      if ( dat.resto$trend_with_age[sel.resto] == "Decreasing" & dat.test$trend_with_disease[i] == "Increasing") { dat.test$decResto_incT2D[i] <- 1 }
      
      # no trend with T2D
      if ( dat.resto$trend_with_age[sel.resto] == "Increasing" & dat.test$trend_with_disease[i] == "No trend") { dat.test$incResto_notT2D[i] <- 1 }
      if ( dat.resto$trend_with_age[sel.resto] == "Decreasing" & dat.test$trend_with_disease[i] == "No trend") { dat.test$decResto_notT2D[i] <- 1 }
      
    }
    
    # now capture any/all ecosystem scenarios
    if ( dat.test$trend_with_disease[i] == "Increasing") { dat.test$anyResto_incT2D[i] <- 1 }
    if ( dat.test$trend_with_disease[i] == "Decreasing") { dat.test$anyResto_decT2D[i] <- 1 }
    if ( dat.test$trend_with_disease[i] == "No trend") { dat.test$anyResto_notT2D[i] <- 1 }
  
  }
  print(paste0("Completed i ",i))
  
} # END i for each compound


summary(dat.test$hl_effect_wilcox)
plot(dat.test$W_statistic, dat.test$p_val)


sel.low <- which(dat.test$data_for_this_cpd == "low data") # 219

length(unique(dat.cpd.collate$cpd_id)) # 1873
length(unique(dat.cpd.collate$cpd_id)) - length(sel.low) # 1654

sel.nonNA <- which(!is.na(dat.test$p_val)) # 1654 applicable tests


saveRDS(dat.test, file = "dat.test-cpp3d-indiv-Compounds-Forslund-SWE-T2D-Hostremoval--RESTO-SIGBH-SELECT-T2D-vs-NORM--Wilcox-v6.RDS")


sel.sig <- which(dat.test$p_val <= 0.05) # 304


summary(dat.test$p_val)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.00080 0.08939 0.20201 0.23099 0.35671 0.87469     219 

hist(dat.test$p_val)

names(dat.test)


sum(dat.test[ ,c("incResto_incT2D", 
                 "decResto_decT2D", 
                 "incResto_decT2D", 
                 "decResto_incT2D", 
                 "incResto_notT2D", 
                 "decResto_notT2D", 
                 "notResto_incT2D", 
                 "notResto_decT2D")], na.rm = TRUE) # 1654

sum(dat.test[ ,c("incResto_incT2D", 
                 "decResto_decT2D", 
                 "incResto_decT2D", 
                 "decResto_incT2D", 
                 #"incResto_notT2D", 
                 #"decResto_notT2D", 
                 "notResto_incT2D", 
                 "notResto_decT2D")], na.rm = TRUE) # 304 - this number had trends in T2D and across restoration

sum(dat.test[ ,c(#"incResto_incT2D", 
                 #"decResto_decT2D", 
                 #"incResto_decT2D", 
                 #"decResto_incT2D", 
                 #"incResto_notT2D", 
                 #"decResto_notT2D", 
                 "notResto_incT2D", 
                 "notResto_decT2D"
                 )], na.rm = TRUE) # 0

sum(dat.test[ ,c("incResto_incT2D", 
                 "decResto_decT2D", 
                 "incResto_decT2D", 
                 "decResto_incT2D" #, 
                 #"incResto_notT2D", 
                 #"decResto_notT2D", 
                 #"notResto_incT2D", 
                 #"notResto_decT2D"
                 )], na.rm = TRUE) # 304

sum(dat.test[ ,c("incResto_incT2D", 
                 #"decResto_decT2D", 
                 #"incResto_decT2D", 
                 "decResto_incT2D" #, 
                 #"incResto_notT2D", 
                 #"decResto_notT2D", 
                 #"notResto_incT2D", 
                 #"notResto_decT2D"
)], na.rm = TRUE) # 139

sum(dat.test[ ,c(#"incResto_incT2D", 
                 "decResto_decT2D", 
                 "incResto_decT2D" #, 
                 #"decResto_incT2D" #, 
                 #"incResto_notT2D", 
                 #"decResto_notT2D", 
                 #"notResto_incT2D", 
                 #"notResto_decT2D"
)], na.rm = TRUE) # 165




# extract sig results

sel.sig <- which(dat.test$p_val <= 0.05) # 304

dat.test.sig <- dat.test[sel.sig, ]

dat.test.sig$minuslog10_p_val <- -log10(dat.test.sig$p_val)

plot(x = dat.test.sig$W_statistic , y =dat.test.sig$minuslog10_p_val , xlab="W statistic", ylab="-log10(P-value)")

dev.print(tiff, filename = paste0(workdir,"/plots/","cpp3d-indiv-compound-VolcanoPlot-P-values--RESTO-SIGBH-SELECT-SWE-T2D-Hostremoval-vs-NORM--Wilcox-v6.tiff"),
          width = 12, height = 14, units = "cm", res=600, compression = "lzw",type="cairo" )


# join compound info??

#length(which(dat.test.sig$sigBH=="sig")) # 2122 1953


dat.test.sig$cpd_names <- NA
dat.test.sig$cpd_forms <- NA

dat.test.sig$OC_x <- NA
dat.test.sig$HC_y <- NA
dat.test.sig$NC_z <- NA

dat.test.sig$trend_group <- NA


for (i in 1:dim(dat.test.sig)[1]) {
  #i<-1
  this_cpd <- dat.test.sig$cpd[i]
  
  sel.cpd <- which(df.comp$id == this_cpd)
  
  dat.test.sig$cpd_names[i] <- df.comp$name[sel.cpd]
  dat.test.sig$cpd_forms[i] <- df.comp$form[sel.cpd]
  
  dat.test.sig$OC_x[i] <- df.comp$OC_ratio[sel.cpd]
  dat.test.sig$HC_y[i] <- df.comp$HC_ratio[sel.cpd]
  dat.test.sig$NC_z[i] <- df.comp$NC_ratio[sel.cpd]
  
  if (!is.na(dat.test.sig$incResto_incT2D[i]) & dat.test.sig$incResto_incT2D[i] == 1) {dat.test.sig$trend_group[i] <- "Increasing in T2D (increased exposure in quality ecosystems)" } # (associated with ecosystem quality)
  if (!is.na(dat.test.sig$decResto_decT2D[i]) & dat.test.sig$decResto_decT2D[i] == 1) {dat.test.sig$trend_group[i] <- "Decreasing in T2D (reduced exposure in quality ecosystems)" } # (associated with ecosystem disturbance)
  if (!is.na(dat.test.sig$incResto_decT2D[i]) & dat.test.sig$incResto_decT2D[i] == 1) {dat.test.sig$trend_group[i] <- "Decreasing in T2D (reduced exposure in degraded ecosystems)" } # (associated with ecosystem quality)
  if (!is.na(dat.test.sig$decResto_incT2D[i]) & dat.test.sig$decResto_incT2D[i] == 1) {dat.test.sig$trend_group[i] <- "Increasing in T2D (increased exposure in degraded ecosystems)" } # (associated with ecosystem disturbance)
  
  print(paste0("completed ",i))
}

write.table(x = dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D-Hostremoval--RESTO-SIGBH-SELECT-T2D-vs-NORM--Wilcox-v6.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D-Hostremoval--RESTO-SIGBH-SELECT-T2D-vs-NORM--Wilcox-v6.rds")

#dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D-Hostremoval--RESTO-SIGBH-SELECT-T2D-vs-NORM--Wilcox-v6.rds")


## plot as Increasing or Decreasing?? in vK space

names(dat.test.sig)
# [1] "cpd"                "data_for_this_cpd"  "median_t2d"         "median_normal"      "alt"                "p_val"              "W_statistic"        "hl_effect_wilcox"  
# [9] "trend_with_disease" "incResto_incT2D"    "decResto_decT2D"    "incResto_decT2D"    "decResto_incT2D"    "incResto_notT2D"    "decResto_notT2D"    "notResto_incT2D"   
# [17] "notResto_decT2D"    "notResto_notT2D"    "anyResto_incT2D"    "anyResto_decT2D"    "anyResto_notT2D"    "minuslog10_p_val"   "cpd_names"          "cpd_forms"         
# [25] "OC_x"               "HC_y"               "NC_z"               "trend_group"

unique( dat.test.sig$trend_group )
# [1] "Decreasing in T2D (reduced exposure in quality ecosystems)"    "Increasing in T2D (increased exposure in degraded ecosystems)"
# [3] "Decreasing in T2D (reduced exposure in degraded ecosystems)"   "Increasing in T2D (increased exposure in quality ecosystems)" 


dat.test.sig$trend_group <- factor(dat.test.sig$trend_group,
                                   levels = c("Decreasing in T2D (reduced exposure in quality ecosystems)",
                                              "Decreasing in T2D (reduced exposure in degraded ecosystems)",
                                              "Increasing in T2D (increased exposure in degraded ecosystems)",
                                              "Increasing in T2D (increased exposure in quality ecosystems)" ),
                                   ordered = TRUE)

col.trend_group <- c("Decreasing in T2D (reduced exposure in quality ecosystems)" = "#fdb462", 
                     "Decreasing in T2D (reduced exposure in degraded ecosystems)" = "#e31a1c",
                     "Increasing in T2D (increased exposure in degraded ecosystems)" = "#33a02c" ,
                     "Increasing in T2D (increased exposure in quality ecosystems)" = "#1f78b4" )

p <- #ggplot(data = filter(dat.test.sig, sigBH == "sig")) +
  ggplot(data = dat.test.sig) +
  coord_equal()+
  ggtitle("Microbiota compound processing potential\n- Type 2 Diabetes case study\n(Restoration-significant compounds)")+
  #xlim(0,3.4)+ ylim(0,4.1)+
  
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + # 
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 1, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  
  xlab("O:C molar ratio")+ ylab("H:C molar ratio")+
  guides(color = guide_legend(title = "Trend with disease\nin functional\ncapacity (%)\nallocated to\ncompounds", nrow = 4))+
  
  geom_mark_rect(data= vkgrouprect, aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ 
  
  annotate(geom="text", x= 0+0.01, y= 2.3+0.02, label = "Lipid", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0.2+0.01, y= 2.2+0.02, label = "Protein", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  #annotate(geom="text", x= 0.52+0.01, y= 2.2+0.02, label = "Amino\nsugar", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0.52+0.01, y= 2.2+0.02, label = "Amino\nsugar", hjust=0, vjust=0, size = 2 , col="#737373", lineheight = 0.8) + # top-left
  annotate(geom="text", x= 0.7+0.01, y= 2.4+0.02, label = "Carbohydrate", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0+0.01, y= 0.5-0.02, label = "Condensed\naromatics", hjust=0, vjust=1, size = 2 , col="#737373", lineheight = 0.8) + # bottom-left
  annotate(geom="text", x= 0.25+0.01, y= 0.75-0.02, label = "Lignin", hjust=0, vjust=1, size = 2 , col="#737373" ) + # bottom-left
  annotate(geom="text", x= 0.67+0.01, y= 0.53-0.02, label = "Tannin", hjust=0, vjust=1, size = 2 , col="#737373") + # bottom-left
  
  theme_bw()+
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.key.spacing.y = unit(0, units = "line"),
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p


hist(dat.test.sig$NC_z)


min(dat.test.sig$NC_z[ dat.test.sig$NC_z > 0 ], na.rm = TRUE) # 0.01449275

dim(dat.test.sig) # 304 28
sel.ok <- which(!is.na(dat.test.sig$NC_z) ) # qty 286

dat.test.sig$z_layer <- NA

subsel <- which(dat.test.sig$NC_z[sel.ok] == 0) # 113
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C = 0"

subsel <- which(dat.test.sig$NC_z[sel.ok] > 0 & dat.test.sig$NC_z[sel.ok] <= 0.2 ) # 105
# Rivas-Ubach et al 2018 Table 1 highlights N/C breaks at 0.126, 0.2, 0.5, 0.7
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C >0 to 0.2"

subsel <- which(dat.test.sig$NC_z[sel.ok] > 0.2) # 68
max(dat.test.sig$NC_z[sel.ok]) # 1
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C >0.2 to 1"

unique(dat.test.sig$z_layer[sel.ok]) # "N:C = 0"       "N:C >0 to 0.2" "N:C >0.2 to 1"

dat.test.sig$z_layer <- factor(dat.test.sig$z_layer, levels = c("N:C = 0",
                                                                "N:C >0 to 0.2",
                                                                "N:C >0.2 to 1"), ordered = TRUE)

saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D-Hostremoval--RESTO-SIGBH-SELECT-T2D-vs-NORM--Wilcox-v6.rds")

dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D-Hostremoval--RESTO-SIGBH-SELECT-T2D-vs-NORM--Wilcox-v6.rds")

dim(dat.test.sig) #  304 29
head(dat.test.sig)
dim(dat.test.sig[ which(dat.test.sig$p_val < 0.05), ]) # 304 29
sel <- which(!is.na(dat.test.sig$OC_x) ) # 286
unique(dat.test.sig$trend_group)
# [1] Decreasing in T2D (reduced exposure in quality ecosystems)    Increasing in T2D (increased exposure in degraded ecosystems)
# [3] Decreasing in T2D (reduced exposure in degraded ecosystems)   Increasing in T2D (increased exposure in quality ecosystems) 
# 4 Levels: Decreasing in T2D (reduced exposure in quality ecosystems) < ... < Increasing in T2D (increased exposure in quality ecosystems)

sel <- which(dat.test.sig$trend_group == "Decreasing in T2D (reduced exposure in degraded ecosystems)") # 70
sel <- which(dat.test.sig$trend_group == "Decreasing in T2D (reduced exposure in degraded ecosystems)" & !is.na(dat.test.sig$OC_x) ) # 66
sel <- which(dat.test.sig$trend_group == "Increasing in T2D (increased exposure in degraded ecosystems)") # 89
sel <- which(dat.test.sig$trend_group == "Increasing in T2D (increased exposure in degraded ecosystems)" & !is.na(dat.test.sig$OC_x) ) # 81

sel <- which(dat.test.sig$trend_group == "Decreasing in T2D (reduced exposure in quality ecosystems)") # 95
sel <- which(dat.test.sig$trend_group == "Decreasing in T2D (reduced exposure in quality ecosystems)" & !is.na(dat.test.sig$OC_x) ) # 90
sel <- which(dat.test.sig$trend_group == "Increasing in T2D (increased exposure in quality ecosystems)") # 50
sel <- which(dat.test.sig$trend_group == "Increasing in T2D (increased exposure in quality ecosystems)" & !is.na(dat.test.sig$OC_x) ) # 49




levels(dat.test.sig$trend_group)


p <- ggplot(data = dat.test.sig[sel.ok, ]) +
  coord_equal()+
  
  #ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study\n(considering only compounds with significant CPP variation in ecosystem restoration)")+
  #ggtitle("Compounds with trending CPP in T2D and soil-ecosystem quality)")+
  
  #xlim(0,3.4)+ ylim(0,4.1)+
  xlim(0,2.6)+ ylim(0,3.1)+
  
  xlab("O:C ratio")+ ylab("H:C ratio")+
  
  facet_wrap(facets = vars(z_layer))+
  
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, group = label), color="#969696", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, group = label), color="#969696", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  #geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0.2 to 2" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0.2 to 1" ), aes(x = OC_x, y = HC_y, group = label), color="#969696", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  
  # mark out groups of interest
  geom_polygon(data = data.frame(x = c(0.96, 1.11, 0.90, 0.75), 
                                 y = c(2.1, 2.05, 1.63, 1.68),
                                 z_layer = rep("N:C = 0", times = 4)), aes(x = x, y = y), fill = NA, color = "green" , alpha = 0.5) +
  geom_text(data = data.frame(x = 1.12, y = 2.1, z_layer = "N:C = 0"), aes(x=x, y=y), label = "Sugars", size = 3.5, color = "darkgreen", hjust=0 )+
  geom_polygon(data = data.frame(x = c(0.28, 0.55, 0.55, 0.20), 
                                 y = c(1.95, 1.85, 1.6, 1.8),
                                 z_layer = rep("N:C >0 to 0.2", times = 4)), aes(x = x, y = y), fill = NA, color = "red" , alpha = 0.5) +
  geom_text(data = data.frame(x = 0.6, y = 1.73, z_layer = "N:C >0 to 0.2"), aes(x=x, y=y), label = "BCFA-ACPs", size = 3.5, color = "darkred", hjust = 0 )+
  
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 )+ # 
  #guides(color = guide_legend(title = "Trend with disease in functional capacity\n(%) allocated to compounds")) +
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 1, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  #guides(color = guide_legend(title = "Trend with disease\nin functional capacity(%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", ncol = 1 )) +
  #guides(color = guide_legend(title = "Trend with disease\nin functional capacity (%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", nrow = 4 ))+ # nrow = 2, ncol = 2
  guides(color = guide_legend(title = "Trend with disease in functional\ncapacity (%) allocated to\ncompounds (and corresponding\ntrend within ecosystems)", nrow = 4 ))+ # nrow = 2, ncol = 2
  
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T,  size = 2.75 , col="#252525" , lineheight = 0.8)+ # size = 2.2 
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2.75 , col="#252525" , lineheight = 0.8)+
  #geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0.2 to 2" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2 , col= "#737373" , lineheight = 0.8)+
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0.2 to 1" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2.75 , col="#252525" , lineheight = 0.8)+
  
  theme_bw()+
  theme(
    
    plot.margin = unit(c(3,3,3,3), unit = "pt"), # t r b l
    #legend.spacing = unit(c(-5,0,0,0), unit = "pt"), # t r b l
    legend.margin = margin(t= -5,r=0,b=0,l=0, unit = "pt"), # t r b l
    
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    
    #legend.key.spacing.y = unit(0, units = "line"),
    legend.key.spacing.y = unit(-0.3, units = "line"),
    
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    #legend.text = element_text(size = rel(0.9), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1),
    strip.background = element_rect(fill = "transparent")
    
  )

p

dev.print(tiff, file = paste0(workdir,"/plots/","cpp3d-Compounds-indiv-vKSpace-TrendGroups-with-Disease--RESTO-SIGBH-SELECT-T2D-SWE-Hostremoval-vs-NORM--Wilcox-axis-lim-v6.tiff"), width = 20, height = 12, units = "cm", res=600, compression="lzw",type="cairo")




## ZOOM-IN
# ZOOM1

p <- ggplot(data = filter( dat.test.sig[sel.ok, ], z_layer == "N:C >0 to 0.2" ) ) +
  #coord_equal()+
  #ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study\n(Focused on compounds significant BH in ecosystem restoration)")+
  
  xlim(0.2, 0.6)+ ylim(1.5,2.0)+ geom_text_repel(aes(x= OC_x, y = HC_y, label = cpd), size = 1.75)+
  #xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + # 
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 1, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  
  xlab("O:C molar ratio")+ ylab("H:C molar ratio")+
  guides(color = guide_legend(title = "Trend with disease\nin functional capacity (%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", nrow = 4 ))+ # nrow = 2, ncol = 2
  
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.key.spacing.y = unit(0, units = "line"),
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    #legend.text = element_text(size = rel(0.9), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","ZOOM1-cpp3d-Compounds-indiv-vKSpace-Trend-with-Disease--RESTO-SIGBH-SELECT-T2D-SWE-Hostremoval-vs-NORM--Wilcox-v6.tiff"), width = 20, height = 16, units = "cm", res=600, compression="lzw",type="cairo")


# ZOOM2

p <- ggplot(data = filter( dat.test.sig[sel.ok, ], z_layer == "N:C = 0" ) ) +
  #coord_equal()+
  #ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study\n(Focused on compounds significant BH in ecosystem restoration)")+
  
  xlim(0.8, 1.1)+ ylim(1.6,2.1)+ geom_text_repel(aes(x= OC_x, y = HC_y, label = cpd), size = 1.75)+
  #xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + # 
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 1, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  
  xlab("O:C molar ratio")+ ylab("H:C molar ratio")+
  guides(color = guide_legend(title = "Trend with disease\nin functional capacity (%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", nrow = 4 ))+ # nrow = 2, ncol = 2
  
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.key.spacing.y = unit(0, units = "line"),
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    #legend.text = element_text(size = rel(0.9), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","ZOOM2-3d-Compounds-indiv-vKSpace-Trend-with-Disease--RESTO-SIGBH-SELECT-SWE-T2D-Hostremoval-vs-NORM--Wilcox-v6.tiff"), width = 20, height = 16, units = "cm", res=600, compression="lzw",type="cairo")



#-------------------------

#### ZOOM plots and examine example cpp3d regions? - w/ Host-removal 
#    Zoom1 - decreasing X-ACPs  (=BCFA-ACPs)
#    Zoom2 - increasing carbs
#-------------------------

dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D-Hostremoval--RESTO-SIGBH-SELECT-T2D-vs-NORM--Wilcox-v6.rds")

sel.ok <- which(!is.na(dat.test.sig$NC_z) ) # qty 286


# build vk layer for each z-layer facet AND age group
for (i in 1:length(levels(dat.test.sig$z_layer))) {
  #i<-1
  temp <- vkgrouprect
  temp$z_layer <- levels(dat.test.sig$z_layer)[i]
  if (i==1) { keep <- temp }
  if (i>1) { keep <- rbind(keep, temp)}
  print(paste0("completed ",i))
}
vkgrouprect.facets <- keep
rm(keep)
str( vkgrouprect.facets )
vkgrouprect.facets$z_layer <- factor( vkgrouprect.facets$z_layer )
vkgrouprect.facets


p <- ggplot(data = dat.test.sig[sel.ok, ]) +
  coord_equal()+
  ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study\n(Focused on compounds significant BH in ecosystem restoration)")+
  #xlim(0,3.4)+ ylim(0,4.1)+
  #xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + # 
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 1, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  
  xlab("O:C molar ratio")+ ylab("H:C molar ratio")+
  guides(color = guide_legend(title = "Trend with disease\nin functional capacity (%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", nrow = 4 ))+ # nrow = 2, ncol = 2
  
  facet_wrap(facets = vars(z_layer))+
  
  #geom_mark_rect(data= vkgrouprect, aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ 
  
  geom_mark_rect(data= filter(vkgrouprect.facets, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets, z_layer == "N:C >0.2 to 1" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  
  annotate(geom="text", x= 0+0.01, y= 2.3+0.02, label = "Lipid", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0.2+0.01, y= 2.2+0.02, label = "Protein", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  #annotate(geom="text", x= 0.52+0.01, y= 2.2+0.02, label = "Amino\nsugar", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0.52+0.01, y= 2.2+0.02, label = "Amino\nsugar", hjust=0, vjust=0, size = 2 , col="#737373", lineheight = 0.8) + # top-left
  annotate(geom="text", x= 0.7+0.01, y= 2.4+0.02, label = "Carbohydrate", hjust=0, vjust=0, size = 2 , col="#737373") + # top-left
  annotate(geom="text", x= 0+0.01, y= 0.5-0.02, label = "Condensed\naromatics", hjust=0, vjust=1, size = 2 , col="#737373", lineheight = 0.8) + # bottom-left
  annotate(geom="text", x= 0.25+0.01, y= 0.75-0.02, label = "Lignin", hjust=0, vjust=1, size = 2 , col="#737373" ) + # bottom-left
  annotate(geom="text", x= 0.67+0.01, y= 0.53-0.02, label = "Tannin", hjust=0, vjust=1, size = 2 , col="#737373") + # bottom-left
  
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.key.spacing.y = unit(0, units = "line"),
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    #legend.text = element_text(size = rel(0.9), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p



## ZOOM-IN - Proteins? N:C >0 to 0.2
# ZOOM1

p <- ggplot(data = filter( dat.test.sig[sel.ok, ], z_layer == "N:C >0 to 0.2" ) ) + # sel.ok[subsel.diab.prot1]. [-subsel.diab.prot2]
  #coord_equal()+
  #ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study\n(Focused on compounds significant BH in ecosystem restoration)")+
  
  xlim(0.2, 0.6)+ ylim(1.5,2.0)+ geom_text_repel(aes(x= OC_x, y = HC_y, label = cpd), size = 1.75)+
  #xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + # 
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 2, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  
  xlab("O:C ratio")+ ylab("H:C ratio")+
  guides(color = guide_legend(title = "Trend with disease\nin functional capacity (%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", nrow = 4 ))+ # nrow = 2, ncol = 2
  
  annotate(geom = "label_npc", npcx = 0.02, npcy = 0.97, label = "N:C >0 to 0.2", size = 3.5)+ # geom = "text_npc"

theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.key.spacing.y = unit(0, units = "line"),
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    #legend.text = element_text(size = rel(0.9), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","ZOOM1-BCFA-ACPs-cpp3d-vKSpace-Trend-with-Disease--RESTO-SIGBH-SELECT-SWE-T2D-Hostremoval-vs-NORM--Wilcox-v6.tiff"), width = 20, height = 16, units = "cm", res=600, compression="lzw",type="cairo")


# check "Proteins" line1 & line2 ??

# examples
sel.diab.prot1 <- which(dat.test.sig$cpd %in% c("cpd11511","cpd11528","cpd11524","cpd11549"))
dat.test.sig[sel.diab.prot1, ]
# cpd data_for_this_cpd median_t2d median_normal  alt      p_val W_statistic hl_effect_wilcox trend_with_disease incResto_incT2D
# 1523 cpd11511              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA
# 1527 cpd11524              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA
# 1529 cpd11528              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA
# 1537 cpd11549              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA
# decResto_decT2D incResto_decT2D decResto_incT2D incResto_notT2D decResto_notT2D notResto_incT2D notResto_decT2D notResto_notT2D
# 1523              NA               1              NA              NA              NA              NA              NA              NA
# 1527              NA               1              NA              NA              NA              NA              NA              NA
# 1529              NA               1              NA              NA              NA              NA              NA              NA
# 1537              NA               1              NA              NA              NA              NA              NA              NA
# anyResto_incT2D anyResto_decT2D anyResto_notT2D minuslog10_p_val                cpd_names     cpd_forms      OC_x     HC_y       NC_z
# 1523              NA               1              NA         1.633983 10-methyl-dodecanoyl-ACP C24H45N2O8PRS 0.3333333 1.875000 0.08333333
# 1527              NA               1              NA         1.633983    5-methyl-hexanoyl-ACP C18H33N2O8PRS 0.4444444 1.833333 0.11111111
# 1529              NA               1              NA         1.633983    7-methyl-octanoyl-ACP C20H37N2O8PRS 0.4000000 1.850000 0.10000000
# 1537              NA               1              NA         1.633983   4-methyl-pentanoyl-ACP C17H31N2O8PRS 0.4705882 1.823529 0.11764706
# trend_group       z_layer
# 1523 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
# 1527 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
# 1529 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
# 1537 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
dat.test.sig$HC_y[sel.diab.prot1]/dat.test.sig$OC_x[sel.diab.prot1]
(1.833333 - 1.875000)/(0.4444444 - 0.3333333) # -0.375

# estimate of line: y = -0.375x + 2

sel.diab.prot1 <- which(dat.test.sig$HC_y == (-0.375*dat.test.sig$OC_x + 2) & dat.test.sig$z_layer == "N:C >0 to 0.2" & dat.test.sig$trend_with_disease == "Decreasing")
# qty 10
# use rounding to 4 d.p
subsel.diab.prot1 <- which(round(dat.test.sig$HC_y[sel.ok], 4) == round(-0.375*dat.test.sig$OC_x[sel.ok] + 2, 4) & dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" & dat.test.sig$trend_with_disease[sel.ok] == "Decreasing")
length(subsel.diab.prot1) # 12

dat.test.sig[sel.ok[subsel.diab.prot1], ]
#           cpd data_for_this_cpd median_t2d median_normal  alt      p_val W_statistic hl_effect_wilcox trend_with_disease incResto_incT2D decResto_decT2D incResto_decT2D decResto_incT2D
# 1517 cpd11499              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA              NA               1              NA
# 1519 cpd11503              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA              NA               1              NA
# 1521 cpd11507              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA              NA               1              NA
# 1523 cpd11511              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA              NA               1              NA
# 1527 cpd11524              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA              NA               1              NA
# 1529 cpd11528              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA              NA               1              NA
# 1531 cpd11532              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA              NA               1              NA
# 1533 cpd11536              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA              NA               1              NA
# 1537 cpd11549              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA              NA               1              NA
# 1539 cpd11553              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA              NA               1              NA
# 1541 cpd11557              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA              NA               1              NA
# 1543 cpd11561              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA              NA               1              NA
# incResto_notT2D decResto_notT2D notResto_incT2D notResto_decT2D notResto_notT2D anyResto_incT2D anyResto_decT2D anyResto_notT2D minuslog10_p_val                cpd_names
# 1517              NA              NA              NA              NA              NA              NA               1              NA         1.633983    4-methyl-hexanoyl-ACP
# 1519              NA              NA              NA              NA              NA              NA               1              NA         1.633983    6-methyl-octanoyl-ACP
# 1521              NA              NA              NA              NA              NA              NA               1              NA         1.633983    8-methyl-decanoyl-ACP
# 1523              NA              NA              NA              NA              NA              NA               1              NA         1.633983 10-methyl-dodecanoyl-ACP
# 1527              NA              NA              NA              NA              NA              NA               1              NA         1.633983    5-methyl-hexanoyl-ACP
# 1529              NA              NA              NA              NA              NA              NA               1              NA         1.633983    7-methyl-octanoyl-ACP
# 1531              NA              NA              NA              NA              NA              NA               1              NA         1.633983    9-methyl-decanoyl-ACP
# 1533              NA              NA              NA              NA              NA              NA               1              NA         1.633983 11-methyl-dodecanoyl-ACP
# 1537              NA              NA              NA              NA              NA              NA               1              NA         1.633983   4-methyl-pentanoyl-ACP
# 1539              NA              NA              NA              NA              NA              NA               1              NA         1.633983   6-methyl-heptanoyl-ACP
# 1541              NA              NA              NA              NA              NA              NA               1              NA         1.633983    8-methyl-nonanoyl-ACP
# 1543              NA              NA              NA              NA              NA              NA               1              NA         1.633983 10-methyl-undecanoyl-ACP
# cpd_forms      OC_x     HC_y       NC_z                                                 trend_group       z_layer
# 1517 C18H33N2O8PRS 0.4444444 1.833333 0.11111111 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
# 1519 C20H37N2O8PRS 0.4000000 1.850000 0.10000000 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
# 1521 C22H41N2O8PRS 0.3636364 1.863636 0.09090909 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
# 1523 C24H45N2O8PRS 0.3333333 1.875000 0.08333333 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
# 1527 C18H33N2O8PRS 0.4444444 1.833333 0.11111111 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
# 1529 C20H37N2O8PRS 0.4000000 1.850000 0.10000000 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
# 1531 C22H41N2O8PRS 0.3636364 1.863636 0.09090909 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
# 1533 C24H45N2O8PRS 0.3333333 1.875000 0.08333333 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
# 1537 C17H31N2O8PRS 0.4705882 1.823529 0.11764706 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
# 1539 C19H35N2O8PRS 0.4210526 1.842105 0.10526316 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
# 1541 C21H39N2O8PRS 0.3809524 1.857143 0.09523810 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
# 1543 C23H43N2O8PRS 0.3478261 1.869565 0.08695652 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2

sel.diab.prot2 <- which(dat.test.sig$cpd %in% c("cpd11572","cpd11531","cpd11473","cpd11548"))
dat.test.sig[sel.diab.prot2, ]
# cpd data_for_this_cpd median_t2d median_normal  alt      p_val W_statistic hl_effect_wilcox trend_with_disease incResto_incT2D decResto_decT2D incResto_decT2D decResto_incT2D
# 1515 cpd11473              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA              NA               1              NA
# 1530 cpd11531              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA              NA               1              NA
# 1536 cpd11548              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA              NA               1              NA
# 1546 cpd11572              <NA>  -4.264883     -4.141466 less 0.02322828         519       -0.1169812         Decreasing              NA              NA               1              NA
# incResto_notT2D decResto_notT2D notResto_incT2D notResto_decT2D notResto_notT2D anyResto_incT2D anyResto_decT2D anyResto_notT2D minuslog10_p_val                       cpd_names
# 1515              NA              NA              NA              NA              NA              NA               1              NA         1.633983             (2E)-Hexenoyl-[acp]
# 1530              NA              NA              NA              NA              NA              NA               1              NA         1.633983  9-methyl-trans-dec-2-enoyl-ACP
# 1536              NA              NA              NA              NA              NA              NA               1              NA         1.633983 4-methyl-trans-pent-2-enoyl-ACP
# 1546              NA              NA              NA              NA              NA              NA               1              NA         1.633983       trans-Octodec-2-enoyl-ACP
# cpd_forms      OC_x     HC_y       NC_z                                                 trend_group       z_layer
# 1515 C17H29N2O8PRS 0.4705882 1.705882 0.11764706 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
# 1530 C22H39N2O8PRS 0.3636364 1.772727 0.09090909 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
# 1536 C17H29N2O8PRS 0.4705882 1.705882 0.11764706 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2
# 1546 C29H53N2O8PRS 0.2758621 1.827586 0.06896552 Decreasing in T2D (reduced exposure in degraded ecosystems) N:C >0 to 0.2

# estimate of line 2: y = -0.625x + 2.0
# round to 4 d.p.

subsel.diab.prot2 <- which(round(dat.test.sig$HC_y[sel.ok], 4) == round(-0.625*dat.test.sig$OC_x[sel.ok] + 2, 4) & dat.test.sig$z_layer[sel.ok] == "N:C >0 to 0.2" & dat.test.sig$trend_with_disease[sel.ok] == "Decreasing")
length(subsel.diab.prot2) # 23


dat.test.sig[ sel.ok[unique(c(subsel.diab.prot1, subsel.diab.prot2))], ]
# cpd data_for_this_cpd      p_val median_t2d median_Normal W_statistic diff_median_perc diff_mean_perc  alt trend_with_disease incResto_incT2D decResto_decT2D
# 1641 cpd11499              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1643 cpd11503              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1645 cpd11507              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1647 cpd11511              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1651 cpd11524              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1653 cpd11528              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1655 cpd11532              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1657 cpd11536              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1661 cpd11549              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1663 cpd11553              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1665 cpd11557              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1667 cpd11561              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1636 cpd11475              <NA> 0.01851005  -4.252574     -4.068426         510    -2.952102e-05  -4.215527e-05 less         Decreasing              NA              NA
# 1637 cpd11465              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1638 cpd11469              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1639 cpd11473              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1640 cpd11498              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1642 cpd11502              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1644 cpd11506              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1646 cpd11510              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1648 cpd11514              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1649 cpd11518              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1650 cpd11523              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1652 cpd11527              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1654 cpd11531              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1656 cpd11535              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1658 cpd11539              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1659 cpd11543              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1660 cpd11548              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1662 cpd11552              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1664 cpd11556              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1666 cpd11560              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1668 cpd11564              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1669 cpd11568              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
# 1670 cpd11572              <NA> 0.02049546  -4.253046     -4.151108         514    -1.477306e-05  -2.056975e-05 less         Decreasing              NA              NA
#      incResto_decT2D decResto_incT2D incResto_notT2D decResto_notT2D minuslog10_p_val                             cpd_names     cpd_forms      OC_x     HC_y       NC_z mass
# 1641               1              NA              NA              NA         1.688342                 4-methyl-hexanoyl-ACP C18H33N2O8PRS 0.4444444 1.833333 0.11111111  469
# 1643               1              NA              NA              NA         1.688342                 6-methyl-octanoyl-ACP C20H37N2O8PRS 0.4000000 1.850000 0.10000000  497
# 1645               1              NA              NA              NA         1.688342                 8-methyl-decanoyl-ACP C22H41N2O8PRS 0.3636364 1.863636 0.09090909  525
# 1647               1              NA              NA              NA         1.688342              10-methyl-dodecanoyl-ACP C24H45N2O8PRS 0.3333333 1.875000 0.08333333  553
# 1651               1              NA              NA              NA         1.688342                 5-methyl-hexanoyl-ACP C18H33N2O8PRS 0.4444444 1.833333 0.11111111  469
# 1653               1              NA              NA              NA         1.688342                 7-methyl-octanoyl-ACP C20H37N2O8PRS 0.4000000 1.850000 0.10000000  497
# 1655               1              NA              NA              NA         1.688342                 9-methyl-decanoyl-ACP C22H41N2O8PRS 0.3636364 1.863636 0.09090909  525
# 1657               1              NA              NA              NA         1.688342              11-methyl-dodecanoyl-ACP C24H45N2O8PRS 0.3333333 1.875000 0.08333333  553
# 1661               1              NA              NA              NA         1.688342                4-methyl-pentanoyl-ACP C17H31N2O8PRS 0.4705882 1.823529 0.11764706  455
# 1663               1              NA              NA              NA         1.688342                6-methyl-heptanoyl-ACP C19H35N2O8PRS 0.4210526 1.842105 0.10526316  483
# 1665               1              NA              NA              NA         1.688342                 8-methyl-nonanoyl-ACP C21H39N2O8PRS 0.3809524 1.857143 0.09523810  511
# 1667               1              NA              NA              NA         1.688342              10-methyl-undecanoyl-ACP C23H43N2O8PRS 0.3478261 1.869565 0.08695652  539
# 1636               1              NA              NA              NA         1.732592                   (2E)-Decenoyl-[acp] C21H37N2O8PRS 0.3809524 1.761905 0.09523810  509
# 1637               1              NA              NA              NA         1.688342    But-2-enoyl-[acyl-carrier protein] C15H25N2O8PRS 0.5333333 1.666667 0.13333333  425
# 1638               1              NA              NA              NA         1.688342                 (2E)-Dodecenoyl-[acp] C23H41N2O8PRS 0.3478261 1.782609 0.08695652  537
# 1639               1              NA              NA              NA         1.688342                   (2E)-Hexenoyl-[acp] C17H29N2O8PRS 0.4705882 1.705882 0.11764706  453
# 1640               1              NA              NA              NA         1.688342        4-methyl-trans-hex-2-enoyl-ACP C18H31N2O8PRS 0.4444444 1.722222 0.11111111  467
# 1642               1              NA              NA              NA         1.688342        6-methyl-trans-oct-2-enoyl-ACP C20H35N2O8PRS 0.4000000 1.750000 0.10000000  495
# 1644               1              NA              NA              NA         1.688342        8-methyl-trans-dec-2-enoyl-ACP C22H39N2O8PRS 0.3636364 1.772727 0.09090909  523
# 1646               1              NA              NA              NA         1.688342     10-methyl-trans-dodec-2-enoyl-ACP C24H43N2O8PRS 0.3333333 1.791667 0.08333333  551
# 1648               1              NA              NA              NA         1.688342 12-methyl-trans-tetra-dec-2-enoyl-ACP C26H47N2O8PRS 0.3076923 1.807692 0.07692308  579
# 1649               1              NA              NA              NA         1.688342  14-methyl-trans-hexa-dec-2-enoyl-ACP C28H51N2O8PRS 0.2857143 1.821429 0.07142857  607
# 1650               1              NA              NA              NA         1.688342        5-methyl-trans-hex-2-enoyl-ACP C18H31N2O8PRS 0.4444444 1.722222 0.11111111  467
# 1652               1              NA              NA              NA         1.688342        7-methyl-trans-oct-2-enoyl-ACP C20H35N2O8PRS 0.4000000 1.750000 0.10000000  495
# 1654               1              NA              NA              NA         1.688342        9-methyl-trans-dec-2-enoyl-ACP C22H39N2O8PRS 0.3636364 1.772727 0.09090909  523
# 1656               1              NA              NA              NA         1.688342     11-methyl-trans-dodec-2-enoyl-ACP C24H43N2O8PRS 0.3333333 1.791667 0.08333333  551
# 1658               1              NA              NA              NA         1.688342 13-methyl-trans-tetra-dec-2-enoyl-ACP C26H47N2O8PRS 0.3076923 1.807692 0.07692308  579
# 1659               1              NA              NA              NA         1.688342  15-methyl-trans-hexa-dec-2-enoyl-ACP C28H51N2O8PRS 0.2857143 1.821429 0.07142857  607
# 1660               1              NA              NA              NA         1.688342       4-methyl-trans-pent-2-enoyl-ACP C17H29N2O8PRS 0.4705882 1.705882 0.11764706  453
# 1662               1              NA              NA              NA         1.688342       6-methyl-trans-hept-2-enoyl-ACP C19H33N2O8PRS 0.4210526 1.736842 0.10526316  481
# 1664               1              NA              NA              NA         1.688342        8-methyl-trans-non-2-enoyl-ACP C21H37N2O8PRS 0.3809524 1.761905 0.09523810  509
# 1666               1              NA              NA              NA         1.688342     10-methyl-trans-undec-2-enoyl-ACP C23H41N2O8PRS 0.3478261 1.782609 0.08695652  537
# 1668               1              NA              NA              NA         1.688342    12-methyl-trans-tridec-2-enoyl-ACP C25H45N2O8PRS 0.3200000 1.800000 0.08000000  565
# 1669               1              NA              NA              NA         1.688342  14-methyl-trans-pentadec-2-enoyl-ACP C27H49N2O8PRS 0.2962963 1.814815 0.07407407  593
# 1670               1              NA              NA              NA         1.688342             trans-Octodec-2-enoyl-ACP C29H53N2O8PRS 0.2758621 1.827586 0.06896552  621
#                                                       trend_group       z_layer
# 1641 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1643 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1645 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1647 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1651 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1653 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1655 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1657 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1661 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1663 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1665 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1667 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1636 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1637 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1638 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1639 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1640 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1642 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1644 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1646 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1648 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1649 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1650 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1652 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1654 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1656 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1658 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1659 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1660 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1662 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1664 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1666 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1668 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1669 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2
# 1670 Decreasing in T2D (reduced exposure in disturbed ecosystems) N:C >0 to 0.2

temp <- dat.test.sig[ sel.ok[unique(c(subsel.diab.prot1, subsel.diab.prot2))], ]
temp$cpd_names <- gsub(pattern = ",", replacement = ";", x = temp$cpd_names)

write.csv(x = temp, file = "Zoom1-Forslund-SWE-T2D-Hostremoval-DECREASING-protein-v6.csv", quote = FALSE, row.names = FALSE )

t2d.zoom1.decreasing.proteins <- dat.test.sig$cpd[ sel.ok[unique(c(subsel.diab.prot1, subsel.diab.prot2)) ]]

length(t2d.zoom1.decreasing.proteins) # 35


## Zoom1 - a) plot in restoration ; b) plot in T2D


# a) in restoration

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds" )
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 15

# Identified Proteins
df <- dat.cpd.collate
sel <- which(df$cpd_id %in% t2d.zoom1.decreasing.proteins) # 525
df <- df[sel, ]

str(df)
# 'data.frame':	525 obs. of  7 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)

  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))

  print(paste0("completed ",i))
}

str(res)
unique(res$group) # "22" "31" "UM" "12" "6"
res$group <- factor(res$group, levels = c("6","12", "22", "31", "UM"), ordered = TRUE)
unique(res$group_label) # "22 yr"   "31 yr"   "Unmined" "12 yr"   "6 yr"
res$group_label <- factor(res$group_label, levels = c("6 yr","12 yr", "22 yr", "31 yr", "Unmined"), ordered = TRUE)
res$age_vec <- as.integer(res$group)

# Kendall Tau correlation
ktcor<- cor.test(x = res$age_vec, y = res$sum_rel_abun, method = "kendall")
# Warning message:
#   In cor.test.default(x = res$age_vec, y = res$sum_rel_abun, method = "kendall") :
#   Cannot compute exact p-value with ties
ktcor
# Kendall's rank correlation tau
# 
# data:  res$age_vec and res$sum_rel_abun
# z = 3.6459, p-value = 0.0002665
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.7406561 

pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("Kendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)


str(res)

p <- ggplot(data = res, aes(x = age_vec, y = sum_rel_abun) )+
  geom_point()+
  #xlab("Reveg age (years)")+ ylab("CPP of indicator group 1a - \U03A3 rel abun (%)")+
  xlab("Reveg age (years)")+ ylab("\U03A3 CPP rel abun BCFA-ACPs (%)")+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_x_continuous(labels= c("6","12", "22", "31", "UM"))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "left", npcy = "top", label = "Forest restoration", size = 3.5)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    )

p

grid.text(label = "(a)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
#dev.print(tiff, file = paste0(workdir,"/plots/","Zoom1-a-Indicator-group1a-BCFA-ACPs-Trend-with-Age-v4.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom1-a-Indicator-group1a-BCFA-ACPs-Trend-with-Age-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



# b) in T2D

dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-v6.rds")
str(dat.cpd.collate)
# 'data.frame':	533824 obs. of  5 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 76

sel <- which(df$cpd_id %in% t2d.zoom1.decreasing.proteins)
df <- df[sel, ]
length(unique(df$cpd_id)) # 35
35*76 # 2660

str(df)
# 'data.frame':	2660 obs. of  5 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
# 'data.frame':	76 obs. of  3 variables:

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

# Wilcoxon-Mann-Whitney Test
x <- res$sum_rel_abun[ res$group_label == "T2D met-" ] # n = 33
y <- res$sum_rel_abun[ res$group_label == "Normal" ] # n = 43

wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %

pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Wilcoxon-Mann-Whitney ~T2D\nW = ",round(wmw.test$statistic,3),"; ",pval)

str(res)

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  #geom_point(alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("\U03A3 CPP rel abun BCFA-ACPs (%)")+ # ylab("CPP of indicator group 1a - \U03A3 rel abun (%)")+
  #geom_smooth(method = "lm")+
  theme_bw()+
  #scale_x_continuous(labels= c("Normal","IGT", "T2D met+", "T2D met-"))+
  
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "right", npcy = "top", label = "T2D case study", size = 3.5)+
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(b)", x = unit(0.04, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
#dev.print(tiff, file = paste0(workdir,"/plots/","Zoom1-b-Indicator-group1a-BCFA-ACPs-Trend-with-SWE-T2D-v6.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom1-b-Indicator-group1a-BCFA-ACPs-Trend-with-SWE-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




# # # # # #
# # # # # #

# ZOOM2

p <- ggplot(data = filter( dat.test.sig[sel.ok, ], z_layer == "N:C = 0" ) ) + # [subsel.diab.carbline]
  #coord_equal()+
  #ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study\n(Focused on compounds significant BH in ecosystem restoration)")+
  
  xlim(0.825, 1.05)+ ylim(1.65,2.1)+ geom_text_repel(aes(x= OC_x, y = HC_y, label = cpd), size = 1.75)+
  #xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + # 
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_group), size = 2, alpha = 0.35 ) + # 
  scale_color_manual(values = col.trend_group)+
  
  annotate(geom = "label_npc", npcx = 0.02, npcy = 0.97, label = "N:C = 0", size = 3.5)+ # geom = "text_npc"
  
  xlab("O:C molar ratio")+ ylab("H:C molar ratio")+
  guides(color = guide_legend(title = "Trend with disease\nin functional capacity (%)\nallocated to compounds\n(and corresponding trend\nwithin ecosystems)", nrow = 4 ))+ # nrow = 2, ncol = 2
  
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.key.spacing.y = unit(0, units = "line"),
    legend.text = element_text(size = rel(0.8), lineheight = 0.8) ,
    #legend.text = element_text(size = rel(0.9), lineheight = 0.8) ,
    title = element_text(size = rel(0.9), lineheight = 1)
    #strip.background = element_rect(fill = "transparent")
  )

p
dev.print(tiff, file = paste0(workdir,"/plots/","ZOOM2-cpp3d-Compounds-indiv-vKSpace-Trend-with-Disease--RESTO-SIGBH-SELECT-SWE-T2D-Hostremoval-vs-NORM--Wilcox-v6.tiff"), width = 20, height = 16, units = "cm", res=600, compression="lzw",type="cairo")



sel.diab.carbs <- which(dat.test.sig$cpd %in% c("cpd00082","cpd00076","cpd03200","cpd01399"))

dat.test.sig[sel.diab.carbs, ]
# cpd data_for_this_cpd median_t2d median_normal     alt       p_val W_statistic hl_effect_wilcox trend_with_disease incResto_incT2D decResto_decT2D
# 495 cpd03200              <NA> -1.5762293    -1.6369490 greater 0.042850471         874       0.04832176         Increasing              NA              NA
# 640 cpd00076              <NA> -1.1834958    -1.2575402 greater 0.034831736         883       0.05897750         Increasing              NA              NA
# 696 cpd00082              <NA> -0.8133804    -0.8440522 greater 0.013953670         919       0.05276007         Increasing              NA              NA
# 820 cpd01399              <NA> -1.2443415    -1.3123971 greater 0.007659197         940       0.05036158         Increasing               1              NA
# incResto_decT2D decResto_incT2D incResto_notT2D decResto_notT2D notResto_incT2D notResto_decT2D notResto_notT2D anyResto_incT2D anyResto_decT2D
# 495              NA               1              NA              NA              NA              NA              NA               1              NA
# 640              NA               1              NA              NA              NA              NA              NA               1              NA
# 696              NA               1              NA              NA              NA              NA              NA               1              NA
# 820              NA              NA              NA              NA              NA              NA              NA               1              NA
# anyResto_notT2D minuslog10_p_val     cpd_names cpd_forms      OC_x     HC_y NC_z                                                   trend_group z_layer
# 495              NA         1.368044 Manninotriose C18H32O16 0.8888889 1.777778    0 Increasing in T2D (increased exposure in degraded ecosystems) N:C = 0
# 640              NA         1.458025       Sucrose C12H22O11 0.9166667 1.833333    0 Increasing in T2D (increased exposure in degraded ecosystems) N:C = 0
# 696              NA         1.855312    D-Fructose   C6H12O6 1.0000000 2.000000    0 Increasing in T2D (increased exposure in degraded ecosystems) N:C = 0
# 820              NA         2.115817 Maltotetraose C24H42O21 0.8750000 1.750000    0  Increasing in T2D (increased exposure in quality ecosystems) N:C = 0


# estimate of line: y = 2x
# (round to 4 d.p.?)
subsel.diab.carbline <- which( round(dat.test.sig$HC_y[sel.ok],4) == round(2*dat.test.sig$OC_x[sel.ok],4) & dat.test.sig$OC_x[sel.ok] > 0.85 & dat.test.sig$OC_x[sel.ok] <= 1 & dat.test.sig$HC_y[sel.ok] > 1.7 & dat.test.sig$HC_y[sel.ok] <= 2 & dat.test.sig$NC_z[sel.ok] == 0 & dat.test.sig$trend_with_disease[sel.ok] == "Increasing")
length(subsel.diab.carbline) # 8

dat.test.sig[ sel.ok[subsel.diab.carbline], ]
# cpd data_for_this_cpd median_t2d median_normal     alt       p_val W_statistic hl_effect_wilcox trend_with_disease incResto_incT2D decResto_decT2D
# 107 cpd00224              <NA> -0.8707599    -0.9006397 greater 0.012151624         924       0.04413430         Increasing              NA              NA
# 494 cpd03198              <NA> -1.1422869    -1.1960000 greater 0.006013673         948       0.05371314         Increasing              NA              NA
# 495 cpd03200              <NA> -1.5762293    -1.6369490 greater 0.042850471         874       0.04832176         Increasing              NA              NA
# 632 cpd00259              <NA> -1.5817833    -1.6374173 greater 0.020250152         905       0.04767979         Increasing              NA              NA
# 640 cpd00076              <NA> -1.1834958    -1.2575402 greater 0.034831736         883       0.05897750         Increasing              NA              NA
# 692 cpd01074              <NA> -2.4026334    -2.4481669 greater 0.003177410         968       0.05699492         Increasing               1              NA
# 696 cpd00082              <NA> -0.8133804    -0.8440522 greater 0.013953670         919       0.05276007         Increasing              NA              NA
# 820 cpd01399              <NA> -1.2443415    -1.3123971 greater 0.007659197         940       0.05036158         Increasing               1              NA
# incResto_decT2D decResto_incT2D incResto_notT2D decResto_notT2D notResto_incT2D notResto_decT2D notResto_notT2D anyResto_incT2D anyResto_decT2D
# 107              NA               1              NA              NA              NA              NA              NA               1              NA
# 494              NA               1              NA              NA              NA              NA              NA               1              NA
# 495              NA               1              NA              NA              NA              NA              NA               1              NA
# 632              NA               1              NA              NA              NA              NA              NA               1              NA
# 640              NA               1              NA              NA              NA              NA              NA               1              NA
# 692              NA              NA              NA              NA              NA              NA              NA               1              NA
# 696              NA               1              NA              NA              NA              NA              NA               1              NA
# 820              NA              NA              NA              NA              NA              NA              NA               1              NA
# anyResto_notT2D minuslog10_p_val     cpd_names cpd_forms      OC_x     HC_y NC_z                                                   trend_group z_layer
# 107              NA         1.915366   L-Arabinose   C5H10O5 1.0000000 2.000000    0 Increasing in T2D (increased exposure in degraded ecosystems) N:C = 0
# 494              NA         2.220860     Melibiose C12H22O11 0.9166667 1.833333    0 Increasing in T2D (increased exposure in degraded ecosystems) N:C = 0
# 495              NA         1.368044 Manninotriose C18H32O16 0.8888889 1.777778    0 Increasing in T2D (increased exposure in degraded ecosystems) N:C = 0
# 632              NA         1.693572    D-Lyxulose   C5H10O5 1.0000000 2.000000    0 Increasing in T2D (increased exposure in degraded ecosystems) N:C = 0
# 640              NA         1.458025       Sucrose C12H22O11 0.9166667 1.833333    0 Increasing in T2D (increased exposure in degraded ecosystems) N:C = 0
# 692              NA         2.497927      Nigerose C12H22O11 0.9166667 1.833333    0  Increasing in T2D (increased exposure in quality ecosystems) N:C = 0
# 696              NA         1.855312    D-Fructose   C6H12O6 1.0000000 2.000000    0 Increasing in T2D (increased exposure in degraded ecosystems) N:C = 0
# 820              NA         2.115817 Maltotetraose C24H42O21 0.8750000 1.750000    0  Increasing in T2D (increased exposure in quality ecosystems) N:C = 0

write.csv(x = dat.test.sig[ sel.ok[subsel.diab.carbline], ], file = "Zoom2-Forslund-SWE-T2D-Hostremoval-Increasing-Carbs-v6.csv", quote = FALSE, row.names = FALSE )

t2d.zoom2.increasing.carbs <- dat.test.sig$cpd[ sel.ok[subsel.diab.carbline] ]

sort( dat.test.sig[ sel.ok[subsel.diab.carbline], "cpd_names" ])
# "D-Fructose"    "D-Lyxulose"    "L-Arabinose"   "Maltotetraose" "Manninotriose" "Melibiose"     "Nigerose"      "Sucrose"    

# adjusted so ONLY THOSE in trend group Increasing in T2D (increased exposure in disturbed ecosystems)
t2d.zoom2B.increasing.carbs <- c("cpd00224", "cpd03198", "cpd03200", "cpd00259", "cpd00076", "cpd00082")

dat.test.sig$cpd_names[ which(dat.test.sig$cpd %in% t2d.zoom2B.increasing.carbs) ]
# "L-Arabinose"   "Melibiose"     "Manninotriose" "D-Lyxulose"    "Sucrose"       "D-Fructose"  


## Zoom2 - a) plot in restoration ; b) plot in T2D


# a) in restoration

dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 15

# Identified Proteins
df <- dat.cpd.collate

sel <- which(df$cpd_id %in% t2d.zoom2B.increasing.carbs) # 90
df <- df[sel, ]

str(df)
# 'data.frame':	90 obs. of  7 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
unique(res$group) # "22" "31" "UM" "12" "6"
res$group <- factor(res$group, levels = c("6","12", "22", "31", "UM"), ordered = TRUE)
unique(res$group_label) # "22 yr"   "31 yr"   "Unmined" "12 yr"   "6 yr"
res$group_label <- factor(res$group_label, levels = c("6 yr","12 yr", "22 yr", "31 yr", "Unmined"), ordered = TRUE)
res$age_vec <- as.integer(res$group)

# Kendall Tau correlation
ktcor<- cor.test(x = res$age_vec, y = res$sum_rel_abun, method = "kendall")
# Warning message:
#   In cor.test.default(x = res$age_vec, y = res$sum_rel_abun, method = "kendall") :
#   Cannot compute exact p-value with ties
ktcor
# Kendall's rank correlation tau
# 
# data:  res$age_vec and res$sum_rel_abun
# z = -4.1522, p-value = 3.292e-05
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.843525 

pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("Kendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

length(t2d.zoom2B.increasing.carbs) # 6

str(res)

p <- ggplot(data = res, aes(x = age_vec, y = sum_rel_abun) )+
  geom_point()+
  xlab("Reveg age (years)")+ ylab("\U03A3 CPP rel abun Sugars (%)")+ # ylab("CPP of indicator group 2 - \U03A3 rel abun (%)")+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_x_continuous(labels= c("6","12", "22", "31", "UM"))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "left", npcy = "top", label = "Forest restoration", size = 3.5)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(a)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
#dev.print(tiff, file = paste0(workdir,"/plots/","Zoom2-a-Indicator-group2-SUGARS-Trend-with-Age-v4.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom2-a-Indicator-group2-SUGARS-Trend-with-Age-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


# b) in T2D

dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-v6.rds")
str(dat.cpd.collate)
# 'data.frame':	533824 obs. of  5 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 76

sel <- which(df$cpd_id %in% t2d.zoom2B.increasing.carbs)
df <- df[sel, ]
length(unique(df$cpd_id)) # 6

str(df)
# 'data.frame':	456 obs. of  5 variables:

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
# 'data.frame':	76 obs. of  3 variables:

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)


# Wilcoxon-Mann-Whitney Test

x <- res$sum_rel_abun[ res$group_label == "T2D met-" ] # n = 33
y <- res$sum_rel_abun[ res$group_label == "Normal" ] # n = 43

wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun % #  "less"

pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Wilcoxon-Mann-Whitney ~T2D\nW = ",round(wmw.test$statistic,3),"; ",pval)

str(res)

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  #geom_point(alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("\U03A3 CPP rel abun Sugars (%)")+ # ylab("CPP of indicator group 2 - \U03A3 rel abun (%)")+
  #geom_smooth(method = "lm")+
  theme_bw()+
  #scale_x_continuous(labels= c("Normal","IGT", "T2D met+", "T2D met-"))+
  
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "right", npcy = "top", label = "T2D case study", size = 3.5)+
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

p

grid.text(label = "(b)", x = unit(0.04, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
#dev.print(tiff, file = paste0(workdir,"/plots/","Zoom2-b-Indicator-group2-SUGARS-Trend-with-SWE-T2D-Hostremoval-v6.tiff"), width = 12, height = 11, units = "cm", res=600, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","Zoom2-b-Indicator-group2-SUGARS-Trend-with-SWE-T2D-Hostremoval-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

#-------------------------

#### Cross-check compounds with functions ??
#    in example of decreasing BCFA-ACPs
#-------------------------

df.out <-readRDS( "df.out--tidy-compounds_indiv-cpp3d--forslund-t2d-swe-hostremoval-pawsey.RDS" )

str(df.out)
# 'data.frame':	507391 obs. of  9 variables:
# $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_2" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" ...
# $ rxn_id             : chr  "rxn25279" "rxn25279" "rxn25279" "rxn03060" ...
# $ cpd_id             : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ cpd_name           : chr  "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" "2-methyl-trans-aconitate" "2-Methylcitrate" ...
# $ cpd_form           : chr  "C7H7O7" "H2O" "C7H5O6" "C7H7O7" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.3333 0.3333 0.3333 0.0556 0.0556 ...

head(df.out)
# superfocus_fxn f                                                                         f__in   rxn_id   cpd_id                                     cpd_name cpd_form
# 2          fxn_2 1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 3          fxn_2 1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd00001                                          H2O      H2O
# 4          fxn_2 1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd25681                     2-methyl-trans-aconitate   C7H5O6
# 5          fxn_3 1                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79) rxn03060 cpd01501                              2-Methylcitrate   C7H7O7
# 6          fxn_3 1                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79) rxn03060 cpd00001                                          H2O      H2O
# 7          fxn_3 1                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79) rxn03060 cpd02597                        cis-2-Methylaconitate   C7H5O6
# cpd_molar_prop cpd_molar_prop_norm
# 2              1          0.33333333
# 3              1          0.33333333
# 4              1          0.33333333
# 5              1          0.05555556
# 6              1          0.05555556
# 7              1          0.05555556


# look for these compounds across all T2D functions?

t2d.zoom1.decreasing.proteins
# [1] "cpd11499" "cpd11503" "cpd11507" "cpd11511" "cpd11524" "cpd11528" "cpd11532" "cpd11536" "cpd11549" "cpd11553" "cpd11557" "cpd11561" "cpd11475" "cpd11465" "cpd11469"
# [16] "cpd11473" "cpd11498" "cpd11502" "cpd11506" "cpd11510" "cpd11514" "cpd11518" "cpd11523" "cpd11527" "cpd11531" "cpd11535" "cpd11539" "cpd11543" "cpd11548" "cpd11552"
# [31] "cpd11556" "cpd11560" "cpd11564" "cpd11568" "cpd11572"


df.abun <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-v6.rds")
str(df.abun)
# 'data.frame':	533824 obs. of  5 variables:
# $ cpd_id      : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 4.982389 0 0 0.000188 ...
# $ log10_abun  : num  -8.187 0.697 -8.187 -8.187 -3.725 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

unique(df.abun$group_label)
# [1] T2D met- Normal  
# Levels: T2D met- < Normal

# now consider only compounds of interest
sel <- which(df.abun$cpd_id %in% t2d.zoom1.decreasing.proteins)
df.abun <- df.abun[sel, ]
str(df.abun)
# 'data.frame':	2660 obs. of  5 variables:
# $ cpd_id      : chr  "cpd11475" "cpd11465" "cpd11469" "cpd11473" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  3.63e-06 3.63e-06 3.63e-06 3.63e-06 3.63e-06 ...
# $ log10_abun  : num  -5.44 -5.44 -5.44 -5.44 -5.44 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

table(df.abun$group_label)
# T2D met-   Normal 
# 1155     1505

p <- ggplot(df.abun,aes(x = cpd_rel_abun, fill = group_label))+
  theme_bw()+
  guides(fill = guide_legend(title = "Diagnosis"))+
  geom_histogram(alpha = 0.6)+
  xlab("\U03A3 CPP (%)")+ ylab("Count")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p
dev.print(tiff, filename = paste0(workdir,"/plots/","Histogram-Decreasing-proteins-T2D-vs-Normal-v6.tiff"),
          width = 10, height = 8, units = "cm", res=600, compression = "lzw",type="cairo" )


# subset function - compound info to compounds of interest
head( df.out )
# superfocus_fxn f                                                                         f__in   rxn_id   cpd_id                                     cpd_name cpd_form
# 2          fxn_2 1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 3          fxn_2 1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd00001                                          H2O      H2O
# 4          fxn_2 1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd25681                     2-methyl-trans-aconitate   C7H5O6
# 5          fxn_3 1                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79) rxn03060 cpd01501                              2-Methylcitrate   C7H7O7
# 6          fxn_3 1                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79) rxn03060 cpd00001                                          H2O      H2O
# 7          fxn_3 1                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79) rxn03060 cpd02597                        cis-2-Methylaconitate   C7H5O6
# cpd_molar_prop cpd_molar_prop_norm
# 2              1          0.33333333
# 3              1          0.33333333
# 4              1          0.33333333
# 5              1          0.05555556
# 6              1          0.05555556
# 7              1          0.05555556

sel <- which(df.out$cpd_id %in% t2d.zoom1.decreasing.proteins) # 72
df.out <- df.out[sel, ]

unique(df.out$superfocus_fxn)
# "fxn_9321" "fxn_9354" "fxn_9355" "fxn_9429"



unique(df.out$f) # 2 1

unique(df.out$f__in)
# [1] "Trans-2-decenoyl-[acyl-carrier-protein] isomerase (EC 5.3.3.14)" "Enoyl-.acyl-carrier-protein. reductase .NADH. .EC 1.3.1.9."     
# [3] "Enoyl-.acyl-carrier-protein. reductase .NADPH. .EC 1.3.1.10." 

unique(df.out$rxn_id)
# [1] "rxn07455" "rxn05322" "rxn05324" "rxn05326" "rxn05327" "rxn05433" "rxn05434" "rxn05435" "rxn05436" "rxn05437" "rxn05438" "rxn05439" "rxn05440" "rxn05441" "rxn05442"
# [16] "rxn05443" "rxn05444" "rxn05445" "rxn05446" "rxn05447" "rxn05448" "rxn05449" "rxn05450" "rxn05464" "rxn05353" "rxn05355" "rxn05356" "rxn05357" "rxn05362" "rxn05366"
# [31] "rxn05370" "rxn05374" "rxn05378" "rxn05382" "rxn05387" "rxn05391" "rxn05395" "rxn05399" "rxn05403" "rxn05407" "rxn05412" "rxn05416" "rxn05420" "rxn05424" "rxn05428"
# [46] "rxn05432" "rxn05463"


phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval.RDS")

df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 17912     4

sel <- which(row.names(df.tax) %in% c("fxn_9321", "fxn_9354", "fxn_9355", "fxn_9429"))
df.tax[sel, ]
# subsys_L1   subsys_L2                         subsys_L3
# fxn_9321 Fatty Acids, Lipids, and Isoprenoids Fatty acids     Fatty Acid Biosynthesis FASII
# fxn_9354 Fatty Acids, Lipids, and Isoprenoids Fatty acids     Fatty Acid Biosynthesis FASII
# fxn_9355 Fatty Acids, Lipids, and Isoprenoids Fatty acids     Fatty Acid Biosynthesis FASII
# fxn_9429 Fatty Acids, Lipids, and Isoprenoids Fatty acids Unsaturated Fatty Acid Metabolism
# fxn
# fxn_9321 3-hydroxyacyl-[acyl-carrier-protein]_dehydratase,_FabA_form_(EC_4.2.1.59)_@_Trans-2-decenoyl-[acyl-carrier-protein]_isomerase_(EC_5.3.3.14)
# fxn_9354                                                                                  Enoyl-[acyl-carrier-protein]_reductase_[NADH]_(EC_1.3.1.9)
# fxn_9355                                                                                Enoyl-[acyl-carrier-protein]_reductase_[NADPH]_(EC_1.3.1.10)
# fxn_9429 3-hydroxyacyl-[acyl-carrier-protein]_dehydratase,_FabA_form_(EC_4.2.1.59)_@_Trans-2-decenoyl-[acyl-carrier-protein]_isomerase_(EC_5.3.3.14)


# linked to these enzymes: 
# EC 5.3.3.14 
# EC 1.3.1.9
# EC 1.3.1.10


# EC 5.3.3.14 - https://www.genome.jp/dbget-bin/www_bget?ec:5.3.3.14
# While the enzyme from Escherichia coli is highly specific for the 10-carbon enoyl-ACP, the enzyme from 
# Streptococcus pneumoniae can also use the 12-carbon enoyl-ACP as substrate in vitro but not 14- or 16-carbon enoyl-ACPs [3]. 
# ACP can be replaced by either CoA or N-acetylcysteamine thioesters. The cis-3-enoyl product is required to form 
# unsaturated fatty acids, such as palmitoleic acid and cis-vaccenic acid, in dissociated (or type II) fatty-acid biosynthesis.

# Orthology
# K18474 - trans-2-decenoyl-[acyl-carrier protein] isomerase
# - https://www.genome.jp/entry/K18474
# Reference	
# PMID:12237320
# Authors	
# Marrakchi H, Choi KH, Rock CO
# Title	
# A new mechanism for anaerobic unsaturated fatty acid formation in Streptococcus pneumoniae.
# Journal	
# J Biol Chem 277:44809-16 (2002)
# DOI:10.1074/jbc.M208920200

# K01716 - 3-hydroxyacyl-[acyl-carrier protein] dehydratase / trans-2-decenoyl-[acyl-carrier protein] isomerase
# - https://www.genome.jp/entry/K01716
# Reference	
# PMID:2832401
# Authors	
# Cronan JE Jr, Li WB, Coleman R, Narasimhan M, de Mendoza D, Schwab JM
# Title	
# Derived amino acid sequence and identification of active site residues of Escherichia coli beta-hydroxydecanoyl thioester dehydrase.
# Journal	
# J Biol Chem 263:4641-6 (1988)
# DOI:10.1016/S0021-9258(18)68830-1


# EC 1.3.1.9 - https://www.genome.jp/dbget-bin/www_bget?enzyme+1.3.1.9
# The enzyme catalyses an essential step in fatty acid biosynthesis, the reduction of the 2,3-double bond 
# in enoyl-acyl-[acyl-carrier-protein] derivatives of the elongating fatty acid moiety. The enzyme from the 
# bacterium Escherichia coli accepts substrates with carbon chain length from 4 to 18 [3]. 
# The FAS-I enzyme from the bacterium Mycobacterium tuberculosis prefers substrates with carbon chain length from 12 to 24 carbons.


# EC 1.3.1.10 - https://www.genome.jp/dbget-bin/www_bget?ec:1.3.1.10
# One of the activities of EC 2.3.1.86, fatty-acyl-CoA synthase system, an enzyme found in yeasts (Ascomycota and Basidiomycota). 
# Catalyses the reduction of enoyl-acyl-[acyl-carrier protein] derivatives of carbon chain length from 4 to 16. 
# The yeast enzyme is Si-specific with respect to NADP+. cf. EC 1.3.1.39, enoyl-[acyl-carrier-protein] reductase (NADPH, Re-specific) and 
# EC 1.3.1.104, enoyl-[acyl-carrier-protein] reductase (NADPH), which describes enzymes whose stereo-specificity towards NADPH is not known. See also EC 1.3.1.9, enoyl-[acyl-carrier-protein] reductase (NADH).



### Collate compounds of interest and append related functions

df.out <-readRDS( "df.out--tidy-compounds_indiv-cpp3d--forslund-t2d-swe-hostremoval-pawsey.RDS" )

str(df.out)
# 'data.frame':	507391 obs. of  9 variables:

phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval.RDS")
df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 17912     4


zoom_results <- list()

## Zoom 1a
temp <- read.csv(file = "Zoom1-Forslund-SWE-T2D-DECREASING-protein.csv", header = TRUE)
temp$zoom_group <- "1a"
temp$zoom_group_description <- "Decreasing BCFA-ACPs"

dim(temp) # 35 27

names(temp)
# [1] "cpd"                    "data_for_this_cpd"      "p_val"                  "median_t2d"             "median_Normal"          "W_statistic"            "diff_median_perc"      
# [8] "diff_mean_perc"         "alt"                    "trend_with_disease"     "incResto_incT2D"        "decResto_decT2D"        "incResto_decT2D"        "decResto_incT2D"       
# [15] "incResto_notT2D"        "decResto_notT2D"        "minuslog10_p_val"       "cpd_names"              "cpd_forms"              "OC_x"                   "HC_y"                  
# [22] "NC_z"                   "mass"                   "trend_group"            "z_layer"                "zoom_group"             "zoom_group_description"

# add functions?
temp$fxns <- NA
temp$superfocus_fxns <- NA
temp$subsys_L3 <- NA

for (i in 1:dim(temp)[1]) {
  #i<-1
  sel <- which(df.out$cpd_id == temp$cpd[i]) # 2
  ( list_fxns <- unique( df.out$superfocus_fxn[sel] ) ) # "fxn_9950"  "fxn_9951"
  sel <- which(row.names(df.tax) %in% list_fxns) # 2
  temp$fxns[i] <- paste0(list_fxns, collapse = "  |  ")
  temp$superfocus_fxns[i] <- paste0( unique(df.tax$fxn[sel]), collapse = "  |  ")
  temp$subsys_L3[i] <- paste0( unique(df.tax$subsys_L3[sel]), collapse = "  |  ")
  print(paste0("completed ",i))
}
temp$superfocus_fxns <- gsub(pattern = ",", replacement = ";", x = temp$superfocus_fxns)
zoom_results[["1a. Decreasing BCFA-ACPs"]] <- temp
write.csv(x = temp, file = "Zoom1-Forslund-SWE-T2D-Hostremoval-DECREASING-BCFA-ACPs--with-functions-v6.csv", quote = FALSE, row.names = FALSE )




## Zoom2-Forslund-SWE-T2D-Increasing-Carbs.csv
temp <- read.csv(file = "Zoom2-Forslund-SWE-T2D-Hostremoval-Increasing-Carbs-v6.csv", header = TRUE)
temp$zoom_group <- "2"
temp$zoom_group_description <- "Increasing carbs"

dim(temp) # 8 31

# add functions?
temp$fxns <- NA
temp$superfocus_fxns <- NA
temp$subsys_L3 <- NA

for (i in 1:dim(temp)[1]) {
  #i<-1
  sel <- which(df.out$cpd_id == temp$cpd[i]) # 2
  ( list_fxns <- unique( df.out$superfocus_fxn[sel] ) ) # "fxn_9950"  "fxn_9951"
  sel <- which(row.names(df.tax) %in% list_fxns) # 2
  temp$fxns[i] <- paste0(list_fxns, collapse = "  |  ")
  temp$superfocus_fxns[i] <- paste0( unique(df.tax$fxn[sel]), collapse = "  |  ")
  temp$subsys_L3[i] <- paste0( unique(df.tax$subsys_L3[sel]), collapse = "  |  ")
  print(paste0("completed ",i))
}

temp$superfocus_fxns <- gsub(pattern = ",", replacement = ";", x = temp$superfocus_fxns)
zoom_results[["2. Increasing carbs"]] <- temp
write.csv(x = temp, file = "Zoom2-Forslund-SWE-T2D-Hostremoval-Increasing-Carbs--with-functions-v6.csv", quote = FALSE, row.names = FALSE )


#-------------------------


##########################
##########################
##########################
##########################


#### Test for selected compounds in restoration (Sun & Badgley) and t2d (Forslund-SWE-T2D - w/ Host-removal)?
#-------------------------

lookfor <- as.data.frame( read_excel(path = "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/compounds-to-look-for.xlsx", sheet = 1, range = "A1:D44", col_names = TRUE) )
class(lookfor) # data.frame

lookfor$cpd_name
# [1] "cholecystokinin"                               "secretin"                                      "somatostatin"                                 
# [4] "motilin"                                       "ghrelin"                                       "glucagon-like peptide 1"                      
# [7] "glucose-dependent insulinotropic peptide"      "insulin-like peptide 5"                        "peptide YY"                                   
# [10] "gastrin"                                       "serotonin"                                     "neurotensin"                                  
# [13] "growth differentiation factor 15"              "fibroblast growth factor 19"                   "guanylin"                                     
# [16] "uroguanylin"                                   "oxyntomodulin"                                 "acetate"                                      
# [19] "propionate"                                    "butyrate"                                      "ammonia"                                      
# [22] "indole"                                        "trimethylamine N-oxide"                        "trimethylamine"                               
# [25] "p-cresyl sulfate"                              "indoxyl sulfate"                               "hydrogen sulfide"                             
# [28] "methane"                                       "menaquinone"                                   "cellulose"                                    
# [31] "xylan"                                         "pectin"                                        "amylopectin"                                  
# [34] "amylose"                                       "Carbohydrate response element binding protein" "glucosamine"                                  
# [37] "muramic acid"                                  "galactose"                                     "mannose"                                      
# [40] "arabinose"                                     "xylose"                                        "Carbon dioxide"                               
# [43] "adenosine triphosphate" 

## compound-level rel abundances in restoration and T2D

dat.cpd.res <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")
str(dat.cpd.res)
# 'data.frame':	125550 obs. of  7 variables:
# $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "20C" "20C" "20C" "20C" ...
# $ cpd_rel_abun: num  0.000435 0.022738 0.000456 5.135411 0.015705 ...
# $ log10_abun  : num  -3.362 -1.643 -3.341 0.711 -1.804 ...
# $ group       : Ord.factor w/ 5 levels "6"<"12"<"22"<..: 3 3 3 3 3 3 3 3 3 3 ...
# $ group_label : Ord.factor w/ 5 levels "6 yr"<"12 yr"<..: 3 3 3 3 3 3 3 3 3 3 ...
# $ ord_group   : num  3 3 3 3 3 3 3 3 3 3 ...

dat.cpd.t2d <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-v6.rds")
str(dat.cpd.t2d)
# 'data.frame':	533824 obs. of  5 variables:
# $ cpd_id      : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 4.982389 0 0 0.000188 ...
# $ log10_abun  : num  -8.187 0.697 -8.187 -8.187 -3.725 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

age_vec
# 6 12 22 31 UM 
# 1  2  3  4  5 

# T2D Met- and Normal only
unique(dat.cpd.t2d$group_label)
# [1] T2D met- Normal  
# Levels: T2D met- < Normal


length( unique(dat.cpd.t2d$cpd_id) ) # 7024
length( unique(dat.cpd.t2d$sample) ) # 76


# compound info is in here
names(df.comp)
# [1] "id"        "abbrev"    "name"      "form"      "OC_ratio"  "HC_ratio"  "NC_ratio"  "PC_ratio"  "NP_ratio"  "O_count"   "N_count"   "P_count"   "S_count"   "mass"     
# [15] "SC_ratio"  "MgC_ratio" "ZnC_ratio" "KC_ratio"  "CaC_ratio" "MnC_ratio" "FeC_ratio" "CoC_ratio" "CuC_ratio" "MoC_ratio"
head(df.comp$abbrev) # "h2o"   "atp"   "nad"   "nadh"  "nadph" "nadp"
head(df.comp$name) # "H2O"   "ATP"   "NAD"   "NADH"  "NADPH" "NADP" 
head(df.comp$form) # "H2O"           "C10H13N5O13P3" "C21H26N7O14P2" "C21H27N7O14P2" "C21H26N7O17P3" "C21H25N7O17P3"


## "Carbon dioxide"                               
sel.cpd <- which(df.comp$name == "CO2")
this_var <- "CO2"

df.comp[sel.cpd, ]
#          id abbrev name
# 11 cpd00011    co2  CO2

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00011")
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 4.1522, p-value = 3.292e-05
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#      tau 
# 0.843525
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00011")
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "cholecystokinin" 
sel.cpd <- which(df.comp$name == "cholecystokinin")
sel.cpd <- grep(pattern = "*holecystokinin", x = df.comp$name)
sel.cpd <- grep(pattern = "CCK|cck", x = df.comp$abbrev)
df.comp[sel.cpd, ]
#             id abbrev   name
# 13205 cpd14835  CCK-8  CCK-8
# 13224 cpd14854 CCK-33 CCK-33
this_var <- "CCK-8"
this_var <- "CCK-33"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14835") # empty
sel <- which(dat.cpd.res$cpd_id == "cpd14854") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14835") # empty
sel <- which(dat.cpd.t2d$cpd_id == "cpd14854")



## "secretin" (KEGG Compound C13523: Secretin; Vitrum)

sel.cpd <- which(df.comp$name == "secretin") # empty
sel.cpd <- grep(pattern = "*ecretin|*itrum", x = df.comp$name)
#sel.cpd <- grep(pattern = "CCK|cck", x = df.comp$abbrev)
df.comp[sel.cpd, ]
#             id abbrev   name             form
# 12827 cpd14242 Vitrum Vitrum C128H219N44O38R2
this_var <- "Secretin (Vitrum)"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14242") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14242") # empty



## "somatostatin"                                 

sel.cpd <- which(df.comp$name == "Somatostatin") # empty
sel.cpd <- grep(pattern = "*omatostatin", x = df.comp$name)
df.comp[sel.cpd, ]
#             id         abbrev           name
# 13113 cpd14743 Somatostatin-2 Somatostatin-2
# 13114 cpd14744 Somatostatin-1 Somatostatin-1
this_var <- "Somatostatin"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14743") # empty
sel <- which(dat.cpd.res$cpd_id == "cpd14744") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14743") # empty
sel <- which(dat.cpd.t2d$cpd_id == "cpd14744") # empty




## "motilin" 
sel.cpd <- which(df.comp$name == "Motilin") # ok
#sel.cpd <- grep(pattern = "*otilin", x = df.comp$name)
df.comp[sel.cpd, ]
#             id  abbrev    name
# 13138 cpd14768 Motilin Motilin
this_var <- "Motilin"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14768") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14768") # empty


## "ghrelin"

sel.cpd <- which(df.comp$name == "Ghrelin") # ok
#sel.cpd <- grep(pattern = "*hrelin", x = df.comp$name)
df.comp[sel.cpd, ]
#             id  abbrev    name
# 13117 cpd14747 Ghrelin Ghrelin
this_var <- "Ghrelin"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14747") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14747") # empty



## "glucagon-like peptide 1" ( or GLP-1;GLP1)                
sel.cpd <- which(df.comp$name == "GLP-1") # ok
sel.cpd <- grep(pattern = "GLP|*lucagon-like peptide 1", x = df.comp$name)
df.comp[sel.cpd, ]
#             id                  abbrev                    name form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio
# 13139 cpd14769 Glucagon-like peptide 1 Glucagon-like peptide 1
this_var <- "GLP-1"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14769") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14769") # empty



## "glucose-dependent insulinotropic peptide"
sel.cpd <- which(df.comp$name == "Glucose-dependent insulinotropic peptide") # empty
sel.cpd <- grep(pattern = "GIP|*lucose-dependent insulinotropic peptide", x = df.comp$name) # ok - GIP
df.comp[sel.cpd, ]
#             id abbrev name
# 13012 cpd14636    GIP  GIP
this_var <- "GIP"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14636") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14636") # empty



## "insulin-like peptide 5"
sel.cpd <- which(df.comp$name == "Glucose-dependent insulinotropic peptide") # empty
sel.cpd <- grep(pattern = "INSl5|*nsulin-like peptide 5", x = df.comp$name) # ok
df.comp[sel.cpd, ]
#             id                 abbrev                   name
# 13268 cpd14898 Insulin-like peptide 5 Insulin-like peptide 5
this_var <- "INSL5"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14898") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14898") # empty


## "peptide YY"                                   
sel.cpd <- which(df.comp$name == "Peptide YY") # ok
sel.cpd <- grep(pattern = "PYY|*eptide YY", x = df.comp$name) # ok
df.comp[sel.cpd, ]
#             id     abbrev       name
# 13209 cpd14839 Peptide YY Peptide YY
this_var <- "PYY"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14839") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14839") # empty


## "gastrin"  
sel.cpd <- which(df.comp$name == "Gastrin") # ok
#sel.cpd <- grep(pattern = "|*", x = df.comp$name) # ok
df.comp[sel.cpd, ]
#             id  abbrev    name
# 13204 cpd14834 Gastrin Gastrin
this_var <- "Gastrin"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14834") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14834") # empty


## "serotonin" 
sel.cpd <- which(df.comp$name == "Serotonin") # ok
#sel.cpd <- grep(pattern = "*erotonin", x = df.comp$name) # ok
df.comp[sel.cpd, ]
#           id    abbrev      name      form
# 568 cpd00579 Serotonin Serotonin C10H13N2O
this_var <- "Serotonin"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00579") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# 
# data:  x and y
# z = -1.8229, p-value = 0.06831
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# -0.370328
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00579") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "neurotensin"                                  
sel.cpd <- which(df.comp$name == "Neurotensin") # ok
#sel.cpd <- grep(pattern = "|*", x = df.comp$name) # ok
df.comp[sel.cpd, ]
#             id      abbrev        name
# 11780 cpd11960 Neurotensin Neurotensin
this_var <- ""
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd11960") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd11960") # empty


## "growth differentiation factor 15" 
sel.cpd <- which(df.comp$name == "Growth differentiation factor 15") # empty
sel.cpd <- grep(pattern = "GDF15|GDF-15", x = df.comp$name) # empty
df.comp[sel.cpd, ]
this_var <- ""
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "") # empty


## "fibroblast growth factor 19"
sel.cpd <- which(df.comp$name == "Fibroblast growth factor 19") # empty
sel.cpd <- grep(pattern = "FGF|*ibroblast growth", x = df.comp$name) # empty
df.comp[sel.cpd, ]
this_var <- ""
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "") # empty


## "guanylin"                                     
sel.cpd <- which(df.comp$name == "Guanylin") # ok
#sel.cpd <- grep(pattern = "|*", x = df.comp$name) # empty
df.comp[sel.cpd, ]
#             id   abbrev     name
# 13099 cpd14729 Guanylin Guanylin
this_var <- ""
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14729") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14729") # empty


## "uroguanylin"
sel.cpd <- which(df.comp$name == "Uroguanylin") # ok
#sel.cpd <- grep(pattern = "|*roguanylin", x = df.comp$name) # empty
df.comp[sel.cpd, ]
#             id      abbrev        name
# 13098 cpd14728 Uroguanylin Uroguanylin
this_var <- ""
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd14728") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd14728") # empty


## "oxyntomodulin"
sel.cpd <- which(df.comp$name == "Oxyntomodulin") # ok
#sel.cpd <- grep(pattern = "|*", x = df.comp$name) # empty
df.comp[sel.cpd, ]
#             id        abbrev          name
# 16751 cpd19466 Oxyntomodulin Oxyntomodulin
this_var <- ""
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd19466") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd19466") # empty


## "acetate"
sel.cpd <- which(df.comp$name == "Acetate") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$name) # ok
df.comp[sel.cpd, ]
#          id abbrev    name   form
# 29 cpd00029     ac Acetate C2H3O2
this_var <- "Acetate"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00029") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# 
# data:  x and y
# z = 3.5446, p-value = 0.0003932
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.7200823
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00029") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "propionate"
sel.cpd <- which(df.comp$name == "Propionate") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$name) # ok
df.comp[sel.cpd, ]
#           id abbrev       name   form
# 140 cpd00141    ppa Propionate C3H5O2
this_var <- "Propionate"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00141") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# 
# data:  x and y
# z = 2.1268, p-value = 0.03344
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.4320494 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00141") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "butyrate"
sel.cpd <- which(df.comp$name == "Butyrate") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$name) # ok
df.comp[sel.cpd, ]
#           id abbrev     name   form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio MnC_ratio
# 210 cpd00211    but Butyrate C4H7O2
this_var <- "Butyrate"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00211") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# 
# data:  x and y
# z = 3.2408, p-value = 0.001192
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#      tau 
# 0.658361
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00211") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "ammonia"
#sel.cpd <- which(df.comp$name == "Ammonia") # empty
#sel.cpd <- grep(pattern = "NH3|*mmonia", x = df.comp$name) # ok
sel.cpd <- which(df.comp$name == "NH3") # ok
df.comp[sel.cpd, ]
#          id abbrev name
# 13 cpd00013    nh4  NH3
this_var <- "Ammonia"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00013") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# 
# data:  x and y
# z = -1.7217, p-value = 0.08513
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#        tau 
# -0.3497543
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)


# Kruskal-Wallis test
kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  cpd_rel_abun by ord_group
# Kruskal-Wallis chi-squared = 8.9, df = 4, p-value = 0.06365


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00013") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "indole" (KEGG C00463: 2,3-Benzopyrrole)
#sel.cpd <- which(df.comp$name == "indole") # empty
#sel.cpd <- which(df.comp$name == "benzopyrrole") # empty
sel.cpd <- which(df.comp$abbrev == "indole" & df.comp$form == "C8H7N") # ok
#sel.cpd <- grep(pattern = "Indole", x = df.comp$abbrev) # empty
#sel.cpd <- grep(pattern = "*ndole", x = df.comp$name) # empty
df.comp[sel.cpd, ]
#           id abbrev  name  form 
# 356 cpd00359 indole indol C8H7N

this_var <- "Indole"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00359") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 0.10127, p-value = 0.9193
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.02057378
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# Kruskal-Wallis test
kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  cpd_rel_abun by ord_group
# Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00359") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "trimethylamine N-oxide"
#sel.cpd <- which(df.comp$name == "TMAO") # empty
sel.cpd <- which(df.comp$form == "C3H9NO") # ok
#sel.cpd <- grep(pattern = "*rimethylamine", x = df.comp$name) # empty
df.comp[sel.cpd, ]
#           id abbrev     name   form
# 797 cpd00811   tmao (CH3)3NO C3H9NO

this_var <- "TMAO"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00811") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# 
# data:  x and y
# z = 2.8357, p-value = 0.004573
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.5760658
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00811") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
#x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = log10_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10 CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "trimethylamine"                               
sel.cpd <- which(df.comp$abbrev == "tma") # ok
#sel.cpd <- which(df.comp$form == "C3H9NO") # ok
#sel.cpd <- grep(pattern = "*rimethylamine", x = df.comp$name) # empty
df.comp[sel.cpd, ]
#           id abbrev    name   form 
# 437 cpd00441    tma (CH3)3N C3H10N

this_var <- "TMA"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00441") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# 
# data:  x and y
# z = 2.7344, p-value = 0.006249
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.5554921
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00441") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
#x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = log10_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10 CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "p-cresyl sulfate" - instead use precursor: p-Cresol
sel.cpd <- which(df.comp$abbrev == "p-Cresol") # ok
#sel.cpd <- which(df.comp$form == "") # ok
#sel.cpd <- grep(pattern = "*resol", x = df.comp$abbrev) # empty
df.comp[sel.cpd, ]
# id   abbrev     name  form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass
# 1022 cpd01042 p-Cresol p-Cresol C7H8O

this_var <- "p-Cresol"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd01042") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 0.20255, p-value = 0.8395
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.04114756
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd01042") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "hydrogen sulfide"                             
sel.cpd <- which(df.comp$name == "H2S") # ok
#sel.cpd <- which(df.comp$form == "") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$abbrev) # empty
df.comp[sel.cpd, ]
#           id abbrev name
# 237 cpd00239    h2s  H2S

this_var <- "H2S"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00239") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -3.0382, p-value = 0.00238
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.6172134
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00239") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "methane"
sel.cpd <- which(df.comp$name == "Methane") # ok
#sel.cpd <- which(df.comp$form == "CH4") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$abbrev) # empty
df.comp[sel.cpd, ]
#            id abbrev    name form
# 1004 cpd01024  metha Methane  CH4

this_var <- "CH4"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd01024") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 2.3293, p-value = 0.01984
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.4731969 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd01024") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
#x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = log10_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10 CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "menaquinone"
sel.cpd <- which(df.comp$abbrev == "Menaquinones") # ok
#sel.cpd <- which(df.comp$form == "CH4") # ok
#sel.cpd <- grep(pattern = "*enaquinones", x = df.comp$abbrev) # multiple
#sel.cpd <- grep(pattern = "C16H16O2*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id       abbrev         name      form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count
# 24771 cpd27501 Menaquinones Menaquinones C16H15O2R

this_var <- "Menaquinones"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd27501") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 1.9242, p-value = 0.05433
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.3909018
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# Kruskal-Wallis test
kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  cpd_rel_abun by ord_group
# Kruskal-Wallis chi-squared = 7.6333, df = 4, p-value = 0.106


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd27501") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "cellulose"                                    
sel.cpd <- which(df.comp$name == "Cellulose") # ok
#sel.cpd <- which(df.comp$form == "CH4") # ok
#sel.cpd <- grep(pattern = "*ellulose", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "C16H16O2*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id    abbrev      name      form
# 11571 cpd11746 Cellulose Cellulose C6H10O5R2

this_var <- "Cellulose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd11746") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -0.50637, p-value = 0.6126
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.1028689 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd11746") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "xylan"
sel.cpd <- which(df.comp$name == "Xylan") # ok
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "xylan", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
# id                                    abbrev                                      name        form  OC_ratio  HC_ratio
# 10924 cpd11070 3beta-Hydroxylanostane-7,11-dione acetate 3beta-Hydroxylanostane-7,11-dione acetate    C32H52O4 0.1250000 1.6250000
# 11789 cpd11970                              Arabinoxylan                              Arabinoxylan  C10H16O8R2 0.8000000 1.6000000
# 12068 cpd12254                     Glucuronoarabinoxylan                     Glucuronoarabinoxylan        null       NaN       NaN
# 12217 cpd12408                            Glucuronoxylan                            Glucuronoxylan C21H31O18R2 0.8571429 1.4761905
# 12354 cpd12550                  4-O-methylglucuronoxylan                  4-O-methylglucuronoxylan C22H33O18R2 0.8181818 1.5000000
# 17230 cpd19945                     11-Deoxylandomycinone                     11-Deoxylandomycinone    C19H14O5 0.2631579 0.7368421
# 19578 cpd22295                               Acetylxylan                               Acetylxylan C27H39O21R3 0.7777778 1.4444444
# 19670 cpd22387                              Arabinoxylan                              Arabinoxylan C35H56O29R2 0.8285714 1.6000000
# 24444 cpd27172    Glucuronoarabinoxylan-Oligosaccharides    Glucuronoarabinoxylan-Oligosaccharides C27H41O23R2 0.8518519 1.5185185
# 24445 cpd27173                    Glucuronoarabinoxylans                    Glucuronoarabinoxylans C52H81O43R2 0.8269231 1.5576923
# 24447 cpd27175           Glucuronoxylan-Oligosaccharides           Glucuronoxylan-Oligosaccharides C22H33O19R2 0.8636364 1.5000000
# 24448 cpd27176                           Glucuronoxylans                           Glucuronoxylans C37H57O31R2 0.8378378 1.5405405
# 25918 cpd28650                      (1, 4-beta-D-xylan)n                      (1, 4-beta-D-xylan)n 

sel.cpd <- which(df.comp$name == "Arabinoxylan") # ok
df.comp[sel.cpd, ]
#             id       abbrev         name        form
# 11789 cpd11970 Arabinoxylan Arabinoxylan  C10H16O8R2 - NOT PRESENT IN DATA
# 19670 cpd22387 Arabinoxylan Arabinoxylan C35H56O29R2


this_var <- "Arabinoxylan"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd22387") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -2.228, p-value = 0.02588
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.4526232
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd22387") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## "Hemicellulose"  
# Hemicelluloses include xyloglucans, xylans, mannans and glucomannans
sel.cpd <- which(df.comp$name == "Hemicellulose") # ok
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*ellulose", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id        abbrev          name form
# 27115 cpd29869 Hemicellulose Hemicellulose null
this_var <- "Hemicellulose"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd29869") # empty - Hemicellulose
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd29869") # empty





## "pectin"
sel.cpd <- which(df.comp$name == "Pectin" & df.comp$form == "C26H34O24R2") # ok
#sel.cpd <- which(df.comp$form == "") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id abbrev   name        form
# 11434 cpd11601 pectin Pectin C26H34O24R2
this_var <- "Pectin"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd11601") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -1.9242, p-value = 0.05433
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.3909018 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd11601") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
#x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = log10_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10 CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "amylopectin"                                  
sel.cpd <- which(df.comp$name == "Amylopectin") # ok
#sel.cpd <- which(df.comp$form == "") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#           id      abbrev        name      form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count
# 262 cpd00265 Amylopectin Amylopectin C30H52O26
this_var <- "Amylopectin"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00265") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 1.6204, p-value = 0.1052
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.3291805 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00265") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "amylose"
sel.cpd <- which(df.comp$name == "Amylose") # ok
#sel.cpd <- which(df.comp$form == "") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id abbrev    name      form
# 11560 cpd11735 14glun Amylose C6H10O5R2
this_var <- "Amylose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd11735") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 1.2153, p-value = 0.2243
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.2468854
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd11735") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "Glycogen"

## Glycogen (animal starch) - Glycogen is a multibranched polysaccharide of glucose that serves as a form of energy storage in animals,[2] fungi, and bacteria
sel.cpd <- which(df.comp$name == "Glycogen") # 
df.comp[sel.cpd, ]
#           id   abbrev     name      form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 154 cpd00155 glycogen Glycogen C24H42O21
this_var <- "Glycogen"


# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00155") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 2.9369, p-value = 0.003315
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.5966396 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00155") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "glucosamine" - instead use Chitin because glucosamine derives from Chitin
# D-Glucosamine is made naturally in the form of glucosamine-6-phosphate
#sel.cpd <- which(df.comp$name == "Glucosamine") # ok - but not in samples
#sel.cpd <- grep(pattern = "D-glucosamine-6-phosphate", x = df.comp$name) # D-Glucosamine is made naturally in the form of glucosamine-6-phosphate
sel.cpd <- grep(pattern = "*D-glucosamine", x = df.comp$name)
sel.cpd <- which(df.comp$name == "alpha-D-glucosamine 6-phosphate") 
# "alpha-D-glucosamine" "alpha-D-glucosamine 6-phosphate"
#sel.cpd <- grep(pattern = "*mino-2-deoxy-glucose", x = df.comp$name)
#sel.cpd <- which(df.comp$name == "2-amino-2-deoxy-glucose") # ok - but not in samples
#sel.cpd <- which(df.comp$form == "C6H13NO5") # ok - but not in samples

unique(df.comp$name[sel.cpd])

df.comp[sel.cpd, ]
#             id                          abbrev                            name      form 
# 21178 cpd23898 alpha-D-glucosamine 6-phosphate alpha-D-glucosamine 6-phosphate C6H13NO8P


this_var <- "D-glucosamine 6-phosphate"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd23898") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -2.5318, p-value = 0.01135
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.5143445
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd23898") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "chitin"
sel.cpd <- which(df.comp$name == "Chitin") #
#sel.cpd <- which(df.comp$form == "") # ok
#sel.cpd <- grep(pattern = "*lucosamine", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id abbrev   name       form 
# 11508 cpd11683 chitin Chitin C8H13NO5R2
this_var <- "Chitin"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd11683") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -3.8484, p-value = 0.0001189
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.7818036 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt
# # Kruskal-Wallis rank sum test
# # data:  cpd_rel_abun by ord_group
# # Kruskal-Wallis chi-squared = 7.9, df = 4, p-value = 0.09531


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd11683") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")







## "galactose"
sel.cpd <- which(df.comp$name == "Galactose") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
df.comp[sel.cpd, ]
#           id  abbrev      name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 108  cpd00108     gal Galactose C6H12O6        1        2        0        0      NaN       6       0       0       0  180        0         0         0        0         0
# 1090 cpd01112 cbs_337 Galactose C6H12O6
this_var <- "Galactose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00108") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -3.0382, p-value = 0.00238
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.6172134
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt



p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00108") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "mannose"                                      

sel.cpd <- which(df.comp$name == "D-Mannose") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*annose$", x = df.comp$name) # multiple - use this D-mannose
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#           id abbrev      name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 138 cpd00138    man D-Mannose C6H12O6
this_var <- "D-Mannose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00138") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -2.9369, p-value = 0.003315
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.5966396
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt



p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00138") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "arabinose"

sel.cpd <- which(df.comp$name == "L-Arabinose") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*rabinose$", x = df.comp$name) # multiple - use this D-Arabinose
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#           id abbrev        name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio
# 223 cpd00224 arab-L L-Arabinose C5H10O5
this_var <- "L-Arabinose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00224") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -4.1522, p-value = 3.292e-05
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.843525
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00224") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "xylose"

sel.cpd <- which(df.comp$name == "Xylose") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*Xylose$", x = df.comp$name) # multiple - use this D-Arabinose
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id             abbrev               name        form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count
# 153   cpd00154              xyl-D             Xylose     C5H10O5 1.0000000 2.000000        0        0      NaN       5       0
# 1047  cpd01068           L-Xylose           L-Xylose     C5H10O5 1.0000000 2.000000        0        0      NaN       5       0
# 1397  cpd01422            cll_392      beta-D-Xylose     C5H10O5 1.0000000 2.000000        0        0      NaN       5       0
# 1460  cpd01487     alpha-D-Xylose     alpha-D-Xylose     C5H10O5 1.0000000 2.000000        0        0      NaN       5       0
# 24104 cpd26831           D-Xylose           D-Xylose     C5H10O5 1.0000000 2.000000        0        0      NaN       5       0
# 25629 cpd28361 Xyloglucans-Xylose Xyloglucans-Xylose C39H64O33R2 0.8461538 1.641026        0        0      NaN      33       0
# 26016 cpd28748         (+)-Xylose         (+)-Xylose     C5H10O5 
this_var <- "Xylose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00154") # Xylose
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -3.7471, p-value = 0.0001789
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.7612299 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt


p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00154") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## glucose
sel.cpd <- which(df.comp$name == "D-Glucose") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*lucose$", x = df.comp$name) # multiple - use this D-Arabinose
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
# id    abbrev      name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 27    cpd00027     glc-D D-Glucose C6H12O6        1        2        0        0      NaN       6       0       0       0  180        0         0         0        0         0
# 24094 cpd26821 D-Glucose D-Glucose C6H12O6
this_var <- "D-Glucose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00027") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -3.9497, p-value = 7.825e-05
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.8023774
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt



p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00027") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## Fructose
sel.cpd <- which(df.comp$name == "D-Fructose") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*ructose$", x = df.comp$name) # multiple - use this D-Arabinose
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#          id abbrev       name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 82 cpd00082    fru D-Fructose C6H12O6
this_var <- "D-Fructose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00082") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -3.7471, p-value = 0.0001789
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.7612299 
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt



p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00082") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## Sucrose
sel.cpd <- which(df.comp$name == "Sucrose") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*ucrose$", x = df.comp$name) # multiple - use this D-Arabinose
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#          id abbrev    name      form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 76 cpd00076   sucr Sucrose C12H22O11
this_var <- "Sucrose"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd00076") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = -3.8484, p-value = 0.0001189
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.7818036
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt



p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00076") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## leptin
sel.cpd <- which(df.comp$name == "Leptin") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*eptin$", x = df.comp$name) # multiple - use this D-Arabinose
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id abbrev   name form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count  mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio MnC_ratio
# 16743 cpd19458 Leptin Leptin
this_var <- "Leptin"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd19458") # empty
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd19458") # empty




## lignin

sel.cpd <- which(df.comp$name == "Lignin") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*ucrose$", x = df.comp$name) # multiple - use this D-Arabinose
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id abbrev   name form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio MnC_ratio
# 12548 cpd12745 Lignin Lignin null
this_var <- "Lignin"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd12745") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 3.9497, p-value = 7.825e-05
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.8023774
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt



p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd12745") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")





## Amylum (plant starch)

sel.cpd <- which(df.comp$name == "Starch") # 
#sel.cpd <- which(df.comp$name == "") # ok but none in samples
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*$", x = df.comp$name) # 
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id abbrev   name        form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count    mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 11482 cpd11657 starch Starch C12H20O10R2 0.8333333 1.666667        0        0      NaN      10       0       0       0 9.9e+02        0         0         0        0         0
# 25462 cpd28193 Starch Starch        null
this_var <- "Starch"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd11657") # 15
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]

# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  x and y
# z = 0.91147, p-value = 0.3621
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.185164
pval <- ifelse(test = ktcor$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(ktcor$p.value,3)) )
test_result <- paste0("RESTORATION\nKendall's tau ~Age\nTau = ",round(ktcor$estimate,3),"; ",pval)

# # Kruskal-Wallis test
# kt <- kruskal.test( cpd_rel_abun ~ ord_group, dat.cpd.res[sel, ]) # Kruskal Wallis test
# kt



p <- ggplot(data = data.frame(age_vec = x, value = y) , aes(x=age_vec, y = value))+ # log10abun
  ggtitle( this_var)+
  geom_point()+
  geom_point(shape = 1)+
  #geom_smooth(method="lm")+
  geom_smooth(method="loess")+
  theme_bw()+
  #xlab("Reveg age (years)")+ ylab("log10 \U03A3 CPP (%)/ vK volume")+
  xlab("Reveg age (years)")+ ylab("CPP rel abun (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  scale_x_continuous(labels= names(age_vec))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Kendall-Tau-Age-vs-",gsub(pattern="/",replacement="-", x = this_var),"-sunbad-resto-v2.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd11657") # 76
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 43



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




### fructan

sel.cpd <- which(df.comp$name == "Fructan") # empty
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*ructan", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#            id               abbrev                 name form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count
# 24363 cpd27091             Fructans             Fructans null      NaN      NaN      NaN      NaN      NaN       0       0       0
# 26031 cpd28763 (2,1-beta-D-fructan) (2,1-beta-D-fructan) null      NaN      NaN      NaN      NaN      NaN       0       0       0
# 26589 cpd29322   2,6-beta-D-fructan   2,6-beta-D-fructan null
this_var <- ""
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd27091") # empty
sel <- which(dat.cpd.res$cpd_id == "cpd28763") # empty
sel <- which(dat.cpd.res$cpd_id == "cpd29322") # empty




### inulin

                       
sel.cpd <- which(df.comp$name == "Inulin") # ok
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "*ructosyl", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id abbrev   name       form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass
# 11435 cpd11602 Inulin Inulin  C24H42O21 0.8750000 1.750000        0        0      NaN      21       0       0       0  342
# 24584 cpd27312 Inulin Inulin C18H31O15R
this_var <- "Inulin"
# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd11602") # empty
sel <- which(dat.cpd.res$cpd_id == "cpd27312") # empty




### galactan

                            
sel.cpd <- which(df.comp$name == "Galactan") # ok
#sel.cpd <- which(df.comp$form == "") # ok
#sel.cpd <- grep(pattern = "*", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
#             id   abbrev     name        form
# 12578 cpd12777 Galactan Galactan C12H20O11R2
this_var <- "Galactan"

# RESTORATION
sel <- which(dat.cpd.res$cpd_id == "cpd12777") # empty


#-------------------------

#### sample summary stats - restoration (Sun & Badgley) and t2d (Forslund-SWE-T2D - w/ Host-removal)?
#-------------------------

## Restoration

dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")

str(dat.cpd.collate)
# 'data.frame':	125550 obs. of  7 variables:
# $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "20C" "20C" "20C" "20C" ...
# $ cpd_rel_abun: num  0.000435 0.022738 0.000456 5.135411 0.015705 ...
# $ log10_abun  : num  -3.362 -1.643 -3.341 0.711 -1.804 ...
# $ group       : Ord.factor w/ 5 levels "6"<"12"<"22"<..: 3 3 3 3 3 3 3 3 3 3 ...
# $ group_label : Ord.factor w/ 5 levels "6 yr"<"12 yr"<..: 3 3 3 3 3 3 3 3 3 3 ...
# $ ord_group   : num  3 3 3 3 3 3 3 3 3 3 ...

length( unique(dat.cpd.collate$cpd_id) ) # 8370
length( unique(dat.cpd.collate$cpd_id[ dat.cpd.collate$cpd_rel_abun > 0] ) ) # 8370
length( unique(dat.cpd.collate$sample) ) # 15

data_in <- dat.cpd.collate

head(data_in)
# cpd_id sample cpd_rel_abun log10_abun group group_label ord_group
# 1 cpd25681    20C 0.0004345842 -3.3619261    22       22 yr         3
# 2 cpd02597    20C 0.0227382574 -1.6432428    22       22 yr         3
# 3 cpd24620    20C 0.0004564208 -3.3406346    22       22 yr         3
# 4 cpd00001    20C 5.1354106008  0.7105752    22       22 yr         3
# 5 cpd01501    20C 0.0157048650 -1.8039658    22       22 yr         3
# 6 cpd00851    20C 0.0121620567 -1.9149930    22       22 yr         3

dim(data_in) # 125550      7
8370*15 # 125550

unique_samps <- unique(data_in$sample)

no_compounds <- numeric(length = length(unique_samps))

for (i in 1:length(unique_samps)) {
  #i<-1
  this_samp <- unique_samps[i]
  sel <- which(data_in$sample == this_samp)
  
  values <- data_in$cpd_rel_abun[sel]
  values <- values[values > 0]
  
  no_compounds[i] <- length( values )
  print(paste0("completed ",i))
}

mean(no_compounds) # 7776.867
sd(no_compounds) # 92.37878

phy <- readRDS("phy-phyloseq-object-sunbad-resto.RDS")
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 30125 taxa and 15 samples ]
# sample_data() Sample Data:       [ 15 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 30125 taxa by 4 taxonomic ranks ]

head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 15
NonZeroFxns

mean(NonZeroFxns) # 20327.93
sd(NonZeroFxns) # 767.4159
#rm(NonZeroFxns)

# check
unique_samps <- sample_names(phy)

unique_samps[1] # "20C"
a <- prune_samples(unique_samps[1], phy)
min(taxa_sums(a)) # 0
# prune taxa that have zero sequence reads
a <- prune_taxa(taxa = taxa_sums(a) > 0, x = a)
a

unique_samps[15] # "30C"
a <- prune_samples(unique_samps[15], phy)
min(taxa_sums(a)) # 0
# prune taxa that have zero sequence reads
a <- prune_taxa(taxa = taxa_sums(a) > 0, x = a)
a






## T2D

dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-v6.rds")

str(dat.cpd.collate)
# 'data.frame':	533824 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 4.982389 0 0 0.000188 ...
# $ log10_abun  : num  -8.187 0.697 -8.187 -8.187 -3.725 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

length( unique(dat.cpd.collate$cpd_id) ) # 7024
length( unique(dat.cpd.collate$cpd_id[ dat.cpd.collate$cpd_rel_abun > 0] ) ) # 7024
length( unique(dat.cpd.collate$sample) ) # 76

data_in <- dat.cpd.collate

head(data_in)
#     cpd_id    sample cpd_rel_abun log10_abun group_label
# 1 cpd24620 ERR260139 0.0000000000 -8.1872270    T2D met-
# 2 cpd00001 ERR260139 4.9823891271  0.6974376    T2D met-
# 3 cpd25681 ERR260139 0.0000000000 -8.1872270    T2D met-
# 4 cpd01501 ERR260139 0.0000000000 -8.1872270    T2D met-
# 5 cpd02597 ERR260139 0.0001882679 -3.7252238    T2D met-
# 6 cpd00851 ERR260139 0.0012152908 -2.9153198    T2D met-

dim(data_in) # 533824      5
7024*76 # 533824

unique_samps <- unique(data_in$sample)

no_compounds <- numeric(length = length(unique_samps))

for (i in 1:length(unique_samps)) {
  #i<-1
  this_samp <- unique_samps[i]
  sel <- which(data_in$sample == this_samp)
  
  values <- data_in$cpd_rel_abun[sel]
  values <- values[values > 0]
  
  no_compounds[i] <- length( values )
  print(paste0("completed ",i))
}

mean(no_compounds) # 5287.329
sd(no_compounds) # 584.9439

phy <- readRDS("phy.cpp-cleaned-Forslund-SWE-T2D-Hostremoval-v6.RDS")
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7024 taxa and 76 samples ]
# sample_data() Sample Data:       [ 76 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 7024 taxa by 1 taxonomic ranks ]
summary( sample_sums(phy) )
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 59.86   68.92   71.10   70.52   72.54   75.30 
sd( sample_sums(phy) )
# 2.881688

# Functions
phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval.RDS")
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 17912 taxa and 76 samples ]
# sample_data() Sample Data:       [ 76 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 17912 taxa by 4 taxonomic ranks ]

head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 76
NonZeroFxns

mean(NonZeroFxns) # 7833.829
sd(NonZeroFxns) # 1699.389
#rm(NonZeroFxns)

# check
unique_samps[1] # "ERR260139"
a <- prune_samples(unique_samps[1], phy)
min(taxa_sums(a)) # o
# prune taxa that have zero sequence reads
a <- prune_taxa(taxa = taxa_sums(a) > 0, x = a)
a

unique_samps[76] # "ERR275252"
a <- prune_samples(unique_samps[76], phy)
min(taxa_sums(a)) # o
# prune taxa that have zero sequence reads
a <- prune_taxa(taxa = taxa_sums(a) > 0, x = a)
a

# check ages of subjects?
df.samp <- readRDS("df.samp.with-t2dclass-age-Forslund-SWE-T2D.RDS")
head(df.samp)

sel <- which(row.names(df.samp) %in% unique_samps)
summary(df.samp$age[sel])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 68.96   69.98   70.34   70.44   71.10   71.84
table(df.samp$group_new[sel])
# T2D met neg T2D met pos         IGT      Normal 
# 33           0           0          43


#-------------------------

##########################
##########################
##########################
##########################


#### Forslund-CHN-T2D - Compare SWE-T2D to Validation T2D dataset in Chinese (CHN) cohort?
#    select > 50 year olds
#    prepare for data download from SRA
#-------------------------

# Forslund et al 2015 states data source: Raw nucleotide data can be found for all samples used
# in the study in the Sequence Read Archive (accession numbers: SRA045646 and SRA050230, CHN samples)

# Forslund et al subject - diagnosis info

subjects <- read_excel(path = "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/test_CHN_data/Samples-and-diagnosis-SItableS1.xlsx", sheet = 1, range = "A1:C785")
subjects <- as.data.frame(subjects)
str(subjects)

sel <- which(subjects$`Country subset` == "CHN") # 256
subjects <- subjects[sel, ]
str(subjects)
# 'data.frame':	256 obs. of  3 variables:
# $ Sample        : chr  "BGI-06A" "BGI-15A" "BGI-17A" "BGI-27A" ...
# $ Country subset: chr  "CHN" "CHN" "CHN" "CHN" ...
# $ Status        : chr  "ND CTRL" "ND CTRL" "ND CTRL" "ND CTRL" ...
table(subjects$Status)
# ND CTRL T2D metformin- T2D metformin+ 
# 185             56             15 

subjects$`Sample Name` <- paste0("bgi-",subjects$Sample)

length(unique(subjects$`Sample Name`)) # 256


sradat <- read_excel(path = "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/test_CHN_data/Metadata-SRA045646-SraRunTable.xlsx", sheet = 1, range = "A1:AD371")
sradat <- as.data.frame(sradat)
str(sradat)
# 'data.frame':	370 obs. of  30 variables:
# $ Run                     : chr  "SRR341581" "SRR341582" "SRR341583" "SRR341584" ...
# $ actual_read_length (run): num  148 148 148 148 148 148 148 148 148 148 ...
# $ Age                     : num  59 43 46 25 60 71 62 54 51 29 ...
# $ Assay Type              : chr  "WGS" "WGS" "WGS" "WGS" ...
# $ AvgSpotLen              : num  148 148 148 148 148 148 148 148 148 148 ...
# $ Bases                   : num  2.46e+09 1.51e+09 2.00e+09 2.17e+09 2.41e+09 ...
# $ BioProject              : chr  "PRJNA422434" "PRJNA422434" "PRJNA422434" "PRJNA422434" ...
# $ BioSample               : chr  "SAMN00715131" "SAMN00715132" "SAMN00715133" "SAMN00715134" ...
# $ Bytes                   : num  1.17e+09 7.49e+08 9.14e+08 9.86e+08 1.20e+09 ...
# $ center_name (exp)       : chr  "BGI" "BGI" "BGI" "BGI" ...
# $ Center Name             : chr  "BGI" "BGI" "BGI" "BGI" ...
# $ Consent                 : chr  "public" "public" "public" "public" ...
# $ DATASTORE filetype      : chr  "fastq,sra" "sra,fastq" "sra,fastq" "sra,fastq" ...
# $ DATASTORE provider      : chr  "s3" "s3" "s3" "s3" ...
# $ DATASTORE region        : chr  "s3.us-east-1" "s3.us-east-1" "s3.us-east-1" "s3.us-east-1" ...
# $ Experiment              : chr  "SRX095662" "SRX095663" "SRX095664" "SRX095665" ...
# $ gender                  : chr  "female" "female" "female" "female" ...
# $ Instrument              : chr  "Illumina Genome Analyzer II" "Illumina Genome Analyzer II" "Illumina Genome Analyzer II" "Illumina Genome Analyzer II" ...
# $ Library Name            : chr  "HGMlijMCFDFAAPEI" "HGMlijMCYDFAAPEI" "HGMlijMBTDFAAPEI" "HGMlijMBXDFAAPEI" ...
# $ LibraryLayout           : chr  "PAIRED" "PAIRED" "PAIRED" "PAIRED" ...
# $ LibrarySelection        : chr  "RANDOM" "RANDOM" "RANDOM" "RANDOM" ...
# $ LibrarySource           : chr  "METAGENOMIC" "METAGENOMIC" "METAGENOMIC" "METAGENOMIC" ...
# $ NATION                  : chr  "China" "China" "China" "China" ...
# $ Organism                : chr  "human gut metagenome" "human gut metagenome" "human gut metagenome" "human gut metagenome" ...
# $ Platform                : chr  "ILLUMINA" "ILLUMINA" "ILLUMINA" "ILLUMINA" ...
# $ ReleaseDate             : chr  "2012-09-05T00:00:00Z" "2012-09-05T00:00:00Z" "2012-09-05T00:00:00Z" "2012-09-05T00:00:00Z" ...
# $ run (run)               : chr  "FC615J5AAXX" "FC61B1KAAXX" "FC61B1DAAXX" "FC61AVBAAXX" ...
# $ Sample Name             : chr  "bgi-DLF001" "bgi-DLF002" "bgi-DLF003" "bgi-DLF004" ...
# $ SRA Study               : chr  "SRP008047" "SRP008047" "SRP008047" "SRP008047" ...
# $ total_bases (run)       : num  2.46e+09 1.51e+09 2.00e+09 2.17e+09 2.41e+09 ...

length(unique(sradat$`Sample Name`)) # 370

unique(sradat$BioProject) # "PRJNA422434"

sel <- which(sradat$`Sample Name` %in% subjects$`Sample Name`) # 256
length(unique(sradat$`Sample Name`[sel])) # 256

sradat.select <- sradat[sel, ]
# order?
identical( sradat.select$`Sample Name`, subjects$`Sample name` ) # FALSE

length(unique(sradat.select$`Sample Name`)) # 256

sradat.select$Diagnosis <- NA

for (i in 1:dim(sradat.select)[1]) {
  #i<-146
  this_samp <- sradat.select$`Sample Name`[i]
  sel.sub <- which(subjects$`Sample Name` == this_samp)
  
  sradat.select$Diagnosis[i] <- subjects$Status[sel.sub]
  
  print(paste0("completed ",i))

}

temp <- sradat.select

# # choose females only
# 
# sel <- which(sradat.select$gender == "female")
# 
# sradat.select.f <- sradat.select[sel, ]
# dim(sradat.select.f) # 115 31
# table(sradat.select.f$Diagnosis)
# # ND CTRL T2D metformin- T2D metformin+ 
# #      90             21              4
# 
# summary( filter(sradat.select.f, Diagnosis == "T2D metformin-" )$Age );hist( filter(sradat.select.f, Diagnosis == "T2D metformin-" )$Age )
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 25.00   52.00   59.00   56.67   62.00   71.00
# 
# summary( filter(sradat.select.f, Diagnosis == "ND CTRL" )$Age );hist( filter(sradat.select.f, Diagnosis == "ND CTRL" )$Age )
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 21.00   34.25   44.00   43.07   51.75   67.00 

# Choose >= 50 in ND CTRL T2D metformin-


sel <- which(sradat.select$Age > 50 & sradat.select$Diagnosis %in% c("ND CTRL", "T2D metformin-")) # 82

sradat.select2 <- sradat.select[sel, ]
dim(sradat.select2) # 82 31
table(sradat.select2$Diagnosis)
# ND CTRL T2D metformin- 
#   52             30 

chn_runs <- sradat.select2$Run
length(chn_runs) # 82

# file for SRA download

write.csv(x = paste0("fastq-dump --outdir /scratch/user/lidd0026/forslund-t2d-chn/ft2d_1_meta_raw --skip-technical --readids --read-filter pass --dumpbase --clip --split-3 ",sort(chn_runs)), file = "forslund-t2d-CHN-testset-sra-runs-download.txt", quote = FALSE, row.names = FALSE)
write.csv(x = paste0("fastq-dump --outdir /cluster/jobs/lidd0026/forslund-t2d-chn/ft2d_1_meta_raw --skip-technical --readids --read-filter pass --dumpbase --clip --split-3 ",sort(chn_runs)), file = "forslund-t2d-CHN-testset-sra-runs-download-cluster.txt", quote = FALSE, row.names = FALSE)


#-------------------------

#### Forslund CHN-T2D - w/ Host-removal - compare no reads w/ vs w/out host-removal
#-------------------------

r1.clean <- read.table(file = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn/4_compare_read_counts/clean_r1_sequence_counts.txt", header = FALSE, sep = " ")
head(r1.clean)
names(r1.clean) <- c("sample","n_reads_cleaned")

r1.nonhost <- read.table(file = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn/4_compare_read_counts/nonhost_r1_sequence_counts.txt", header = FALSE, sep = " ")
head(r1.nonhost)
names(r1.nonhost) <- c("sample","n_reads_nonhost")

identical(r1.clean$sample, r1.nonhost$sample) # TRUE

r1.chn <- r1.clean
r1.chn$n_reads_nonhost <- r1.nonhost$n_reads_nonhost
r1.chn$percent_host_removed <- 100*(r1.chn$n_reads_cleaned - r1.chn$n_reads_nonhost)/r1.chn$n_reads_cleaned
hist(r1.chn$percent_host_removed);summary(r1.chn$percent_host_removed)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000000 0.0002930 0.0008462 0.0043117 0.0023843 0.1132800

r1.chn$n_reads_host_removed <- r1.chn$n_reads_cleaned - r1.chn$n_reads_nonhost

median(r1.chn$percent_host_removed) # 0.0008462543

# For the Chinese T2D cohort study samples contained a median (and interquartile range, IQR): 0.0008 % (0.0003 - 0.002 %) of reads that were classified as human sequences.

hist(r1.chn$n_reads_cleaned);summary(r1.chn$n_reads_cleaned)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 758416 10224587 12700891 14145936 19187333 28613123 

hist(r1.chn$n_reads_nonhost);summary(r1.chn$n_reads_nonhost)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 758416 10224381 12700702 14145399 19183919 28613098 

# For the Chinese T2D cohort study data we used median 12700702 (IQR: 10224381 - 19183919) cleaned non-host sequences

hist(r1.chn$n_reads_host_removed);summary(r1.chn$n_reads_host_removed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0    34.0   114.0   536.6   375.5 13216.0

# median 114 (IQR 34 - 376)

#-------------------------

#### Forslund-CHN-T2D - w/ Host-removal - read in superfocus - fxn potential outputs
#-------------------------

sampid <- sradat.select2$Run

superfocus_out_dir <- "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn/3_fxn_superfocus_copy_hostremoval"

list.dirs(superfocus_out_dir)
head( list.dirs(superfocus_out_dir) )

# # don't keep 1st two 
# ( results_dirs <- list.dirs(superfocus_out_dir)[-c(1,2)] )

# # don't keep 1st directory
( results_dirs <- list.dirs(superfocus_out_dir)[-c(1)] )

head(results_dirs)
# [1] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR341581"
# [2] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR341585"
# [3] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR341586"
# [4] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR341587"
# [5] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR341588"
# [6] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR341589"

names(results_dirs) <- gsub(pattern = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn/3_fxn_superfocus_copy_hostremoval/superfocus_out_", replacement = "", x = results_dirs)
head(results_dirs)
# SRR341581 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR341581" 
# SRR341585 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR341585" 
# SRR341586 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR341586" 
# SRR341587 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR341587" 
# SRR341588 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR341588" 
# SRR341589 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR341589"

length(results_dirs) # 82

sel <- which(names(results_dirs) %in% sampid) # qty 82
#results_dirs <- results_dirs[sel]

length( which(names(results_dirs) %in% sampid)) # 82

# check identical order
identical(sampid, names(results_dirs)) # TRUE
length(results_dirs) # 82


# In this data one Run corresponds to a single Sample_ID !!!

# collate results into a long-format table

sfx.long <- data.frame(sampleID=NA, subsys_L1=NA, subsys_L2=NA, subsys_L3=NA,fxn=NA,percent_abun=NA)

for (i in 1:length(sampid)) {
  #i<-1
  this_samp <- sampid[i]
  sel.folder <- grep(pattern = this_samp, x = results_dirs)
  this_folder <- results_dirs[sel.folder]
  
  #tab1 <- read_excel(path = paste0(this_folder,"/output_all_levels_and_function.xlsx"), skip = 4, col_names = TRUE)
  
  tab <- read.csv(file = paste0(this_folder,"/output_all_levels_and_function.xls"), sep = "\t", skip = 4 )
  # names(tab)
  # [1] "Subsystem.Level.1"                                                                         
  # [2] "Subsystem.Level.2"                                                                         
  # [3] "Subsystem.Level.3"                                                                         
  # [4] "Function"                                                                                  
  # [5] "X.scratch.pawsey1216.cliddicoat.ft2d_chn.2b_clean_hostremoval.SRR341581_non_host.1.fastq"  
  # [6] "X.scratch.pawsey1216.cliddicoat.ft2d_chn.2b_clean_hostremoval.SRR341581_non_host.1.fastq.."
  
  
  # [1] "Subsystem.Level.1"
  # [2] "Subsystem.Level.2"
  # [3] "Subsystem.Level.3"
  # [4] "Function"
  # [5] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq"
  # [6] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq"
  # [7] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq"
  # [8] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq"
  # [9] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq.." # this is %
  # [10] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq.." # this is %
  # [11] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq.." # this is %
  # [12] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq.." # this is %
  
  
  tab$sampid <- this_samp
  names(tab)
  
  #tab <- tab[,c(7,1,2,3,4,6)]
  
  # last column is sampid
  # take average of percentages
  
  #sel.col.percent <- grep(pattern = "R1.good.fastq..$", x = names(tab))
  sel.col.percent <- grep(pattern = "_non_host.1.fastq..$", x = names(tab))
  if (length(sel.col.percent)>1) {
    tab$percent_abun <- apply(X = tab[ ,sel.col.percent], MARGIN = 1, FUN = mean )
  } else {
    tab$percent_abun <- tab[ ,sel.col.percent]
  }
  
  # sum(tab$percent_abun) # 100
  # mean(tab$percent_abun) # 0.004338583
  
  names(sfx.long) # "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"    "percent_abun"
  # names(tab)
  # [1] "Subsystem.Level.1"
  # [2] "Subsystem.Level.2"
  # [3] "Subsystem.Level.3"
  # [4] "Function"
  # ...
  # [13] "sampid"
  # [14] "percent_abun"
  
  tab <- tab[ ,c("sampid","Subsystem.Level.1","Subsystem.Level.2","Subsystem.Level.3","Function","percent_abun")]
  names(tab) <- names(sfx.long)
  
  sfx.long <- rbind(sfx.long,tab)
  
  print(paste0("completed ",i," - sample ID: ",sampid[i]))
}


head(sfx.long)
# remove empty 1st row
sfx.long <- sfx.long[-1, ]
dim(sfx.long) # 709268      6
head(sfx.long)
#    sampleID                   subsys_L1                    subsys_L2                           subsys_L3
# 2 SRR341581 Amino Acids and Derivatives                            - Creatine and Creatinine Degradation
# 3 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine                Glycine Biosynthesis
# 4 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine      Glycine and Serine Utilization
# 5 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine      Glycine and Serine Utilization
# 6 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine      Glycine and Serine Utilization
# 7 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine             Glycine cleavage system
# fxn percent_abun
# 2                                                              Creatinine_amidohydrolase_(EC_3.5.2.10)  0.013708019
# 3                                                           L-threonine_3-dehydrogenase_(EC_1.1.1.103)  0.027416038
# 4                                                     D-3-phosphoglycerate_dehydrogenase_(EC_1.1.1.95)  0.019648161
# 5 L-serine_dehydratase,_beta_subunit_(EC_4.3.1.17)_/_L-serine_dehydratase,_alpha_subunit_(EC_4.3.1.17)  0.009595613
# 6                                                                   L-serine_dehydratase_(EC_4.3.1.17)  0.009595613
# 7                                                                   L-serine_dehydratase_(EC_4.3.1.17)  0.002741604


sfx.long$full_fxn_tax <- paste0(sfx.long$subsys_L1,"___", sfx.long$subsys_L2,"___", sfx.long$subsys_L3,"___", sfx.long$fxn)


## translate from long to wide format

names(sfx.long)
# "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"          "percent_abun" "full_fxn_tax"

sfx.wide <- dcast(sfx.long, formula = full_fxn_tax ~ sampleID, value.var = "percent_abun")
dim(sfx.wide) # 19352    83

sel.na <- which(is.na(sfx.wide),arr.ind = TRUE)
sfx.wide[sel.na] <- 0

# function taxonomy
full_fxn_names <- sfx.wide$full_fxn_tax

length(full_fxn_names) # 19352
length(unique(full_fxn_names)) # 19352

names(full_fxn_names) <- paste0("fxn_",c(1:length(full_fxn_names)))
head(full_fxn_names)
# fxn_1 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylaconitate_cis-trans_isomerase" 
# fxn_2 
# "Amino Acids and Derivatives___-___Amino acid racemase___4-hydroxyproline_epimerase_(EC_5.1.1.8)" 
# fxn_3 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)" 
# fxn_4 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic" 
# fxn_5 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_catabolic" 
# fxn_6 
# "Amino Acids and Derivatives___-___Amino acid racemase___Amino_acid_racemase_RacX" 


tax.fxn <- separate(sfx.wide, full_fxn_tax, c("subsys_L1", "subsys_L2", "subsys_L3", "fxn"), sep= "___", remove=TRUE)
# remove sample ids
tax.fxn <- tax.fxn[ ,-which(names(tax.fxn) %in% sampid)]

row.names(tax.fxn) <- names(full_fxn_names)


head(sfx.wide)

names(sfx.wide)
# [1] "full_fxn_tax" "SRR341581"    "SRR341585"    "SRR341586"    "SRR341587"    "SRR341588"    "SRR341589"    "SRR341599"    "SRR341600"    "SRR341601"    "SRR341602"   
# [12] "SRR341604"    "SRR341606"    "SRR341609"    "SRR341636"    "SRR341645"    "SRR341646"    "SRR341652"    "SRR341654"    "SRR341655"    "SRR341657"    "SRR341660"   
# [23] "SRR341661"    "SRR341663"    "SRR341664"    "SRR341665"    "SRR341669"    "SRR341670"    "SRR341673"    "SRR341674"    "SRR341675"    "SRR341676"    "SRR341681"   
# [34] "SRR341684"    "SRR341687"    "SRR341693"    "SRR341713"    "SRR413575"    "SRR413576"    "SRR413578"    "SRR413579"    "SRR413580"    "SRR413581"    "SRR413582"   
# [45] "SRR413584"    "SRR413585"    "SRR413587"    "SRR413592"    "SRR413593"    "SRR413594"    "SRR413597"    "SRR413598"    "SRR413599"    "SRR413600"    "SRR413601"   
# [56] "SRR413603"    "SRR413605"    "SRR413606"    "SRR413607"    "SRR413608"    "SRR413610"    "SRR413613"    "SRR413614"    "SRR413615"    "SRR413616"    "SRR413617"   
# [67] "SRR413618"    "SRR413619"    "SRR413620"    "SRR413621"    "SRR413623"    "SRR413625"    "SRR413626"    "SRR413634"    "SRR413637"    "SRR413642"    "SRR413652"   
# [78] "SRR413660"    "SRR413661"    "SRR413670"    "SRR413688"    "SRR413758"    "SRR413768"  

#names(sfx.wide) <- gsub(pattern = "-", replacement = "_", x = names(sfx.wide))

identical(as.character(full_fxn_names), sfx.wide$full_fxn_tax) # TRUE

row.names(sfx.wide) <- names(full_fxn_names)
sfx.wide <- sfx.wide[ ,-1]

names(sfx.wide)


head(sampid)
# "SRR341581" "SRR341585" "SRR341586" "SRR341587" "SRR341588" "SRR341589"

length(sampid) # 82

names(sampid) # NULL - in this case there is NOT an alternative sample name being used

# check alignment of sample IDs and sample names
identical(names(sfx.wide) , as.character(sampid)) # TRUE
# identical(names(sfx.wide) , as.character(gsub(pattern = "-",replacement = "_",x = sampid))) # FALSE
# length( names(sfx.wide) %in% as.character(gsub(pattern = "-",replacement = "_",x = sampid)) ) # 113 - i.e. matching but order different

#NOT RUN THIS TIME
#names(sfx.wide) <- names(sampid)


names(tax.fxn) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn"
dim(tax.fxn) # 19352     4

length(unique(tax.fxn$subsys_L1)) # 35
length(unique(tax.fxn$subsys_L2)) # 187
length(unique(tax.fxn$subsys_L3)) # 1105
length(unique(tax.fxn$fxn)) # 10248

#-------------------------

#### Forslund-CHN-T2D - w/ Host-removal - get into Phyloseq object
#-------------------------

# sfx.wide - is equiv to OTU table

# tax.fxn - is equiv to TAX table

# meta - is equiv to sample table

## Create 'taxonomyTable'
#  tax_table - Works on any character matrix. 
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
tax.m <- as.matrix( tax.fxn )
dim(tax.m) # 19352     4

TAX <- tax_table( tax.m )


## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
otu.m <- as.matrix( sfx.wide )
dim(otu.m)
# 19352    82

OTU <- otu_table(otu.m, taxa_are_rows = TRUE)


## Create a phyloseq object, merging OTU & TAX tables
phy = phyloseq(OTU, TAX)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 19352 taxa and 82 samples ]
# tax_table()   Taxonomy Table:    [ 19352 taxa by 4 taxonomic ranks ]

sample_names(phy)
# [1] "SRR341581" "SRR341585" "SRR341586" "SRR341587" "SRR341588" "SRR341589" "SRR341599" "SRR341600" "SRR341601" "SRR341602" "SRR341604" "SRR341606" "SRR341609"
# [14] "SRR341636" "SRR341645" "SRR341646" "SRR341652" "SRR341654" "SRR341655" "SRR341657" "SRR341660" "SRR341661" "SRR341663" "SRR341664" "SRR341665" "SRR341669"
# [27] "SRR341670" "SRR341673" "SRR341674" "SRR341675" "SRR341676" "SRR341681" "SRR341684" "SRR341687" "SRR341693" "SRR341713" "SRR413575" "SRR413576" "SRR413578"
# [40] "SRR413579" "SRR413580" "SRR413581" "SRR413582" "SRR413584" "SRR413585" "SRR413587" "SRR413592" "SRR413593" "SRR413594" "SRR413597" "SRR413598" "SRR413599"
# [53] "SRR413600" "SRR413601" "SRR413603" "SRR413605" "SRR413606" "SRR413607" "SRR413608" "SRR413610" "SRR413613" "SRR413614" "SRR413615" "SRR413616" "SRR413617"
# [66] "SRR413618" "SRR413619" "SRR413620" "SRR413621" "SRR413623" "SRR413625" "SRR413626" "SRR413634" "SRR413637" "SRR413642" "SRR413652" "SRR413660" "SRR413661"
# [79] "SRR413670" "SRR413688" "SRR413758" "SRR413768"

### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

samp <- sradat.select2

dim(samp) # 82 31

head(row.names(samp)) # "1" "5" "6" "7" "8" "9"

row.names(samp) <- samp$Run

identical(row.names(samp), sample_names(phy)) # TRUE

SAMP <- sample_data(samp)



### Combine SAMPDATA into phyloseq object
phy <- merge_phyloseq(phy, SAMP)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 19352 taxa and 82 samples ]
# sample_data() Sample Data:       [ 82 samples by 31 sample variables ]
# tax_table()   Taxonomy Table:    [ 19352 taxa by 4 taxonomic ranks ]

head(taxa_names(phy))
# "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"

head(phy@tax_table)
# Taxonomy Table:     [6 taxa by 4 taxonomic ranks]:
#   subsys_L1                     subsys_L2 subsys_L3             fxn                                            
# fxn_1 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "2-methylaconitate_cis-trans_isomerase"        
# fxn_2 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "4-hydroxyproline_epimerase_(EC_5.1.1.8)"      
# fxn_3 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)"                
# fxn_4 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic"
# fxn_5 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_catabolic"   
# fxn_6 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Amino_acid_racemase_RacX"    


getwd()  # "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R"

saveRDS(object = phy, file = "phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal.RDS")


str(samp)
# 'data.frame':	82 obs. of  31 variables:
table( samp$gender )
# female   male 
#     41     41
sel <- which(samp$Diagnosis == "T2D metformin-")
table( samp$gender[sel] )
# female   male 
#     17     13
summary( samp$Age[ which(samp$Diagnosis == "T2D metformin-" & samp$gender == "female")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.00   56.00   61.00   60.29   63.00   71.00
length( samp$Age[ which(samp$Diagnosis == "T2D metformin-" & samp$gender == "female")] )
# [1] 17
summary( samp$Age[ which(samp$Diagnosis == "T2D metformin-" & samp$gender == "male")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.00   53.00   58.00   61.15   70.00   75.00 
length( samp$Age[ which(samp$Diagnosis == "T2D metformin-" & samp$gender == "male")] )
# [1] 13

sel <- which(samp$Diagnosis == "ND CTRL")
table( samp$gender[sel] )
# female   male 
#     24     28 
summary( samp$Age[ which(samp$Diagnosis == "ND CTRL" & samp$gender == "female")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.00   53.00   56.00   56.58   60.00   67.00 
length( samp$Age[ which(samp$Diagnosis == "ND CTRL" & samp$gender == "female")] )
# [1] 24
summary( samp$Age[ which(samp$Diagnosis == "ND CTRL" & samp$gender == "male")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 52.00   53.75   56.00   58.71   62.25   74.00 
length( samp$Age[ which(samp$Diagnosis == "ND CTRL" & samp$gender == "male")] )
# [1] 28

# T2D met- (total n = 30 total; females n = 17, ages 51-71; males n = 13, ages 51-75)
# Normal (total n = 52 total; females n = 24, ages 51-67; males n = 28, ages 52-74)


# get stats??
head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 82
NonZeroFxns

mean(NonZeroFxns) # 8649.61
sd(NonZeroFxns) # 4213.538


#-------------------------

#### Forslund-CHN-T2D - w/ Host removal - COPY of R code to run CPP steps on HPC
#    1) build reaction search - get reactions and compounds
#    2) get cpd rel abun per sample
#    3) collate compounds for each sample
#-------------------------

# # # # # # # # # # # # #
# #
# # R script for cpp3d
# # - build reaction search in parallel - get_reactions & compounds
# # - get cpd rel abun per sample
# # - collate_compounds
# #
# # Forslund T2D - CHN - Host-removal
# # Craig Liddicoat - Flinders University
# # Running on Pawsey Setonix
# # # # # # # # # # # # #
# 
# # Add a new path
# .libPaths(c("/software/projects/pawsey1216/cliddicoat/setonix/2024.05/r/4.4.1",
#             "/software/projects/pawsey1216/cliddicoat/setonix/2024.05/r/4.3", .libPaths()))
# 
# R.Version()
# 
# # load packages
# #library(readxl); packageVersion("readxl")
# library(parallel); packageVersion("parallel")
# library(doParallel); packageVersion("doParallel")
# library(dplyr); packageVersion("dplyr")
# library(stringr); packageVersion("stringr")
# library(phyloseq); packageVersion("phyloseq") # '1.44.0'
# 
# message("\n# establish folders and input files")
# 
# message("\nworkdir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/5_rcode_cpp'")
# workdir <- "/scratch/pawsey1216/cliddicoat/ft2d_chn/5_rcode_cpp"
# message("\nsetwd(workdir)")
# setwd(workdir)
# message("\ntemp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/5_rcode_cpp/working'")
# temp_dir <- "/scratch/pawsey1216/cliddicoat/ft2d_chn/5_rcode_cpp/working"
# 
# message("\nthis_study <- '-forslund-t2d-chn-hostremoval-pawsey'")
# this_study <- "-forslund-t2d-chn-hostremoval-pawsey"
# message("\nphy <- readRDS('phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal.RDS')")
# phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal.RDS")
# 
# 
# subsys.lut <- readRDS("subsys.lut.RDS")
# rxns.lut <- readRDS("rxns.lut.RDS")
# rxn_pathways.lut <- readRDS("rxn_pathways.lut.RDS")
# compounds.lut <- readRDS("compounds.lut.RDS")
# 
# 
# 
# message("\n### 1) build reaction search in parallel - get_reactions & compounds")
# message("\n# # # # # # # # # #")
# message("\ndf.tax <- as.data.frame(phy@tax_table)")
# df.tax <- as.data.frame(phy@tax_table)
# message("\nhead(row.names(df.tax))")
# head(row.names(df.tax))
# message("\ndim(df.tax)")
# dim(df.tax)
# 
# 
# get_rxns_and_compounds_indiv <- function( df.tax, subsys.lut, rxns.lut, rxn_pathways.lut ) {
#   
#   rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
#   rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
#   
#   sub1 <- df.tax$subsys_L1[i]
#   sub2 <- df.tax$subsys_L2[i]
#   sub3 <- df.tax$subsys_L3[i]
#   
#   fxn.temp <- df.tax$fxn[i]
#   fxn.superfocus.rowlabel <- row.names(df.tax)[i]
#   
#   # store results corresponding to each Superfocus row
#   fxn.list <- list()
#   fxn.list[[ fxn.superfocus.rowlabel  ]] <- list()
#   
#   # check for multiple functions/reactions?
#   flag1 <- grepl(pattern = "_/_|/", x = fxn.temp)
#   flag2 <- grepl(pattern = "_@_", x = fxn.temp)
#   if (!any(flag1,flag2)==TRUE) {
#     # no multiples
#     fxns <- fxn.temp
#   } else if (flag1==TRUE) {
#     fxns <- unlist( strsplit(fxn.temp, split = "_/_") )  ###### WHAT ABOUT SPLIT FOR "/" WITHOUT UNDERSCORES ??
#   } else {
#     fxns <- unlist( strsplit(fxn.temp, split = "_@_") )
#   }
#   # remove underscores
#   ( fxns <- gsub(pattern = "_", replacement = " ", x = fxns) )
#   
#   # process each fxn & store attributes
#   df.fxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=1:length(fxns),`f__in`=fxns, matching_method=NA, rxns=NA)
#   
#   # Identify '/' separators with no '_'  ??
#   
#   for (f in 1:length(fxns)) {  # this accounts for multiple functions/reactions reported in Superfocus outputs
#     #f<-1
#     #f<-2
#     f.in <- fxns[f]
#     
#     # these concatenated expressions will be used to look for exact match using hierarchy in ModelSEED Subsystem table
#     full_hier_target <- paste0(sub1,"__",sub2,"__",sub3,"__",f.in)
#     full_hier_list <- paste0(subsys.lut$Class,"__",subsys.lut$Subclass,"__",gsub("_"," ",subsys.lut$Name),"__",subsys.lut$Role)
#     
#     ## data cleaning
#     
#     # trim off '_#' and '_##' tags
#     trim_nchar <- str_locate(string = f.in, pattern = " # | ## ")[1]
#     if (!is.na(trim_nchar) & length(trim_nchar)==1) {
#       f.in <- substring(text = f.in , first = 1, last = trim_nchar-1)
#     }
#     
#     # Eliminate unwanted parsing of regular expressions: '[', ']','***', '(', ')'
#     f.in <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\} ", replacement ="." , x = f.in) # used later
#     
#     #rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
#     #rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
#     
#     full_hier_target <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_target)
#     full_hier_list <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_list)
#     
#     sel.rx <- grep(pattern = full_hier_target, x = full_hier_list)
#     
#     ## ALTERNATIVE #1 == FULL HIERACHICAL MATCH
#     if (length(sel.rx)>=1) {
#       df.fxns$matching_method[f] <- "Exact hierachy match"
#       df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
#       
#     } else if (str_detect(string = fxns[f], pattern = " \\(EC ")) {  ## ALTERNATIVE #2 == MATCHING ECs
#       # search by EC id if present
#       
#       f.in <- fxns[f] # this goes back to string with brackets for EC
#       ## LOOK FOR MULTIPLE ECs ?
#       
#       how_many_ECs <- str_count(string = f.in, pattern = "\\(EC.*?\\)")
#       
#       ECs <- as.character( str_extract_all(string = f.in, pattern = "\\(EC.*?\\)", simplify = TRUE) )
#       #class(ECs)
#       ECs <- gsub(pattern = "\\(EC |\\)", replacement = "", x = ECs)
#       ECs.collapse <- paste0(ECs, collapse = "|")
#       
#       sel.rx <- which(rxns.lut$ec_numbers == ECs.collapse)
#       
#       if (length(how_many_ECs)==0 | length(ECs)==0) {
#         # there was a glitch, database typo, or some error in identifying the EC number
#         df.fxns$matching_method[f] <- "No match found"
#         df.fxns$rxns[f] <- NA
#         
#       } else if (length(sel.rx)>=1) {
#         # combined EC hits identified
#         df.fxns$matching_method[f] <- "EC number"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if (length(which(rxns.lut$ec_numbers %in% ECs)) >=1) {
#         # treat EC hits individually
#         sel.rx <- which(rxns.lut$ec_numbers %in% ECs) # look 1st where ECs are exact matches for EC numbers in Reactions lookup table
#         
#         df.fxns$matching_method[f] <- "EC number"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if (length(grep(pattern = ECs, x = rxns.lut$ec_numbers)) >=1) {
#         # this allows EC to be part of a combination of EC numbers that are listed in Reactions lookup table
#         sel.rx <- grep(pattern = ECs, x = rxns.lut$ec_numbers)
#         
#         df.fxns$matching_method[f] <- "EC number"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else {
#         # it had an EC number but couldn't find a match in the EC numbers listed in Reaction lookup table
#         df.fxns$matching_method[f] <- "No match found"
#         df.fxns$rxns[f] <- NA
#         
#       }
#       # END EC matching
#       
#       
#     } else {  ## ALTERNATIVE 3 == FXN NAME MATCHING
#       ## otherwise attempt to match function name - a) first look for exact matches   ########## then b) closest match above a threshold
#       # 1. 'reactions' table by name: rxns.lut$name
#       # 2. 'reactions' table by aliases: rxns.lut$aliases
#       # 3. 'Model SEED Subsystems' table by Role: subsys.lut$Role
#       # 4. 'Unique_ModelSEED_Reaction_Pathways' table by External ID: rxn_pathways.lut$External_rxn_name
#       
#       if ( length( grep(pattern = f.in, x = rxns.lut$name) )>=1 ) {
#         # 1a - exact match - rxns.lut$name
#         sel.rx <- grep(pattern = f.in, x = rxns.lut$name)
#         #rxns.lut$name[sel.rx]
#         df.fxns$matching_method[f] <- "Matched Reactions name"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if ( length( grep(pattern = f.in, x = rxns.lut$aliases) )>=1 ) {
#         # 2a - exact match - rxns.lut$aliases
#         sel.rx <- grep(pattern = f.in, x = rxns.lut$aliases)
#         #rxns.lut$aliases[sel.rx]
#         #rxns.lut$name[sel.rx]
#         
#         df.fxns$matching_method[f] <- "Matched Reactions aliases"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if ( length( grep(pattern = f.in, x = subsys.lut$Role) )>=1 ) {
#         # 3a - exact match - subsys.lut$Role
#         sel.rx <- grep(pattern = f.in, x = subsys.lut$Role)
#         #subsys.lut$Role[sel.rx]
#         #subsys.lut$Reaction[sel.rx]
#         
#         df.fxns$matching_method[f] <- "Matched Subsytem role"
#         df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
#         
#       } else if ( length( grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name) )>=1 ) {
#         # 4a - exact match - rxn_pathways.lut$External_rxn_name
#         sel.rx <- grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name)
#         
#         df.fxns$matching_method[f] <- "Matched ModelSEED Reaction pathways"
#         df.fxns$rxns[f] <- paste0( unique(rxn_pathways.lut$rxn_id[sel.rx]), collapse = ";")
#         
#         
#       } else {
#         df.fxns$matching_method[f] <- "No match found"
#         df.fxns$rxns[f] <- NA
#         
#       }
#       
#       ## DON'T RUN PARTIAL MATCHING AT THIS STAGE
#       
#       
#     } # END function - reaction search
#     
#     #fxn.list[[ fxn.superfocus.rowlabel  ]][[ f ]][[ "fxns" ]] <- df.fxns
#     
#     print(paste0("completed fxn ", f))
#     
#     
#     ## now investigate these reactions ...
#     # Reactions lookup table: 
#     # - "equation": Definition of reaction expressed using compound IDs and after protonation
#     # Compounds lookup table:
#     # - "formula": Standard chemical format (using Hill system) in protonated form to match reported charge
#     #df.fxns
#     
#     
#     #if (df.fxns$matching_method == "No match found") {
#     if (df.fxns$rxns[f] == "" | is.na(df.fxns$rxns[f])) {
#       
#       df.Rxns <- NA
#       df.Compounds <- NA
#       
#     } else { # reaction(s) were identified
#       
#       # consider reactions for this f.in only (possibly > 1 f.in per Superfocus row)
#       f.in.rxns <- unique(unlist(str_split(string = df.fxns$rxns[f], pattern = ";")))
#       
#       df.Rxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel, f=f, f__in=f.in,rxn_id= f.in.rxns,
#                             rxn_name=NA, rxn_eqn=NA, rxn_defn=NA,compds=NA,compd_coef=NA, chem_formx=NA )
#       
#       for (r in 1:dim(df.Rxns)[1]) {
#         #r<-1
#         #this_rxn <- "rxn00004"
#         this_rxn <- df.Rxns$rxn_id[r]
#         sel <- which(rxns.lut$id == this_rxn)
#         ( df.Rxns$rxn_name[r] <- rxns.lut$name[sel] )
#         ( df.Rxns$rxn_eqn[r] <- rxns.lut$equation[sel] )
#         ( df.Rxns$rxn_defn[r] <- rxns.lut$definition[sel] )
#         
#         # extract compound info
#         
#         #df.Rxns$rxn_eqn[r]
#         #[1] "(1) cpd00010[0] + (1) cpd29672[0] <=> (1) cpd00045[0] + (1) cpd11493[0]"
#         #[1] "(45) cpd00144[0] + (45) cpd00175[0] <=> (45) cpd00014[0] + (45) cpd00091[0] + (1) cpd15634[0]"
#         
#         ( compds.idx <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd")[[1]][,"start"] )
#         # 5 23 43 61
#         # 6 25 46 65 83
#         
#         ( compds <- as.character( str_extract_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd.....", simplify = TRUE) ) )
#         # "cpd00010" "cpd29672" "cpd00045" "cpd11493"
#         
#         if (length(compds)>=1) {
#           
#           df.Rxns$compds[r] <- paste0(compds, collapse = ";")
#           
#           ## get compound coefficients?
#           start_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\(")[[1]][,"start"]
#           end_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\)")[[1]][,"start"]
#           ( compd.coeff <- as.numeric( substring(text = df.Rxns$rxn_eqn[r], first = start_brackets+1, last = end_brackets-1)) )
#           
#           df.Rxns$compd_coef[r] <- paste0(compd.coeff, collapse = ";")
#           
#           # get formulas of compounds
#           
#           formx <-filter(compounds.lut, id %in% compds )
#           row.names(formx) <- formx$id
#           ( formx.char <- formx[compds, ]$formula )
#           # "C21H32N7O16P3S" "HOR"            "C10H11N5O10P2"  "C11H22N2O7PRS" 
#           # "C15H19N2O18P2"      "C17H25N3O17P2"      "C9H12N2O12P2"       "C9H11N2O9P"         "C630H945N45O630P45"
#           # "C7H7O7" "H2O"    "C7H5O6"
#           df.Rxns$chem_formx[r] <- paste0(formx.char, collapse = ";")
#           
#           ( compd.names <- formx[compds, ]$name )
#           # "2-methyl-trans-aconitate" "cis-2-Methylaconitate"
#           
#           temp.df.Compounds <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=f, f__in=f.in,rxn_id= f.in.rxns[r], 
#                                           cpd_id=compds, cpd_name=compd.names, cpd_form=formx.char, cpd_molar_prop=compd.coeff #, 
#                                           #OC_x=OC_ratio, HC_y=HC_ratio , NC_z=NC_ratio 
#           )
#           
#         } else {
#           # No specified reaction equation or chemical formula info
#           df.Rxns$compds[r] <- NA
#           df.Rxns$compd_coef[r] <- NA
#           df.Rxns$chem_formx[r] <- NA
#           
#           temp.df.Compounds <- NA
#           
#         }
#         
#         if (r==1) { df.Compounds <- temp.df.Compounds }
#         
#         if (r>1 & is.data.frame(df.Compounds) & is.data.frame(temp.df.Compounds)) { df.Compounds <- rbind(df.Compounds, temp.df.Compounds) }
#         
#         # clean up - if there are additional reactions?
#         temp.df.Compounds <- NA
#         
#       } # END loop for r - rxn_id's per f/f.in
#       
#     } # END else loop when reactions identified
#     
#     # store results corresponding to each sub-reaction of each Superfocus row
#     fxn.list[[ fxn.superfocus.rowlabel  ]][[ "fxns" ]] <- df.fxns
#     
#     if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]] <- list() } # set this only once
#     fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]][[ f ]] <- df.Rxns
#     
#     if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]] <- list() } # set this only once
#     fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]][[ f ]] <- df.Compounds
#     
#     
#   } # END loop - f in 1:length(fxns)) - to account for multiple functions/reactions reported in each row of Superfocus outputs
#   
#   saveRDS(object = fxn.list, file = paste0(temp_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds") )
#   
# } # END function to be run in parallel for each superfocus row
# 
# 
# # # # # # # # # # # # # # # # # # #
# 
# no_forks <- 8
# 
# # this makes clusters on Unix-like system (may need to use other alternative for Windows)
# cl<-makeForkCluster(nnodes = no_forks)      # no of nodes will depend on your HPC facility
# registerDoParallel(cl)
# 
# foreach(i=1:dim(df.tax)[1] , .packages=c('stringr', 'dplyr')) %dopar%  #
#   get_rxns_and_compounds_indiv( df.tax=df.tax, subsys.lut=subsys.lut, rxns.lut=rxns.lut, rxn_pathways.lut=rxn_pathways.lut )
# 
# stopCluster(cl)
# 
# 
# message("\n## assemble results")
# 
# message("\n(num_results_files <- dim(df.tax)[1])")
# (num_results_files <- dim(df.tax)[1])
# 
# # assemble all compound data outputs
# # start with blank row
# 
# df.out <- data.frame(superfocus_fxn=NA, f=NA, f__in=NA, rxn_id=NA, cpd_id=NA, cpd_name=NA, cpd_form=NA, cpd_molar_prop=NA )
# 
# for (i in 1:num_results_files) {
#   fxn.superfocus.rowlabel <- row.names(df.tax)[i]
#   temp <- readRDS(paste0(temp_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))
#   
#   f_no <- length( temp[[1]][["compounds"]] )
#   
#   for (f in 1:f_no) {
#     #f<-2
#     # only add non-NA results
#     if (is.data.frame( temp[[1]][["compounds"]][[f]] )) {
#       
#       df.temp <- temp[[1]][["compounds"]][[f]]
#       ok <- complete.cases(df.temp)
#       df.temp <- df.temp[ which(ok==TRUE), ] # updated version will include some compounds with vK coordinates that are NA. vK coordinates are considered later
#       df.out <- rbind(df.out,df.temp)
#     }
#   }
#   print(paste0("added df ",i," of ",num_results_files ))
#   
# }
# 
# 
# message("\nstr(df.out)")
# str(df.out)
# 
# 
# saveRDS(object = df.out, file = paste0("df.out--get_rxns_and_compounds_indiv-",this_study,".RDS"))
# 
# # remove NA first row
# message("\nhead(df.out)")
# head(df.out)
# 
# df.out <- df.out[-1, ]
# 
# message("\ndim(df.out)")
# dim(df.out)
# 
# 
# message("\n## normalise molar_prop to cpd_relabun so total of 1 per superfocus function")
# 
# df.out$cpd_molar_prop_norm <- NA
# 
# message("\nlength(unique(df.out$superfocus_fxn))")
# length(unique(df.out$superfocus_fxn))
# 
# message("\nphy")
# phy
# 
# message("\n% of functions represented - with compound information")
# 100*(length(unique(df.out$superfocus_fxn)) / ntaxa(phy))
# 
# 
# fxns_found <- unique(df.out$superfocus_fxn)
# 
# for (k in 1:length(fxns_found)) {
#   #k<-1
#   this_fxn <- fxns_found[k]
#   sel <- which(df.out$superfocus_fxn == this_fxn)
#   
#   sum_molar_prop <- sum( df.out$cpd_molar_prop[sel], na.rm = TRUE)
#   # calculate 
#   
#   df.out$cpd_molar_prop_norm[sel] <- df.out$cpd_molar_prop[sel]/sum_molar_prop
#   
#   print(paste0("completed ",k))
#   
# }
# 
# message("\nsum(df.out$cpd_molar_prop_norm)")
# sum(df.out$cpd_molar_prop_norm)
# 
# message("\nsample_sums(phy)")
# sample_sums(phy)
# 
# message("\ngetwd()")
# getwd()
# 
# saveRDS(object = df.out, file = paste0("df.out--tidy-compounds_indiv-cpp3d-",this_study,".RDS"))
# 
# 
# 
# message("\n### 2) get cpd rel abun per sample")
# message("\n# # # # # # # # # #")
# 
# 
# df.OTU <- as.data.frame( phy@otu_table ) # this is Superfocus functional relative abundance data represented in phyloseq OTU abundance table
# message("\ndim(df.OTU)")
# dim(df.OTU)
# 
# 
# get_cpd_relabun_per_sample <- function(phy_in, dat.cpd) {
#   
#   this_samp <- sample_names(phy_in)[i]
#   df.OTU <- as.data.frame( phy_in@otu_table[ ,this_samp] )
#   
#   dat.cpd$sample <- this_samp
#   
#   dat.cpd$cpd_rel_abun_norm <- NA
#   
#   fxns_all <- row.names(df.OTU)
#   
#   for (k in 1:length(fxns_all)) {
#     #k<-1
#     this_fxn <- fxns_all[k]
#     sel <- which(dat.cpd$superfocus_fxn == this_fxn)
#     
#     if (length(sel)>=1) {
#       dat.cpd$cpd_rel_abun_norm[sel] <- df.OTU[this_fxn, ]*dat.cpd$cpd_molar_prop_norm[sel]
#       
#     }
#   } # END rel abun values for all relevant functions added
#   
#   saveRDS(object = dat.cpd, file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
#   
# } # END
# 
# 
# no_forks <- 8
# 
# # this makes clusters on Unix-like system
# cl<-makeForkCluster(nnodes = no_forks)      # no of nodes will depend on your HPC facility
# registerDoParallel(cl)
# 
# foreach(i=1: length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
#   get_cpd_relabun_per_sample( phy_in = phy, dat.cpd = df.out)
# 
# stopCluster(cl)
# 
# 
# message("\n## assemble results")
# 
# # output 1
# i<-1
# this_samp <- sample_names(phy)[i]
# dat <- readRDS( file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
# head(dat)
# 
# for ( i in 2:length(sample_names(phy)) ) {
#   this_samp <- sample_names(phy)[i]
#   temp <- readRDS( file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
#   dat <- rbind(dat, temp)
#   print(paste0("completed ",i))
# }
# 
# 
# saveRDS(object = dat, file = paste0("dat.cpd-long-all-samps-cpp3d-",this_study,".rds") )
# 
# rm(temp)
# 
# message("\nstr(dat)")
# str(dat)
# 
# message("\nsum(dat$cpd_rel_abun_norm)")
# sum(dat$cpd_rel_abun_norm)
# 
# message("\naverage functional relative abundance per sample")
# message("\nsum(dat$cpd_rel_abun_norm)/nsamples(phy)")
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# 
# message("\nnames(dat)")
# names(dat)
# 
# message("\nlength(unique(dat$cpd_id))")
# length(unique(dat$cpd_id))
# 
# 
# 
# 
# message("\n### 3) collate_compounds within each sample")
# message("\n# # # # # # # # # #")
# 
# 
# unique_cpd <- unique(dat$cpd_id)
# samp_names <- sample_names(phy)
# 
# 
# collate_compounds <- function(dat.cpd, unique_cpd, samp) {
#   #i<-1
#   #samp = samp_names[i]
#   #dat.cpd = dat[which(dat$sample == samp_names[i]), ]
#   
#   this_samp <- samp
#   
#   cpd_data <- data.frame(cpd_id = unique_cpd, sample=this_samp, cpd_rel_abun=NA)
#   
#   for (c in 1:length(unique_cpd)) {
#     #c<-1
#     this_cpd <- unique_cpd[c]
#     sel.cpd <- which(dat.cpd$cpd_id == this_cpd)
#     
#     if (length(sel.cpd) >=1) {
#       cpd_data$cpd_rel_abun[c] <- sum(dat.cpd$cpd_rel_abun_norm[sel.cpd])
#     }
#     
#   } # END all compounds
#   
#   saveRDS(object = cpd_data, file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
#   
# } # END
# 
# 
# 
# no_forks <- 4
# 
# # this makes clusters on Unix-like system
# cl<-makeForkCluster(nnodes = no_forks)   # no of nodes will depend on your HPC facility
# registerDoParallel(cl)
# 
# foreach(i=1:length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
#   collate_compounds(dat.cpd = dat[which(dat$sample == samp_names[i]), ], unique_cpd = unique_cpd, samp = samp_names[i])
# 
# stopCluster(cl)
# 
# 
# message("\n## assemble results")
# 
# # output 1
# i<-1
# this_samp <- sample_names(phy)[i]
# dat.cpd.collate <- readRDS( file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
# head(dat.cpd.collate)
# 
# for ( i in 2:length(sample_names(phy)) ) {
#   this_samp <- sample_names(phy)[i]
#   temp <- readRDS( file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
#   
#   dat.cpd.collate <- rbind(dat.cpd.collate, temp)
#   
#   print(paste0("completed ",i))
# }
# 
# 
# message("\nstr(dat.cpd.collate)")
# str(dat.cpd.collate)
# 
# message("\nsum(dat.cpd.collate$cpd_rel_abun)")
# sum(dat.cpd.collate$cpd_rel_abun)
# 
# message("\nsum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))")
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# 
# saveRDS(object = dat.cpd.collate, file = paste0("dat.cpd.collate-all-samps-cpp3d-",this_study,".rds" ))
# 
# # END


#-------------------------

#### Forslund-CHN-T2D - w/ Host-removal - COPY of OUTOUTS from R code after running CPP steps on HPC
#-------------------------

# $platform
# [1] "x86_64-pc-linux-gnu"
# 
# $arch
# [1] "x86_64"
# 
# $os
# [1] "linux-gnu"
# 
# $system
# [1] "x86_64, linux-gnu"
# 
# $status
# [1] ""
# 
# $major
# [1] "4"
# 
# $minor
# [1] "4.1"
# 
# $year
# [1] "2024"
# 
# $month
# [1] "06"
# 
# $day
# [1] "14"
# 
# $`svn rev`
# [1] "86737"
# 
# $language
# [1] "R"
# 
# $version.string
# [1] "R version 4.4.1 (2024-06-14)"
# 
# $nickname
# [1] "Race for Your Life"
# 
# [1] ‘4.4.1’
# Loading required package: foreach
# Loading required package: iterators
# [1] ‘1.0.17’
# 
# Attaching package: ‘dplyr’
# 
# The following objects are masked from ‘package:stats’:
#   
#   filter, lag
# 
# The following objects are masked from ‘package:base’:
#   
#   intersect, setdiff, setequal, union
# 
# [1] ‘1.1.4’
# [1] ‘1.5.2’
# [1] ‘1.46.0’
# 
# # establish folders and input files
# 
# workdir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/5_rcode_cpp'
# 
# setwd(workdir)
# 
# temp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/5_rcode_cpp/working'
# 
# this_study <- '-forslund-t2d-chn-hostremoval-pawsey'
# 
# phy <- readRDS('phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal.RDS')
# 
# ### 1) build reaction search in parallel - get_reactions & compounds
# 
# # # # # # # # # # #
# 
# df.tax <- as.data.frame(phy@tax_table)
# 
# head(row.names(df.tax))
# [1] "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"
# 
# dim(df.tax)
# [1] 19352     4
# 
# 
# 
# ## assemble results
# 
# (num_results_files <- dim(df.tax)[1])
# [1] 19352
# [1] "added df 1 of 19352"
# [1] "added df 2 of 19352"
# etc...
# [1] "added df 19350 of 19352"
# [1] "added df 19351 of 19352"
# [1] "added df 19352 of 19352"
# 
# str(df.out)
# 'data.frame':	554563 obs. of  8 variables:
#   $ superfocus_fxn: chr  NA "fxn_2" "fxn_2" "fxn_3" ...
# $ f             : int  NA 1 1 1 1 1 1 1 1 1 ...
# $ f__in         : chr  NA "4-hydroxyproline epimerase (EC 5.1.1.8)" "4-hydroxyproline epimerase (EC 5.1.1.8)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id        : chr  NA "rxn02360" "rxn02360" "rxn00283" ...
# $ cpd_id        : chr  NA "cpd00851" "cpd02175" "cpd00035" ...
# $ cpd_name      : chr  NA "trans-4-Hydroxy-L-proline" "cis-4-Hydroxy-D-proline" "L-Alanine" ...
# $ cpd_form      : chr  NA "C5H9NO3" "C5H9NO3" "C3H7NO2" ...
# $ cpd_molar_prop: num  NA 1 1 1 1 1 1 1 1 1 ...
# 
# head(df.out)
# superfocus_fxn  f                                   f__in   rxn_id   cpd_id
# 1           <NA> NA                                    <NA>     <NA>     <NA>
#   2          fxn_2  1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd00851
# 3          fxn_2  1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd02175
# 4          fxn_3  1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00035
# 5          fxn_3  1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00117
# 6          fxn_3  1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00035
# cpd_name cpd_form cpd_molar_prop
# 1                      <NA>     <NA>             NA
# 2 trans-4-Hydroxy-L-proline  C5H9NO3              1
# 3   cis-4-Hydroxy-D-proline  C5H9NO3              1
# 4                 L-Alanine  C3H7NO2              1
# 5                 D-Alanine  C3H7NO2              1
# 6                 L-Alanine  C3H7NO2              1
# 
# dim(df.out)
# [1] 554562      8
# 
# ## normalise molar_prop to cpd_relabun so total of 1 per superfocus function
# 
# length(unique(df.out$superfocus_fxn))
# [1] 10577
# 
# phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 19352 taxa and 82 samples ]
# sample_data() Sample Data:       [ 82 samples by 31 sample variables ]
# tax_table()   Taxonomy Table:    [ 19352 taxa by 4 taxonomic ranks ]
# 
# % of functions represented - with compound information
# [1] 54.65585
# [1] "completed 1"
# [1] "completed 2"
# [1] "completed 3"
# [1] "completed 4"
# [1] "completed 5"
# etc ...
# [1] "completed 10575"
# [1] "completed 10576"
# [1] "completed 10577"
# 
# sum(df.out$cpd_molar_prop_norm)
# [1] 10577
# 
# sample_sums(phy)
# SRR341581 SRR341585 SRR341586 SRR341587 SRR341588 SRR341589 SRR341599 SRR341600 
# 100       100       100       100       100       100       100       100 
# SRR341601 SRR341602 SRR341604 SRR341606 SRR341609 SRR341636 SRR341645 SRR341646 
# 100       100       100       100       100       100       100       100 
# SRR341652 SRR341654 SRR341655 SRR341657 SRR341660 SRR341661 SRR341663 SRR341664 
# 100       100       100       100       100       100       100       100 
# SRR341665 SRR341669 SRR341670 SRR341673 SRR341674 SRR341675 SRR341676 SRR341681 
# 100       100       100       100       100       100       100       100 
# SRR341684 SRR341687 SRR341693 SRR341713 SRR413575 SRR413576 SRR413578 SRR413579 
# 100       100       100       100       100       100       100       100 
# SRR413580 SRR413581 SRR413582 SRR413584 SRR413585 SRR413587 SRR413592 SRR413593 
# 100       100       100       100       100       100       100       100 
# SRR413594 SRR413597 SRR413598 SRR413599 SRR413600 SRR413601 SRR413603 SRR413605 
# 100       100       100       100       100       100       100       100 
# SRR413606 SRR413607 SRR413608 SRR413610 SRR413613 SRR413614 SRR413615 SRR413616 
# 100       100       100       100       100       100       100       100 
# SRR413617 SRR413618 SRR413619 SRR413620 SRR413621 SRR413623 SRR413625 SRR413626 
# 100       100       100       100       100       100       100       100 
# SRR413634 SRR413637 SRR413642 SRR413652 SRR413660 SRR413661 SRR413670 SRR413688 
# 100       100       100       100       100       100       100       100 
# SRR413758 SRR413768 
# 100       100 
# 
# getwd()
# [1] "/scratch/pawsey1216/cliddicoat/ft2d_chn/5_rcode_cpp"
# 
# ### 2) get cpd rel abun per sample
# 
# # # # # # # # # # #
# 
# dim(df.OTU)
# [1] 19352    82
# 
# 
# 
# 
# ## assemble results
# superfocus_fxn f                                   f__in   rxn_id   cpd_id
# 2          fxn_2 1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd00851
# 3          fxn_2 1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd02175
# 4          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00035
# 5          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00117
# 6          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00035
# 7          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00117
# cpd_name cpd_form cpd_molar_prop cpd_molar_prop_norm
# 2 trans-4-Hydroxy-L-proline  C5H9NO3              1           0.5000000
# 3   cis-4-Hydroxy-D-proline  C5H9NO3              1           0.5000000
# 4                 L-Alanine  C3H7NO2              1           0.1666667
# 5                 D-Alanine  C3H7NO2              1           0.1666667
# 6                 L-Alanine  C3H7NO2              1           0.1666667
# 7                 D-Alanine  C3H7NO2              1           0.1666667
# sample cpd_rel_abun_norm
# 2 SRR341581                 0
# 3 SRR341581                 0
# 4 SRR341581                 0
# 5 SRR341581                 0
# 6 SRR341581                 0
# 7 SRR341581                 0
# [1] "completed 2"
# [1] "completed 3"
# [1] "completed 4"
# etc ...
# [1] "completed 80"
# [1] "completed 81"
# [1] "completed 82"
# 
# str(dat)
# 'data.frame':	45474084 obs. of  11 variables:
#   $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_3" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "4-hydroxyproline epimerase (EC 5.1.1.8)" "4-hydroxyproline epimerase (EC 5.1.1.8)" "Alanine racemase (EC 5.1.1.1)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id             : chr  "rxn02360" "rxn02360" "rxn00283" "rxn00283" ...
# $ cpd_id             : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ cpd_name           : chr  "trans-4-Hydroxy-L-proline" "cis-4-Hydroxy-D-proline" "L-Alanine" "D-Alanine" ...
# $ cpd_form           : chr  "C5H9NO3" "C5H9NO3" "C3H7NO2" "C3H7NO2" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.167 0.167 0.167 ...
# $ sample             : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun_norm  : num  0 0 0 0 0 0 0 0 0 0 ...
# 
# sum(dat$cpd_rel_abun_norm)
# [1] 5388.071
# 
# average functional relative abundance per sample
# 
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# [1] 65.70818
# 
# names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"              
# [4] "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm"
# [10] "sample"              "cpd_rel_abun_norm"  
# 
# length(unique(dat$cpd_id))
# [1] 7213
# 
# ### 3) collate_compounds within each sample
# 
# # # # # # # # # # #
# 
# 
# ## assemble results
# cpd_id    sample cpd_rel_abun
# 1 cpd00851 SRR341581   0.00000000
# 2 cpd02175 SRR341581   0.00000000
# 3 cpd00035 SRR341581   0.12458586
# 4 cpd00117 SRR341581   0.08746449
# 5 cpd00051 SRR341581   0.05825908
# 6 cpd00586 SRR341581   0.00000000
# 
# 
# 
# str(dat.cpd.collate)
# 'data.frame':	591466 obs. of  3 variables:
# $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.1246 0.0875 0.0583 ...
# 
# sum(dat.cpd.collate$cpd_rel_abun)
# [1] 5388.071
# 
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# [1] 65.70818


#-------------------------

#### Forslund CHN-T2D - w/ Host-removal - continue CPP analysis
#-------------------------

# copy output file from HPC
dat.cpd.collate <- readRDS("/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn/run_cpp_r_on_hpc_gcb_r1_chn/dat.cpd.collate-all-samps-cpp3d--forslund-t2d-chn-hostremoval-pawsey.rds")

str(dat.cpd.collate)
# 'data.frame':	591466 obs. of  3 variables:
# $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.1246 0.0875 0.0583 ...

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000000  0.000000  0.000077  0.009110  0.000978 11.770873 

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf    -Inf  -4.111    -Inf  -3.010   1.071 


# log10 abun
dat.cpd.collate$log10_abun <- dat.cpd.collate$cpd_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(dat.cpd.collate$log10_abun == 0) # qty 155997
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(dat.cpd.collate$log10_abun[ -subsel.zero ])
  dat.cpd.collate$log10_abun[ subsel.zero ] <- zero_replace
}
dat.cpd.collate$log10_abun <- log10(dat.cpd.collate$log10_abun)

hist(dat.cpd.collate$log10_abun); summary( dat.cpd.collate$log10_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -8.446  -8.446  -4.111  -4.940  -3.010   1.071 


# make group variable from sample name

dat.cpd.collate$group <- NA

# from above
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 19352 taxa and 82 samples ]
# sample_data() Sample Data:       [ 82 samples by 31 sample variables ]
# tax_table()   Taxonomy Table:    [ 19352 taxa by 4 taxonomic ranks ]

identical( phy@sam_data$Run , samp$Run ) # TRUE
identical( sample_names(phy), samp$Run ) # TRUE
unique(samp$Diagnosis)
# "T2D metformin-" "ND CTRL"   

samp$group_new <- factor(samp$Diagnosis, 
                         levels = c("T2D metformin-", "ND CTRL"),
                         labels = c("T2D met-", "Normal"),
                         ordered = TRUE )

#for (i in 1:length(sample_names(phy))) {
for (i in 1:length( samp$Run )) {
  #i<-1
  #this_samp <- sample_names(phy)[i]
  this_samp <- samp$Run[i]
  sel <- which(dat.cpd.collate$sample == this_samp)
  #dat.cpd.collate$group[sel] <- phy@sam_data$age[i]
  dat.cpd.collate$group[sel] <- as.character( samp$group_new[i] )
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group) # "T2D met-" "Normal"  
dat.cpd.collate$group <- factor(dat.cpd.collate$group, levels = c("T2D met-", "Normal"), ordered = TRUE)

head(dat.cpd.collate)

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval.rds")


str(dat.cpd.collate)
# 'data.frame':	591466 obs. of  5 variables:
# $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.1246 0.0875 0.0583 ...
# $ log10_abun  : num  -8.446 -8.446 -0.905 -1.058 -1.235 ...
# $ group       : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...


length( unique(dat.cpd.collate$cpd_id) ) # 7213
7213*82 # 591466


## CPP stats ?

data_in <- dat.cpd.collate

head(data_in)
# cpd_id    sample cpd_rel_abun log10_abun    group
# 1 cpd00851 SRR341581   0.00000000 -8.4459524 T2D met-
# 2 cpd02175 SRR341581   0.00000000 -8.4459524 T2D met-
# 3 cpd00035 SRR341581   0.12458586 -0.9045313 T2D met-
# 4 cpd00117 SRR341581   0.08746449 -1.0581682 T2D met-
# 5 cpd00051 SRR341581   0.05825908 -1.2346364 T2D met-
# 6 cpd00586 SRR341581   0.00000000 -8.4459524 T2D met-

dim(data_in) # 591466      5

unique_samps <- unique(data_in$sample)

no_compounds <- numeric(length = length(unique_samps))
sample_sum_relabun <- numeric(length = length(unique_samps))

for (i in 1:length(unique_samps)) {
  #i<-1
  this_samp <- unique_samps[i]
  sel <- which(data_in$sample == this_samp)
  
  values <- data_in$cpd_rel_abun[sel]
  values <- values[values > 0]
  
  no_compounds[i] <- length( values )
  sample_sum_relabun[i] <- sum(values)
  print(paste0("completed ",i))
}

mean(no_compounds) # 5310.598
sd(no_compounds) #  1670.226

mean(sample_sum_relabun) # 65.70818
sd(sample_sum_relabun) # 4.138492

length(unique(data_in$cpd_id)) # 7213

#-------------------------

#### Forslund-CHN-T2D - w/ Host-removal - Test for selected compounds / CPP in T2D?
#-------------------------

## compound-level rel abundances in T2D

dat.cpd.t2d <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval.rds")
str(dat.cpd.t2d)
# 'data.frame':	591466 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.1246 0.0875 0.0583 ...
# $ log10_abun  : num  -8.446 -8.446 -0.905 -1.058 -1.235 ...
# $ group       : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

dat.cpd.t2d$group_label <- factor(dat.cpd.t2d$group, levels = c("T2D met-", "Normal"),
                                  ordered = TRUE)

length( unique(dat.cpd.t2d$cpd_id) ) # 7213
length( unique(dat.cpd.t2d$sample) ) # 82


# compound info is in here
names(df.comp)
# [1] "id"        "abbrev"    "name"      "form"      "OC_ratio"  "HC_ratio"  "NC_ratio"  "PC_ratio"  "NP_ratio"  "O_count"   "N_count"   "P_count"   "S_count"   "mass"     
# [15] "SC_ratio"  "MgC_ratio" "ZnC_ratio" "KC_ratio"  "CaC_ratio" "MnC_ratio" "FeC_ratio" "CoC_ratio" "CuC_ratio" "MoC_ratio"
head(df.comp$abbrev) # "h2o"   "atp"   "nad"   "nadh"  "nadph" "nadp"
head(df.comp$name) # "H2O"   "ATP"   "NAD"   "NADH"  "NADPH" "NADP" 
head(df.comp$form) # "H2O"           "C10H13N5O13P3" "C21H26N7O14P2" "C21H27N7O14P2" "C21H26N7O17P3" "C21H25N7O17P3"


## Fructose
sel.cpd <- which(df.comp$name == "D-Fructose") # 
df.comp[sel.cpd, ]
#          id abbrev       name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 82 cpd00082    fru D-Fructose C6H12O6
this_var <- "D-Fructose"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00082") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] #
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] #


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## Sucrose
sel.cpd <- which(df.comp$name == "Sucrose") # 
df.comp[sel.cpd, ]
#          id abbrev    name      form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 76 cpd00076   sucr Sucrose C12H22O11
this_var <- "Sucrose"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00076") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun % # "greater"
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "Arabinose"

sel.cpd <- which(df.comp$name == "L-Arabinose") # 
df.comp[sel.cpd, ]
#           id abbrev        name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio
# 223 cpd00224 arab-L L-Arabinose C5H10O5
this_var <- "L-Arabinose"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00224") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n =


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "ammonia"
#sel.cpd <- which(df.comp$name == "Ammonia") # empty
#sel.cpd <- grep(pattern = "NH3|*mmonia", x = df.comp$name) # ok
sel.cpd <- which(df.comp$name == "NH3") # ok
df.comp[sel.cpd, ]
#          id abbrev name
# 13 cpd00013    nh4  NH3
this_var <- "Ammonia"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00013") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## glucose
sel.cpd <- which(df.comp$name == "D-Glucose") # 
df.comp[sel.cpd, ]
#             id    abbrev      name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 27    cpd00027     glc-D D-Glucose C6H12O6        1        2        0        0      NaN       6       0       0       0  180        0         0         0        0         0
# 24094 cpd26821 D-Glucose D-Glucose C6H12O6
this_var <- "D-Glucose"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00027") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n =

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "Galactose"
sel.cpd <- which(df.comp$name == "Galactose") # 
df.comp[sel.cpd, ]
#            id  abbrev      name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 108  cpd00108     gal Galactose C6H12O6        1        2        0        0      NaN       6       0       0       0  180        0         0         0        0         0
# 1090 cpd01112 cbs_337 Galactose C6H12O6
this_var <- "Galactose"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00108") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n =


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "xylose"
sel.cpd <- which(df.comp$name == "Xylose") # 
df.comp[sel.cpd, ]
#             id             abbrev               name        form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count
# 153   cpd00154              xyl-D             Xylose     C5H10O5 1.0000000 2.000000        0        0      NaN       5       0
this_var <- "Xylose"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00154") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "mannose"                                      
sel.cpd <- which(df.comp$name == "D-Mannose") # 
df.comp[sel.cpd, ]
#           id abbrev      name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 138 cpd00138    man D-Mannose C6H12O6
this_var <- "D-Mannose"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00138") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n =

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "glucosamine"
# D-Glucosamine is made naturally in the form of glucosamine-6-phosphate
#sel.cpd <- which(df.comp$name == "Glucosamine") # ok - but not in samples
#sel.cpd <- grep(pattern = "D-glucosamine-6-phosphate", x = df.comp$name) # D-Glucosamine is made naturally in the form of glucosamine-6-phosphate
sel.cpd <- grep(pattern = "*D-glucosamine", x = df.comp$name)
sel.cpd <- which(df.comp$name == "alpha-D-glucosamine 6-phosphate") 
# "alpha-D-glucosamine" "alpha-D-glucosamine 6-phosphate"
#sel.cpd <- grep(pattern = "*mino-2-deoxy-glucose", x = df.comp$name)
#sel.cpd <- which(df.comp$name == "2-amino-2-deoxy-glucose") # ok - but not in samples
#sel.cpd <- which(df.comp$form == "C6H13NO5") # ok - but not in samples

unique(df.comp$name[sel.cpd])
df.comp[sel.cpd, ]
#             id                          abbrev                            name      form 
# 21178 cpd23898 alpha-D-glucosamine 6-phosphate alpha-D-glucosamine 6-phosphate C6H13NO8P

this_var <- "D-glucosamine 6-phosphate"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd23898") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n =

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "chitin"
sel.cpd <- which(df.comp$name == "Chitin") #
df.comp[sel.cpd, ]
#             id abbrev   name       form 
# 11508 cpd11683 chitin Chitin C8H13NO5R2
this_var <- "Chitin"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd11683") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n =


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "xylan"
sel.cpd <- which(df.comp$name == "Xylan") # ok
#sel.cpd <- which(df.comp$form == "") # ok
sel.cpd <- grep(pattern = "xylan", x = df.comp$name) # multiple
#sel.cpd <- grep(pattern = "*", x = df.comp$form) # 
df.comp[sel.cpd, ]
# id                                    abbrev                                      name        form  OC_ratio  HC_ratio
# 10924 cpd11070 3beta-Hydroxylanostane-7,11-dione acetate 3beta-Hydroxylanostane-7,11-dione acetate    C32H52O4 0.1250000 1.6250000
# 11789 cpd11970                              Arabinoxylan                              Arabinoxylan  C10H16O8R2 0.8000000 1.6000000
# 12068 cpd12254                     Glucuronoarabinoxylan                     Glucuronoarabinoxylan        null       NaN       NaN
# 12217 cpd12408                            Glucuronoxylan                            Glucuronoxylan C21H31O18R2 0.8571429 1.4761905
# 12354 cpd12550                  4-O-methylglucuronoxylan                  4-O-methylglucuronoxylan C22H33O18R2 0.8181818 1.5000000
# 17230 cpd19945                     11-Deoxylandomycinone                     11-Deoxylandomycinone    C19H14O5 0.2631579 0.7368421
# 19578 cpd22295                               Acetylxylan                               Acetylxylan C27H39O21R3 0.7777778 1.4444444
# 19670 cpd22387                              Arabinoxylan                              Arabinoxylan C35H56O29R2 0.8285714 1.6000000
# 24444 cpd27172    Glucuronoarabinoxylan-Oligosaccharides    Glucuronoarabinoxylan-Oligosaccharides C27H41O23R2 0.8518519 1.5185185
# 24445 cpd27173                    Glucuronoarabinoxylans                    Glucuronoarabinoxylans C52H81O43R2 0.8269231 1.5576923
# 24447 cpd27175           Glucuronoxylan-Oligosaccharides           Glucuronoxylan-Oligosaccharides C22H33O19R2 0.8636364 1.5000000
# 24448 cpd27176                           Glucuronoxylans                           Glucuronoxylans C37H57O31R2 0.8378378 1.5405405
# 25918 cpd28650                      (1, 4-beta-D-xylan)n                      (1, 4-beta-D-xylan)n 

sel.cpd <- which(df.comp$name == "Arabinoxylan") # ok
df.comp[sel.cpd, ]
#             id       abbrev         name        form
# 11789 cpd11970 Arabinoxylan Arabinoxylan  C10H16O8R2 - NOT PRESENT IN DATA
# 19670 cpd22387 Arabinoxylan Arabinoxylan C35H56O29R2

this_var <- "Arabinoxylan"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd22387") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n =

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "hydrogen sulfide"                             
sel.cpd <- which(df.comp$name == "H2S") # ok
df.comp[sel.cpd, ]
#           id abbrev name
# 237 cpd00239    h2s  H2S

this_var <- "H2S"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00239") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n =


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## Starch
sel.cpd <- which(df.comp$name == "Starch") # 
df.comp[sel.cpd, ]
#             id abbrev   name        form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count    mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 11482 cpd11657 starch Starch C12H20O10R2 0.8333333 1.666667        0        0      NaN      10       0       0       0 9.9e+02        0         0         0        0         0
# 25462 cpd28193 Starch Starch        null
this_var <- "Starch"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd11657") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n =


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "amylopectin"                                  
sel.cpd <- which(df.comp$name == "Amylopectin") # ok
df.comp[sel.cpd, ]
#           id      abbrev        name      form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count
# 262 cpd00265 Amylopectin Amylopectin C30H52O26
this_var <- "Amylopectin"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00265") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n =


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "amylose"
sel.cpd <- which(df.comp$name == "Amylose") # ok
df.comp[sel.cpd, ]
#             id abbrev    name      form
# 11560 cpd11735 14glun Amylose C6H10O5R2
this_var <- "Amylose"


# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd11735") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n =


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## Glycogen (animal starch) - Glycogen is a multibranched polysaccharide of glucose that serves as a form of energy storage in animals,[2] fungi, and bacteria
sel.cpd <- which(df.comp$name == "Glycogen") # 
df.comp[sel.cpd, ]
#           id   abbrev     name      form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 154 cpd00155 glycogen Glycogen C24H42O21
this_var <- "Glycogen"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00155") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "acetate"
sel.cpd <- which(df.comp$name == "Acetate") # ok
df.comp[sel.cpd, ]
#          id abbrev    name   form
# 29 cpd00029     ac Acetate C2H3O2
this_var <- "Acetate"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00029") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "propionate"
sel.cpd <- which(df.comp$name == "Propionate") # ok
df.comp[sel.cpd, ]
#           id abbrev       name   form
# 140 cpd00141    ppa Propionate C3H5O2
this_var <- "Propionate"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00141") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "butyrate"
sel.cpd <- which(df.comp$name == "Butyrate") # ok
df.comp[sel.cpd, ]
#           id abbrev     name   form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio MnC_ratio
# 210 cpd00211    but Butyrate C4H7O2
this_var <- "Butyrate"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00211") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "methane"
sel.cpd <- which(df.comp$name == "Methane") # ok
df.comp[sel.cpd, ]
#            id abbrev    name form
# 1004 cpd01024  metha Methane  CH4

this_var <- "CH4"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd01024") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
#x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 52



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = log10_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10 CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "Carbon dioxide"                               
sel.cpd <- which(df.comp$name == "CO2")
this_var <- "CO2"

df.comp[sel.cpd, ]
#          id abbrev name
# 11 cpd00011    co2  CO2


# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00011")
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## lignin
sel.cpd <- which(df.comp$name == "Lignin") # 
df.comp[sel.cpd, ]
#             id abbrev   name form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio MnC_ratio
# 12548 cpd12745 Lignin Lignin null
this_var <- "Lignin"


# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd12745") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n =
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n =



# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")





## "trimethylamine N-oxide"
#sel.cpd <- which(df.comp$name == "TMAO") # empty
sel.cpd <- which(df.comp$form == "C3H9NO") # ok
df.comp[sel.cpd, ]
#           id abbrev     name   form
# 797 cpd00811   tmao (CH3)3NO C3H9NO

this_var <- "TMAO"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00811") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
#x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n =
x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 52


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = log10_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10 CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "trimethylamine"                               
sel.cpd <- which(df.comp$abbrev == "tma") # ok
df.comp[sel.cpd, ]
#           id abbrev    name   form 
# 437 cpd00441    tma (CH3)3N C3H10N

this_var <- "TMA"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00441") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
#x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n =
x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 52


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = log10_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10 CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "pectin"
sel.cpd <- which(df.comp$name == "Pectin" & df.comp$form == "C26H34O24R2") # ok
df.comp[sel.cpd, ]
#             id abbrev   name        form
# 11434 cpd11601 pectin Pectin C26H34O24R2
this_var <- "Pectin"


# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd11601") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
#x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
#y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n = 52


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = log10_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10 CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "p-cresyl sulfate" - instead use precursor: p-Cresol
sel.cpd <- which(df.comp$abbrev == "p-Cresol") # ok
df.comp[sel.cpd, ]
# id   abbrev     name  form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass
# 1022 cpd01042 p-Cresol p-Cresol C7H8O

this_var <- "p-Cresol"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd01042") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n =


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "menaquinone"
sel.cpd <- which(df.comp$abbrev == "Menaquinones") # ok
df.comp[sel.cpd, ]
#             id       abbrev         name      form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count
# 24771 cpd27501 Menaquinones Menaquinones C16H15O2R

this_var <- "Menaquinones"


# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd27501") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n =


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "serotonin" 
sel.cpd <- which(df.comp$name == "Serotonin") # ok
df.comp[sel.cpd, ]
#           id    abbrev      name      form
# 568 cpd00579 Serotonin Serotonin C10H13N2O
this_var <- "Serotonin"


# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00579") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "cellulose"                                    
sel.cpd <- which(df.comp$name == "Cellulose") # ok
df.comp[sel.cpd, ]
#             id    abbrev      name      form
# 11571 cpd11746 Cellulose Cellulose C6H10O5R2

this_var <- "Cellulose"


# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd11746") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
#x <- dat.cpd.t2d$log10_abun[sel[subsel.t2d]] # n = 
#y <- dat.cpd.t2d$log10_abun[sel[subsel.norm]] # n =


# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## "indole" (KEGG C00463: 2,3-Benzopyrrole)
#sel.cpd <- which(df.comp$name == "indole") # empty
sel.cpd <- which(df.comp$abbrev == "indole" & df.comp$form == "C8H7N") # ok
df.comp[sel.cpd, ]
#           id abbrev  name  form 
# 356 cpd00359 indole indol C8H7N

this_var <- "Indole"


# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00359") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("T2D vs NORMAL\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)


p <- ggplot(data = dat.cpd.t2d[sel, ] , aes(x=group_label, y = cpd_rel_abun))+ # log10_abun. cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
#grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-T2D-v6.tiff"), width = 7, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


#-------------------------

#### Forslund-CHN-T2D - w/ Host-removal - test compound groupings - ACPs & sugars?
#-------------------------

dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval.rds")

str(dat.cpd.collate)
# 'data.frame':	591466 obs. of  5 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 82
length(unique(df$cpd_id)) # 7213
82*7213 # 591466

sel <- which(df$cpd_id %in% t2d.zoom1.decreasing.proteins)
df <- df[sel, ]
length(unique(df$cpd_id)) # 35
35*82 # 2870

str(df)
# 'data.frame':	2870 obs. of  5 variables:

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
# 'data.frame':	82 obs. of  4 variables:

# unique(res$group) # "T2D met neg" "Normal"
# res$group <- factor(res$group, levels = c("T2D met neg", "Normal"), ordered = TRUE)
unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)


# Wilcoxon-Mann-Whitney Test

# x <- res$sum_rel_abun[ res$group == "T2D met neg" ] # 
# y <- res$sum_rel_abun[ res$group == "Normal" ] #
x <- res$sum_rel_abun[ res$group_label == "T2D met-" ] # n = 30
y <- res$sum_rel_abun[ res$group_label == "Normal" ] # n = 52

wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Wilcoxon-Mann-Whitney ~T2D\nW = ",round(wmw.test$statistic,3),"; ",pval)

str(res)

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  #geom_point(alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("\U03A3 CPP rel abun BCFA-ACPs (%)")+ # ylab("CPP of indicator group 1 - \U03A3 rel abun (%)")+
  ylim(-0.00005,0.014)+
  
  #geom_smooth(method = "lm")+
  theme_bw()+
  #scale_x_continuous(labels= c("Normal","IGT", "T2D met+", "T2D met-"))+
  
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "right", npcy = "top", label = "T2D case study", size = 3.5)+
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = rel(0.9))
  )

p

#grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
#grid.text(label = "(a)", x = unit(0.12, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
grid.text(label = "(b)", x = unit(0.12, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Zoom1-Indicator-group1-BCFA-ACPs-Trend-with-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")





## Increasing carbs

dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval.rds")

str(dat.cpd.collate)
# 'data.frame':	591466 obs. of  5 variables:

df <- dat.cpd.collate
length(unique(df$sample)) # 82
length(unique(df$cpd_id)) # 7213
82*7213 # 591466

t2d.zoom2B.increasing.carbs # "cpd00224" "cpd03198" "cpd03200" "cpd00259" "cpd00076" "cpd00082"

sel <- which(df$cpd_id %in% t2d.zoom2B.increasing.carbs)
df <- df[sel, ]
length(unique(df$cpd_id)) # 6
6*82 # 492

str(df)
# 'data.frame':	492 obs. of  5 variables:

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
# 'data.frame':	82 obs. of  4 variables:

# unique(res$group) # "T2D met neg" "Normal"
# res$group <- factor(res$group, levels = c("T2D met neg", "Normal"), ordered = TRUE)
unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)


# Wilcoxon-Mann-Whitney Test

# x <- res$sum_rel_abun[ res$group == "T2D met neg" ] # 
# y <- res$sum_rel_abun[ res$group == "Normal" ] #
x <- res$sum_rel_abun[ res$group_label == "T2D met-" ] # n = 30
y <- res$sum_rel_abun[ res$group_label == "Normal" ] # n = 52


wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun % #  "less"
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Wilcoxon-Mann-Whitney ~T2D\nW = ",round(wmw.test$statistic,3),"; ",pval)

str(res)

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  #geom_point(alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("\U03A3 CPP rel abun Sugars (%)")+ # ylab("CPP of indicator group 2 - \U03A3 rel abun (%)")+
  #geom_smooth(method = "lm")+
  theme_bw()+
  #scale_x_continuous(labels= c("Normal","IGT", "T2D met+", "T2D met-"))+
  
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3)+
  #annotate(geom="text_npc", npcx = "right", npcy = "top", label = "T2D case study", size = 3.5)+
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = rel(0.9))
  )

p

#grid.text(label = "(b)", x = unit(0.1, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CHN-Zoom2-b-Indicator-group2-SUGARS-Trend-with-T2D-v6.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


#-------------------------

#### Forslund-CHN-T2D - w/ Host-removal - CPP as phyloseq object
#    PCoA using CPP vs Functions
#-------------------------

dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval.rds")
data_in <- dat.cpd.collate
str(data_in)
# 'data.frame':	591466 obs. of  5 variables:

length( unique(data_in$cpd_id) ) # 7213
length( unique(data_in$cpd_id[data_in$cpd_rel_abun > 0]) ) # 7213
length( unique(data_in$sample) ) # 82
7213*82 # 591466


### get data into phyloseq object ...

head(data_in)
# cpd_id    sample cpd_rel_abun log10_abun    group
# 1 cpd00851 SRR341581   0.00000000 -8.4459524 T2D met-
# 2 cpd02175 SRR341581   0.00000000 -8.4459524 T2D met-
# 3 cpd00035 SRR341581   0.12458586 -0.9045313 T2D met-
# 4 cpd00117 SRR341581   0.08746449 -1.0581682 T2D met-
# 5 cpd00051 SRR341581   0.05825908 -1.2346364 T2D met-
# 6 cpd00586 SRR341581   0.00000000 -8.4459524 T2D met-

data_in$group_label <- data_in$group

df.wide <- dcast(data_in, formula = sample + group_label ~ cpd_id , value.var = "cpd_rel_abun" )

df.wide[1:5, 1:10]
#      sample group_label  cpd00001 cpd00002  cpd00003  cpd00004  cpd00005  cpd00006   cpd00007  cpd00008
# 1 SRR341581    T2D met- 11.372788 1.357728 0.1731380 0.1704617 0.4159819 0.4159819 0.04165925 0.3838933
# 2 SRR341585    T2D met-  8.329547 2.390563 0.3058945 0.2863093 0.4428879 0.4428879 0.04030459 0.4304293
# 3 SRR341586    T2D met-  7.761999 2.405028 0.3363914 0.3363914 0.1032110 0.1032110 0.10193680 1.5889146
# 4 SRR341587    T2D met- 10.048373 1.072346 0.2858257 0.2708795 0.3111261 0.3111261 0.05880801 0.3543000
# 5 SRR341588    T2D met- 10.454767 1.424459 0.2693916 0.2547486 0.4114752 0.4114752 0.13160432 0.5565338

unique(paste0(df.wide$sample,"--",df.wide$group_label))
# [1] "SRR341581--T2D met-" "SRR341585--T2D met-" "SRR341586--T2D met-" "SRR341587--T2D met-" "SRR341588--T2D met-" "SRR341589--T2D met-" "SRR341599--T2D met-"
# [8] "SRR341600--T2D met-" "SRR341601--T2D met-" "SRR341602--T2D met-" "SRR341604--T2D met-" "SRR341606--T2D met-" "SRR341609--T2D met-" "SRR341636--Normal"  
# [15] "SRR341645--Normal"   "SRR341646--Normal"   "SRR341652--Normal"   "SRR341654--T2D met-" "SRR341655--T2D met-" "SRR341657--T2D met-" "SRR341660--T2D met-"
# [22] "SRR341661--T2D met-" "SRR341663--T2D met-" "SRR341664--T2D met-" "SRR341665--T2D met-" "SRR341669--T2D met-" "SRR341670--T2D met-" "SRR341673--T2D met-"
# [29] "SRR341674--T2D met-" "SRR341675--T2D met-" "SRR341676--T2D met-" "SRR341681--T2D met-" "SRR341684--T2D met-" "SRR341687--T2D met-" "SRR341693--Normal"  
# [36] "SRR341713--Normal"   "SRR413575--Normal"   "SRR413576--Normal"   "SRR413578--Normal"   "SRR413579--Normal"   "SRR413580--Normal"   "SRR413581--Normal"  
# [43] "SRR413582--Normal"   "SRR413584--Normal"   "SRR413585--Normal"   "SRR413587--Normal"   "SRR413592--Normal"   "SRR413593--Normal"   "SRR413594--Normal"  
# [50] "SRR413597--Normal"   "SRR413598--Normal"   "SRR413599--Normal"   "SRR413600--Normal"   "SRR413601--Normal"   "SRR413603--Normal"   "SRR413605--Normal"  
# [57] "SRR413606--Normal"   "SRR413607--Normal"   "SRR413608--Normal"   "SRR413610--Normal"   "SRR413613--Normal"   "SRR413614--Normal"   "SRR413615--Normal"  
# [64] "SRR413616--Normal"   "SRR413617--Normal"   "SRR413618--Normal"   "SRR413619--Normal"   "SRR413620--Normal"   "SRR413621--Normal"   "SRR413623--Normal"  
# [71] "SRR413625--Normal"   "SRR413626--Normal"   "SRR413634--Normal"   "SRR413637--Normal"   "SRR413642--Normal"   "SRR413652--Normal"   "SRR413660--Normal"  
# [78] "SRR413661--Normal"   "SRR413670--Normal"   "SRR413688--Normal"   "SRR413758--Normal"   "SRR413768--Normal"

# save group variable
samp <- df.wide[ ,1:2]
row.names(samp) <- samp$sample

# transpose
df.wide <- t(df.wide[ ,-2]) # minus 'group' column

head(df.wide)

samp_names <- df.wide[1, ]
tax_names <- row.names(df.wide[-1, ])
head(tax_names) # "cpd00001" "cpd00002" "cpd00003" "cpd00004" "cpd00005" "cpd00006"
otu.df <- df.wide[-1, ] # remove sample labels in 1st row
# this is necessary to create numeric matrix

colnames(otu.df) <- samp_names

# convert OTU table to matrix
class(otu.df) # "matrix" "array"
#otu.df <- as.matrix(otu.df)

# convert to numeric matrix
# https://stackoverflow.com/questions/20791877/convert-character-matrix-into-numeric-matrix
otu.df <- apply(otu.df, 2, as.numeric)

rownames(otu.df) # NULL
dim(otu.df) #  7213   82
rownames(otu.df) <- tax_names

## Create 'otuTable'
#  otu_table - Works on any numeric matrix.
#  You must also specify if the species are rows or columns
OTU <- otu_table(otu.df, taxa_are_rows = TRUE)


# # convert Taxonomy table to matrix

tax <- data.frame(cpd_id = tax_names)
row.names(tax) <- tax_names

tax <- as.matrix(tax)

identical( row.names(otu.df), row.names(tax) ) # TRUE


## Create 'taxonomyTable'
#  tax_table - Works on any character matrix.
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
TAX <- tax_table(tax)


## Create a phyloseq object, merging OTU & TAX tables
phy.cpp = phyloseq(OTU, TAX)
phy.cpp
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7213 taxa and 82 samples ]
# tax_table()   Taxonomy Table:    [ 7213 taxa by 1 taxonomic ranks ]


sample_names(phy.cpp)
#  [1] "SRR341581" "SRR341585" "SRR341586" "SRR341587" "SRR341588" ... etc.

#identical(sample_names(phy.cpp), samp$sample) # TRUE

identical(sample_names(phy.cpp), sradat.select2$Run) # TRUE

row.names(sradat.select2) <- sradat.select2$Run
samp <- sradat.select2

# row.names need to match sample_names() from phyloseq object
#row.names(samp) <- samp$sample
#identical(row.names(samp), samp$sample) # TRUE



### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

SAMP <- sample_data(samp)


### Combine SAMPDATA into phyloseq object
phy.cpp <- merge_phyloseq(phy.cpp, SAMP)
phy.cpp
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7213 taxa and 82 samples ]
# sample_data() Sample Data:       [ 82 samples by 31 sample variables ]
# tax_table()   Taxonomy Table:    [ 7213 taxa by 1 taxonomic ranks ]

phy.cpp@sam_data

min(taxa_sums(phy.cpp)) #  3.667808e-08

saveRDS(object = phy.cpp, file = "phy.cpp-cleaned-Forslund-CHN-T2D-Hostremoval-v6.RDS")


phy_in <- phy.cpp

sum(sample_sums(phy_in)) # 5388.071
sample_sums(phy_in)
# SRR341581 SRR341585 SRR341586 SRR341587 SRR341588 SRR341589 SRR341599 SRR341600 SRR341601 SRR341602 SRR341604 SRR341606 SRR341609 SRR341636 SRR341645 SRR341646 SRR341652
# 68.96271  77.40765  65.61562  62.06336  67.40098  76.13943  71.18610  68.96971  73.23230  67.43343  69.95068  71.88263  67.84404  74.73696  68.02291  61.31068  73.47749
# SRR341654 SRR341655 SRR341657 SRR341660 SRR341661 SRR341663 SRR341664 SRR341665 SRR341669 SRR341670 SRR341673 SRR341674 SRR341675 SRR341676 SRR341681 SRR341684 SRR341687
# 58.11243  67.58085  68.85250  69.02749  64.71532  61.12007  62.63412  65.97605  61.27748  62.74963  70.46170  61.74862  62.39391  63.33472  59.83734  61.27360  64.34223
# SRR341693 SRR341713 SRR413575 SRR413576 SRR413578 SRR413579 SRR413580 SRR413581 SRR413582 SRR413584 SRR413585 SRR413587 SRR413592 SRR413593 SRR413594 SRR413597 SRR413598
# 61.70671  67.25532  65.48353  68.70297  60.33908  64.00873  65.55010  68.52288  66.74919  63.26478  60.29283  62.98827  61.04347  66.95152  63.23313  63.00998  68.96300
# SRR413599 SRR413600 SRR413601 SRR413603 SRR413605 SRR413606 SRR413607 SRR413608 SRR413610 SRR413613 SRR413614 SRR413615 SRR413616 SRR413617 SRR413618 SRR413619 SRR413620
# 62.19227  61.93779  61.11408  61.32069  62.37917  65.33822  63.60747  61.10597  69.97433  66.35465  64.08552  66.59683  65.92939  65.90732  61.28665  66.18646  63.05098
# SRR413621 SRR413623 SRR413625 SRR413626 SRR413634 SRR413637 SRR413642 SRR413652 SRR413660 SRR413661 SRR413670 SRR413688 SRR413758 SRR413768
# 63.30039  65.02099  67.77512  62.54607  61.94472  66.19041  72.03254  70.68306  68.77104  70.36588  59.88877  70.38936  62.26137  70.42075

summary( sample_sums(phy_in) )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 58.09   62.29   65.52   65.71   68.75   77.16 

sd( sample_sums(phy_in) )
# 4.138492

max(taxa_sums(phy_in)) # 541.9709


phy_in@sam_data$group_label <- factor(phy_in@sam_data$Diagnosis,
                                      levels = c("T2D metformin-", "ND CTRL"),
                                      labels = c("T2D met-", "Normal"),
                                      ordered = TRUE
                                      )

phy_in@sam_data$sex <- factor(phy_in@sam_data$gender,
                                      levels = c("female", "male"),
                                      labels = c("Female", "Male"),
                                      ordered = TRUE
)

# don't rarefy - already in form of relative abundance %


table(phy_in@sam_data$group_label)
# T2D met-   Normal
# 30       52





## ordination plot
## PCoA + Bray-Curtis

set.seed(123)
ord <- ordinate(phy_in, "PCoA", "bray")


ord

unique(phy_in@sam_data$group_label)
# [1] T2D met- Normal
# Levels: T2D met- < Normal


p <- plot_ordination(phy_in, ord, type="samples", color="group_label")
p

p$labels$x # "Axis.1   [66.1%]"
x_lab <- "PCo1 (66.1%)"

p$labels$y # "Axis.2   [10.5%]"
y_lab <- "PCo2 (10.5%)"

66.1 + 10.5 # 76.6


#temp <- r1.ps
p_df <- p$data


cols.group <- c("Normal" = "#7fbf7b",
                "T2D met-" = "#af8dc3")

shape.sex <- c("Female" = 1,
                "Male" = 2)


p <- #plot_ordination(temp, ord, type="samples", color="group") +
  #ggplot(data = p_df, aes(x = NMDS1, y = NMDS2, color = group))+
  ggplot(data = p_df, aes(x = Axis.1, y = Axis.2, color = group_label, shape = sex))+
  theme_bw()+
  geom_point()+

  stat_ellipse(linetype = "dashed")+

  xlab(x_lab) + ylab(y_lab)+
  scale_color_manual(values = cols.group, name = "Diagnosis") +
  scale_shape_manual(values = shape.sex, name = "Sex") +
  #annotate(geom="text", x= -1.2, y= 1.3, label = paste0("Stress = ",round(ord$stress,digits=4)),size = 3, hjust=0, vjust=1) +
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = "T2D: CPP", size = 3 )+
  theme(
    panel.grid.major = element_blank(),

    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.85)),
    legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = -8, unit = "pt"), # randomForest::margin() interferes !!

    panel.grid.minor = element_blank())
p

#grid.text(label = "(d)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )

dev.print(tiff, file = paste0(workdir,"/plots/","PCoA-Cpp3d-Forslund-CHN-T2D-Hostremoval-v6.tiff"), width = 9.5, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## PERMANOVA

# Calculate bray curtis distance matrix
set.seed(123)
bray <- phyloseq::distance(phy_in, method = "bray")
sampledf <- data.frame(sample_data(phy_in))
str(sampledf)

names(phy_in@sam_data)
# [1] "Run"                      "actual_read_length (run)" "Age"                      "Assay Type"               "AvgSpotLen"               "Bases"
# [7] "BioProject"               "BioSample"                "Bytes"                    "center_name (exp)"        "Center Name"              "Consent"
# [13] "DATASTORE filetype"       "DATASTORE provider"       "DATASTORE region"         "Experiment"               "gender"                   "Instrument"
# [19] "Library Name"             "LibraryLayout"            "LibrarySelection"         "LibrarySource"            "NATION"                   "Organism"
# [25] "Platform"                 "ReleaseDate"              "run (run)"                "Sample Name"              "SRA Study"                "total_bases (run)"
# [31] "Diagnosis"                "group_label"              "sex"

# Adonis test
set.seed(123)
adonis2(bray ~ group_label + sex , data = sampledf)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ group_label + sex, data = sampledf)
#             Df SumOfSqs      R2      F Pr(>F)    
# group_label  1  0.30636 0.10964 10.216  0.001 ***
# sex          1  0.11876 0.04250  3.960  0.016 *  
# Residual    79  2.36914 0.84786                  
# Total       81  2.79426 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



beta <- betadisper(bray, sampledf$group_label)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
# Groups     1 0.19198 0.191981 15.794    999  0.001 ***
# Residuals 80 0.97240 0.012155                         
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#-------------------------

#### Wilcoxon-Mann-Whitney Tests - Forslund-CHN-T2D - w/ Host-removal - CPP differences for ALL compounds
#    Test sig diff Normal vs T2D - consider all p <= 0.05
#-------------------------

dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval.rds")

str(dat.cpd.collate)
# 'data.frame':	591466 obs. of  5 variables:
# $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.1246 0.0875 0.0583 ...
# $ log10_abun  : num  -8.446 -8.446 -0.905 -1.058 -1.235 ...
# $ group       : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

data_in <- dat.cpd.collate

unique(data_in$group)
# [1] T2D met- Normal  
# Levels: T2D met- < Normal

dat.test <- data.frame(cpd = unique(data_in$cpd_id), data_for_this_cpd=NA , 
                       median_t2d = NA, median_normal = NA,  
                       alt = NA,
                       p_val = NA, W_statistic = NA, hl_effect_wilcox = NA,
                       trend_with_disease = NA
)

for (i in 1:dim(dat.test)[1]) {
  #i<-1
  this_cpd <- dat.test$cpd[i]
  sel <- which(data_in$cpd_id == this_cpd)
  # prepare data
  df = data.frame(group = as.character( data_in$group[sel]), 
                  value = data_in$log10_abun[sel], 
                  value_perc = data_in$cpd_rel_abun[sel] )
  
  x <- df$value[df$group == "T2D met-"]
  y <- df$value[df$group == "Normal"]
  
  if ( length(which(x == min(x) )) > 0.5*length(x) & length(which(y == min(y) )) > 0.5*length(y) ) {
    # low data: if at least one dataset does not have at least 50% non-zero cases
    dat.test$data_for_this_cpd[i] <- "low data"
    
  } else {
    
    # # test for homogeneity of variances
    # var.test(x,y) # e.g., j<-100: F = 2.3477, num df = 32, denom df = 42, p-value = 0.009872
    
    dat.test$median_t2d[i] <- median( x ) # log10 abun  
    dat.test$median_normal[i] <- median( y ) # log10 abun
    
    alt <- NA
    if (median( x ) > median( y ) ) {
      alt <- "greater"
    } else {
      alt <- "less"
    }
    dat.test$alt[i] <- alt
    
    # Wilcoxon-Mann-Whitney Test
    wmw.test <- wilcox.test(x, y, alternative = alt, paired = FALSE, conf.int = TRUE) # based on log10 abun
    
    dat.test$p_val[i] <- wmw.test$p.value
    dat.test$W_statistic[i] <- wmw.test$statistic
    
    # https://search.r-project.org/CRAN/refmans/DescTools/html/HodgesLehmann.html
    # https://aakinshin.net/posts/r-hodges-lehmann-problems
    
    dat.test$hl_effect_wilcox[i] <- wmw.test$estimate
    #dat.test$hl_effect_hlfxn[i] <- hl(x, y)
    
    if (!(is.na(dat.test$p_val[i])|is.na(dat.test$W_statistic[i]))) {
      if (dat.test$p_val[i] <= 0.05 & alt == "greater") { dat.test$trend_with_disease[i] <- "Increasing" }
      else if (dat.test$p_val[i] <= 0.05 & alt == "less") { dat.test$trend_with_disease[i] <- "Decreasing" }
      else { dat.test$trend_with_disease[i] <- "No trend" }
    }
  }
  print(paste0("Completed ",i))
}


sel.low <- which(dat.test$data_for_this_cpd == "low data") # 978

length(unique(dat.cpd.collate$cpd_id)) # 7213
length(unique(dat.cpd.collate$cpd_id)) - length(sel.low) # 6235

sel.nonNA <- which(!is.na(dat.test$p_val)) # 6235 applicable tests

head(dat.test)
str(dat.test)

dat.test %>% arrange(., p_val) %>% slice(1:50)

write.csv(x = dat.test, file = "dat.test-cpp3d-indiv-Compounds-Forslund-CHN-T2D-over50s-Hostremoval--All-Compounds-T2D-VS-NORM--Wilcox-v6.csv")

saveRDS(dat.test, file = "dat.test-cpp3d-indiv-Compounds-Forslund-CHN-T2D-over50s-Hostremoval--All-Compounds-T2D-VS-NORM--Wilcox-v6.RDS")


# extract sig results ... (no p-adjustment)

sel.sig <- which(dat.test$p_val <= 0.05) # 4292

dat.test.sig <- dat.test[sel.sig, ]

dat.test.sig$minuslog10_p_val <- -log10(dat.test.sig$p_val)

plot(x = dat.test.sig$W_statistic , y =dat.test.sig$minuslog10_p_val , xlab="W statistic", ylab="-log10(P-value)")

dev.print(tiff, filename = paste0(workdir,"/plots/","3d-indiv-compound-VolcanoPlot-P-values--All-Compounds-CHN-T2D-VS-NORM--Wilcox-v4.tiff"),
          width = 12, height = 14, units = "cm", res=600, compression = "lzw",type="cairo" )


# join compound info??

dat.test.sig$cpd_names <- NA
dat.test.sig$cpd_forms <- NA

dat.test.sig$OC_x <- NA
dat.test.sig$HC_y <- NA
dat.test.sig$NC_z <- NA


for (i in 1:dim(dat.test.sig)[1]) {
  #i<-1
  this_cpd <- dat.test.sig$cpd[i]
  
  sel.cpd <- which(df.comp$id == this_cpd)
  
  dat.test.sig$cpd_names[i] <- df.comp$name[sel.cpd]
  dat.test.sig$cpd_forms[i] <- df.comp$form[sel.cpd]
  
  dat.test.sig$OC_x[i] <- df.comp$OC_ratio[sel.cpd]
  dat.test.sig$HC_y[i] <- df.comp$HC_ratio[sel.cpd]
  dat.test.sig$NC_z[i] <- df.comp$NC_ratio[sel.cpd]
 
  print(paste0("completed ",i))
}

write.table(x = dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-Forslund-CHN-T2D-Hostremoval-All-Compounds-T2D-VS-NORM--Wilcox-v6.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-Forslund-CHN-T2D-Hostremoval-All-Compounds-T2D-VS-NORM--Wilcox-v6.rds")


hist(dat.test.sig$NC_z); summary(dat.test.sig$NC_z)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.0000  0.0377  0.1150  0.2000  1.0000     328


dim(dat.test.sig) # 4292 15
sel.ok <- which(!is.na(dat.test.sig$NC_z) ) # qty 3964


dat.test.sig$z_layer <- NA

subsel <- which(dat.test.sig$NC_z[sel.ok] == 0) # 1896
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C = 0"

subsel <- which(dat.test.sig$NC_z[sel.ok] > 0 & dat.test.sig$NC_z[sel.ok] <= 0.2 ) # 1180
# Rivas-Ubach et al 2018 Table 1 highlights N/C breaks at 0.126, 0.2, 0.5, 0.7
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C >0 to 0.2"

subsel <- which(dat.test.sig$NC_z[sel.ok] > 0.2) # 888
max(dat.test.sig$NC_z[sel.ok]) # 1
dat.test.sig$z_layer[sel.ok [subsel]] <- "N:C >0.2 to 1"

unique(dat.test.sig$z_layer[sel.ok]) # "N:C >0 to 0.2" "N:C >0.2 to 1" "N:C = 0"   

dat.test.sig$z_layer <- factor(dat.test.sig$z_layer, levels = c("N:C = 0",
                                                                "N:C >0 to 0.2",
                                                                "N:C >0.2 to 1"), ordered = TRUE)

saveRDS(object =  dat.test.sig, file = "dat.test.sig-cpp3d-indiv-Compounds-Forslund-CHN-T2D-Hostremoval-All-Compounds-T2D-VS-NORM--Wilcox-v6.rds")


#dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-Forslund-CHN-T2D-Hostremoval-All-Compounds-T2D-VS-NORM--Wilcox-v6.rds")

dim(dat.test.sig) #  4292   16
head(dat.test.sig)
dim(dat.test.sig[ which(dat.test.sig$p_val <= 0.05), ]) # 4292   16
sel <- which(!is.na(dat.test.sig$OC_x) ) # 3964
sel <- which(dat.test.sig$trend_with_disease == "Decreasing") # 3984
sel <- which(dat.test.sig$trend_with_disease == "Decreasing" & !is.na(dat.test.sig$OC_x) ) # 3695
sel <- which(dat.test.sig$trend_with_disease == "Increasing") # 308
sel <- which(dat.test.sig$trend_with_disease == "Increasing" & !is.na(dat.test.sig$OC_x) ) # 269

head( sort(dat.test.sig$p_val) ) # 4.243639e-09 1.063086e-07 1.084533e-07 1.084533e-07 1.084533e-07 1.084533e-07

## plot as Increasing or Decreasing?? in vK space
## Use adjusted compound classes (adapted from Wu 2018, D'Andrilli, Rivas-Ubach 2018, and Minor et al 2015)

# # Zones and labels as above
# vkgrouprect.facets2 <- read.table(file = "cpp3d-compound-classes.tsv", header = TRUE, sep = "\t" )
# vkgrouprect.facets2.labels <- read.table(file = "cpp3d-compound-classes-labels.tsv", header = TRUE, sep = "\t" )




p <- ggplot(data = dat.test.sig[sel.ok, ]) +
  coord_equal()+
  #ggtitle("Compound processing potential of microbiota - Type 2 Diabetes case study")+
  #xlim(0,3.4)+ ylim(0,4.1)+
  xlim(0,2.6)+ ylim(0,3.1)+
  #geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 ) + #
  xlab("O:C ratio")+ ylab("H:C ratio")+
  
  guides(color = guide_legend(title = "Trend with disease in functional capacity\n(%) allocated to compounds"))+
  
  facet_wrap(facets = vars(z_layer))+
  
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, group = label), color="#969696", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, group = label), color="#969696", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  #geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0.2 to 2" ), aes(x = OC_x, y = HC_y, group = label), color="grey", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  geom_mark_rect(data= filter(vkgrouprect.facets2, z_layer == "N:C >0.2 to 1" ), aes(x = OC_x, y = HC_y, group = label), color="#969696", expand = unit(0, "mm"),radius = unit(0, "mm")  )+ # color="#737373",
  
  geom_point(aes(x = OC_x, y = HC_y, color = trend_with_disease), size = 1, alpha = 0.35 )+ # 
  
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C = 0" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T,  size = 2.75 , col="#252525" , lineheight = 0.8)+
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0 to 0.2" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2.75 , col="#252525" , lineheight = 0.8)+
  #geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0.2 to 2" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2 , col="#737373" , lineheight = 0.8)+
  geom_text(data = filter(vkgrouprect.facets2.labels, z_layer == "N:C >0.2 to 1" ), aes(x = OC_x, y = HC_y, label = label, hjust = hjust, vjust = vjust ), parse = T, size = 2.75 , col="#252525" , lineheight = 0.8)+
  
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_text(size = rel(1)),
    legend.text = element_text(size = rel(0.9)) ,
    title = element_text(size = rel(0.9), lineheight = 1),
    strip.background = element_rect(fill = "transparent")
  )

p
# dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-Trend-with-Disease----All-Compounds-SELECT-T2D-VS-NORM--Wilcox-v2c-axis-large.tiff"), width = 20, height = 12, units = "cm", res=500, compression="lzw",type="cairo")
# dev.print(tiff, file = paste0(workdir,"/plots/","3d-Compounds-indiv-vKSpace-Trend-with-Disease----All-Compounds-SELECT-T2D-VS-NORM--Wilcox-v2c-axis-lim.tiff"), width = 22, height = 12, units = "cm", res=500, compression="lzw",type="cairo")

dev.print(tiff, file = paste0(workdir,"/plots/","noGgtitle-3d-Compounds-indiv-vKSpace-Trend-with-Disease----All-Compounds-CHN-Hostremoval-T2D-VS-NORM--Wilcox-v6-axis-lim.tiff"), width = 20, height = 12, units = "cm", res=500, compression="lzw",type="cairo")

#-------------------------


##########################
##########################
##########################
##########################

#### Merged results - w/ Host-removal
#### Tidy plot - Combined plot - consistent CPP trends across Restoration, T2D (SWE) and T2D (CHN) - w/ Host-removal
#-------------------------

# repeat statistical tests for tidy plotting
# store results in this format for facet_grid and annotation of significance results using 
sig_symb
# function(x) {
#   out <- character(length = length(x))
#   for (i in seq_along(x)) {
#     if (is.na(x[i])) {
#       out[i] <- NA
#     } else if (x[i] < 0.001) {
#       out[i] <- "***"
#     } else if (x[i] >= 0.001 & x[i] < 0.01) {
#       out[i] <- "**"
#     } else if (x[i] >= 0.01 & x[i] <= 0.05) {
#       out[i] <- "*"
#     } else {
#       out[i] <- "ns"
#     }
#   } # END loop
#   return(out)
# }
# <bytecode: 0x7fe35d1592e8>

annot.df2 <- data.frame(
  study = rep(c("Restoration","T2D (SWE)", "T2D (CHN)"),each = 5),
  compounds = rep(c("D-Fructose", "L-Arabinose", "Sugars", "Lignin", "BCFA-ACPs"), times =3), # \U03A3 = Sum
  y_val = c(1,2,3,4,5,6,7,8,9,10,1,2,3,4,5),
  x_label = rep(1.5, times = 15),
  x_seg_start = rep(1, times = 15),
  x_seg_end = rep(2, times = 15),
  #label = c( sig_symb( <p-values> ), ... ), # populate later
  stringsAsFactors = FALSE)
unique(annot.df2$compounds)
annot.df2$compounds <- factor(annot.df2$compounds,
                            levels = c("D-Fructose", "L-Arabinose", "Sugars", "Lignin", "BCFA-ACPs"),
                            ordered = TRUE)
annot.df2$study <- factor(annot.df2$study, levels = c("Restoration","T2D (SWE)", "T2D (CHN)"), ordered = TRUE)

annot.df2[ ,c(1,2)]
#          study   compounds
# 1  Restoration  D-Fructose
# 2  Restoration L-Arabinose
# 3  Restoration      Sugars
# 4  Restoration      Lignin
# 5  Restoration   BCFA-ACPs
# 6    T2D (SWE)  D-Fructose
# 7    T2D (SWE) L-Arabinose
# 8    T2D (SWE)      Sugars
# 9    T2D (SWE)      Lignin
# 10   T2D (SWE)   BCFA-ACPs
# 11   T2D (CHN)  D-Fructose
# 12   T2D (CHN) L-Arabinose
# 13   T2D (CHN)      Sugars
# 14   T2D (CHN)      Lignin
# 15   T2D (CHN)   BCFA-ACPs

# store test results in order 1-15
# x,y format is suitable for x-y scatterplot trends
x_values <- list()
y_values <- list()
p_values <- list()
stat_values <- list()
estimate_values <- list()

# x1,x2 - y1,y2 format needed for comparing two different groups (Wilcoxon-Mann-Whitney Tests)
x1_values <- list()
x2_values <- list()
y1_values <- list()
y2_values <- list()

#          study      compounds
## 1  Restoration     D-Fructose

dat.cpd.res <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-sunbad-resto.rds")
age_vec
# 6 12 22 31 UM 
# 1  2  3  4  5 

## Fructose
sel.cpd <- which(df.comp$name == "D-Fructose") # 
df.comp[sel.cpd, ] 
# id abbrev       name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 82 cpd00082    fru D-Fructose C6H12O6
sel <- which(dat.cpd.res$cpd_id == "cpd00082")
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]
# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
#x_values[[1]] <- x
x_values[[1]] <- dat.cpd.res$group_label[sel]
y_values[[1]] <- y
p_values[[1]] <- ktcor$p.value
stat_values[[1]] <- ktcor$statistic
estimate_values[[1]] <- ktcor$estimate

## 2  Restoration    L-Arabinose

sel.cpd <- which(df.comp$name == "L-Arabinose") # 
df.comp[sel.cpd, ]
#           id abbrev        name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio
# 223 cpd00224 arab-L L-Arabinose C5H10O5
sel <- which(dat.cpd.res$cpd_id == "cpd00224") # 223
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]
# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
#x_values[[2]] <- x
x_values[[2]] <- dat.cpd.res$group_label[sel]
y_values[[2]] <- y
p_values[[2]] <- ktcor$p.value
stat_values[[2]] <- ktcor$statistic
estimate_values[[2]] <- ktcor$estimate


## 3  Restoration Sugars
t2d.zoom2B.increasing.carbs # "cpd00224" "cpd03198" "cpd03200" "cpd00259" "cpd00076" "cpd00082"

df <- dat.cpd.res
sel <- which(df$cpd_id %in% t2d.zoom2B.increasing.carbs) # 90
df <- df[sel, ]
res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )
for (i in 1:length(unique(df$sample))) {
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  print(paste0("completed ",i))
}
str(res)
unique(res$group) # "22" "31" "UM" "12" "6"
res$group <- factor(res$group, levels = c("6","12", "22", "31", "UM"), ordered = TRUE)
unique(res$group_label) # "22 yr"   "31 yr"   "Unmined" "12 yr"   "6 yr"
res$group_label <- factor(res$group_label, levels = c("6 yr","12 yr", "22 yr", "31 yr", "Unmined"), ordered = TRUE)
res$age_vec <- as.integer(res$group)
# Kendall Tau correlation
ktcor<- cor.test(x = res$age_vec, y = res$sum_rel_abun, method = "kendall")
ktcor
x_values[[3]] <- res$group_label
y_values[[3]] <- res$sum_rel_abun
p_values[[3]] <- ktcor$p.value
stat_values[[3]] <- ktcor$statistic
estimate_values[[3]] <- ktcor$estimate


## 4  Restoration         Lignin

sel.cpd <- which(df.comp$name == "Lignin") # 
df.comp[sel.cpd, ]
#             id abbrev   name form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio MnC_ratio
# 12548 cpd12745 Lignin Lignin null
sel <- which(dat.cpd.res$cpd_id == "cpd12745")
x <- dat.cpd.res$ord_group[sel]
y <- dat.cpd.res$cpd_rel_abun[sel]
# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor
#x_values[[4]] <- x
x_values[[4]] <- dat.cpd.res$group_label[sel]
y_values[[4]] <- y
p_values[[4]] <- ktcor$p.value
stat_values[[4]] <- ktcor$statistic
estimate_values[[4]] <- ktcor$estimate


## 5  Restoration BCFA-ACPs

df <- dat.cpd.res
sel <- which(df$cpd_id %in% t2d.zoom1.decreasing.proteins) # 525
df <- df[sel, ]
res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )
for (i in 1:length(unique(df$sample))) {
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  print(paste0("completed ",i))
}
str(res)
unique(res$group) # "22" "31" "UM" "12" "6"
res$group <- factor(res$group, levels = c("6","12", "22", "31", "UM"), ordered = TRUE)
unique(res$group_label) # "22 yr"   "31 yr"   "Unmined" "12 yr"   "6 yr"
res$group_label <- factor(res$group_label, levels = c("6 yr","12 yr", "22 yr", "31 yr", "Unmined"), ordered = TRUE)
res$age_vec <- as.integer(res$group)
# Kendall Tau correlation
ktcor<- cor.test(x = res$age_vec, y = res$sum_rel_abun, method = "kendall")
ktcor
x_values[[5]] <- res$group_label
y_values[[5]] <- res$sum_rel_abun
p_values[[5]] <- ktcor$p.value
stat_values[[5]] <- ktcor$statistic
estimate_values[[5]] <- ktcor$estimate


## 6    T2D (SWE)     D-Fructose
dat.cpd.t2d <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-v6.rds")
# select only Normal and T2D
unique(dat.cpd.t2d$group_label)
# [1] T2D met- Normal  
# Levels: T2D met- < Normal

sel.cpd <- which(df.comp$name == "D-Fructose") # 
df.comp[sel.cpd, ]
#          id abbrev       name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 82 cpd00082    fru D-Fructose C6H12O6

sel <- which(dat.cpd.t2d$cpd_id == "cpd00082") # 82
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
wmw.test

x_values[[6]] <- x
y_values[[6]] <- y
p_values[[6]] <- wmw.test$p.value
stat_values[[6]] <- wmw.test$statistic # W

x1_values[[6]] <- rep( "T2D met-", times = length(x_values[[6]]) )
x2_values[[6]] <- rep( "Normal", times = length(y_values[[6]]) )
y1_values[[6]] <- x_values[[6]] # T2D
y2_values[[6]] <- y_values[[6]] # Normal



## 7    T2D (SWE)    L-Arabinose

sel.cpd <- which(df.comp$name == "L-Arabinose") # 
df.comp[sel.cpd, ]
#           id abbrev        name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio
# 223 cpd00224 arab-L L-Arabinose C5H10O5
this_var <- "L-Arabinose"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00224")
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
wmw.test

x_values[[7]] <- x # T2D CPP% values
y_values[[7]] <- y # Normal CPP% values
p_values[[7]] <- wmw.test$p.value
stat_values[[7]] <- wmw.test$statistic # W

x1_values[[7]] <- rep( "T2D met-", times = length(x_values[[7]]) )
x2_values[[7]] <- rep( "Normal", times = length(y_values[[7]]) )
y1_values[[7]] <- x_values[[7]] # T2D
y2_values[[7]] <- y_values[[7]] # Normal


## 8    T2D (SWE) Sugars

df <- dat.cpd.t2d
sel <- which(df$cpd_id %in% t2d.zoom2B.increasing.carbs)
df <- df[sel, ]
str(df) # 'data.frame':	456 obs. of  7 variables:
length(unique(df$cpd_id)) # 6
length(unique(df$sample)) # 76
res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )
for (i in 1:length(unique(df$sample))) {
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  print(paste0("completed ",i))
}
str(res) # 'data.frame':	76 obs. of  4 variables:
unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)
# Wilcoxon-Mann-Whitney Test
x <- res$sum_rel_abun[ res$group_label == "T2D met-" ] # n = 33
y <- res$sum_rel_abun[ res$group_label == "Normal" ] # n = 43
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun % #  "less"
wmw.test

x_values[[8]] <- x # T2D CPP% values
y_values[[8]] <- y # Normal CPP% values
p_values[[8]] <- wmw.test$p.value
stat_values[[8]] <- wmw.test$statistic # W

x1_values[[8]] <- rep( "T2D met-", times = length(x_values[[8]]) )
x2_values[[8]] <- rep( "Normal", times = length(y_values[[8]]) )
y1_values[[8]] <- x_values[[8]] # T2D
y2_values[[8]] <- y_values[[8]] # Normal


## 9    T2D (SWE)         Lignin

sel.cpd <- which(df.comp$name == "Lignin") # 
df.comp[sel.cpd, ]
#             id abbrev   name form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio MnC_ratio
# 12548 cpd12745 Lignin Lignin null
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd12745")
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 33
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 43
# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
wmw.test

x_values[[9]] <- x
y_values[[9]] <- y
p_values[[9]] <- wmw.test$p.value
stat_values[[9]] <- wmw.test$statistic # W

x1_values[[9]] <- rep( "T2D met-", times = length(x_values[[9]]) )
x2_values[[9]] <- rep( "Normal", times = length(y_values[[9]]) )
y1_values[[9]] <- x_values[[9]] # T2D
y2_values[[9]] <- y_values[[9]] # Normal


## 10   T2D (SWE) BCFA-ACPs

df <- dat.cpd.t2d
sel <- which(df$cpd_id %in% t2d.zoom1.decreasing.proteins)
df <- df[sel, ]
str(df) # 'data.frame':	2660 obs. of  7 variables:
length(unique(df$cpd_id)) # 35
length(unique(df$sample)) # 76
res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )
for (i in 1:length(unique(df$sample))) {
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  print(paste0("completed ",i))
}
str(res) # 'data.frame':	76 obs. of  4 variables:
unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)
# Wilcoxon-Mann-Whitney Test
x <- res$sum_rel_abun[ res$group_label == "T2D met-" ] # n = 33
y <- res$sum_rel_abun[ res$group_label == "Normal" ] # n = 43
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun % #  "greater"
wmw.test

x_values[[10]] <- x # T2D CPP% values
y_values[[10]] <- y # Normal CPP% values
p_values[[10]] <- wmw.test$p.value
stat_values[[10]] <- wmw.test$statistic # W

x1_values[[10]] <- rep( "T2D met-", times = length(x_values[[10]]) )
x2_values[[10]] <- rep( "Normal", times = length(y_values[[10]]) )
y1_values[[10]] <- x_values[[10]] # T2D
y2_values[[10]] <- y_values[[10]] # Normal




## 11   T2D (CHN)     D-Fructose

dat.cpd.t2d <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval.rds")
str(dat.cpd.t2d) # 'data.frame':	591466 obs. of  5 variables:
dat.cpd.t2d$group_label <- factor(dat.cpd.t2d$group, levels = c("T2D met-", "Normal"), ordered = TRUE)

## Fructose
sel.cpd <- which(df.comp$name == "D-Fructose") # 
df.comp[sel.cpd, ]
#          id abbrev       name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 82 cpd00082    fru D-Fructose C6H12O6

sel <- which(dat.cpd.t2d$cpd_id == "cpd00082")
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
wmw.test

x_values[[11]] <- x # T2D CPP% values
y_values[[11]] <- y # Normal CPP% values
p_values[[11]] <- wmw.test$p.value
stat_values[[11]] <- wmw.test$statistic # W

x1_values[[11]] <- rep( "T2D met-", times = length(x_values[[11]]) )
x2_values[[11]] <- rep( "Normal", times = length(y_values[[11]]) )
y1_values[[11]] <- x_values[[11]] # T2D
y2_values[[11]] <- y_values[[11]] # Normal


## 12   T2D (CHN)    L-Arabinose

sel.cpd <- which(df.comp$name == "L-Arabinose") # 
df.comp[sel.cpd, ]
#           id abbrev        name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio
# 223 cpd00224 arab-L L-Arabinose C5H10O5
this_var <- "L-Arabinose"

# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd00224")
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
wmw.test

x_values[[12]] <- x # T2D CPP% values
y_values[[12]] <- y # Normal CPP% values
p_values[[12]] <- wmw.test$p.value
stat_values[[12]] <- wmw.test$statistic # W

x1_values[[12]] <- rep( "T2D met-", times = length(x_values[[12]]) )
x2_values[[12]] <- rep( "Normal", times = length(y_values[[12]]) )
y1_values[[12]] <- x_values[[12]] # T2D
y2_values[[12]] <- y_values[[12]] # Normal



## 13   T2D (CHN) Sugars

df <- dat.cpd.t2d
sel <- which(df$cpd_id %in% t2d.zoom2B.increasing.carbs)
df <- df[sel, ]
str(df) # 'data.frame':	492 obs. of  6 variables:
length(unique(df$cpd_id)) # 6
length(unique(df$sample)) # 82
res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )
for (i in 1:length(unique(df$sample))) {
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  print(paste0("completed ",i))
}
str(res) # 'data.frame':	82 obs. of  4 variables:
unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)
# Wilcoxon-Mann-Whitney Test
x <- res$sum_rel_abun[ res$group_label == "T2D met-" ] # n = 30
y <- res$sum_rel_abun[ res$group_label == "Normal" ] # n = 52
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun % #  "less"
wmw.test

x_values[[13]] <- x # T2D CPP% values
y_values[[13]] <- y # Normal CPP% values
p_values[[13]] <- wmw.test$p.value
stat_values[[13]] <- wmw.test$statistic # W

x1_values[[13]] <- rep( "T2D met-", times = length(x_values[[13]]) )
x2_values[[13]] <- rep( "Normal", times = length(y_values[[13]]) )
y1_values[[13]] <- x_values[[13]] # T2D
y2_values[[13]] <- y_values[[13]] # Normal


## 14   T2D (CHN)         Lignin

sel.cpd <- which(df.comp$name == "Lignin") # 
df.comp[sel.cpd, ]
#             id abbrev   name form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio MnC_ratio
# 12548 cpd12745 Lignin Lignin null
# T2D
sel <- which(dat.cpd.t2d$cpd_id == "cpd12745")
subsel.t2d <- which(dat.cpd.t2d$group_label[sel] == "T2D met-")
subsel.norm <- which(dat.cpd.t2d$group_label[sel] == "Normal")
dat.cpd.t2d[sel ,]
dat.cpd.t2d[sel[subsel.t2d] ,]
dat.cpd.t2d[sel[subsel.norm] ,]
x <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.t2d]] # n = 30
y <- dat.cpd.t2d$cpd_rel_abun[sel[subsel.norm]] # n = 52
# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
wmw.test

x_values[[14]] <- x
y_values[[14]] <- y
p_values[[14]] <- wmw.test$p.value
stat_values[[14]] <- wmw.test$statistic # W

x1_values[[14]] <- rep( "T2D met-", times = length(x_values[[14]]) )
x2_values[[14]] <- rep( "Normal", times = length(y_values[[14]]) )
y1_values[[14]] <- x_values[[14]] # T2D
y2_values[[14]] <- y_values[[14]] # Normal


## 15   T2D (CHN) BCFA-ACPs

df <- dat.cpd.t2d
sel <- which(df$cpd_id %in% t2d.zoom1.decreasing.proteins)
df <- df[sel, ]
str(df) # 'data.frame':	2870 obs. of  6 variables:
length(unique(df$cpd_id)) # 35
length(unique(df$sample)) # 82
res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA, group_label = NA )
for (i in 1:length(unique(df$sample))) {
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  print(paste0("completed ",i))
}
str(res) # 'data.frame':	82 obs. of  4 variables:
unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)
# Wilcoxon-Mann-Whitney Test
x <- res$sum_rel_abun[ res$group_label == "T2D met-" ] # n = 30
y <- res$sum_rel_abun[ res$group_label == "Normal" ] # n = 52
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun % #  "greater"
wmw.test

x_values[[15]] <- x # T2D CPP% values
y_values[[15]] <- y # Normal CPP% values
p_values[[15]] <- wmw.test$p.value
stat_values[[15]] <- wmw.test$statistic # W

x1_values[[15]] <- rep( "T2D met-", times = length(x_values[[15]]) )
x2_values[[15]] <- rep( "Normal", times = length(y_values[[15]]) )
y1_values[[15]] <- x_values[[15]] # T2D
y2_values[[15]] <- y_values[[15]] # Normal

# plot with significance annotation
annot.df2 # note this matches order of list for saved p-values

# study   compounds y_val x_label x_seg_start x_seg_end
# 1  Restoration  D-Fructose     1     1.5           1         2
# 2  Restoration L-Arabinose     2     1.5           1         2
# 3  Restoration      Sugars     3     1.5           1         2
# 4  Restoration      Lignin     4     1.5           1         2
# 5  Restoration   BCFA-ACPs     5     1.5           1         2
# 6    T2D (SWE)  D-Fructose     6     1.5           1         2
# 7    T2D (SWE) L-Arabinose     7     1.5           1         2
# 8    T2D (SWE)      Sugars     8     1.5           1         2
# 9    T2D (SWE)      Lignin     9     1.5           1         2
# 10   T2D (SWE)   BCFA-ACPs    10     1.5           1         2
# 11   T2D (CHN)  D-Fructose     1     1.5           1         2
# 12   T2D (CHN) L-Arabinose     2     1.5           1         2
# 13   T2D (CHN)      Sugars     3     1.5           1         2
# 14   T2D (CHN)      Lignin     4     1.5           1         2
# 15   T2D (CHN)   BCFA-ACPs     5     1.5           1         2

annot.df2$p_val <- c(p_values[[1]], p_values[[2]], p_values[[3]], p_values[[4]], p_values[[5]],
                     p_values[[6]], p_values[[7]], p_values[[8]], p_values[[9]], p_values[[10]],
                     p_values[[11]], p_values[[12]], p_values[[13]], p_values[[14]], p_values[[15]]
                     )

annot.df2$sig_symb <- sig_symb( annot.df2$p_val )
annot.df2$sig_symb
# "***" "***" "***" "***" "***" "*"   "*"   "**"  "*"   "*"   "**"  "***" "***" "*"   "***"


annot.df2$estimate <- c(estimate_values[[1]],estimate_values[[2]],estimate_values[[3]],estimate_values[[4]],estimate_values[[5]],
                        NA,NA,NA,NA,NA,
                        NA,NA,NA,NA,NA)
annot.df2$text_npcx <- c("right","right","right","left","left",
                         NA,NA,NA,"left",NA,
                         NA,NA,NA,"left",NA)
annot.df2$text_npcy <- c("top","top","top","top","top",
                         NA,NA,NA,"top",NA,
                         NA,NA,NA,"top",NA)
annot.df2$ktcor_text <- c(paste0("Kendall's tau ~Age\nTau = ",round(annot.df2$estimate[1],3)," ",annot.df2$sig_symb[1]),
                          paste0("Kendall's tau ~Age\nTau = ",round(annot.df2$estimate[2],3)," ",annot.df2$sig_symb[2]),
                          paste0("Kendall's tau ~Age\nTau = ",round(annot.df2$estimate[3],3)," ",annot.df2$sig_symb[3]),
                          paste0("Kendall's tau ~Age\nTau = ",round(annot.df2$estimate[4],3)," ",annot.df2$sig_symb[4]),
                          paste0("Kendall's tau ~Age\nTau = ",round(annot.df2$estimate[5],3)," ",annot.df2$sig_symb[5]),
                          NA,NA,NA,NA,NA,
                          NA,NA,NA,NA,NA
                         )

annot.df2$compounds
# [1] D-Fructose  L-Arabinose Sugars      Lignin      BCFA-ACPs   D-Fructose  L-Arabinose Sugars      Lignin      BCFA-ACPs   D-Fructose  L-Arabinose
# [13] Sugars      Lignin      BCFA-ACPs  
# Levels: D-Fructose < L-Arabinose < Sugars < Lignin < BCFA-ACPs

# plot data
pdat <- data.frame(
  x_class = c( as.character(x_values[[1]]), as.character(x_values[[2]]), as.character(x_values[[3]]), as.character(x_values[[4]]), as.character(x_values[[5]]),
               x1_values[[6]], x2_values[[6]], x1_values[[7]], x2_values[[7]], x1_values[[8]], x2_values[[8]], x1_values[[9]], x2_values[[9]], x1_values[[10]], x2_values[[10]],
               x1_values[[11]], x2_values[[11]], x1_values[[12]], x2_values[[12]], x1_values[[13]], x2_values[[13]], x1_values[[14]], x2_values[[14]], x1_values[[15]], x2_values[[15]]
               ),
  
  y_cpp = c(y_values[[1]], y_values[[2]], y_values[[3]], y_values[[4]], y_values[[5]],
             y1_values[[6]], y2_values[[6]], y1_values[[7]], y2_values[[7]], y1_values[[8]], y2_values[[8]], y1_values[[9]], y2_values[[9]], y1_values[[10]], y2_values[[10]],
             y1_values[[11]], y2_values[[11]], y1_values[[12]], y2_values[[12]], y1_values[[13]], y2_values[[13]], y1_values[[14]], y2_values[[14]], y1_values[[15]], y2_values[[15]]
             ),
  
  study = c( rep( "Restoration", times = length( c( x_values[[1]], x_values[[2]], x_values[[3]], x_values[[4]], x_values[[5]]) ) ), # 75 = 15*5
             rep( "T2D (SWE)", times = length( c( x1_values[[6]], x2_values[[6]], x1_values[[7]], x2_values[[7]], x1_values[[8]], x2_values[[8]], x1_values[[9]], x2_values[[9]], x1_values[[10]], x2_values[[10]] ) )  ), # 380 = 76*5
             rep( "T2D (CHN)", times = length( c( x1_values[[11]], x2_values[[11]], x1_values[[12]], x2_values[[12]], x1_values[[13]], x2_values[[13]], x1_values[[14]], x2_values[[14]], x1_values[[15]], x2_values[[15]]) ) ) # 410 = 82*5
             ),
  compounds = c( rep( "D-Fructose", times = length(x_values[[1]]) ), # 15
                 rep( "L-Arabinose", times = length(x_values[[2]]) ), # 15
                 rep( "Sugars", times = length(x_values[[3]]) ), # 15
                 rep( "Lignin", times = length(x_values[[4]]) ), # 15
                 rep( "BCFA-ACPs", times = length(x_values[[5]]) ), # 15
                 
                 rep( "D-Fructose", times = length( c( x1_values[[6]], x2_values[[6]] ) ) ), # 76
                 rep( "L-Arabinose", times = length( c( x1_values[[7]], x2_values[[7]] ) ) ), # 76
                 rep( "Sugars", times = length( c( x1_values[[8]], x2_values[[8]] ) ) ), # 76
                 rep( "Lignin", times = length( c( x1_values[[9]], x2_values[[9]] ) ) ), # 76
                 rep( "BCFA-ACPs", times = length( c( x1_values[[10]], x2_values[[10]] ) ) ), # 76
                 
                 rep( "D-Fructose", times = length( c( x1_values[[11]], x2_values[[11]] ) ) ), # 76
                 rep( "L-Arabinose", times = length( c( x1_values[[12]], x2_values[[12]] ) ) ), # 76
                 rep( "Sugars", times = length( c( x1_values[[13]], x2_values[[13]] ) ) ), # 76
                 rep( "Lignin", times = length( c( x1_values[[14]], x2_values[[14]] ) ) ), # 76
                 rep( "BCFA-ACPs", times = length( c( x1_values[[15]], x2_values[[15]] ) ) ) # 76
                 )
  )

unique( pdat$x_class )
# "22 yr"    "31 yr"    "Unmined"  "12 yr"    "6 yr"     "T2D met-" "Normal" 
pdat$x_class <- factor(pdat$x_class, levels = c("6 yr", "12 yr", "22 yr", "31 yr", "Unmined", "T2D met-", "Normal"), ordered = TRUE )
pdat$x_class2 <- factor(pdat$x_class, 
                        levels = c("6 yr", "12 yr", "22 yr", "31 yr", "Unmined", "T2D met-", "Normal"),
                        labels = c("6 yr", "12 yr", "22 yr", "31 yr", "UM", "T2D met-", "Normal"),
                        ordered = TRUE )
pdat$x_class_integer <- as.integer( pdat$x_class2 )

unique(pdat$study) # "Restoration" "T2D (SWE)"   "T2D (CHN)"  
pdat$study <- factor(pdat$study, levels = c("Restoration", "T2D (SWE)", "T2D (CHN)"), ordered = TRUE)

unique(pdat$compounds) # "D-Fructose"  "L-Arabinose" "Sugars"      "Lignin"      "BCFA-ACPs" 
pdat$compounds <- factor(pdat$compounds, levels = c("D-Fructose", "L-Arabinose", "Sugars", "Lignin", "BCFA-ACPs"), ordered = TRUE)

# include log10 transform for Lignin:
pdat$y_cpp_withlog10_lignin <- pdat$y_cpp
sel <- which(pdat$compounds == "Lignin") # qty 173
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(pdat$y_cpp_withlog10_lignin[sel] == 0) # qty 63
#if (length(subsel.zero) > 0) {
zero_replace <- 0.5*min( pdat$y_cpp_withlog10_lignin[sel [ -subsel.zero ]] )
pdat$y_cpp_withlog10_lignin[sel [ subsel.zero ]] <- zero_replace
#}
pdat$y_cpp_withlog10_lignin[sel] <- log10(pdat$y_cpp_withlog10_lignin[sel])

pdat[sel, ]

annot.df2 
# adjust p_val and sig_symb using log10 data for Lignin 4, 9, 14
annot.df2$p_val_withlog10_lignin <- annot.df2$p_val
annot.df2[c(4,9,14), ]
# study compounds y_val x_label x_seg_start x_seg_end        p_val sig_symb  estimate text_npcx text_npcy                          ktcor_text
# 4  Restoration    Lignin     4     1.5           1         2 7.825459e-05      *** 0.8023774      left       top Kendall's tau ~Age\nTau = 0.802 ***
# 9    T2D (SWE)    Lignin     9     1.5           1         2 2.687130e-02        *        NA      left       top                                <NA>
# 14   T2D (CHN)    Lignin     4     1.5           1         2 2.805392e-02        *        NA      left       top                                <NA>
#    p_val_withlog10_lignin
# 4            7.825459e-05
# 9            2.687130e-02
# 14           2.805392e-02

head(pdat)
#   x_class      y_cpp       study  compounds x_class2 x_class_integer y_cpp_withlog10_lignin
# 1   22 yr 0.02919476 Restoration D-Fructose    22 yr               3             0.02919476
# 2   31 yr 0.02061704 Restoration D-Fructose    31 yr               4             0.02061704
# 3   31 yr 0.02424720 Restoration D-Fructose    31 yr               4             0.02424720
# 4 Unmined 0.02005817 Restoration D-Fructose       UM               5             0.02005817
# 5   12 yr 0.02661002 Restoration D-Fructose    12 yr               2             0.02661002
# 6    6 yr 0.03193183 Restoration D-Fructose     6 yr               1             0.03193183

## 4  Restoration    Lignin
sel <- which(pdat$study == "Restoration" & pdat$compounds == "Lignin") # 15
x <- as.integer( pdat$x_class[sel] )
y <- pdat$y_cpp_withlog10_lignin[sel]
# Kendall Tau correlation
ktcor<- cor.test(x = x, y = y, method = "kendall")
ktcor$p.value # 7.825459e-05
annot.df2$p_val_withlog10_lignin[4] # 7.825459e-05 = same
annot.df2$p_val[4] # 7.825459e-05

## 9    T2D (SWE)    Lignin 
sel <- which(pdat$study == "T2D (SWE)" & pdat$compounds == "Lignin") # 76
subsel.t2d <- which(pdat$x_class[sel] == "T2D met-") # 33
subsel.norm <- which(pdat$x_class[sel] == "Normal") # 43
x <- pdat$y_cpp_withlog10_lignin[sel[subsel.t2d]]
y <- pdat$y_cpp_withlog10_lignin[sel[subsel.norm]]
# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
wmw.test$p.value # 0.0268713
annot.df2$p_val_withlog10_lignin[9] # 0.0268713 = same
annot.df2$p_val[9] # 0.0268713 = same

## 14   T2D (CHN)    Lignin
sel <- which(pdat$study == "T2D (CHN)" & pdat$compounds == "Lignin") # 82
subsel.t2d <- which(pdat$x_class[sel] == "T2D met-") # 30
subsel.norm <- which(pdat$x_class[sel] == "Normal") # 52
x <- pdat$y_cpp_withlog10_lignin[sel[subsel.t2d]]
y <- pdat$y_cpp_withlog10_lignin[sel[subsel.norm]]
# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
wmw.test$p.value # 0.02805392
annot.df2$p_val_withlog10_lignin[14] # 0.02805392 = same
annot.df2$p_val[14] # 0.02805392


p <- ggplot(data = pdat, aes(x=x_class, y = y_cpp))+
  geom_violin(alpha = 0.3)+
  geom_boxplot(width = 0.2, outlier.shape = NA )+ # width = 0.2, alpha = 0.3 outlier.shape = 1
  
  #facet_grid(rows = vars(compounds), cols = vars(study), scales = "free")+ # , scales = "free_y"
  facet_grid(rows = vars(compounds), cols = vars(study), scales = "free", space = "free_x")+ # , scales = "free_y"
  
  #geom_text(data= annot.df, mapping = aes(x = x_label, y = y_val, label = label), size = 5, vjust = 0.3)+ #  label = "**", size = 6   vjust = "bottom"
  #geom_segment(data=annot.df, mapping = aes(x = x_seg_start, y = y_val, xend = x_seg_end, yend = y_val))+
  
  xlab(NULL)+ ylab("CPP rel abun (%)")+
  theme_classic()+
  theme(#panel.grid.major = element_blank(), 
    strip.background = element_rect(fill="white", linetype = "blank"),
    strip.text = element_text(size = rel(0.92)), # margin=margin(t = 0,r = 0,b = 0,l = 0,"pt")
    axis.text.x = element_text(size = rel(1.1)) #,
  )
p

# split this in parts
# 1)
p <- #ggplot(data = filter(pdat, study == "Restoration"), aes(x=x_class_integer, y = y_cpp_withlog10_lignin))+
  ggplot(data = filter(pdat, study == "Restoration" & compounds %in% c("Sugars", "Lignin", "BCFA-ACPs")), aes(x=x_class_integer, y = y_cpp_withlog10_lignin))+
  
  #geom_violin(alpha = 0.3)+
  #geom_boxplot(width = 0.2, outlier.shape = NA )+ # width = 0.2, alpha = 0.3 outlier.shape = 1
  
  geom_point(shape = 1)+
  geom_smooth(method="loess")+
  
  #facet_grid(rows = vars(compounds), cols = vars(study), scales = "free")+ # , scales = "free_y"
  facet_grid(rows = vars(compounds), cols = vars(study), scales = "free", space = "free_x")+ # , scales = "free_y"
  
  scale_x_continuous(labels= levels(pdat$x_class2)[1:5])+
  
  #annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  geom_text_npc(data= filter(annot.df2, study == "Restoration" & compounds %in% c("Sugars", "Lignin", "BCFA-ACPs")), mapping = aes(npcx = text_npcx, npcy = text_npcy, label = ktcor_text), size = 3 , lineheight = 0.85 )+
  
  #geom_text(data= annot.df, mapping = aes(x = x_label, y = y_val, label = label), size = 5, vjust = 0.3)+ #  label = "**", size = 6   vjust = "bottom"
  #geom_segment(data=annot.df, mapping = aes(x = x_seg_start, y = y_val, xend = x_seg_end, yend = y_val))+
  
  xlab(NULL)+ ylab("CPP rel abun (%)")+
  theme_classic()+
  theme(
    strip.background = element_rect(fill="white", linetype = "blank"),
    strip.text = element_text(size = rel(0.92)), # margin=margin(t = 0,r = 0,b = 0,l = 0,"pt")
    
    strip.text.y.right = element_blank(),
    axis.text.x = element_text(size = rel(1.1)) #,
  )
p
#dev.print(tiff, file = paste0(workdir,"/plots/","Part1-Sunbad-SWE-CHN-T2D-consistent-CPP-trends-v4.tiff"), width = 8, height = 16, units = "cm", res=600, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","Part1-Sunbad-SWE-CHN-T2D-consistent-CPP-trends-Sugar-Lignin-BCFA-v4.tiff"), width = 8, height = 10, units = "cm", res=600, compression="lzw",type="cairo")

# 2)

# adjust height of significance annotation
annot.df2$y_val <- c(NA, NA, NA, NA, NA,
                     0.42, 0.34, 1.0, -3.8, 0.007,
                     0.56, 0.56, 1.44, -3.8, 0.011)

annot.df2$x_seg_start = rep(1.1, times = 15)  # rep(1, times = 15)
annot.df2$x_seg_end = rep(1.9, times = 15)   # rep(2, times = 15)

p <- #ggplot(data = filter(pdat, study %in% c("T2D (SWE)", "T2D (CHN)") ), aes(x=x_class2, y = y_cpp_withlog10_lignin))+
  ggplot(data = filter(pdat, study %in% c("T2D (SWE)", "T2D (CHN)") & compounds %in% c("Sugars", "Lignin", "BCFA-ACPs") ), aes(x=x_class2, y = y_cpp_withlog10_lignin))+
  
  geom_violin(alpha = 0.3)+
  geom_boxplot(width = 0.2, outlier.shape = NA )+ # width = 0.2, alpha = 0.3 outlier.shape = 1
  
  #facet_grid(rows = vars(compounds), cols = vars(study), scales = "free")+ # , scales = "free_y"
  facet_grid(rows = vars(compounds), cols = vars(study), scales = "free", space = "free_x")+ # , scales = "free_y"
  
  geom_text(data= filter(annot.df2, study %in% c("T2D (SWE)", "T2D (CHN)") & compounds %in% c("Sugars", "Lignin", "BCFA-ACPs") ), mapping = aes(x = x_label, y = y_val, label = sig_symb), size = 5, vjust = 0.3)+ #  label = "**", size = 6   vjust = "bottom"
  geom_segment(data=filter(annot.df2, study %in% c("T2D (SWE)", "T2D (CHN)") & compounds %in% c("Sugars", "Lignin", "BCFA-ACPs") ), mapping = aes(x = x_seg_start, y = y_val, xend = x_seg_end, yend = y_val))+
  
  #geom_text_npc(data= filter(annot.df2, study %in% c("T2D (SWE)", "T2D (CHN)") ), mapping = aes(npcx = text_npcx, npcy = text_npcy, label = ktcor_text), size = 3 , lineheight = 0.85 )+
  #geom_text_npc(data= filter(annot.df2, study %in% c("T2D (SWE)", "T2D (CHN)") & compounds == "Lignin" ), mapping = aes(npcx = text_npcx, npcy = text_npcy, label = ktcor_text), size = 3 , lineheight = 0.85 , vjust = 1)+ # nudge_y = 1 not working
  geom_text(data= filter(annot.df2, study %in% c("T2D (SWE)", "T2D (CHN)") & compounds == "Lignin" ), x = 0.5, y = -2.7, label = "Log10-transformed", size = 2.8 , vjust = 1, hjust = 0 )+ # nudge_y = 1 not working
  
  #xlab(NULL)+ ylab("CPP rel abun (%)")+
  xlab(NULL)+ ylab(NULL)+
  
  theme_classic()+
  theme(
    strip.background = element_rect(fill="white", linetype = "blank"),
    strip.text = element_text(size = rel(0.92)), # margin=margin(t = 0,r = 0,b = 0,l = 0,"pt")
    axis.text.x = element_text(size = rel(1.1)) #,
  )
p
#dev.print(tiff, file = paste0(workdir,"/plots/","Part2-Sunbad-SWE-CHN-T2D-consistent-CPP-trends-v4.tiff"), width = 12, height = 16, units = "cm", res=600, compression="lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","Part2-Sunbad-SWE-CHN-T2D-consistent-CPP-trends-Sugar-Lignin-BCFA-v4.tiff"), width = 12, height = 10, units = "cm", res=600, compression="lzw",type="cairo")



## Print panels individually
# study == "Restoration", "T2D (SWE)", "T2D (CHN)"  & compounds %in% c("Sugars", "Lignin", "BCFA-ACPs")
 
# a)
p <- #ggplot(data = filter(pdat, study == "Restoration"), aes(x=x_class_integer, y = y_cpp_withlog10_lignin))+
  ggplot(data = filter(pdat, study == "Restoration" & compounds == "Sugars"), aes(x=x_class_integer, y = y_cpp_withlog10_lignin))+
  geom_point(shape = 1)+
  geom_smooth(method="loess")+
  scale_x_continuous(labels= levels(pdat$x_class2)[1:5])+
  geom_text_npc(data= filter(annot.df2, study == "Restoration" & compounds == "Sugars"), mapping = aes(npcx = text_npcx, npcy = text_npcy, label = ktcor_text), size = 3 , lineheight = 0.85 )+
  xlab(NULL)+ ylab(NULL)+ # + ylab("CPP rel abun (%)")+
  theme_classic()+
  theme( axis.text.x = element_text(size = rel(1.1)) )
p
#grid.text(label = "(a)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Part-a-consistent-CPP-trends-Sugar-Lignin-BCFA-v4b.tiff"), width = 8, height = 3.8, units = "cm", res=600, compression="lzw",type="cairo")


# d)
p <- #ggplot(data = filter(pdat, study == "Restoration"), aes(x=x_class_integer, y = y_cpp_withlog10_lignin))+
  ggplot(data = filter(pdat, study == "Restoration" & compounds == "Lignin"), aes(x=x_class_integer, y = y_cpp_withlog10_lignin))+
  geom_point(shape = 1)+
  geom_smooth(method="loess")+
  scale_x_continuous(labels= levels(pdat$x_class2)[1:5])+
  geom_text_npc(data= filter(annot.df2, study == "Restoration" & compounds == "Lignin"), mapping = aes(label = ktcor_text), npcx = "right", npcy = "bottom", size = 3 , lineheight = 0.85 )+ #mapping = aes(npcx = text_npcx, npcy = text_npcy, label = ktcor_text), size = 3 , lineheight = 0.85 )+
  xlab(NULL)+ ylab(NULL)+ # + ylab("CPP rel abun (%)")+
  theme_classic()+
  theme( axis.text.x = element_text(size = rel(1.1)) )
p
dev.print(tiff, file = paste0(workdir,"/plots/","Part-d-consistent-CPP-trends-Sugar-Lignin-BCFA-v4b.tiff"), width = 8, height = 3.8, units = "cm", res=600, compression="lzw",type="cairo")


# g)
p <- #ggplot(data = filter(pdat, study == "Restoration"), aes(x=x_class_integer, y = y_cpp_withlog10_lignin))+
  ggplot(data = filter(pdat, study == "Restoration" & compounds == "BCFA-ACPs"), aes(x=x_class_integer, y = y_cpp_withlog10_lignin))+
  geom_point(shape = 1)+
  geom_smooth(method="loess")+
  scale_x_continuous(labels= levels(pdat$x_class2)[1:5])+
  geom_text_npc(data= filter(annot.df2, study == "Restoration" & compounds == "BCFA-ACPs"), mapping = aes(label = ktcor_text), npcx = "right", npcy = "bottom", size = 3 , lineheight = 0.85 )+ #mapping = aes(npcx = text_npcx, npcy = text_npcy, label = ktcor_text), size = 3 , lineheight = 0.85 )+
  xlab(NULL)+ ylab(NULL)+ # + ylab("CPP rel abun (%)")+
  theme_classic()+
  theme( axis.text.x = element_text(size = rel(1.1)) )
p
dev.print(tiff, file = paste0(workdir,"/plots/","Part-g-consistent-CPP-trends-Sugar-Lignin-BCFA-v4b.tiff"), width = 8, height = 3.8, units = "cm", res=600, compression="lzw",type="cairo")


# b)
# adjust height of significance annotation
annot.df2$y_val <- c(NA, NA, NA, NA, NA,            # Restoration
                     0.42, 0.34, 1.0, -3.8, 0.007,  # T2D SWE, Lignin is -3.8
                     0.56, 0.56, 1.44, -3.8, 0.011) # T2D CHN, Lignin is -3.8

p <- ggplot(data = filter(pdat, study == "T2D (SWE)" & compounds == "Sugars" ), aes(x=x_class2, y = y_cpp_withlog10_lignin))+
  ylim(0.25,1.6)+
  geom_violin(alpha = 0.3)+
  geom_boxplot(width = 0.2, outlier.shape = NA )+ # width = 0.2, alpha = 0.3 outlier.shape = 1
  geom_text(data= filter(annot.df2, study == "T2D (SWE)" & compounds == "Sugars" ), mapping = aes(x = x_label, y = y_val, label = sig_symb), size = 5, vjust = 0.3)+ #  label = "**", size = 6   vjust = "bottom"
  geom_segment(data=filter(annot.df2, study == "T2D (SWE)" & compounds == "Sugars" ), mapping = aes(x = x_seg_start, y = y_val, xend = x_seg_end, yend = y_val))+
  xlab(NULL)+ ylab(NULL)+
  theme_classic()+
  theme( axis.text.x = element_text(size = rel(1.1)) )
p
dev.print(tiff, file = paste0(workdir,"/plots/","Part-b-consistent-CPP-trends-Sugar-Lignin-BCFA-v6.tiff"), width = 5.8, height = 3.8, units = "cm", res=600, compression="lzw",type="cairo")

# c)
# adjust height of significance annotation
annot.df2$y_val <- c(NA, NA, NA, NA, NA,            # Restoration
                     0.42, 0.34, 1.0, -3.8, 0.007,  # T2D SWE, Lignin is -3.8
                     0.56, 0.56, 1.55, -3.8, 0.011) # T2D CHN, Lignin is -3.8

p <- ggplot(data = filter(pdat, study == "T2D (CHN)" & compounds == "Sugars" ), aes(x=x_class2, y = y_cpp_withlog10_lignin))+
  ylim(0.25,1.6)+
  geom_violin(alpha = 0.3)+
  geom_boxplot(width = 0.2, outlier.shape = NA )+ # width = 0.2, alpha = 0.3 outlier.shape = 1
  geom_text(data= filter(annot.df2, study == "T2D (CHN)" & compounds == "Sugars" ), mapping = aes(x = x_label, y = y_val, label = sig_symb), size = 5, vjust = 0.3)+ #  label = "**", size = 6   vjust = "bottom"
  geom_segment(data=filter(annot.df2, study == "T2D (CHN)" & compounds == "Sugars" ), mapping = aes(x = x_seg_start, y = y_val, xend = x_seg_end, yend = y_val))+
  xlab(NULL)+ ylab(NULL)+
  theme_classic()+
  theme( axis.text.x = element_text(size = rel(1.1)) )
p
dev.print(tiff, file = paste0(workdir,"/plots/","Part-c-consistent-CPP-trends-Sugar-Lignin-BCFA-v6.tiff"), width = 5.8, height = 3.8, units = "cm", res=600, compression="lzw",type="cairo")


# e)
# adjust height of significance annotation
annot.df2$y_val <- c(NA, NA, NA, NA, NA,            # Restoration
                     0.42, 0.34, 1.0, -2.9, 0.007,  # T2D SWE, Lignin is -3.8
                     0.56, 0.56, 1.44, -2.9, 0.011) # T2D CHN, Lignin is -3.8

p <- ggplot(data = filter(pdat, study == "T2D (SWE)" & compounds == "Lignin" ), aes(x=x_class2, y = y_cpp_withlog10_lignin))+
  ylim(-6.3,-2.7)+
  geom_violin(alpha = 0.3)+
  geom_boxplot(width = 0.2, outlier.shape = NA )+ # width = 0.2, alpha = 0.3 outlier.shape = 1
  geom_text(data= filter(annot.df2, study == "T2D (SWE)" & compounds == "Lignin" ), mapping = aes(x = x_label, y = y_val, label = sig_symb), size = 5, vjust = 0.3)+ #  label = "**", size = 6   vjust = "bottom"
  geom_segment(data=filter(annot.df2, study == "T2D (SWE)" & compounds == "Lignin" ), mapping = aes(x = x_seg_start, y = y_val, xend = x_seg_end, yend = y_val))+
  xlab(NULL)+ ylab(NULL)+
  theme_classic()+
  theme( axis.text.x = element_text(size = rel(1.1)) )
p
dev.print(tiff, file = paste0(workdir,"/plots/","Part-e-consistent-CPP-trends-Sugar-Lignin-BCFA-v6.tiff"), width = 5.8, height = 3.8, units = "cm", res=600, compression="lzw",type="cairo")

# f)

p <- ggplot(data = filter(pdat, study == "T2D (CHN)" & compounds == "Lignin" ), aes(x=x_class2, y = y_cpp_withlog10_lignin))+
  ylim(-6.3,-2.7)+
  geom_violin(alpha = 0.3)+
  geom_boxplot(width = 0.2, outlier.shape = NA )+ # width = 0.2, alpha = 0.3 outlier.shape = 1
  geom_text(data= filter(annot.df2, study == "T2D (CHN)" & compounds == "Lignin" ), mapping = aes(x = x_label, y = y_val, label = sig_symb), size = 5, vjust = 0.3)+ #  label = "**", size = 6   vjust = "bottom"
  geom_segment(data=filter(annot.df2, study == "T2D (CHN)" & compounds == "Lignin" ), mapping = aes(x = x_seg_start, y = y_val, xend = x_seg_end, yend = y_val))+
  xlab(NULL)+ ylab(NULL)+
  theme_classic()+
  theme( axis.text.x = element_text(size = rel(1.1)) )
p
dev.print(tiff, file = paste0(workdir,"/plots/","Part-f-consistent-CPP-trends-Sugar-Lignin-BCFA-v6.tiff"), width = 5.8, height = 3.8, units = "cm", res=600, compression="lzw",type="cairo")



# h)
# adjust height of significance annotation
annot.df2$y_val <- c(NA, NA, NA, NA, NA,            # Restoration
                     0.42, 0.34, 1.0, -2.9, 0.0065,  # T2D SWE, Lignin is -3.8
                     0.56, 0.56, 1.44, -2.9, 0.0128) # T2D CHN, Lignin is -3.8

p <- ggplot(data = filter(pdat, study == "T2D (SWE)" & compounds == "BCFA-ACPs" ), aes(x=x_class2, y = y_cpp_withlog10_lignin))+
  ylim(0,0.0133)+
  geom_violin(alpha = 0.3)+
  geom_boxplot(width = 0.2, outlier.shape = NA )+ # width = 0.2, alpha = 0.3 outlier.shape = 1
  geom_text(data= filter(annot.df2, study == "T2D (SWE)" & compounds == "BCFA-ACPs" ), mapping = aes(x = x_label, y = y_val, label = sig_symb), size = 5, vjust = 0.3)+ #  label = "**", size = 6   vjust = "bottom"
  geom_segment(data=filter(annot.df2, study == "T2D (SWE)" & compounds == "BCFA-ACPs" ), mapping = aes(x = x_seg_start, y = y_val, xend = x_seg_end, yend = y_val))+
  xlab(NULL)+ ylab(NULL)+
  theme_classic()+
  theme( axis.text.x = element_text(size = rel(1.1)) )
p
dev.print(tiff, file = paste0(workdir,"/plots/","Part-h-consistent-CPP-trends-Sugar-Lignin-BCFA-v6.tiff"), width = 5.8, height = 3.8, units = "cm", res=600, compression="lzw",type="cairo")

# i)

p <- ggplot(data = filter(pdat, study == "T2D (CHN)" & compounds == "BCFA-ACPs" ), aes(x=x_class2, y = y_cpp_withlog10_lignin))+
  ylim(0,0.0133)+
  geom_violin(alpha = 0.3)+
  geom_boxplot(width = 0.2, outlier.shape = NA )+ # width = 0.2, alpha = 0.3 outlier.shape = 1
  geom_text(data= filter(annot.df2, study == "T2D (CHN)" & compounds == "BCFA-ACPs" ), mapping = aes(x = x_label, y = y_val, label = sig_symb), size = 5, vjust = 0.3)+ #  label = "**", size = 6   vjust = "bottom"
  geom_segment(data=filter(annot.df2, study == "T2D (CHN)" & compounds == "BCFA-ACPs" ), mapping = aes(x = x_seg_start, y = y_val, xend = x_seg_end, yend = y_val))+
  xlab(NULL)+ ylab(NULL)+
  theme_classic()+
  theme( axis.text.x = element_text(size = rel(1.1)) )
p
dev.print(tiff, file = paste0(workdir,"/plots/","Part-i-consistent-CPP-trends-Sugar-Lignin-BCFA-v6.tiff"), width = 5.8, height = 3.8, units = "cm", res=600, compression="lzw",type="cairo")


#-------------------------


##########################
##########################
##########################
##########################

#### Add function info to full tables of compounds with trending CPP
#-------------------------

# get sig results tables for each case study
# trim to keep key columns
# join function info - relevant to each case study
# export as table, ordered by P-values


## Sunbad resto
dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-sunbad-resto.rds")
head(dat.test.sig)
sel.sigBH <- which( dat.test.sig$sigBH == "sig") # 2122
resto.sig.cpds <- unique(dat.test.sig$cpd[sel.sigBH]) # qty 2122

# link compounds to functions
df.out <-readRDS( "df.out--tidy-compounds_indiv--cpp3d-sunbad-resto.RDS" )

str(df.out)
# 'data.frame':	1154547 obs. of  9 variables:
# $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_3" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylaconitate isomerase" "2-methylaconitate isomerase" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" ...
# $ rxn_id             : chr  "rxn25278" "rxn25278" "rxn25279" "rxn25279" ...
# $ cpd_id             : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ cpd_name           : chr  "2-methyl-trans-aconitate" "cis-2-Methylaconitate" "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" ...
# $ cpd_form           : chr  "C7H5O6" "C7H5O6" "C7H7O7" "H2O" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.333 0.333 0.333 ...

phy <- readRDS("phy-phyloseq-object-sunbad-resto.RDS")
df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 30125     4

min(taxa_sums(phy)) # 2.521458e-06


## tidy sig test results table and add function info
temp <- dat.test.sig
# only include sigBH results
sel.sigBH <- which( temp$sigBH == "sig") # 2122
temp <- temp[sel.sigBH, ]

head(temp)
names(temp)
# [1] "cpd"               "data_for_this_cpd" "p_val"             "kendall_tau"       "trend_with_age"    "sigBH"             "minuslog10_p_val" 
# [8] "cpd_names"         "cpd_forms"         "OC_x"              "HC_y"              "NC_z"              "z_layer"           "mass"       

temp <- temp[ ,c(
  "cpd", "p_val", "kendall_tau", "sigBH", 
  "trend_with_age", "minuslog10_p_val",  "cpd_names", "cpd_forms",
  "OC_x", "HC_y", "NC_z", "z_layer")]

names(temp) <- c(
  "cpd", "p_val_kendall", "kendall_tau", "sigBH", 
  "trend_with_age", "minuslog10_p_val_kendall",  "cpd_names", "cpd_forms",
  "OC_x", "HC_y", "NC_z", "z_layer")

dim(temp) #  2122 12

# add functions
temp$fxns <- NA
#temp$superfocus_fxns <- NA
temp$sub_fxns <- NA
temp$subsys_L3 <- NA
temp$no_linked_superfocus_fxns <- NA
temp$no_linked_subsys_L3 <- NA

for (i in 1:dim(temp)[1]) {
  #i<-1
  sel <- which(df.out$cpd_id == temp$cpd[i]) # 2
  ( list_fxns <- unique( df.out$superfocus_fxn[sel] ) ) # "fxn_9950"  "fxn_9951"
  temp$fxns[i] <- paste0(list_fxns, collapse = "  |  ")
  temp$sub_fxns[i] <- paste0( unique(df.out$f__in[sel]), collapse = "  |  ")
  
  sel <- which(row.names(df.tax) %in% list_fxns) # 2
  
  #temp$superfocus_fxns[i] <- paste0( unique(df.tax$fxn[sel]), collapse = "  |  ")
  temp$subsys_L3[i] <- paste0( unique(df.tax$subsys_L3[sel]), collapse = "  |  ")
  
  temp$no_linked_superfocus_fxns[i] <- length( unique(df.tax$fxn[sel]) )
  temp$no_linked_subsys_L3[i] <- length( unique(df.tax$subsys_L3[sel]) )
  
  print(paste0("completed ",i))
}
#temp$superfocus_fxns <- gsub(pattern = ",", replacement = ";", x = temp$superfocus_fxns)
temp$sub_fxns <- gsub(pattern = ",", replacement = ";", x = temp$sub_fxns)

head(temp)



saveRDS(object = temp, file = "temp-Sig-results-Restoration-data-linked-fxns-subsyL3-v6.RDS")
temp <- readRDS("temp-Sig-results-Restoration-data-linked-fxns-subsyL3-v6.RDS")
sel <- which(temp$no_linked_superfocus_fxns < 10 & temp$no_linked_subsys_L3 > 400 )
temp$cpd[sel] # "cpd12559" "cpd29932" "cpd30045" "cpd30050" "cpd30057" - these include 'hypothetical_protein'

dim(phy@tax_table) #  30125     4

colnames(phy@tax_table) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn" 
length(unique(phy@tax_table[ , "subsys_L1" ])) #
length(unique(phy@tax_table[ , "subsys_L2" ])) #
length(unique(phy@tax_table[ , "subsys_L3" ])) #
length(unique(phy@tax_table[ , "fxn" ])) #

write.table(x = temp, file = "Sig-results-table-Restoration-compounds-with-functions-v6.tsv", sep = "\t", quote = FALSE, row.names = FALSE )
write.table(x = temp[ order(temp$p_val_kendall, decreasing = FALSE) , ], file = "Sig-results-table-Restoration-compounds-with-functions-ordered-v6.tsv", sep = "\t", quote = FALSE, row.names = FALSE )



## T2D-SWE

dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-Forslund-SWE-T2D-Hostremoval--All-Compounds-SELECT-T2D-vs-NORM--Wilcox-v6.rds")
str(dat.test.sig)

# link compounds to functions
df.out <-readRDS( "df.out--tidy-compounds_indiv-cpp3d--forslund-t2d-swe-hostremoval-pawsey.RDS" )
str(df.out)
# 'data.frame':	507391 obs. of  9 variables:
# $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_2" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" ...
# $ rxn_id             : chr  "rxn25279" "rxn25279" "rxn25279" "rxn03060" ...
# $ cpd_id             : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ cpd_name           : chr  "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" "2-methyl-trans-aconitate" "2-Methylcitrate" ...
# $ cpd_form           : chr  "C7H7O7" "H2O" "C7H5O6" "C7H7O7" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.3333 0.3333 0.3333 0.0556 0.0556 ...

phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval.RDS")
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 17912 taxa and 76 samples ]
# sample_data() Sample Data:       [ 76 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 17912 taxa by 4 taxonomic ranks ]
df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 17912    4
table(phy@sam_data$Status)
# ND CTRL T2D metformin- 
#   43             33


## tidy sig test results table and add function info
temp <- dat.test.sig
head(temp)
names(temp)
# [1] "cpd"                "data_for_this_cpd"  "median_t2d"         "median_normal"      "alt"                "p_val"              "W_statistic"       
# [8] "hl_effect_wilcox"   "trend_with_disease" "minuslog10_p_val"   "cpd_names"          "cpd_forms"          "OC_x"               "HC_y"              
# [15] "NC_z"               "mass"               "class"              "z_layer"    
temp <- temp[ ,c(
  "cpd", "median_t2d", "median_normal", "alt", "p_val", "W_statistic", "hl_effect_wilcox", 
  "trend_with_disease", "minuslog10_p_val", "cpd_names", "cpd_forms",
  "OC_x", "HC_y", "NC_z", "z_layer")]
names(temp) <- c(
  "cpd", "median_log10_cpp_relabun_t2d", "median_log10_cpp_relabun_normal", "alt_wilcox", "p_val_wilcox", "W_statistic", "hl_effect_wilcox", 
  "trend_with_t2d", "minuslog10_p_val_wilcox", "cpd_names", "cpd_forms",
  "OC_x", "HC_y", "NC_z", "z_layer")

dim(temp) #  1269   15


# add functions
temp$fxns <- NA
#temp$superfocus_fxns <- NA
temp$sub_fxns <- NA
temp$subsys_L3 <- NA
temp$no_linked_superfocus_fxns <- NA
temp$no_linked_subsys_L3 <- NA

for (i in 1:dim(temp)[1]) {
  #i<-1
  sel <- which(df.out$cpd_id == temp$cpd[i]) # 2
  ( list_fxns <- unique( df.out$superfocus_fxn[sel] ) ) # "fxn_9950"  "fxn_9951"
  temp$fxns[i] <- paste0(list_fxns, collapse = "  |  ")
  temp$sub_fxns[i] <- paste0( unique(df.out$f__in[sel]), collapse = "  |  ")
  
  sel <- which(row.names(df.tax) %in% list_fxns) # 2
  
  #temp$superfocus_fxns[i] <- paste0( unique(df.tax$fxn[sel]), collapse = "  |  ")
  temp$subsys_L3[i] <- paste0( unique(df.tax$subsys_L3[sel]), collapse = "  |  ")
  
  temp$no_linked_superfocus_fxns[i] <- length( unique(df.tax$fxn[sel]) )
  temp$no_linked_subsys_L3[i] <- length( unique(df.tax$subsys_L3[sel]) )
  
  print(paste0("completed ",i))
}
#temp$superfocus_fxns <- gsub(pattern = ",", replacement = ";", x = temp$superfocus_fxns)
temp$sub_fxns <- gsub(pattern = ",", replacement = ";", x = temp$sub_fxns)

head(temp)


saveRDS(object = temp, file = "temp-Sig-results-T2D-SWE-data-linked-fxns-subsyL3-v6.RDS")

dim(temp) # 1269   20

colnames(phy@tax_table) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn" 
length(unique(phy@tax_table[ , "subsys_L1" ])) #
length(unique(phy@tax_table[ , "subsys_L2" ])) #
length(unique(phy@tax_table[ , "subsys_L3" ])) #
length(unique(phy@tax_table[ , "fxn" ])) # 

write.table(x = temp, file = "Sig-results-table-T2D-SWE-compounds-with-functions-v6.tsv", sep = "\t", quote = FALSE, row.names = FALSE )
write.table(x = temp[ order(temp$p_val_wilcox, decreasing = FALSE) , ], file = "Sig-results-table-T2D-SWE-compounds-with-functions-ordered-v6.tsv", sep = "\t", quote = FALSE, row.names = FALSE )




## T2D-SWE - w/ Host-removal - p <= 0.05 AND trending +/- with soil ecosystem condition - i.e., intersect Restoration and T2D-SWE

#saveRDS(object = temp, file = "temp-Sig-results-T2D-SWE-data-linked-fxns-subsyL3-v6.RDS")
temp <- readRDS("temp-Sig-results-T2D-SWE-data-linked-fxns-subsyL3-v6.RDS")
#saveRDS(object = temp, file = "temp-Sig-results-Restoration-data-linked-fxns-subsyL3-v6.RDS")
temp.resto <- readRDS("temp-Sig-results-Restoration-data-linked-fxns-subsyL3-v6.RDS")

length(resto.sig.cpds) # 2122
identical(temp.resto$cpd, resto.sig.cpds) # TRUE

sel <- which(temp$cpd %in% resto.sig.cpds) # 304
temp <- temp[sel, ]

names(temp)
# [1] "cpd"                             "median_log10_cpp_relabun_t2d"    "median_log10_cpp_relabun_normal" "alt_wilcox"                     
# [5] "p_val_wilcox"                    "W_statistic"                     "hl_effect_wilcox"                "trend_with_t2d"                 
# [9] "minuslog10_p_val_wilcox"         "cpd_names"                       "cpd_forms"                       "OC_x"                           
# [13] "HC_y"                            "NC_z"                            "z_layer"                         "fxns"                           
# [17] "sub_fxns"                        "subsys_L3"                       "no_linked_superfocus_fxns"       "no_linked_subsys_L3" 

plot(temp$no_linked_superfocus_fxns, temp$no_linked_subsys_L3)

# only keep intersecting compounds. Join up data
row.names(temp) <- temp$cpd
row.names(temp.resto) <- temp.resto$cpd
temp.resto <- temp.resto[ row.names(temp) , ]
identical(temp$cpd, temp.resto$cpd) # TRUE

names(temp.resto)
# [1] "cpd"                       "p_val_kendall"             "kendall_tau"               "sigBH"                     "trend_with_age"           
# [6] "minuslog10_p_val_kendall"  "cpd_names"                 "cpd_forms"                 "OC_x"                      "HC_y"                     
# [11] "NC_z"                      "z_layer"                   "fxns"                      "sub_fxns"                  "subsys_L3"                
# [16] "no_linked_superfocus_fxns" "no_linked_subsys_L3"   
names(temp.resto) <- c(
  "cpd"         ,          "p_val_kendall_Restoration"   ,    "kendall_tau_Restoration" ,  "sigBH"   ,  "trend_with_age_Restoration"  ,
  "minuslog10_p_val_kendall_Restoration" , "cpd_names"    ,   "cpd_forms"          ,       "OC_x"       ,           "HC_y"           ,
  "NC_z"        ,    "z_layer"          ,   "fxns"         ,  "sub_fxns_Restoration"    ,   "subsys_L3_Restoration"     ,   
  "no_linked_superfocus_fxns_Restoration", "no_linked_subsys_L3_Restoration"
  )

temp.join <- cbind( temp, temp.resto[ ,c( "p_val_kendall_Restoration", "kendall_tau_Restoration", "trend_with_age_Restoration",    "minuslog10_p_val_kendall_Restoration" ,
                                          "sub_fxns_Restoration", "subsys_L3_Restoration", "no_linked_superfocus_fxns_Restoration", "no_linked_subsys_L3_Restoration"  )])
dim(temp.join)
# 304 28


saveRDS(object = temp.join, file = "temp.join-Sig-results-T2D-SWE-Hostremoval-AND-Restoration-data-linked-fxns-subsyL3-v6.RDS")


write.table(x = temp.join, file = "Sig-results-table-T2D-SWE-Hostremoval-AND-Restoration-compounds-with-functions-v6.tsv", sep = "\t", quote = FALSE, row.names = FALSE )
write.table(x = temp.join[ order(temp.join$p_val_wilcox, decreasing = FALSE) , ], file = "Sig-results-table-T2D-SWE-Hostremoval-AND-Restoration-compounds-with-functions-ordered-v6.tsv", sep = "\t", quote = FALSE, row.names = FALSE )





## T2D-CHN

dat.test.sig <- readRDS("dat.test.sig-cpp3d-indiv-Compounds-Forslund-CHN-T2D-Hostremoval-All-Compounds-T2D-VS-NORM--Wilcox-v6.rds")
str(dat.test.sig)

# link compounds to functions
df.out <-readRDS( "df.out--tidy-compounds_indiv-cpp3d--forslund-t2d-chn-hostremoval-pawsey.RDS" )
str(df.out)
# 'data.frame':	554562 obs. of  9 variables:
# $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_3" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "4-hydroxyproline epimerase (EC 5.1.1.8)" "4-hydroxyproline epimerase (EC 5.1.1.8)" "Alanine racemase (EC 5.1.1.1)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id             : chr  "rxn02360" "rxn02360" "rxn00283" "rxn00283" ...
# $ cpd_id             : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ cpd_name           : chr  "trans-4-Hydroxy-L-proline" "cis-4-Hydroxy-D-proline" "L-Alanine" "D-Alanine" ...
# $ cpd_form           : chr  "C5H9NO3" "C5H9NO3" "C3H7NO2" "C3H7NO2" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.167 0.167 0.167 ...

phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal.RDS")
df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 19352     4

table(phy@sam_data$Diagnosis)
# ND CTRL T2D metformin- 
#   52             30 

min(taxa_sums(phy)) # 1.758995e-06

## tidy sig test results table and add function info
temp <- dat.test.sig
head(temp)
names(temp)
# [1] "cpd"                "data_for_this_cpd"  "median_t2d"         "median_normal"      "alt"                "p_val"              "W_statistic"       
# [8] "hl_effect_wilcox"   "trend_with_disease" "minuslog10_p_val"   "cpd_names"          "cpd_forms"          "OC_x"               "HC_y"              
# [15] "NC_z"               "z_layer"  
temp <- temp[ ,c(
  "cpd", "median_t2d", "median_normal", "alt", "p_val", "W_statistic", "hl_effect_wilcox", 
  "trend_with_disease", "minuslog10_p_val", "cpd_names", "cpd_forms",
  "OC_x", "HC_y", "NC_z", "z_layer")]
names(temp) <- c(
  "cpd", "median_log10_cpp_relabun_t2d", "median_log10_cpp_relabun_normal", "alt_wilcox", "p_val_wilcox", "W_statistic", "hl_effect_wilcox", 
  "trend_with_t2d", "minuslog10_p_val_wilcox", "cpd_names", "cpd_forms",
  "OC_x", "HC_y", "NC_z", "z_layer")

dim(temp) #  4292   15


# add functions
temp$fxns <- NA
#temp$superfocus_fxns <- NA
temp$sub_fxns <- NA
temp$subsys_L3 <- NA
temp$no_linked_superfocus_fxns <- NA
temp$no_linked_subsys_L3 <- NA

for (i in 1:dim(temp)[1]) {
  #i<-1
  sel <- which(df.out$cpd_id == temp$cpd[i]) # 2
  ( list_fxns <- unique( df.out$superfocus_fxn[sel] ) ) # "fxn_9950"  "fxn_9951"
  temp$fxns[i] <- paste0(list_fxns, collapse = "  |  ")
  temp$sub_fxns[i] <- paste0( unique(df.out$f__in[sel]), collapse = "  |  ")
  
  sel <- which(row.names(df.tax) %in% list_fxns) # 2
  
  #temp$superfocus_fxns[i] <- paste0( unique(df.tax$fxn[sel]), collapse = "  |  ")
  temp$subsys_L3[i] <- paste0( unique(df.tax$subsys_L3[sel]), collapse = "  |  ")
  
  temp$no_linked_superfocus_fxns[i] <- length( unique(df.tax$fxn[sel]) )
  temp$no_linked_subsys_L3[i] <- length( unique(df.tax$subsys_L3[sel]) )
  
  print(paste0("completed ",i))
}
#temp$superfocus_fxns <- gsub(pattern = ",", replacement = ";", x = temp$superfocus_fxns)
temp$sub_fxns <- gsub(pattern = ",", replacement = ";", x = temp$sub_fxns)

head(temp)


saveRDS(object = temp, file = "temp-Sig-results-T2D-CHN-data-linked-fxns-subsyL3-v6.RDS")

dim(temp) # 4292  20

colnames(phy@tax_table) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn" 
length(unique(phy@tax_table[ , "subsys_L1" ])) #
length(unique(phy@tax_table[ , "subsys_L2" ])) #
length(unique(phy@tax_table[ , "subsys_L3" ])) #
length(unique(phy@tax_table[ , "fxn" ])) #

write.table(x = temp, file = "Sig-results-table-T2D-CHN-compounds-with-functions-v6.tsv", sep = "\t", quote = FALSE, row.names = FALSE )
write.table(x = temp[ order(temp$p_val_wilcox, decreasing = FALSE) , ], file = "Sig-results-table-T2D-CHN-compounds-with-functions-ordered-v6.tsv", sep = "\t", quote = FALSE, row.names = FALSE )




## Overlap Soil AND T2D-SWE AND T2D-CHN - w/ Host-removal
rm(temp)

temp.join <- readRDS("temp.join-Sig-results-T2D-SWE-Hostremoval-AND-Restoration-data-linked-fxns-subsyL3-v6.RDS")
temp.chn <- readRDS("temp-Sig-results-T2D-CHN-data-linked-fxns-subsyL3-v6.RDS")

names(temp.join)
# [1] "cpd"                                   "median_log10_cpp_relabun_t2d"          "median_log10_cpp_relabun_normal"       "alt_wilcox"                           
# [5] "p_val_wilcox"                          "W_statistic"                           "hl_effect_wilcox"                      "trend_with_t2d"                       
# [9] "minuslog10_p_val_wilcox"               "cpd_names"                             "cpd_forms"                             "OC_x"                                 
# [13] "HC_y"                                  "NC_z"                                  "z_layer"                               "fxns"                                 
# [17] "sub_fxns"                              "subsys_L3"                             "no_linked_superfocus_fxns"             "no_linked_subsys_L3"                  
# [21] "p_val_kendall_Restoration"             "kendall_tau_Restoration"               "trend_with_age_Restoration"            "minuslog10_p_val_kendall_Restoration" 
# [25] "sub_fxns_Restoration"                  "subsys_L3_Restoration"                 "no_linked_superfocus_fxns_Restoration" "no_linked_subsys_L3_Restoration"

names(temp.join) <- c(
  "cpd" , "median_log10_cpp_relabun_t2d_SWE", "median_log10_cpp_relabun_normal_SWE", "alt_wilcox_SWE",
  "p_val_wilcox_SWE", "W_statistic_SWE", "hl_effect_wilcox_SWE", "trend_with_t2d_SWE",
  "minuslog10_p_val_wilcox_SWE", 
  "cpd_names", "cpd_forms",  "OC_x", "HC_y", "NC_z",  "z_layer", 
  "fxns", "sub_fxns_SWE", "subsys_L3_SWE",  "no_linked_superfocus_fxns_SWE", "no_linked_subsys_L3_SWE",
  "p_val_kendall_Restoration", "kendall_tau_Restoration", "trend_with_age_Restoration", "minuslog10_p_val_kendall_Restoration",
  "sub_fxns_Restoration", "subsys_L3_Restoration", "no_linked_superfocus_fxns_Restoration", "no_linked_subsys_L3_Restoration" 
)

sel <- which(names(temp.join) == "fxns") # 16
temp.join <- temp.join[ , -sel]

identical(temp.join$cpd, row.names(temp.join)) # TRUE
class(temp.join$cpd) # character
class(temp.chn$cpd) # character

row.names(temp.join) # cpd ids
row.names(temp.chn) # numbers
length(unique(temp.chn$cpd)) # 4292
row.names(temp.chn) <- temp.chn$cpd

sel1 <- which(temp.join$cpd %in% temp.chn$cpd) # 251
temp.join <- temp.join[ sel1, ]
row.names(temp.join)

# now cut back temp.chn
temp.chn <- temp.chn[ row.names(temp.join), ]

identical(temp.join$cpd, temp.chn$cpd) # TRUE

names(temp.chn)
# [1] "cpd"                             "median_log10_cpp_relabun_t2d"    "median_log10_cpp_relabun_normal" "alt_wilcox"                     
# [5] "p_val_wilcox"                    "W_statistic"                     "hl_effect_wilcox"                "trend_with_t2d"                 
# [9] "minuslog10_p_val_wilcox"         "cpd_names"                       "cpd_forms"                       "OC_x"                           
# [13] "HC_y"                            "NC_z"                            "z_layer"                         "fxns"                           
# [17] "sub_fxns"                        "subsys_L3"                       "no_linked_superfocus_fxns"       "no_linked_subsys_L3"  

# rename
names(temp.chn) <- c(
  "cpd",   "median_log10_cpp_relabun_t2d_CHN", "median_log10_cpp_relabun_normal_CHN", "alt_wilcox_CHN" ,
  "p_val_wilcox_CHN",  "W_statistic_CHN",  "hl_effect_wilcox_CHN",  "trend_with_t2d_CHN",
  "minuslog10_p_val_wilcox_CHN",  "cpd_names", "cpd_forms",  "OC_x",
  "HC_y" , "NC_z", "z_layer",  "fxns",
  "sub_fxns_CHN" , "subsys_L3_CHN", "no_linked_superfocus_fxns_CHN", "no_linked_subsys_L3_CHN"
)

# merge tables
identical(temp.join$cpd, temp.chn$cpd) # TRUE
row.names(temp.join)
row.names(temp.chn)


temp.match <- cbind(temp.join, temp.chn[ ,c(
  "median_log10_cpp_relabun_t2d_CHN", "median_log10_cpp_relabun_normal_CHN", "alt_wilcox_CHN" ,
  "p_val_wilcox_CHN",  "W_statistic_CHN",  "hl_effect_wilcox_CHN",  "trend_with_t2d_CHN",
  "minuslog10_p_val_wilcox_CHN",
  "sub_fxns_CHN" , "subsys_L3_CHN", "no_linked_superfocus_fxns_CHN", "no_linked_subsys_L3_CHN" )])

names(temp.match)
# [1] "cpd"                                   "median_log10_cpp_relabun_t2d_SWE"      "median_log10_cpp_relabun_normal_SWE"   "alt_wilcox_SWE"                       
# [5] "p_val_wilcox_SWE"                      "W_statistic_SWE"                       "hl_effect_wilcox_SWE"                  "trend_with_t2d_SWE"                   
# [9] "minuslog10_p_val_wilcox_SWE"           "cpd_names"                             "cpd_forms"                             "OC_x"                                 
# [13] "HC_y"                                  "NC_z"                                  "z_layer"                               "sub_fxns_SWE"                         
# [17] "subsys_L3_SWE"                         "no_linked_superfocus_fxns_SWE"         "no_linked_subsys_L3_SWE"               "p_val_kendall_Restoration"            
# [21] "kendall_tau_Restoration"               "trend_with_age_Restoration"            "minuslog10_p_val_kendall_Restoration"  "sub_fxns_Restoration"                 
# [25] "subsys_L3_Restoration"                 "no_linked_superfocus_fxns_Restoration" "no_linked_subsys_L3_Restoration"       "median_log10_cpp_relabun_t2d_CHN"     
# [29] "median_log10_cpp_relabun_normal_CHN"   "alt_wilcox_CHN"                        "p_val_wilcox_CHN"                      "W_statistic_CHN"                      
# [33] "hl_effect_wilcox_CHN"                  "trend_with_t2d_CHN"                    "minuslog10_p_val_wilcox_CHN"           "sub_fxns_CHN"                         
# [37] "subsys_L3_CHN"                         "no_linked_superfocus_fxns_CHN"         "no_linked_subsys_L3_CHN"     

temp.match[ ,c("trend_with_t2d_SWE", "trend_with_t2d_CHN" )]

sel <- which(temp.match$trend_with_t2d_CHN == temp.match$trend_with_t2d_SWE) # 155
temp.match <- temp.match[sel, ]

temp.match$trend_with_eco_degradation <- NA
sel <- which(temp.match$trend_with_age_Restoration == "Increasing") # 73
temp.match$trend_with_eco_degradation[sel] <- "Decreasing"
sel <- which(temp.match$trend_with_age_Restoration == "Decreasing") # 82
temp.match$trend_with_eco_degradation[sel] <- "Increasing"

sel <- which(temp.match$trend_with_t2d_SWE == temp.match$trend_with_eco_degradation) # 76 (qty is same w/ Host-removal vs no host-removal)

temp.match.degrad <- temp.match[sel, ]

table(temp.match.degrad$trend_with_eco_degradation)
# Decreasing Increasing 
# 60         16 
table(temp.match.degrad$trend_with_t2d_SWE)
# Decreasing Increasing 
# 60         16 
table(temp.match.degrad$trend_with_t2d_CHN)
# Decreasing Increasing 
# 60         16 


temp.match.degrad[ , c("cpd_names", "trend_with_t2d_SWE", "trend_with_t2d_CHN", "trend_with_eco_degradation", "no_linked_subsys_L3_Restoration", "no_linked_subsys_L3_SWE", "no_linked_subsys_L3_CHN"  )]

saveRDS(temp.match, file = "temp.match-Sig-results-Restoration-T2D-SWE-CHN-Hostremoval-data-linked-fxns-subsyL3-v6.RDS")


temp.match.degrad[ order(temp.match.degrad$kendall_tau_Restoration) , c("cpd_names", "trend_with_t2d_SWE", "trend_with_t2d_CHN", "trend_with_eco_degradation", "no_linked_subsys_L3_Restoration", "no_linked_subsys_L3_SWE", "no_linked_subsys_L3_CHN", "kendall_tau_Restoration"  )]

plot( temp.match.degrad$hl_effect_wilcox_SWE, temp.match.degrad$hl_effect_wilcox_CHN )
plot( temp.match.degrad$minuslog10_p_val_wilcox_SWE, temp.match.degrad$minuslog10_p_val_wilcox_CHN )

temp.match.degrad$minuslog10_p_val_wilcox_T2D_mean <- rowMeans( cbind( temp.match.degrad$minuslog10_p_val_wilcox_SWE, temp.match.degrad$minuslog10_p_val_wilcox_CHN  ) )
temp.match.degrad$no_linked_subsys_L3_T2D_mean <- rowMeans( cbind( temp.match.degrad$no_linked_subsys_L3_SWE, temp.match.degrad$no_linked_subsys_L3_CHN  ) )

saveRDS(temp.match.degrad, file = "temp.match.degrad-Sig-results-EcoDegradation-consistent-T2D-SWE-CHN-Hostremoval-data-linked-fxns-subsyL3-v6.RDS")

write.table(x = temp.match.degrad[ order(temp.match.degrad$kendall_tau_Restoration) , ], file = "Sig-Matching-EcoDegradation-T2D-SWE-CHN-Hostremoval-compounds-with-functions-ordered-by-kendall-tau-restoration-v6.tsv", sep = "\t", quote = FALSE, row.names = TRUE )

write.table( temp.match.degrad[ order(temp.match.degrad$kendall_tau_Restoration) , c("cpd_names", "trend_with_eco_degradation", "trend_with_t2d_SWE", "trend_with_t2d_CHN", 
                                                                                     "no_linked_subsys_L3_Restoration", "no_linked_subsys_L3_SWE", "no_linked_subsys_L3_CHN", "no_linked_subsys_L3_T2D_mean",
                                                                                     "kendall_tau_Restoration", "minuslog10_p_val_wilcox_SWE", "minuslog10_p_val_wilcox_CHN", "minuslog10_p_val_wilcox_T2D_mean"
                                                                                      )] , file = "Sig-Matching-SUMMARY-EcoDegradation-T2D-SWE-CHN-Hostremoval-compounds-ordered-by-kendall-tau-restoration-v6.tsv", sep = "\t", quote = FALSE, row.names = TRUE )


pdat <- temp.match.degrad[  , c("cpd", "cpd_names", "trend_with_eco_degradation", "trend_with_t2d_SWE", "trend_with_t2d_CHN", 
                                "no_linked_subsys_L3_Restoration", "no_linked_subsys_L3_SWE", "no_linked_subsys_L3_CHN", "no_linked_subsys_L3_T2D_mean",
                                "kendall_tau_Restoration", "minuslog10_p_val_wilcox_SWE", "minuslog10_p_val_wilcox_CHN", "minuslog10_p_val_wilcox_T2D_mean")]

pdat$kendall_tau_Restoration

length(which( t2d.zoom1.decreasing.proteins %in% temp.match.degrad$cpd)) # 35
sel <- which( pdat$cpd %in% t2d.zoom1.decreasing.proteins )
pdat[sel, c(2, 9, 10, 13)]
sel <- which( pdat$cpd %in% t2d.zoom2B.increasing.carbs ) # "cpd00224" "cpd03198" "cpd03200" "cpd00259" "cpd00076" "cpd00082"
pdat$cpd_names[sel] # "L-Arabinose" "Melibiose"   "D-Fructose" 
sel <- which( pdat$cpd_names %in% c("L-Arabinose", "Melibiose", "D-Fructose", "6-Phosphosucrose") )
pdat[sel, c(2, 9, 10, 13)]
#                 cpd_names no_linked_subsys_L3_T2D_mean kendall_tau_Restoration minuslog10_p_val_wilcox_T2D_mean
# cpd00224      L-Arabinose                         46.5               -4.152231                         2.496697
# cpd03198        Melibiose                         17.5               -2.734396                         2.578962
# cpd00082       D-Fructose                         42.0               -3.747136                         2.208783
# cpd01693 6-Phosphosucrose                         27.5               -2.734396                         2.184635



p <- ggplot(pdat, aes(x = kendall_tau_Restoration, y = minuslog10_p_val_wilcox_T2D_mean, size = no_linked_subsys_L3_T2D_mean, color = kendall_tau_Restoration))+
  #scale_x_reverse()+
  geom_point()+
  scale_color_viridis(option='turbo')+
  xlab("Kendall's tau correlation with reveg age")+
  #ylab("Mean -log10(P-value) from T2D (SWE, CHN) Wilcox tests")+
  ylab("Mean -log10(P-value) from T2D Wilcoxon tests")+
  
  #annotate(geom="text", x = -4.15, y = 2.57, label = "L-Arabinose", size = 3.25 , hjust = 0, vjust =0 )+
  annotate(geom="text", x = -3.85, y = 2.5, label = "L-Arabinose", size = 3.25 , hjust = 0, vjust =0 )+
  #
  
  #annotate(geom="text", x = -3.75, y = 2.55, label = "D-Fructose", size = 3.25 , hjust = 0, vjust =0 )+
  annotate(geom="text", x = -3.69, y = 2.115, label = "D-Fructose", size = 3.25 , hjust = 0, vjust =1 )+
  annotate("segment",x = -3.75, y = 2.2087, xend = -3.69, yend = 2.12, color = "#252525" )+
  
  #annotate(geom="text", x = -2.73, y = 2.55, label = "6-Phosphosucrose", size = 3.25 , hjust = 0, vjust =0 )+
  annotate(geom="text", x = -2.50, y = 2.175, label = "6-Phosphosucrose", size = 3.25 , hjust = 0, vjust =0 )+
  
  #annotate(geom="text", x = -2.73, y = 2.92, label = "Melibiose", size = 3.25 , hjust = 0, vjust =0 )+
  annotate(geom="text", x = -2.60, y = 2.60, label = "Melibiose", size = 3.25 , hjust = 0, vjust =0 )+
  
  #annotate(geom="text", x = 3.54, y = 3.31, label = "BCFA-ACPs (34 of 35)", size = 3.25 , hjust = 1, vjust =1 )+
  annotate(geom="text", x = 3.53, y = 3.265, label = "BCFA-ACPs (34 of 35)", size = 3.25 , hjust = 1, vjust =1 )+
  
  #annotate(geom="text", x = 3.14, y = 1.99, label = "BCFA-ACPs (1 of 35)", size = 3.25 , hjust = 1, vjust =1 )+
  annotate(geom="text", x = 3.0, y = 1.95, label = "BCFA-ACPs (1 of 35)", size = 3.25 , hjust = 1, vjust =1 )+
  annotate("segment",x = 2.95, y = 1.95, xend = 3.14, yend = 1.978, color = "#252525" )+
  
  #annotate(geom="text", x = 3.95, y = 1.56, label = "Lignin", size = 3.25 , hjust = 1, vjust =1 )+
  annotate(geom="text", x = 3.65, y = 1.57, label = "Lignin", size = 3.25 , hjust = 1, vjust =0 )+
  annotate("segment",x = 3.65, y = 1.57, xend = 3.95, yend = 1.56, color = "#252525" )+
  
  annotate("segment",x = 0, y = 2.43, xend = 4, yend = 2.43, color = "red", arrow=arrow(length = unit(5,"pt")) )+
  annotate(geom="text", x = 2, y = 2.4, label = "Decreased in T2D &\ndegraded ecosystems", size = 3.25 , hjust = 0.5, vjust =1 , color = "red", lineheight = 0.8)+
  
  annotate("segment",x = 0, y = 2.28, xend = -4, yend = 2.28, color = "blue", arrow=arrow(length = unit(5,"pt")) )+
  annotate(geom="text", x = -2, y = 2.4, label = "Increased in T2D &\ndegraded ecosystems", size = 3.25 , hjust = 0.5, vjust =1 , color = "blue", lineheight = 0.8)+
  
  #annotate("segment",x = -4.15, y = 1.45, xend = -3.15, yend = 1.45, color = "#252525", linetype = "dotted" )+
  #annotate(geom="text", x = -4.15, y = 1.43, label = "*", size = 5 , hjust = 0, vjust =0, color = "#252525")+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = "*", size = 5 , hjust = 0, vjust = 1, color = "#252525")+
  
  annotate("segment",x = -4.15, y = 2, xend = -3.15, yend = 2, color = "#252525", linetype = "dotted" )+
  annotate(geom="text", x = -4.15, y = 2, label = "**", size = 5 , hjust = 0, vjust =0, color = "#252525")+
  
  annotate("segment",x = -4.15, y = 3, xend = -3.15, yend = 3, color = "#252525", linetype = "dotted" )+
  annotate(geom="text", x = -4.15, y = 3, label = "***", size = 5 , hjust = 0, vjust =0, color = "#252525" )+
  
  
  guides(color = guide_colorbar(title = "Correlation of\nCPP% with\nrestoration"),
         size = guide_legend(title = "Mean no. linked\nsubsystems (L3)\nin T2D"))+
  theme_classic()

p

dev.print(tiff, file = paste0(workdir,"/plots/","Matching-trends-EcoDegradation-T2D-SWE-CHN-Hostremoval-KendallTauResto-vs-MeanMinuslog10PWilcox-v6.tiff"), width = 14, height = 14, units = "cm", res=600, compression="lzw",type="cairo")


table(temp.match.degrad$trend_with_eco_degradation)
# Decreasing Increasing 
# 60         16 

sel <- which(temp.match.degrad$cpd %in% t2d.zoom1.decreasing.proteins) # qty 35 !
sel <- which(temp.match.degrad$cpd %in% t2d.zoom2.increasing.carbs) # qty 3
temp.match.degrad$cpd_names [sel]
# "L-Arabinose" "Melibiose"   "D-Fructose" 

#-------------------------

#### Track compounds / Enzymes (sub-functions) / Explore SUPER-FOCUS functions associated with BCFA-ACPs?
#-------------------------

sort(t2d.zoom1.decreasing.proteins)
# [1] "cpd11465" "cpd11469" "cpd11473" "cpd11475" "cpd11498" "cpd11499" "cpd11502" "cpd11503" "cpd11506" "cpd11507" "cpd11510" "cpd11511"
# [13] "cpd11514" "cpd11518" "cpd11523" "cpd11524" "cpd11527" "cpd11528" "cpd11531" "cpd11532" "cpd11535" "cpd11536" "cpd11539" "cpd11543"
# [25] "cpd11548" "cpd11549" "cpd11552" "cpd11553" "cpd11556" "cpd11557" "cpd11560" "cpd11561" "cpd11564" "cpd11568" "cpd11572"

blcfa <- t2d.zoom1.decreasing.proteins

## Sunbad resto

# link compounds to functions
df.out <-readRDS( "df.out--tidy-compounds_indiv--cpp3d-sunbad-resto.RDS" )

str(df.out)
# 'data.frame':	1154547 obs. of  9 variables:
# $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_3" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylaconitate isomerase" "2-methylaconitate isomerase" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" ...
# $ rxn_id             : chr  "rxn25278" "rxn25278" "rxn25279" "rxn25279" ...
# $ cpd_id             : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ cpd_name           : chr  "2-methyl-trans-aconitate" "cis-2-Methylaconitate" "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" ...
# $ cpd_form           : chr  "C7H5O6" "C7H5O6" "C7H7O7" "H2O" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.333 0.333 0.333 ...

phy <- readRDS("phy-phyloseq-object-sunbad-resto.RDS")
df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 30125     4
min(taxa_sums(phy)) # 2.521458e-06

length(blcfa) # 35

sel <- which(df.out$cpd_id %in% blcfa) # 131

df.out[sel, ]

# blcfa names
unique(df.out$cpd_name[sel])
# [1] "But-2-enoyl-[acyl-carrier protein]"    "(2E)-Dodecenoyl-[acp]"                 "(2E)-Hexenoyl-[acp]"                   "(2E)-Decenoyl-[acp]"                  
# [5] "4-methyl-trans-hex-2-enoyl-ACP"        "4-methyl-hexanoyl-ACP"                 "6-methyl-trans-oct-2-enoyl-ACP"        "6-methyl-octanoyl-ACP"                
# [9] "8-methyl-trans-dec-2-enoyl-ACP"        "8-methyl-decanoyl-ACP"                 "10-methyl-trans-dodec-2-enoyl-ACP"     "10-methyl-dodecanoyl-ACP"             
# [13] "12-methyl-trans-tetra-dec-2-enoyl-ACP" "14-methyl-trans-hexa-dec-2-enoyl-ACP"  "5-methyl-trans-hex-2-enoyl-ACP"        "5-methyl-hexanoyl-ACP"                
# [17] "7-methyl-trans-oct-2-enoyl-ACP"        "7-methyl-octanoyl-ACP"                 "9-methyl-trans-dec-2-enoyl-ACP"        "9-methyl-decanoyl-ACP"                
# [21] "11-methyl-trans-dodec-2-enoyl-ACP"     "11-methyl-dodecanoyl-ACP"              "13-methyl-trans-tetra-dec-2-enoyl-ACP" "15-methyl-trans-hexa-dec-2-enoyl-ACP" 
# [25] "4-methyl-trans-pent-2-enoyl-ACP"       "4-methyl-pentanoyl-ACP"                "6-methyl-trans-hept-2-enoyl-ACP"       "6-methyl-heptanoyl-ACP"               
# [29] "8-methyl-trans-non-2-enoyl-ACP"        "8-methyl-nonanoyl-ACP"                 "10-methyl-trans-undec-2-enoyl-ACP"     "10-methyl-undecanoyl-ACP"             
# [33] "12-methyl-trans-tridec-2-enoyl-ACP"    "14-methyl-trans-pentadec-2-enoyl-ACP"  "trans-Octodec-2-enoyl-ACP" 


# no of unique reactions
length( unique(df.out$rxn_id[sel]) ) # 59
# no of unique sub-functions 
length( unique(df.out$f__in[sel])) # 5
unique(df.out$f__in[sel])
# [1] "Enoyl-.acyl-carrier-protein. reductase .NADH. .EC 1.3.1.9."      "Trans-2-decenoyl-[acyl-carrier-protein] isomerase (EC 5.3.3.14)"
# [3] "3-oxoacyl-.acyl-carrier-protein. synthase, KASII .EC 2.3.1.41."  "3-oxoacyl-.acyl-carrier-protein. synthase, KASIII .EC 2.3.1.41."
# [5] "Enoyl-.acyl-carrier-protein. reductase .NADPH. .EC 1.3.1.10."
# no of unique superfocus functions
length( unique(df.out$superfocus_fxn[sel])) # 7

(sf_fxn <- sort(unique(df.out$superfocus_fxn[sel])))
# "fxn_15794" "fxn_15811" "fxn_15813" "fxn_15859" "fxn_15860" "fxn_16195" "fxn_8733" 

sel.all <- which(df.out$superfocus_fxn %in% sf_fxn)
length(sel.all) # 722
length(unique(df.out$rxn_id[sel.all])) # 89 total reactions associated with those functions

( cpd_all <- unique(df.out$cpd_id[sel.all]) )
length(cpd_all) # 110

sel.other <- which(!cpd_all %in% blcfa) # 75 other compounds are in the same set of reactions
cpd_other <- cpd_all[sel.other]
length(cpd_other) # 75
length(unique(cpd_other)) # 75

sel.wider <- which(df.out$cpd_id %in% cpd_other)
length(sel.wider) # 357050
( other_cpd_names <- sort( unique(df.out$cpd_name[sel.wider]) ) )
# [1] "(2E)-Hexadecenoyl-[acp]"                 "(2E)-Octadecenoyl-[acp]"                 "(2E)-Octenoyl-[acp]"                    
# [4] "(2E)-Tetradecenoyl-[acp]"                "10-methyl-3-oxo-dodecanoyl-ACP"          "10-methyl-3-oxo-undecanoyl-ACP"         
# [7] "11-methyl-3-oxo-dodecanoyl-ACP"          "12-methyl-3-oxo-tetra-decanoyl-ACP"      "12-methyl-3-oxo-tridecanoyl-ACP"        
# [10] "12-methyl-tetra-decanoyl-ACP"            "12-methyl-tridecanoyl-ACP"               "13-methyl-3-oxo-tetra-decanoyl-ACP"     
# [13] "13-methyl-tetra-decanoyl-ACP"            "14-methyl-3-oxo-hexa-decanoyl-ACP"       "14-methyl-3-oxo-pentadecanoyl-ACP"      
# [16] "14-methyl-hexa-decanoyl-ACP"             "14-methyl-pentadecanoyl-ACP"             "15-methyl-3-oxo-hexa-decanoyl-ACP"      
# [19] "15-methyl-hexa-decanoyl-ACP"             "2-methylbutyryl-ACP"                     "3-hydroxy-cis-D7-tetraecenoyl-ACPs"     
# [22] "3-hydroxy-cis-D9-hexaecenoyl-ACPs"       "3-Hydroxy-octanoyl-ACPs"                 "3-Hydroxyglutaryl-ACP-methyl-ester"     
# [25] "3-hydroxypimeloyl-ACP-methyl-esters"     "3-oxodecanoyl-acp"                       "3-oxododecanoyl-acp"                    
# [28] "3-oxohexadecanoyl-acp"                   "3-Oxohexanoyl-[acp]"                     "3-oxooctanoyl-acp"                      
# [31] "3-Oxooctodecanoyl-ACP"                   "3-oxotetradecanoyl-acp"                  "4-methyl-3-oxo-hexanoyl-ACP"            
# [34] "4-methyl-3-oxo-pentanoyl-ACP"            "5-methyl-3-oxo-hexanoyl-ACP"             "6-methyl-3-oxo-heptanoyl-ACP"           
# [37] "6-methyl-3-oxo-octanoyl-ACP"             "7-methyl-3-oxo-octanoyl-ACP"             "8-methyl-3-oxo-decanoyl-ACP"            
# [40] "8-methyl-3-oxo-nonanoyl-ACP"             "9-methyl-3-oxo-decanoyl-ACP"             "a 3-hydroxy cis delta5-dodecenoyl-[acp]"
# [43] "Acetoacetyl-ACP"                         "Acetyl-ACP"                              "Acetyl-CoA"                             
# [46] "ACP"                                     "Butyryl-ACP"                             "cis-3-Decenoyl-[acyl-carrier protein]"  
# [49] "Cis-delta-3-decenoyl-ACPs"               "CO2"                                     "CoA"                                    
# [52] "Decanoyl-ACP"                            "Dodecanoyl-ACP"                          "Enoylglutaryl-ACP-methyl-esters"        
# [55] "Enoylpimeloyl-ACP-methyl-esters"         "H+"                                      "H2O"                                    
# [58] "hexadecanoyl-acp"                        "Hexanoyl-ACP"                            "isobutyryl-ACP"                         
# [61] "isovaleryl-ACP"                          "Malonyl-acyl-carrierprotein-"            "Myristoyl-ACP"                          
# [64] "NAD"                                     "NADH"                                    "NADP"                                   
# [67] "NADPH"                                   "Octanoyl-ACP"                            "Octodecanoyl-ACP"                       
# [70] "R-3-hydroxystearoyl-ACPs"                "Stearoyl-[acyl-carrier protein]"         "Trans-D2-decenoyl-ACPs"                 
# [73] "Trans-D3-cis-D5-dodecenoyl-ACPs"         "Trans-D3-cis-D7-tetradecenoyl-ACPs"      "Trans-D3-cis-D9-hexadecenoyl-ACPs" 

# only keep other ACPs

other_ACP_names <- c(
  "(2E)-Hexadecenoyl-[acp]"          ,       "(2E)-Octadecenoyl-[acp]"           ,      "(2E)-Octenoyl-[acp]"                  ,  
  "(2E)-Tetradecenoyl-[acp]"        ,        "10-methyl-3-oxo-dodecanoyl-ACP"   ,       "10-methyl-3-oxo-undecanoyl-ACP"      ,   
  "11-methyl-3-oxo-dodecanoyl-ACP"  ,        "12-methyl-3-oxo-tetra-decanoyl-ACP" ,     "12-methyl-3-oxo-tridecanoyl-ACP"     ,   
  "12-methyl-tetra-decanoyl-ACP"    ,        "12-methyl-tridecanoyl-ACP"          ,     "13-methyl-3-oxo-tetra-decanoyl-ACP"  ,   
  "13-methyl-tetra-decanoyl-ACP"    ,        "14-methyl-3-oxo-hexa-decanoyl-ACP"  ,     "14-methyl-3-oxo-pentadecanoyl-ACP"   ,   
  "14-methyl-hexa-decanoyl-ACP"     ,        "14-methyl-pentadecanoyl-ACP"        ,     "15-methyl-3-oxo-hexa-decanoyl-ACP"   ,   
  "15-methyl-hexa-decanoyl-ACP"     ,        "2-methylbutyryl-ACP"                ,     "3-hydroxy-cis-D7-tetraecenoyl-ACPs"  ,   
  "3-hydroxy-cis-D9-hexaecenoyl-ACPs",       "3-Hydroxy-octanoyl-ACPs"            ,     "3-Hydroxyglutaryl-ACP-methyl-ester"  ,   
  "3-hydroxypimeloyl-ACP-methyl-esters",     "3-oxodecanoyl-acp"                  ,     "3-oxododecanoyl-acp"                 ,   
  "3-oxohexadecanoyl-acp"             ,      "3-Oxohexanoyl-[acp]"                ,     "3-oxooctanoyl-acp"                   ,   
  "3-Oxooctodecanoyl-ACP"             ,      "3-oxotetradecanoyl-acp"             ,     "4-methyl-3-oxo-hexanoyl-ACP"         ,   
  "4-methyl-3-oxo-pentanoyl-ACP"      ,      "5-methyl-3-oxo-hexanoyl-ACP"        ,     "6-methyl-3-oxo-heptanoyl-ACP"        ,   
  "6-methyl-3-oxo-octanoyl-ACP"       ,      "7-methyl-3-oxo-octanoyl-ACP"         ,    "8-methyl-3-oxo-decanoyl-ACP"         ,   
  "8-methyl-3-oxo-nonanoyl-ACP"       ,      "9-methyl-3-oxo-decanoyl-ACP"        ,     "a 3-hydroxy cis delta5-dodecenoyl-[acp]",
  "Acetoacetyl-ACP"                   ,      "Acetyl-ACP"                         ,     
  
  # "Acetyl-CoA"                             ,
  # "ACP"                               ,      
  
  "Butyryl-ACP"                        ,     "cis-3-Decenoyl-[acyl-carrier protein]"  ,
  "Cis-delta-3-decenoyl-ACPs"          ,    
  
  # "CO2"                                ,     "CoA"                                   , 
  
  "Decanoyl-ACP"                      ,      "Dodecanoyl-ACP"                     ,     "Enoylglutaryl-ACP-methyl-esters"       , 
  "Enoylpimeloyl-ACP-methyl-esters"   ,    
  
  # "H+"                                 ,     "H2O"                                   , 
  
  "hexadecanoyl-acp"                  ,      "Hexanoyl-ACP"                       ,     "isobutyryl-ACP"                        , 
  "isovaleryl-ACP"                    ,      "Malonyl-acyl-carrierprotein-"       ,     "Myristoyl-ACP"                         , 
  
  # "NAD"                               ,      "NADH"                               ,     "NADP"                                  , 
  # "NADPH"                              ,   
  
  "Octanoyl-ACP"                       ,     "Octodecanoyl-ACP"                      , 
  "R-3-hydroxystearoyl-ACPs"           ,     "Stearoyl-[acyl-carrier protein]"    ,     "Trans-D2-decenoyl-ACPs"                , 
  "Trans-D3-cis-D5-dodecenoyl-ACPs"    ,     "Trans-D3-cis-D7-tetradecenoyl-ACPs" ,     "Trans-D3-cis-D9-hexadecenoyl-ACPs" 
  )

length(other_ACP_names) # 65


# how many wider functions do the other ACP name compounds occur in?
# not including the initial 7 reactions linked to BLCFA?

sel <- which(df.out$cpd_name %in% other_ACP_names ) # 3890
sel <- which(df.out$cpd_name %in% other_ACP_names & !df.out$superfocus_fxn %in% sf_fxn ) # 3667

sf_fxn_other_ACP <- unique(df.out$superfocus_fxn[sel]) # 266
length( unique(df.out$f__in[sel]) ) # 102
length( unique(df.out$rxn_id[sel]) ) # 226
unique(df.out$f__in[sel])
# [1] "1-acyl-sn-glycerol-3-phosphate acyltransferase (EC 2.3.1.51)"                                                            
# [2] "Acyl-CoA:1-acyl-sn-glycerol-3-phosphate acyltransferase (EC 2.3.1.51)"                                                   
# [3] "Beta-ketoadipyl CoA thiolase (EC 2.3.1.-)"                                                                               
# [4] "Phenylitaconyl-CoA hydratase (EC 4.2.1.-)"                                                                               
# [5] "Beta-lysine acetyltransferase (EC 2.3.1.-)"                                                                              
# [6] "3-oxoacyl-[acyl-carrier protein] reductase (EC 1.1.1.100)"                                                               
# [7] "Acetyl-CoA acetyltransferase (EC 2.3.1.9); Beta-ketoadipyl CoA thiolase (EC 2.3.1.-)"                                    
# [8] "Aminotransferase, class II (EC 2.3.1.-)"                                                                                 
# [9] "Octanoate-[acyl-carrier-protein]-protein-N-octanoyltransferase (EC 2.3.1.181)"                                           
# [10] "CDP-4-dehydro-6-deoxy-D-glucose 3-dehydratase (EC 4.2.1.-)"                                                              
# [11] "2-ketobutyrate formate-lyase (EC 2.3.1.-)"                                                                               
# [12] "3-oxoacyl-[acyl-carrier protein] reductase paralog (EC 1.1.1.100) in cluster with unspecified monosaccharide transporter"
# [13] "decarboxylase"                                                                                                           
# [14] "Ribosomal-protein-S18p-alanine acetyltransferase (EC 2.3.1.-)"                                                           
# [15] "Probable poly(beta-D-mannuronate) O-acetylase (EC 2.3.1.-)"                                                              
# [16] "Colanic acid biosynthesis acetyltransferase WcaB (EC 2.3.1.-)"                                                           
# [17] "Colanic acid biosynthesis acetyltransferase WcaF (EC 2.3.1.-)"                                                           
# [18] "UDP-N-acetylglucosamine 4,6-dehydratase (EC 4.2.1.-)"                                                                    
# [19] "Apolipoprotein N-acyltransferase (EC 2.3.1.-) in lipid-linked oligosaccharide synthesis cluster"                         
# [20] "Apolipoprotein N-acyltransferase (EC 2.3.1.-) in lipid-linked oligosaccharide synthesis cluster # truncated"             
# [21] "Phenolpthiocerol synthesis type-I polyketide synthase PpsA (EC 2.3.1.41)"                                                
# [22] "Phenolpthiocerol synthesis type-I polyketide synthase PpsA (EC 2.3.1.41) # fragment"                                     
# [23] "Phenolpthiocerol synthesis type-I polyketide synthase PpsB (EC 2.3.1.41)"                                                
# [24] "Phenolpthiocerol synthesis type-I polyketide synthase PpsB (EC 2.3.1.41) # fragment"                                     
# [25] "Phthiocerol/phthiodiolone dimycocerosyl transferase PapA5 (EC 2.3.1.-)"                                                  
# [26] "3-oxoacyl-.acyl-carrier protein. reductase .EC 1.1.1.100."                                                               
# [27] "3-oxoacyl-[acyl-carrier-protein] synthase, KASII (EC 2.3.1.179)"                                                         
# [28] "Mycolyl transferase 85C (EC 2.3.1.-)"                                                                                    
# [29] "Malonyl CoA-acyl carrier protein transacylase .EC 2.3.1.39."                                                             
# [30] "Lipid A biosynthesis .KDO. 2-.lauroyl.-lipid IVA acyltransferase .EC 2.3.1.-."                                           
# [31] "Lipid A biosynthesis lauroyl acyltransferase .EC 2.3.1.-."                                                               
# [32] "Lipid A biosynthesis lauroyl acyltransferase (EC 2.3.1.241)"                                                             
# [33] "UDP-3-O-[3-hydroxymyristoyl] glucosamine N-acyltransferase (EC 2.3.1.-)"                                                 
# [34] "tRNA pseudouridine 13 synthase (EC 4.2.1.-)"                                                                             
# [35] "Apolipoprotein N-acyltransferase (EC 2.3.1.-)"                                                                           
# [36] "Enoyl-[acyl-carrier-protein] reductase [FMN] (EC 1.3.1.9); Glutamate-1-semialdehyde aminotransferase (EC 5.4.3.8)"       
# [37] "3-oxoacyl-.acyl-carrier protein. reductase"                                                                              
# [38] "Malonyl-[acyl-carrier protein] O-methyltransferase (EC 2.1.1.197)"                                                       
# [39] "Malonyl-[acyl-carrier protein] O-methyltransferase (EC 2.1.1.197) ## BioC"                                               
# [40] "Acyl-[acyl-carrier-protein] synthetase (EC 6.2.1.20)"                                                                    
# [41] "Cyclohexanecarboxylate-CoA ligase AliA (EC 6.2.1.-)"                                                                     
# [42] "O-methyltransferase"                                                                                                     
# [43] "Lipoamide acyltransferase component of branched-chain alpha-keto acid dehydrogenase complex (EC 2.3.1.-)"                
# [44] "Lipoyl synthase (EC 2.8.1.8)"                                                                                            
# [45] "Dihydrolipoamide acetyltransferase component (E2) of acetoin dehydrogenase complex (EC 2.3.1.-)"                         
# [46] "Octanoate-[acyl-carrier-protein]-protein-N-octanoyltransferase, LipM (EC 2.3.1.181)"                                     
# [47] "2-vinyl bacteriochlorophyllide hydratase BchF (EC 4.2.1.-)"                                                              
# [48] "3-oxoacyl-[acyl-carrier-protein] synthase (EC 2.3.1.41)"                                                                 
# [49] "Enoyl-[acyl-carrier-protein] reductase [FMN] (EC 1.3.1.9), inferred for PFA pathway"                                     
# [50] "Probable acyl-ACP desaturase, Stearoyl-ACP desaturase (EC 1.14.19.2)"                                                    
# [51] "Carnitinyl-CoA dehydratase (EC 4.2.1.-)"                                                                                 
# [52] "Enoyl-[acyl-carrier-protein] reductase of FASI (EC 1.3.1.9)"                                                             
# [53] "3-hydroxypalmitoyl-[acyl-carrier-protein] dehydratase of FASI (EC 4.2.1.61)"                                             
# [54] "[Acyl-carrier-protein] palmitoyl transferase of FASI (EC 2.3.1.-)"                                                       
# [55] "3-oxoacyl-[acyl-carrier-protein] reductase of FASI (EC 1.1.1.100)"                                                       
# [56] "3-oxoacyl-[acyl-carrier-protein] synthase of FASI (EC 2.3.1.41)"                                                         
# [57] "Oleoyl-[acyl carrier protein] thioesterase of FASI (EC 3.1.2.14)"                                                        
# [58] "\"3-oxoacyl-[acyl-carrier-protein] synthase, KASIII (EC 2.3.1.41)\""                                                     
# [59] "3-hydroxyacyl-[acyl-carrier-protein] dehydratase, FabA form (EC 4.2.1.59)"                                               
# [60] "3-hydroxyacyl-[acyl-carrier-protein] dehydratase, FabZ form (EC 4.2.1.59)"                                               
# [61] "3-hydroxydecanoyl-[acyl-carrier-protein] dehydratase (EC 4.2.1.59)"                                                      
# [62] "3-oxoacyl-[ACP] reductase (EC 1.1.1.100)"                                                                                
# [63] "3-oxoacyl-[ACP] synthase (EC 2.3.1.41) FabV like"                                                                        
# [64] "3-oxoacyl-[acyl-carrier-protein] synthase (EC 2.3.1.41); Malonyl CoA-acyl carrier protein transacylase (EC 2.3.1.39)"    
# [65] "3-oxoacyl-[acyl-carrier-protein] synthase II (EC 2.3.1.41)"                                                              
# [66] "3-oxoacyl-[acyl-carrier-protein] synthase III (EC 2.3.1.41)"                                                             
# [67] "3-oxoacyl-.acyl-carrier-protein. synthase, KASI .EC 2.3.1.41."                                                           
# [68] "3-oxoacyl-[acyl-carrier-protein] synthase, KASIII (EC 2.3.1.180)"                                                        
# [69] "beta-ketodecanoyl-[acyl-carrier-protein] synthase (EC 2.3.1.207)"                                                        
# [70] "Enoyl-[acyl-carrier-protein] reductase [FMN] (EC 1.3.1.9)"                                                               
# [71] "Malonyl CoA-acyl carrier protein transacylase (EC 2.3.1.39); Enoyl-[acyl-carrier-protein] reductase [FMN] (EC 1.3.1.9)"  
# [72] "3-oxoacyl-[acyl-carrier-protein] synthase, KASI (EC 2.3.1.41)"                                                           
# [73] "Enoyl-[acyl-carrier-protein] reductase [NADH] (EC 1.3.1.9)"                                                              
# [74] "Oleoyl-[acyl carrier protein] thioesterase (EC 3.1.2.14)"                                                                
# [75] "(3R)-hydroxymyristoyl-[ACP] dehydratase (EC 4.2.1.-)"                                                                    
# [76] "3-hydroxydecanoyl-[ACP] dehydratase (EC 4.2.1.59)"                                                                       
# [77] "3-hydroxydecanoyl-[ACP] dehydratase (EC 4.2.1.60)"                                                                       
# [78] "FIG138576: 3-oxoacyl-[ACP] synthase (EC 2.3.1.41)"                                                                       
# [79] "3-oxoacyl-[acyl-carrier-protein] synthase, mitochondrial (EC 2.3.1.41)"                                                  
# [80] "Glycosyl-4,4'-diaponeurosporenoate acyltransferase precursor (EC 2.3.1.-)"                                               
# [81] "Glycerol-3-phosphate acyltransferase (EC 2.3.1.15)"                                                                      
# [82] "1-acyl-sn-glycerol-3-phosphate acyltransferase .EC 2.3.1.51."                                                            
# [83] "Acyl-phosphate:glycerol-3-phosphate O-acyltransferase PlsY"                                                              
# [84] "Glycerol-3-phosphate acyltransferase .EC 2.3.1.15."                                                                      
# [85] "Phosphate:acyl-ACP acyltransferase PlsX"                                                                                 
# [86] "Lysine N-acyltransferase MbtK (EC 2.3.1.-)"                                                                              
# [87] "2-hydroxypenta-2,4-dienoate hydratase (CmtF) (EC 4.2.1.-)"                                                               
# [88] "2-hydroxypenta-2,4-dienoate hydratase (CmtF) (EC 4.2.1.-) # *** Dlit ***"                                                
# [89] "2,3-dihydroxy-2,3-dihydro-p-cumate dehydrogenase(CmtB) (EC 1.1.1.100)"                                                   
# [90] "2,3-dihydroxy-2,3-dihydro-p-cumate dehydrogenase(CmtB) (EC 1.1.1.100) # *** Dlit ***"                                    
# [91] "Acyl-homoserine-lactone synthase LuxI (EC 2.3.1.184)"                                                                    
# [92] "Rhamnolipids biosynthesis 3-oxoacyl-[acyl-carrier-protein] reductase RhlG (EC 1.1.1.100)"                                
# [93] "Fatty acyl-CoA synthetase of sporopollenin biosynthesis (EC 6.2.1.-)"                                                    
# [94] " Acyl-homoserine-lactone synthase BpsI (EC 2.3.1.184)"                                                                   
# [95] " Acyl-homoserine-lactone synthase LuxI (EC 2.3.1.184)"                                                                   
# [96] " Autoinducer synthesis protein LuxI (EC 2.3.1.184)"                                                                      
# [97] " N-acyl-L-homoserine lactone synthetase TraI (EC 2.3.1.184)"                                                             
# [98] "Autoinducer 2 (AI-2) aldolase LsrF (EC 4.2.1.-)"                                                                         
# [99] "L-2,4-diaminobutyric acid acetyltransferase (EC 2.3.1.-)"                                                                
# [100] "L-2,4-diaminobutyric acid acetyltransferase (EC 2.3.1.-) # EctA"                                                         
# [101] "L-ectoine synthase (EC 4.2.1.-)"                                                                                         
# [102] "L-ectoine synthase (EC 4.2.1.-) # EctC"




## T2D-SWE

# link compounds to functions
df.out <-readRDS( "df.out--tidy-compounds_indiv-cpp3d--forslund-t2d-swe-hostremoval-pawsey.RDS" )
str(df.out)
# 'data.frame':	507391 obs. of  9 variables:
# $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_2" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" ...
# $ rxn_id             : chr  "rxn25279" "rxn25279" "rxn25279" "rxn03060" ...
# $ cpd_id             : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ cpd_name           : chr  "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" "2-methyl-trans-aconitate" "2-Methylcitrate" ...
# $ cpd_form           : chr  "C7H7O7" "H2O" "C7H5O6" "C7H7O7" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.3333 0.3333 0.3333 0.0556 0.0556 ...

phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval.RDS")
df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 17912     4
table(phy@sam_data$Status)
# ND CTRL T2D metformin- 
#   43             33 
length(row.names(phy@tax_table)) # 17912
head(row.names(phy@tax_table)) # "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"

dim(df.out) # 507391      9
sel <- which(df.out$superfocus_fxn %in% row.names(phy@tax_table)) # 507391

length(blcfa) # 35

sel <- which(df.out$cpd_id %in% blcfa) # 72

df.out[sel, ]

# blcfa names
sort( unique(df.out$cpd_name[sel]) )
# [1] "(2E)-Decenoyl-[acp]"                   "(2E)-Dodecenoyl-[acp]"                 "(2E)-Hexenoyl-[acp]"                   "10-methyl-dodecanoyl-ACP"             
# [5] "10-methyl-trans-dodec-2-enoyl-ACP"     "10-methyl-trans-undec-2-enoyl-ACP"     "10-methyl-undecanoyl-ACP"              "11-methyl-dodecanoyl-ACP"             
# [9] "11-methyl-trans-dodec-2-enoyl-ACP"     "12-methyl-trans-tetra-dec-2-enoyl-ACP" "12-methyl-trans-tridec-2-enoyl-ACP"    "13-methyl-trans-tetra-dec-2-enoyl-ACP"
# [13] "14-methyl-trans-hexa-dec-2-enoyl-ACP"  "14-methyl-trans-pentadec-2-enoyl-ACP"  "15-methyl-trans-hexa-dec-2-enoyl-ACP"  "4-methyl-hexanoyl-ACP"                
# [17] "4-methyl-pentanoyl-ACP"                "4-methyl-trans-hex-2-enoyl-ACP"        "4-methyl-trans-pent-2-enoyl-ACP"       "5-methyl-hexanoyl-ACP"                
# [21] "5-methyl-trans-hex-2-enoyl-ACP"        "6-methyl-heptanoyl-ACP"                "6-methyl-octanoyl-ACP"                 "6-methyl-trans-hept-2-enoyl-ACP"      
# [25] "6-methyl-trans-oct-2-enoyl-ACP"        "7-methyl-octanoyl-ACP"                 "7-methyl-trans-oct-2-enoyl-ACP"        "8-methyl-decanoyl-ACP"                
# [29] "8-methyl-nonanoyl-ACP"                 "8-methyl-trans-dec-2-enoyl-ACP"        "8-methyl-trans-non-2-enoyl-ACP"        "9-methyl-decanoyl-ACP"                
# [33] "9-methyl-trans-dec-2-enoyl-ACP"        "But-2-enoyl-[acyl-carrier protein]"    "trans-Octodec-2-enoyl-ACP"


# no of unique reactions
length( unique(df.out$rxn_id[sel]) ) # 47
# no of unique sub-functions 
length( unique(df.out$f__in[sel])) # 3
unique(df.out$f__in[sel])
# [1] "Trans-2-decenoyl-[acyl-carrier-protein] isomerase (EC 5.3.3.14)" "Enoyl-.acyl-carrier-protein. reductase .NADH. .EC 1.3.1.9."     
# [3] "Enoyl-.acyl-carrier-protein. reductase .NADPH. .EC 1.3.1.10." 
# no of unique superfocus functions
length( unique(df.out$superfocus_fxn[sel])) # 4

(sf_fxn <- sort(unique(df.out$superfocus_fxn[sel])))

sel.all <- which(df.out$superfocus_fxn %in% sf_fxn)
length(sel.all) # 317
length(unique(df.out$rxn_id[sel.all])) # 62 total reactions associated with those functions

( cpd_all <- unique(df.out$cpd_id[sel.all]) )
length(cpd_all) # 75

sel.other <- which(!cpd_all %in% blcfa) # 40 other compounds are in the same set of reactions
cpd_other <- cpd_all[sel.other]
length(cpd_other) # 40
length(unique(cpd_other)) # 40

sel.wider <- which(df.out$cpd_id %in% cpd_other)
length(sel.wider) # 120462
( other_cpd_names <- sort( unique(df.out$cpd_name[sel.wider]) ) )
# [1] "(2E)-Hexadecenoyl-[acp]"                 "(2E)-Octadecenoyl-[acp]"                 "(2E)-Octenoyl-[acp]"                    
# [4] "(2E)-Tetradecenoyl-[acp]"                "12-methyl-tetra-decanoyl-ACP"            "12-methyl-tridecanoyl-ACP"              
# [7] "13-methyl-tetra-decanoyl-ACP"            "14-methyl-hexa-decanoyl-ACP"             "14-methyl-pentadecanoyl-ACP"            
# [10] "15-methyl-hexa-decanoyl-ACP"             "3-hydroxy-cis-D7-tetraecenoyl-ACPs"      "3-hydroxy-cis-D9-hexaecenoyl-ACPs"      
# [13] "3-Hydroxy-octanoyl-ACPs"                 "3-Hydroxyglutaryl-ACP-methyl-ester"      "3-hydroxypimeloyl-ACP-methyl-esters"    
# [16] "a 3-hydroxy cis delta5-dodecenoyl-[acp]" "Butyryl-ACP"                             "cis-3-Decenoyl-[acyl-carrier protein]"  
# [19] "Cis-delta-3-decenoyl-ACPs"               "Decanoyl-ACP"                            "Dodecanoyl-ACP"                         
# [22] "Enoylglutaryl-ACP-methyl-esters"         "Enoylpimeloyl-ACP-methyl-esters"         "H+"                                     
# [25] "H2O"                                     "hexadecanoyl-acp"                        "Hexanoyl-ACP"                           
# [28] "Myristoyl-ACP"                           "NAD"                                     "NADH"                                   
# [31] "NADP"                                    "NADPH"                                   "Octanoyl-ACP"                           
# [34] "Octodecanoyl-ACP"                        "R-3-hydroxystearoyl-ACPs"                "Stearoyl-[acyl-carrier protein]"        
# [37] "Trans-D2-decenoyl-ACPs"                  "Trans-D3-cis-D5-dodecenoyl-ACPs"         "Trans-D3-cis-D7-tetradecenoyl-ACPs"     
# [40] "Trans-D3-cis-D9-hexadecenoyl-ACPs"

# only keep other ACPs

other_ACP_names <- c(
  
  "(2E)-Hexadecenoyl-[acp]"                , "(2E)-Octadecenoyl-[acp]"          ,       "(2E)-Octenoyl-[acp]"                   , 
  "(2E)-Tetradecenoyl-[acp]"              ,  "12-methyl-tetra-decanoyl-ACP"     ,       "12-methyl-tridecanoyl-ACP"             , 
  "13-methyl-tetra-decanoyl-ACP"          ,  "14-methyl-hexa-decanoyl-ACP"      ,       "14-methyl-pentadecanoyl-ACP"           , 
  "15-methyl-hexa-decanoyl-ACP"           ,  "3-hydroxy-cis-D7-tetraecenoyl-ACPs",      "3-hydroxy-cis-D9-hexaecenoyl-ACPs"     , 
  "3-Hydroxy-octanoyl-ACPs"               ,  "3-Hydroxyglutaryl-ACP-methyl-ester" ,     "3-hydroxypimeloyl-ACP-methyl-esters"    ,
  "a 3-hydroxy cis delta5-dodecenoyl-[acp]", "Butyryl-ACP"                        ,     "cis-3-Decenoyl-[acyl-carrier protein]" , 
  "Cis-delta-3-decenoyl-ACPs"             ,  "Decanoyl-ACP"                       ,     "Dodecanoyl-ACP"                        , 
  "Enoylglutaryl-ACP-methyl-esters"       ,  "Enoylpimeloyl-ACP-methyl-esters"    ,  
  
  # "H+"                                    , 
  # "H2O"                                   , 
  
  "hexadecanoyl-acp"                   ,     "Hexanoyl-ACP"                          , 
  "Myristoyl-ACP"                         ,  
  
  #"NAD"                                ,     "NADH"                                  , 
  # "NADP"                                  ,  "NADPH"                              ,   
  
  "Octanoyl-ACP"                          , 
  "Octodecanoyl-ACP"                      ,  "R-3-hydroxystearoyl-ACPs"           ,     "Stearoyl-[acyl-carrier protein]"       , 
  "Trans-D2-decenoyl-ACPs"                ,  "Trans-D3-cis-D5-dodecenoyl-ACPs"    ,     "Trans-D3-cis-D7-tetradecenoyl-ACPs"    , 
  "Trans-D3-cis-D9-hexadecenoyl-ACPs"

)

length(other_ACP_names) # 34



# how many wider functions do the other ACP name compounds occur in?
# not including the initial 4 superfocus functions linked to BLCFA?

sel <- which(df.out$cpd_name %in% other_ACP_names ) # 704
sel <- which(df.out$cpd_name %in% other_ACP_names & !df.out$superfocus_fxn %in% sf_fxn ) # 634

sf_fxn_other_ACP <- unique(df.out$superfocus_fxn[sel]) # 115
length( unique(df.out$f__in[sel]) ) # 52
length( unique(df.out$rxn_id[sel]) ) # 87
unique(df.out$f__in[sel])
# [1] "Beta-ketoadipyl CoA thiolase (EC 2.3.1.-)"                                                                               
# [2] "Acetyl-CoA acetyltransferase (EC 2.3.1.9); Beta-ketoadipyl CoA thiolase (EC 2.3.1.-)"                                    
# [3] "Octanoate-[acyl-carrier-protein]-protein-N-octanoyltransferase (EC 2.3.1.181)"                                           
# [4] "CDP-4-dehydro-6-deoxy-D-glucose 3-dehydratase (EC 4.2.1.-)"                                                              
# [5] "2-ketobutyrate formate-lyase (EC 2.3.1.-)"                                                                               
# [6] "Ribosomal-protein-S18p-alanine acetyltransferase (EC 2.3.1.-)"                                                           
# [7] "Probable poly(beta-D-mannuronate) O-acetylase (EC 2.3.1.-)"                                                              
# [8] "Colanic acid biosynthesis acetyltransferase WcaB (EC 2.3.1.-)"                                                           
# [9] "Colanic acid biosynthesis acetyltransferase WcaF (EC 2.3.1.-)"                                                           
# [10] "UDP-N-acetylglucosamine 4,6-dehydratase (EC 4.2.1.-)"                                                                    
# [11] "Apolipoprotein N-acyltransferase (EC 2.3.1.-) in lipid-linked oligosaccharide synthesis cluster"                         
# [12] "Lipid A biosynthesis .KDO. 2-.lauroyl.-lipid IVA acyltransferase .EC 2.3.1.-."                                           
# [13] "Lipid A biosynthesis lauroyl acyltransferase .EC 2.3.1.-."                                                               
# [14] "Lipid A biosynthesis lauroyl acyltransferase (EC 2.3.1.241)"                                                             
# [15] "UDP-3-O-[3-hydroxymyristoyl] glucosamine N-acyltransferase (EC 2.3.1.-)"                                                 
# [16] "tRNA pseudouridine 13 synthase (EC 4.2.1.-)"                                                                             
# [17] "1-acyl-sn-glycerol-3-phosphate acyltransferase (EC 2.3.1.51)"                                                            
# [18] "Acyl-CoA:1-acyl-sn-glycerol-3-phosphate acyltransferase (EC 2.3.1.51)"                                                   
# [19] "Apolipoprotein N-acyltransferase (EC 2.3.1.-)"                                                                           
# [20] "3-oxoacyl-[acyl-carrier protein] reductase paralog (EC 1.1.1.100) in cluster with unspecified monosaccharide transporter"
# [21] "Lipoyl synthase (EC 2.8.1.8)"                                                                                            
# [22] "Dihydrolipoamide acetyltransferase component (E2) of acetoin dehydrogenase complex (EC 2.3.1.-)"                         
# [23] "Octanoate-[acyl-carrier-protein]-protein-N-octanoyltransferase, LipM (EC 2.3.1.181)"                                     
# [24] "Carnitinyl-CoA dehydratase (EC 4.2.1.-)"                                                                                 
# [25] "Enoyl-[acyl-carrier-protein] reductase of FASI (EC 1.3.1.9)"                                                             
# [26] "3-hydroxypalmitoyl-[acyl-carrier-protein] dehydratase of FASI (EC 4.2.1.61)"                                             
# [27] "[Acyl-carrier-protein] palmitoyl transferase of FASI (EC 2.3.1.-)"                                                       
# [28] "3-oxoacyl-[acyl-carrier-protein] reductase of FASI (EC 1.1.1.100)"                                                       
# [29] "3-oxoacyl-[acyl-carrier-protein] synthase of FASI (EC 2.3.1.41)"                                                         
# [30] "3-hydroxyacyl-[acyl-carrier-protein] dehydratase, FabA form (EC 4.2.1.59)"                                               
# [31] "3-hydroxyacyl-[acyl-carrier-protein] dehydratase, FabZ form (EC 4.2.1.59)"                                               
# [32] "3-hydroxydecanoyl-[acyl-carrier-protein] dehydratase (EC 4.2.1.59)"                                                      
# [33] "3-oxoacyl-[ACP] reductase (EC 1.1.1.100)"                                                                                
# [34] "3-oxoacyl-[acyl-carrier-protein] synthase (EC 2.3.1.41)"                                                                 
# [35] "3-oxoacyl-[acyl-carrier-protein] synthase (EC 2.3.1.41); Malonyl CoA-acyl carrier protein transacylase (EC 2.3.1.39)"    
# [36] "3-oxoacyl-[acyl-carrier-protein] synthase II (EC 2.3.1.41)"                                                              
# [37] "3-oxoacyl-.acyl-carrier-protein. synthase, KASI .EC 2.3.1.41."                                                           
# [38] "Enoyl-[acyl-carrier-protein] reductase [FMN] (EC 1.3.1.9)"                                                               
# [39] "Malonyl CoA-acyl carrier protein transacylase .EC 2.3.1.39."                                                             
# [40] "Malonyl CoA-acyl carrier protein transacylase (EC 2.3.1.39); Enoyl-[acyl-carrier-protein] reductase [FMN] (EC 1.3.1.9)"  
# [41] "(3R)-hydroxymyristoyl-[ACP] dehydratase (EC 4.2.1.-)"                                                                    
# [42] "3-hydroxydecanoyl-[ACP] dehydratase (EC 4.2.1.59)"                                                                       
# [43] "3-hydroxydecanoyl-[ACP] dehydratase (EC 4.2.1.60)"                                                                       
# [44] "3-oxoacyl-[ACP] synthase (EC 2.3.1.41) FabV like"                                                                        
# [45] "FIG138576: 3-oxoacyl-[ACP] synthase (EC 2.3.1.41)"                                                                       
# [46] "Glycosyl-4,4'-diaponeurosporenoate acyltransferase precursor (EC 2.3.1.-)"                                               
# [47] "1-acyl-sn-glycerol-3-phosphate acyltransferase .EC 2.3.1.51."                                                            
# [48] "Acyl-phosphate:glycerol-3-phosphate O-acyltransferase PlsY"                                                              
# [49] "Glycerol-3-phosphate acyltransferase .EC 2.3.1.15."                                                                      
# [50] "Phosphate:acyl-ACP acyltransferase PlsX"                                                                                 
# [51] "Uncharacterized acetyltransferase YafP (EC 2.3.1.-)"                                                                     
# [52] "Autoinducer 2 (AI-2) aldolase LsrF (EC 4.2.1.-)" 









## T2D-CHN

# link compounds to functions
df.out <-readRDS( "df.out--tidy-compounds_indiv-cpp3d--forslund-t2d-chn-hostremoval-pawsey.RDS" )
str(df.out)
# 'data.frame':	554562 obs. of  9 variables:
# $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_3" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "4-hydroxyproline epimerase (EC 5.1.1.8)" "4-hydroxyproline epimerase (EC 5.1.1.8)" "Alanine racemase (EC 5.1.1.1)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id             : chr  "rxn02360" "rxn02360" "rxn00283" "rxn00283" ...
# $ cpd_id             : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ cpd_name           : chr  "trans-4-Hydroxy-L-proline" "cis-4-Hydroxy-D-proline" "L-Alanine" "D-Alanine" ...
# $ cpd_form           : chr  "C5H9NO3" "C5H9NO3" "C3H7NO2" "C3H7NO2" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.167 0.167 0.167 ...

phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal.RDS")
df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 19352     4

table(phy@sam_data$Diagnosis)
# ND CTRL T2D metformin- 
#   52             30 

min(taxa_sums(phy)) # 1.734989e-06

length(blcfa) # 35

sel <- which(df.out$cpd_id %in% blcfa) # 72

df.out[sel, ]

# blcfa names
sort( unique(df.out$cpd_name[sel]) )
# [1] "(2E)-Decenoyl-[acp]"                   "(2E)-Dodecenoyl-[acp]"                 "(2E)-Hexenoyl-[acp]"                   "10-methyl-dodecanoyl-ACP"             
# [5] "10-methyl-trans-dodec-2-enoyl-ACP"     "10-methyl-trans-undec-2-enoyl-ACP"     "10-methyl-undecanoyl-ACP"              "11-methyl-dodecanoyl-ACP"             
# [9] "11-methyl-trans-dodec-2-enoyl-ACP"     "12-methyl-trans-tetra-dec-2-enoyl-ACP" "12-methyl-trans-tridec-2-enoyl-ACP"    "13-methyl-trans-tetra-dec-2-enoyl-ACP"
# [13] "14-methyl-trans-hexa-dec-2-enoyl-ACP"  "14-methyl-trans-pentadec-2-enoyl-ACP"  "15-methyl-trans-hexa-dec-2-enoyl-ACP"  "4-methyl-hexanoyl-ACP"                
# [17] "4-methyl-pentanoyl-ACP"                "4-methyl-trans-hex-2-enoyl-ACP"        "4-methyl-trans-pent-2-enoyl-ACP"       "5-methyl-hexanoyl-ACP"                
# [21] "5-methyl-trans-hex-2-enoyl-ACP"        "6-methyl-heptanoyl-ACP"                "6-methyl-octanoyl-ACP"                 "6-methyl-trans-hept-2-enoyl-ACP"      
# [25] "6-methyl-trans-oct-2-enoyl-ACP"        "7-methyl-octanoyl-ACP"                 "7-methyl-trans-oct-2-enoyl-ACP"        "8-methyl-decanoyl-ACP"                
# [29] "8-methyl-nonanoyl-ACP"                 "8-methyl-trans-dec-2-enoyl-ACP"        "8-methyl-trans-non-2-enoyl-ACP"        "9-methyl-decanoyl-ACP"                
# [33] "9-methyl-trans-dec-2-enoyl-ACP"        "But-2-enoyl-[acyl-carrier protein]"    "trans-Octodec-2-enoyl-ACP"


# no of unique reactions
length( unique(df.out$rxn_id[sel]) ) # 47
# no of unique sub-functions 
length( unique(df.out$f__in[sel])) # 3
unique(df.out$f__in[sel])
# [1] "Trans-2-decenoyl-[acyl-carrier-protein] isomerase (EC 5.3.3.14)" "Enoyl-.acyl-carrier-protein. reductase .NADH. .EC 1.3.1.9."     
# [3] "Enoyl-.acyl-carrier-protein. reductase .NADPH. .EC 1.3.1.10."
# no of unique superfocus functions
length( unique(df.out$superfocus_fxn[sel])) # 4

(sf_fxn <- sort(unique(df.out$superfocus_fxn[sel])))

sel.all <- which(df.out$superfocus_fxn %in% sf_fxn)
length(sel.all) # 317
length(unique(df.out$rxn_id[sel.all])) # 62 total reactions associated with those functions

( cpd_all <- unique(df.out$cpd_id[sel.all]) )
length(cpd_all) # 75

sel.other <- which(!cpd_all %in% blcfa) # 40 other compounds are in the same set of reactions
cpd_other <- cpd_all[sel.other]
length(cpd_other) # 40
length(unique(cpd_other)) # 40

sel.wider <- which(df.out$cpd_id %in% cpd_other)
length(sel.wider) # 127826
( other_cpd_names <- sort( unique(df.out$cpd_name[sel.wider]) ) )
# [1] "(2E)-Hexadecenoyl-[acp]"                 "(2E)-Octadecenoyl-[acp]"                 "(2E)-Octenoyl-[acp]"                    
# [4] "(2E)-Tetradecenoyl-[acp]"                "12-methyl-tetra-decanoyl-ACP"            "12-methyl-tridecanoyl-ACP"              
# [7] "13-methyl-tetra-decanoyl-ACP"            "14-methyl-hexa-decanoyl-ACP"             "14-methyl-pentadecanoyl-ACP"            
# [10] "15-methyl-hexa-decanoyl-ACP"             "3-hydroxy-cis-D7-tetraecenoyl-ACPs"      "3-hydroxy-cis-D9-hexaecenoyl-ACPs"      
# [13] "3-Hydroxy-octanoyl-ACPs"                 "3-Hydroxyglutaryl-ACP-methyl-ester"      "3-hydroxypimeloyl-ACP-methyl-esters"    
# [16] "a 3-hydroxy cis delta5-dodecenoyl-[acp]" "Butyryl-ACP"                             "cis-3-Decenoyl-[acyl-carrier protein]"  
# [19] "Cis-delta-3-decenoyl-ACPs"               "Decanoyl-ACP"                            "Dodecanoyl-ACP"                         
# [22] "Enoylglutaryl-ACP-methyl-esters"         "Enoylpimeloyl-ACP-methyl-esters"         "H+"                                     
# [25] "H2O"                                     "hexadecanoyl-acp"                        "Hexanoyl-ACP"                           
# [28] "Myristoyl-ACP"                           "NAD"                                     "NADH"                                   
# [31] "NADP"                                    "NADPH"                                   "Octanoyl-ACP"                           
# [34] "Octodecanoyl-ACP"                        "R-3-hydroxystearoyl-ACPs"                "Stearoyl-[acyl-carrier protein]"        
# [37] "Trans-D2-decenoyl-ACPs"                  "Trans-D3-cis-D5-dodecenoyl-ACPs"         "Trans-D3-cis-D7-tetradecenoyl-ACPs"     
# [40] "Trans-D3-cis-D9-hexadecenoyl-ACPs" 

# only keep other ACPs

other_ACP_names <- c(
  
  "(2E)-Hexadecenoyl-[acp]"               ,  "(2E)-Octadecenoyl-[acp]"           ,      "(2E)-Octenoyl-[acp]"                    ,
  "(2E)-Tetradecenoyl-[acp]"              ,  "12-methyl-tetra-decanoyl-ACP"       ,     "12-methyl-tridecanoyl-ACP"             , 
  "13-methyl-tetra-decanoyl-ACP"          ,  "14-methyl-hexa-decanoyl-ACP"         ,    "14-methyl-pentadecanoyl-ACP"           , 
  "15-methyl-hexa-decanoyl-ACP"           ,  "3-hydroxy-cis-D7-tetraecenoyl-ACPs" ,     "3-hydroxy-cis-D9-hexaecenoyl-ACPs"     , 
  "3-Hydroxy-octanoyl-ACPs"               ,  "3-Hydroxyglutaryl-ACP-methyl-ester"  ,    "3-hydroxypimeloyl-ACP-methyl-esters"   , 
  "a 3-hydroxy cis delta5-dodecenoyl-[acp]", "Butyryl-ACP"                        ,     "cis-3-Decenoyl-[acyl-carrier protein]" , 
  "Cis-delta-3-decenoyl-ACPs"           ,    "Decanoyl-ACP"                       ,     "Dodecanoyl-ACP"                        , 
  "Enoylglutaryl-ACP-methyl-esters"     ,    "Enoylpimeloyl-ACP-methyl-esters"    ,   
  
  #"H+"                                    , 
  # "H2O"                                 ,  
  
  "hexadecanoyl-acp"                   ,     "Hexanoyl-ACP"                          , 
  "Myristoyl-ACP"                       ,  
  
  #"NAD"                                ,     "NADH"                                  , 
  # "NADP"                                ,    "NADPH"                              ,   
  
  "Octanoyl-ACP"                           ,
  "Octodecanoyl-ACP"                    ,    "R-3-hydroxystearoyl-ACPs"           ,     "Stearoyl-[acyl-carrier protein]"       , 
  "Trans-D2-decenoyl-ACPs"              ,    "Trans-D3-cis-D5-dodecenoyl-ACPs"    ,     "Trans-D3-cis-D7-tetradecenoyl-ACPs"    , 
  "Trans-D3-cis-D9-hexadecenoyl-ACPs" 
  
)

length(other_ACP_names) # 34


# how many wider functions do the other ACP name compounds occur in?
# not including the initial 4 superfocus functions linked to BLCFA?

sel <- which(df.out$cpd_name %in% other_ACP_names ) # 743
sel <- which(df.out$cpd_name %in% other_ACP_names & !df.out$superfocus_fxn %in% sf_fxn ) # 673

sf_fxn_other_ACP <- unique(df.out$superfocus_fxn[sel]) # 129
length( unique(df.out$f__in[sel]) ) # 54
length( unique(df.out$rxn_id[sel]) ) # 96
unique(df.out$f__in[sel])
# [1] "Beta-ketoadipyl CoA thiolase (EC 2.3.1.-)"                                                                               
# [2] "Acetyl-CoA acetyltransferase (EC 2.3.1.9); Beta-ketoadipyl CoA thiolase (EC 2.3.1.-)"                                    
# [3] "Beta-lysine acetyltransferase (EC 2.3.1.-)"                                                                              
# [4] "Octanoate-[acyl-carrier-protein]-protein-N-octanoyltransferase (EC 2.3.1.181)"                                           
# [5] "CDP-4-dehydro-6-deoxy-D-glucose 3-dehydratase (EC 4.2.1.-)"                                                              
# [6] "2-ketobutyrate formate-lyase (EC 2.3.1.-)"                                                                               
# [7] "Ribosomal-protein-S18p-alanine acetyltransferase (EC 2.3.1.-)"                                                           
# [8] "Probable poly(beta-D-mannuronate) O-acetylase (EC 2.3.1.-)"                                                              
# [9] "Colanic acid biosynthesis acetyltransferase WcaB (EC 2.3.1.-)"                                                           
# [10] "Colanic acid biosynthesis acetyltransferase WcaF (EC 2.3.1.-)"                                                           
# [11] "UDP-N-acetylglucosamine 4,6-dehydratase (EC 4.2.1.-)"                                                                    
# [12] "Apolipoprotein N-acyltransferase (EC 2.3.1.-) in lipid-linked oligosaccharide synthesis cluster"                         
# [13] "Lipid A biosynthesis .KDO. 2-.lauroyl.-lipid IVA acyltransferase .EC 2.3.1.-."                                           
# [14] "Lipid A biosynthesis lauroyl acyltransferase .EC 2.3.1.-."                                                               
# [15] "Lipid A biosynthesis lauroyl acyltransferase (EC 2.3.1.241)"                                                             
# [16] "UDP-3-O-[3-hydroxymyristoyl] glucosamine N-acyltransferase (EC 2.3.1.-)"                                                 
# [17] "tRNA pseudouridine 13 synthase (EC 4.2.1.-)"                                                                             
# [18] "1-acyl-sn-glycerol-3-phosphate acyltransferase (EC 2.3.1.51)"                                                            
# [19] "Acyl-CoA:1-acyl-sn-glycerol-3-phosphate acyltransferase (EC 2.3.1.51)"                                                   
# [20] "Apolipoprotein N-acyltransferase (EC 2.3.1.-)"                                                                           
# [21] "3-oxoacyl-[acyl-carrier protein] reductase paralog (EC 1.1.1.100) in cluster with unspecified monosaccharide transporter"
# [22] "Acyl-[acyl-carrier-protein] synthetase (EC 6.2.1.20)"                                                                    
# [23] "Lipoyl synthase (EC 2.8.1.8)"                                                                                            
# [24] "Dihydrolipoamide acetyltransferase component (E2) of acetoin dehydrogenase complex (EC 2.3.1.-)"                         
# [25] "Carnitinyl-CoA dehydratase (EC 4.2.1.-)"                                                                                 
# [26] "Enoyl-[acyl-carrier-protein] reductase of FASI (EC 1.3.1.9)"                                                             
# [27] "3-hydroxypalmitoyl-[acyl-carrier-protein] dehydratase of FASI (EC 4.2.1.61)"                                             
# [28] "[Acyl-carrier-protein] palmitoyl transferase of FASI (EC 2.3.1.-)"                                                       
# [29] "3-oxoacyl-[acyl-carrier-protein] reductase of FASI (EC 1.1.1.100)"                                                       
# [30] "3-oxoacyl-[acyl-carrier-protein] synthase of FASI (EC 2.3.1.41)"                                                         
# [31] "3-hydroxyacyl-[acyl-carrier-protein] dehydratase, FabA form (EC 4.2.1.59)"                                               
# [32] "3-hydroxyacyl-[acyl-carrier-protein] dehydratase, FabZ form (EC 4.2.1.59)"                                               
# [33] "3-hydroxydecanoyl-[acyl-carrier-protein] dehydratase (EC 4.2.1.59)"                                                      
# [34] "3-oxoacyl-[ACP] reductase (EC 1.1.1.100)"                                                                                
# [35] "3-oxoacyl-[acyl-carrier-protein] synthase (EC 2.3.1.41)"                                                                 
# [36] "3-oxoacyl-[acyl-carrier-protein] synthase (EC 2.3.1.41); Malonyl CoA-acyl carrier protein transacylase (EC 2.3.1.39)"    
# [37] "3-oxoacyl-[acyl-carrier-protein] synthase II (EC 2.3.1.41)"                                                              
# [38] "3-oxoacyl-.acyl-carrier-protein. synthase, KASI .EC 2.3.1.41."                                                           
# [39] "Enoyl-[acyl-carrier-protein] reductase [FMN] (EC 1.3.1.9)"                                                               
# [40] "Malonyl CoA-acyl carrier protein transacylase .EC 2.3.1.39."                                                             
# [41] "Malonyl CoA-acyl carrier protein transacylase (EC 2.3.1.39); Enoyl-[acyl-carrier-protein] reductase [FMN] (EC 1.3.1.9)"  
# [42] "Oleoyl-[acyl carrier protein] thioesterase (EC 3.1.2.14)"                                                                
# [43] "(3R)-hydroxymyristoyl-[ACP] dehydratase (EC 4.2.1.-)"                                                                    
# [44] "3-hydroxydecanoyl-[ACP] dehydratase (EC 4.2.1.59)"                                                                       
# [45] "3-hydroxydecanoyl-[ACP] dehydratase (EC 4.2.1.60)"                                                                       
# [46] "3-oxoacyl-[ACP] synthase (EC 2.3.1.41) FabV like"                                                                        
# [47] "FIG138576: 3-oxoacyl-[ACP] synthase (EC 2.3.1.41)"                                                                       
# [48] "1-acyl-sn-glycerol-3-phosphate acyltransferase .EC 2.3.1.51."                                                            
# [49] "Acyl-phosphate:glycerol-3-phosphate O-acyltransferase PlsY"                                                              
# [50] "Glycerol-3-phosphate acyltransferase .EC 2.3.1.15."                                                                      
# [51] "Phosphate:acyl-ACP acyltransferase PlsX"                                                                                 
# [52] "Uncharacterized acetyltransferase YafP (EC 2.3.1.-)"                                                                     
# [53] " Acyl-homoserine-lactone synthase LuxI (EC 2.3.1.184)"                                                                   
# [54] "Autoinducer 2 (AI-2) aldolase LsrF (EC 4.2.1.-)" 


#-------------------------

#### II. Track compounds / Enzymes (sub-functions) / SUPER-FOCUS functions associated with BCFA-ACPs?
#### Get all ACPs that don't feature in the list of 35 ACP compounds
#-------------------------

sort(t2d.zoom1.decreasing.proteins)
# [1] "cpd11465" "cpd11469" "cpd11473" "cpd11475" "cpd11498" "cpd11499" "cpd11502" "cpd11503" "cpd11506" "cpd11507" "cpd11510" "cpd11511"
# [13] "cpd11514" "cpd11518" "cpd11523" "cpd11524" "cpd11527" "cpd11528" "cpd11531" "cpd11532" "cpd11535" "cpd11536" "cpd11539" "cpd11543"
# [25] "cpd11548" "cpd11549" "cpd11552" "cpd11553" "cpd11556" "cpd11557" "cpd11560" "cpd11561" "cpd11564" "cpd11568" "cpd11572"

acp35 <- t2d.zoom1.decreasing.proteins

## Sunbad resto

# link compounds to functions
df.out <-readRDS( "df.out--tidy-compounds_indiv--cpp3d-sunbad-resto.RDS" )

str(df.out)
# 'data.frame':	1154547 obs. of  9 variables:
# $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_3" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylaconitate isomerase" "2-methylaconitate isomerase" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" ...
# $ rxn_id             : chr  "rxn25278" "rxn25278" "rxn25279" "rxn25279" ...
# $ cpd_id             : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ cpd_name           : chr  "2-methyl-trans-aconitate" "cis-2-Methylaconitate" "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" ...
# $ cpd_form           : chr  "C7H5O6" "C7H5O6" "C7H7O7" "H2O" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.333 0.333 0.333 ...

phy <- readRDS("phy-phyloseq-object-sunbad-resto.RDS")
df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 30125     4
min(taxa_sums(phy)) # 2.521458e-06

length(acp35) # 35

sel1 <- which(df.out$cpd_id %in% acp35) # 131

df.out[sel1, ]

# ACP-35 names
unique(df.out$cpd_name[sel1])
# [1] "But-2-enoyl-[acyl-carrier protein]"    "(2E)-Dodecenoyl-[acp]"                 "(2E)-Hexenoyl-[acp]"                   "(2E)-Decenoyl-[acp]"                  
# [5] "4-methyl-trans-hex-2-enoyl-ACP"        "4-methyl-hexanoyl-ACP"                 "6-methyl-trans-oct-2-enoyl-ACP"        "6-methyl-octanoyl-ACP"                
# [9] "8-methyl-trans-dec-2-enoyl-ACP"        "8-methyl-decanoyl-ACP"                 "10-methyl-trans-dodec-2-enoyl-ACP"     "10-methyl-dodecanoyl-ACP"             
# [13] "12-methyl-trans-tetra-dec-2-enoyl-ACP" "14-methyl-trans-hexa-dec-2-enoyl-ACP"  "5-methyl-trans-hex-2-enoyl-ACP"        "5-methyl-hexanoyl-ACP"                
# [17] "7-methyl-trans-oct-2-enoyl-ACP"        "7-methyl-octanoyl-ACP"                 "9-methyl-trans-dec-2-enoyl-ACP"        "9-methyl-decanoyl-ACP"                
# [21] "11-methyl-trans-dodec-2-enoyl-ACP"     "11-methyl-dodecanoyl-ACP"              "13-methyl-trans-tetra-dec-2-enoyl-ACP" "15-methyl-trans-hexa-dec-2-enoyl-ACP" 
# [25] "4-methyl-trans-pent-2-enoyl-ACP"       "4-methyl-pentanoyl-ACP"                "6-methyl-trans-hept-2-enoyl-ACP"       "6-methyl-heptanoyl-ACP"               
# [29] "8-methyl-trans-non-2-enoyl-ACP"        "8-methyl-nonanoyl-ACP"                 "10-methyl-trans-undec-2-enoyl-ACP"     "10-methyl-undecanoyl-ACP"             
# [33] "12-methyl-trans-tridec-2-enoyl-ACP"    "14-methyl-trans-pentadec-2-enoyl-ACP"  "trans-Octodec-2-enoyl-ACP" 


## what about other ACPs?

#sel2 <- grep(pattern = "-[acyl-carrier protein]|-[acp]|-ACP", x = df.out$cpd_name)
sel2 <- grep(pattern = "acyl-carrier protein|acp|ACP", x = df.out$cpd_name) # qty 20035

subsel <- which( sel2 %in% sel1 ) # 131

length( df.out$cpd_name[sel2[-subsel]] ) # 19904
length( unique(df.out$cpd_name[sel2[-subsel]]) ) # 296

( other_ACP_names <- unique(df.out$cpd_name[sel2[-subsel]]) )
# [1] "Pimelyl-[acyl-carrier protein]"                                      "ACP"                                                                 "Dodecanoyl-ACP"                                                     
# [4] "Myristoyl-ACP"                                                       "Tetradecenoyl-ACP"                                                   "Palmitoyl-ACP"                                                      
# [7] "Hexadecenoyl-ACP"                                                    "Octadecanoyl-ACP"                                                    "Octadecenoyl-ACP"                                                   
# [10] "Acyl-[acyl-carrier protein]"                                         "Saturated-Fatty-Acyl-ACPs"                                           "18-Carbamoyl-3,5,7,9,11,13,15,17-octaoxo-octadecanoyl-[acp]"        
# [13] "18-Carbamoyl-9-hydroxy-3,5,7,11,13,15,17-octaoxo-octadecanoyl-[acp]" "dihomo gamma linolenoyl-2-enoyl [acp]"                               "a dihomo gamma linolenoyl [acp]"                                    
# [16] "eicosatrienoyl-2-enoyl [acp]"                                        "EICOSATRIENOYL-ACP"                                                  "DOCOSAPENTAENOYL-2-ENOYL-ACP"                                       
# [19] "DOCOSAPENTAENOYL-ACP"                                                "Petroselinoyl-ACPs"                                                  "Petrosel-2-enoyl-ACPs"                                              
# [22] "3,5,7,9,11,13,15-Heptaoxo-hexadecanoyl-[acp]"                        "3,5,7,9,11,13,15-Heptaoxo-octadecanoyl-[acp]"                        "3,5,7,9,11,13,15,17,19-Nonaoxo-eicosanoyl-[acp]"                    
# [25] "3,5,7,9,11,13,15,17,19-Nonaoxo-henicosanoyl-[acp]"                   "3,5,7,9,11,13,15,17,19-Nonaoxo-docosanoyl-[acp]"                     "Myristoyl-ACPs"                                                     
# [28] "Lignoceroyl-ACPs"                                                    "3-oxo-cerotoyl-ACPs"                                                 "R-3-hydroxycerotoyl-ACPs"                                           
# [31] "Trans-D2-hexacos-2-enoyl-ACPs"                                       "Cerotoyl-ACPs"                                                       "Acetoacetyl-ACP"                                                    
# [34] "Heptadecanoyl-ACPs"                                                  "Long-Chain-Acyl-ACPs"                                                "Stearoyl-[acyl-carrier protein]"                                    
# [37] "Polyketide-ACP-Proteins"                                             "Hepta-oxo-hexadecanoyl-ACPs"                                         "R-3-hydroxyarachidoyl-ACPs"                                         
# [40] "a trans-eicos-2-enoyl-[acp]"                                         "3-oxo-arachidoyl-ACPs"                                               "Arachidoyl-ACPs"                                                    
# [43] "3-oxo-behenoyl-ACPs"                                                 "a trans-docos-2-enoyl-[acp]"                                         "Behenoyl-ACPs"                                                      
# [46] "3-oxo-lignoceroyl-ACPs"                                              "R-3-hydroxylignoceroyl-ACPs"                                         "a trans-tetracos-2-enoyl-[acp]"                                     
# [49] "3-Hydroxystearoyl-[acp]"                                             "(2E)-Octadecenoyl-[acp]"                                             "(R)-3-Hydroxyacyl-[acyl-carrier protein]"                           
# [52] "trans-2,3-Dehydroacyl-[acyl-carrier protein]"                        "3-hydroxy-dihomo gamma linolenoyl [acp]"                             "3-hydroxy-eicosatrienoyl [acp]"                                     
# [55] "3-HYDROXY-DOCOSAPENTAENOYL-ACP"                                      "R-3-hydroxypetroselinoyl-ACPs"                                       "3-Oxoacyl-[acyl-carrier protein]"                                   
# [58] "3-Oxostearoyl-[acp]"                                                 "3-oxo-cis-dodec-5-enoyl-[acyl-carrier protein]"                      "(R)-3-hydroxy-cis-dodec-5-enoyl-[acyl-carrier protein]"             
# [61] "3-oxo-cis-myristol-7-eoyl-[acyl-carrier protein]"                    "(R)-3-hydroxy-cis-myristol-7-eoyl-[acyl-carrier protein]"            "3-oxo-cis-palm-9-eoyl-[acyl-carrier protein]"                       
# [64] "(R)-3-hydroxy-cis-palm-9-eoyl-[acyl-carrier protein]"                "3-Oxooctadecanoyl-[acyl-carrier protein]"                            "3-oxo-cis-vacc-11-enoyl-[acyl-carrier protein]"                     
# [67] "(R)-3-hydroxy-cis-vacc-11-enoyl-[acyl-carrier protein]"              "3-oxo-cis-D7-tetradecenoyl-ACPs"                                     "3-hydroxy-cis-D7-tetraecenoyl-ACPs"                                 
# [70] "3-oxo-cis-D9-hexadecenoyl-ACPs"                                      "3-hydroxy-cis-D9-hexaecenoyl-ACPs"                                   "3-Ketoglutaryl-ACP-methyl-ester"                                    
# [73] "3-Hydroxyglutaryl-ACP-methyl-ester"                                  "3-Ketopimeloyl-ACP-methyl-esters"                                    "3-hydroxypimeloyl-ACP-methyl-esters"                                
# [76] "3-oxo-dihomo gamma linolenoyl-[acp]"                                 "3-oxo-eicosatrienoyl [acp]"                                          "3-OXO-EICOSAPENTAENOYL-ACP"                                         
# [79] "Beta-3-hydroxybutyryl-ACPs"                                          "3-Hydroxy-octanoyl-ACPs"                                             "3-oxooctanoyl-acp"                                                  
# [82] "R-3-Hydroxypalmitoyl-ACPs"                                           "3-oxohexadecanoyl-acp"                                               "3-oxo-petroselinoyl-ACPs"                                           
# [85] "R-3-hydroxy-cis-vaccenoyl-ACPs"                                      "3-oxo-cis-vaccenoyl-ACPs"                                            "R-3-hydroxystearoyl-ACPs"                                           
# [88] "a 3-oxo-cis-delta5-dodecenoyl-[acp]"                                 "a 3-hydroxy cis delta5-dodecenoyl-[acp]"                             "a cis,cis-delta13,31-3-oxo-C50:2-[acp]"                             
# [91] "a cis,cis-delta13,31-3-hydroxyC50:2-[acp]"                           "a cis-delta19-3-oxo-C38:1-[acp]"                                     "a cis-delta19-3-hydroxyC38:1-[acp]"                                 
# [94] "a cis,cis-delta15,33-3-oxo-C52:2-[acp]"                              "a cis,cis-delta15,33-3-hydroxyC52:2-[acp]"                           "R-3-hydroxybehenoyl-ACPs"                                           
# [97] "a cis-delta11-3-oxo-C30:1-[acp]"                                     "a cis-delta11-3-hydroxyC30:1-[acp]"                                  "a cis,cis-delta5,23-3-oxo-C42:2-[acp]"                              
# [100] "a cis,cis-delta5,23-3-hydroxyC42:2-[acp]"                            "a cis,cis-delta17,35-3-oxo-C54:2-[acp]"                              "a cis,cis-delta17,35-3-hydroxyC54:2-[acp]"                          
# [103] "a cis,cis-delta13,25-3-oxo-C44:2-[acp]"                              "a cis,cis-delta13,25-3-hydroxyC44:2-[acp]"                           "a cis-delta9-3-oxo-C28:1-[acp]"                                     
# [106] "a cis-delta9-3-hydroxyC28:1-[acp]"                                   "a cis,cis-delta9,21-3-oxo-C40:2-[acp]"                               "a cis,cis-delta9,21-3-hydroxyC40:2-[acp]"                           
# [109] "a cis-delta13-3-oxo-C32:1-[acp]"                                     "a cis-delta13-3-hydroxyC32:1-[acp]"                                  "a cis,cis-delta7,19-3-oxo-C38:2-[acp]"                              
# [112] "a cis,cis-delta7,19-3-hydroxyC38:2-[acp]"                            "a cis-delta21-3-oxo-C40:1-[acp]"                                     "a cis-delta21-3-hydroxyC40:1-[acp]"                                 
# [115] "a cis-delta15-3-oxo-C34:1-[acp]"                                     "a cis-delta15-3-hydroxyC34:1-[acp]"                                  "a cis-delta7-3-oxo-C26:1-[acp]"                                     
# [118] "a cis-delta7-3-hydroxyC26:1-[acp]"                                   "a cis,cis-delta7,25-3-oxo-C44:2-[acp]"                               "a cis,cis-delta7,25-3-hydroxyC44:2-[acp]"                           
# [121] "a cis,cis-delta5,17-3-oxo-C36:2-[acp]"                               "a cis,cis-delta5,17-3-hydroxyC36:2[acp]"                             "a cis,cis-delta19,37-3-oxo-C56:2-[acp]"                             
# [124] "a cis,cis-delta19,37-3-hydroxy-C56:2-[acp]"                          "a cis,cis-delta15,27-3-oxo-C46:2-[acp]"                              "a cis,cis-delta15,27-3-hydroxyC46:2-[acp]"                          
# [127] "a cis,cis-delta9,27-3-oxo-C46:2-[acp]"                               "a cis,cis-delta9,27-3-hydroxyC46:2-[acp]"                            "a cis,cis-delta21,39-3-oxo-C58:2-[acp]"                             
# [130] "a cis,cis-delta21,39-3-hydroxyC58:2-[acp]"                           "a cis,cis-delta11,23-3-oxo-C42:2-[acp]"                              "a cis,cis-delta11,23-3-hydroxyC42:2-[acp]"                          
# [133] "a cis,cis-delta17,29-3-oxo-C48:2-[acp]"                              "a cis,cis-delta17,29-3-hydroxyC48:2-[acp]"                           "a cis-delta5-3-oxo-C24:1-[acp]"                                     
# [136] "a cis-delta5-3-hydroxyC24:1-[acp]"                                   "a cis,cis-delta11,29-3-oxo-C48:2-[acp]"                              "a cis,cis-delta11,29-3-hydroxyC48:2-[acp]"                          
# [139] "a cis,cis-delta19,31-3-oxo-C50:2-[acp]"                              "a cis,cis-delta19,31-3-hydroxyC50:2-[acp]"                           "(R)-3-Hydroxydecanoyl-[acyl-carrier protein]"                       
# [142] "Octanoyl-ACP"                                                        "Lipoyl-ACP"                                                          "OCTANOYL-ACP"                                                       
# [145] "GABA-ACP"                                                            "apo-ACP"                                                             "a cis,cis-delta19-37-C56:2-[acp]"                                   
# [148] "a cis-delta19-37-hydroxy-38-methyl-C57:1-[acp]"                      "a cis-methoxy-C59-meroacyl-[acp]"                                    "a trans-delta18-37-hydroxy-38-methyl-C58:1-[acp]"                   
# [151] "a cis,cis-delta21,39-C58:2-[acp]"                                    "a cis-delta21-39-hydroxy-40-methyl-C59:1-[acp]"                      "a cis-keto-C60-meroacyl-[acp]"                                      
# [154] "a trans-delta20-39-hydroxy-40-methyl-C60:1-[acp]"                    "a trans-keto-C61-meroacyl-[acp]"                                     "a cis,cis-delta19,31-C50:2-[acp]"                                   
# [157] "a C52-alpha-meroacyl-[acp]"                                          "cis-dec-3-enoyl-[acyl-carrier protein]"                              "cis-dodec-5-enoyl-[acyl-carrier protein]"                           
# [160] "Iso-C7:0 ACP"                                                        "Iso-C13:0 ACP"                                                       "Iso-C6:0 ACP"                                                       
# [163] "Iso-C14:0 ACP"                                                       "iso-C15:0 ACP"                                                       "Fatty Acid (n-C5:0 ACP)"                                            
# [166] "pentadecanoyl-ACP (n-C15:0ACP)"                                      "Iso-C16:0 ACP"                                                       "Iso-C17:0 ACP"                                                      
# [169] "heptadecanoyl-ACP (n-C17:0ACP)"                                      "heptadecenoyl ACP (C17:1ACP)"                                        "Cis-Delta5-dodecenoyl-ACPs"                                         
# [172] "Cis-Delta7-tetradecenoyl-ACPs"                                       "Glutaryl-ACP-methyl-esters"                                          "Butyryl-ACP"                                                        
# [175] "3-Oxohexanoyl-[acp]"                                                 "Hexanoyl-ACP"                                                        "3-oxodecanoyl-acp"                                                  
# [178] "Decanoyl-ACP"                                                        "3-oxododecanoyl-acp"                                                 "3-oxotetradecanoyl-acp"                                             
# [181] "a cis-delta17-C36:1-[acp]"                                           "Holo-ACP-Synthases"                                                  "a cis-delta19-C38:1-[acp]"                                          
# [184] "a cis-delta9-C28:1-[acp]"                                            "a cis-delta7-C26:1-[acp]"                                            "a cis-delta11-C30:1-[acp]"                                          
# [187] "a cis-docos-3-enoyl-[acp]"                                           "a cis-delta5-C24:1-[acp]"                                            "a cis-delta13-C32:1-[acp]"                                          
# [190] "a cis-delta15-C34:1-[acp]"                                           "a cis-delta17-3-oxo-C36:1-[acp]"                                     "Acetyl-ACP"                                                         
# [193] "All-ACPs"                                                            "a trans-methoxy-C60-meroacyl-[acp]"                                  "D-3-Hydroxyhexanoyl-[acp]"                                          
# [196] "(R)-3-Hydroxybutanoyl-[acyl-carrier protein]"                        "D-3-Hydroxydodecanoyl-[acp]"                                         "(R)-3-Hydroxyoctanoyl-[acyl-carrier protein]"                       
# [199] "4-methyl-3-oxo-hexanoyl-ACP"                                         "4-methyl-3-hydroxy-hexanoyl-ACP"                                     "6-methyl-3-oxo-octanoyl-ACP"                                        
# [202] "6-methyl-3-hydroxy-octanoyl-ACP"                                     "8-methyl-3-oxo-decanoyl-ACP"                                         "8-methyl-3-hydroxy-decanoyl-ACP"                                    
# [205] "10-methyl-3-oxo-dodecanoyl-ACP"                                      "10-methyl-3-hydroxy-dodecanoyl-ACP"                                  "12-methyl-3-oxo-tetra-decanoyl-ACP"                                 
# [208] "12-methyl-3-hydroxy-tetra-decanoyl-ACP"                              "14-methyl-3-oxo-hexa-decanoyl-ACP"                                   "14-methyl-3-hydroxy-hexa-decanoyl-ACP"                              
# [211] "5-methyl-3-oxo-hexanoyl-ACP"                                         "5-methyl-3-hydroxy-hexanoyl-ACP"                                     "7-methyl-3-oxo-octanoyl-ACP"                                        
# [214] "7-methyl-3-hydroxy-octanoyl-ACP"                                     "9-methyl-3-oxo-decanoyl-ACP"                                         "9-methyl-3-hydroxy-decanoyl-ACP"                                    
# [217] "11-methyl-3-oxo-dodecanoyl-ACP"                                      "11-methyl-3-hydroxy-dodecanoyl-ACP"                                  "13-methyl-3-oxo-tetra-decanoyl-ACP"                                 
# [220] "13-methyl-3-hydroxy-tetra-decanoyl-ACP"                              "15-methyl-3-oxo-hexa-decanoyl-ACP"                                   "15-methyl-3-hydroxy-hexa-decanoyl-ACP"                              
# [223] "4-methyl-3-oxo-pentanoyl-ACP"                                        "4-methyl-3-hydroxy-pentanoyl-ACP"                                    "6-methyl-3-oxo-heptanoyl-ACP"                                       
# [226] "6-methyl-3-hydroxy-heptanoyl-ACP"                                    "8-methyl-3-oxo-nonanoyl-ACP"                                         "8-methyl-3-hydroxy-nonanoyl-ACP"                                    
# [229] "10-methyl-3-oxo-undecanoyl-ACP"                                      "10-methyl-3-hydroxy-undecanoyl-ACP"                                  "12-methyl-3-oxo-tridecanoyl-ACP"                                    
# [232] "12-methyl-3-hydroxy-tridecanoyl-ACP"                                 "14-methyl-3-oxo-pentadecanoyl-ACP"                                   "14-methyl-3-hydroxy-pentadecanoyl-ACP"                              
# [235] "3-Oxooctodecanoyl-ACP"                                               "3-Hydroxyoctodecanoyl-ACP"                                           "Palmitoleoyl-ACPs"                                                  
# [238] "a cis,cis-delta11,29-C48:2-[acp]"                                    "a cis,cis-delta13,31-C50:2-[acp]"                                    "a cis,cis-delta7,19-C38:2-[acp]"                                    
# [241] "a cis,cis-delta3,15-C34:2-[acp]"                                     "a cis,cis-delta11,23-C42:2-[acp]"                                    "a cis,cis-delta5,23-C42:2-[acp]"                                    
# [244] "a cis,cis-delta17,35-C54:2-[acp]"                                    "a cis,cis-delta13,25-C44:2-[acp]"                                    "a cis,cis-delta7,25-C44:2-[acp]"                                    
# [247] "a cis,cis-delta9,21-C40:2-[acp]"                                     "a cis,cis-delta15,27-C46:2-[acp]"                                    "a cis,cis-delta3,21-C40:2-[acp]"                                    
# [250] "a cis,cis-delta9,27-C46:2-[acp]"                                     "a cis,cis-delta17,29-C48:2-[acp]"                                    "a cis,cis-delta5,17-C36:2-[acp]"                                    
# [253] "a cis,cis-delta15,33-C52:2-[acp]"                                    "(Z)-hexadec-11-enoyl-ACP"                                            "(Z)-3-oxooctadec-13-enoyl-ACP"                                      
# [256] "(2E)-Tetradecenoyl-[acp]"                                            "(2E)-Octenoyl-[acp]"                                                 "hexadecanoyl-acp"                                                   
# [259] "(2E)-Hexadecenoyl-[acp]"                                             "12-methyl-tetra-decanoyl-ACP"                                        "14-methyl-hexa-decanoyl-ACP"                                        
# [262] "13-methyl-tetra-decanoyl-ACP"                                        "15-methyl-hexa-decanoyl-ACP"                                         "12-methyl-tridecanoyl-ACP"                                          
# [265] "14-methyl-pentadecanoyl-ACP"                                         "Octodecanoyl-ACP"                                                    "2-methylbutyryl-ACP"                                                
# [268] "isovaleryl-ACP"                                                      "isobutyryl-ACP"                                                      "Delta6-hexadecenoyl-ACPs"                                           
# [271] "Delta4-hexadecenoyl-ACPs"                                            "trans-3-cis-5-dodecenoyl-[acyl-carrier protein]"                     "trans-3-cis-7-myristoleoyl-[acyl-carrier protein]"                  
# [274] "trans-3-cis-9-palmitoleoyl-[acyl-carrier protein]"                   "trans-octadec-2-enoyl-[acyl-carrier protein]"                        "trans-3-cis-11-vacceoyl-[acyl-carrier protein]"                     
# [277] "Crotonyl-ACPs"                                                       "Trans-D2-decenoyl-ACPs"                                              "Cis-vaccenoyl-ACPs"                                                 
# [280] "a cis-vaccen-2-enoyl-[acp]"                                          "Oleoyl-ACPs"                                                         "2E-9Z-octadeca-2-9-dienoyl-ACPs"                                    
# [283] "Malonyl-acp-methyl-ester"                                            "Pimelyl-[acyl-carrier protein] methyl ester"                         "pentadecenoyl-ACP (C15:1ACP)"                                       
# [286] "gamma-L-Glutamyl-4-aminobutyryl-ACP"                                 "hexadecenoyl-[acyl-carrier protein]"                                 "2-Hexadecenoyl-[acyl-carrier protein]"                              
# [289] "Octadecynoyl-ACP"                                                    "Trans-D3-cis-D7-tetradecenoyl-ACPs"                                  "Trans-D3-cis-D9-hexadecenoyl-ACPs"                                  
# [292] "Enoylglutaryl-ACP-methyl-esters"                                     "Enoylpimeloyl-ACP-methyl-esters"                                     "Trans-D3-cis-D5-dodecenoyl-ACPs"                                    
# [295] "cis-3-Decenoyl-[acyl-carrier protein]"                               "Cis-delta-3-decenoyl-ACPs"                            


unique(df.out$f__in[sel2[-subsel]]) # 223 enzyme / sub-functions


length(other_ACP_names)

# other mono methyl ACPs?
sel.other_mm_acps <- grep(pattern = "methyl", x = other_ACP_names)
other_ACP_names[sel.other_mm_acps]
# [1] "3-Ketoglutaryl-ACP-methyl-ester"                  "3-Hydroxyglutaryl-ACP-methyl-ester"               "3-Ketopimeloyl-ACP-methyl-esters"                 "3-hydroxypimeloyl-ACP-methyl-esters"             
# [5] "a cis-delta19-37-hydroxy-38-methyl-C57:1-[acp]"   "a trans-delta18-37-hydroxy-38-methyl-C58:1-[acp]" "a cis-delta21-39-hydroxy-40-methyl-C59:1-[acp]"   "a trans-delta20-39-hydroxy-40-methyl-C60:1-[acp]"
# [9] "Glutaryl-ACP-methyl-esters"                       "4-methyl-3-oxo-hexanoyl-ACP"                      "4-methyl-3-hydroxy-hexanoyl-ACP"                  "6-methyl-3-oxo-octanoyl-ACP"                     
# [13] "6-methyl-3-hydroxy-octanoyl-ACP"                  "8-methyl-3-oxo-decanoyl-ACP"                      "8-methyl-3-hydroxy-decanoyl-ACP"                  "10-methyl-3-oxo-dodecanoyl-ACP"                  
# [17] "10-methyl-3-hydroxy-dodecanoyl-ACP"               "12-methyl-3-oxo-tetra-decanoyl-ACP"               "12-methyl-3-hydroxy-tetra-decanoyl-ACP"           "14-methyl-3-oxo-hexa-decanoyl-ACP"               
# [21] "14-methyl-3-hydroxy-hexa-decanoyl-ACP"            "5-methyl-3-oxo-hexanoyl-ACP"                      "5-methyl-3-hydroxy-hexanoyl-ACP"                  "7-methyl-3-oxo-octanoyl-ACP"                     
# [25] "7-methyl-3-hydroxy-octanoyl-ACP"                  "9-methyl-3-oxo-decanoyl-ACP"                      "9-methyl-3-hydroxy-decanoyl-ACP"                  "11-methyl-3-oxo-dodecanoyl-ACP"                  
# [29] "11-methyl-3-hydroxy-dodecanoyl-ACP"               "13-methyl-3-oxo-tetra-decanoyl-ACP"               "13-methyl-3-hydroxy-tetra-decanoyl-ACP"           "15-methyl-3-oxo-hexa-decanoyl-ACP"               
# [33] "15-methyl-3-hydroxy-hexa-decanoyl-ACP"            "4-methyl-3-oxo-pentanoyl-ACP"                     "4-methyl-3-hydroxy-pentanoyl-ACP"                 "6-methyl-3-oxo-heptanoyl-ACP"                    
# [37] "6-methyl-3-hydroxy-heptanoyl-ACP"                 "8-methyl-3-oxo-nonanoyl-ACP"                      "8-methyl-3-hydroxy-nonanoyl-ACP"                  "10-methyl-3-oxo-undecanoyl-ACP"                  
# [41] "10-methyl-3-hydroxy-undecanoyl-ACP"               "12-methyl-3-oxo-tridecanoyl-ACP"                  "12-methyl-3-hydroxy-tridecanoyl-ACP"              "14-methyl-3-oxo-pentadecanoyl-ACP"               
# [45] "14-methyl-3-hydroxy-pentadecanoyl-ACP"            "12-methyl-tetra-decanoyl-ACP"                     "14-methyl-hexa-decanoyl-ACP"                      "13-methyl-tetra-decanoyl-ACP"                    
# [49] "15-methyl-hexa-decanoyl-ACP"                      "12-methyl-tridecanoyl-ACP"                        "14-methyl-pentadecanoyl-ACP"                      "2-methylbutyryl-ACP"                             
# [53] "Malonyl-acp-methyl-ester"                         "Pimelyl-[acyl-carrier protein] methyl ester"      "Enoylglutaryl-ACP-methyl-esters"                  "Enoylpimeloyl-ACP-methyl-esters" 

# contain additional hydroxyl groups, carbonyl (oxo), 

sel.methyl <- grep(pattern = "methyl", x = other_ACP_names) # 56
sel.enoyl <- grep(pattern = "enoyl", x = other_ACP_names) # 58
sel.match <- which(sel.methyl %in% sel.enoyl) # Empty <<< i.e., no methyl unsaturated bcfa !

sel.oxo <- grep(pattern = "-oxo-", x = other_ACP_names) # 60
length(which(sel.methyl %in% sel.oxo)) # 18 <<<
sel.hydroxy <- grep(pattern = "-hydroxy-", x = other_ACP_names) # 33
length(which(sel.methyl %in% sel.hydroxy)) # 22 <<<

sel.esters <- grep(pattern = "ester", x = other_ACP_names) # 9
length(which(sel.methyl %in% sel.esters)) # 9 <<<

sel.anoyl <- grep(pattern = "anoyl", x = other_ACP_names) # 75
length(which(sel.methyl %in% sel.anoyl & !sel.methyl %in% sel.oxo & !sel.methyl %in% sel.hydroxy)) # 6 <<<
subsel <- which(sel.methyl %in% sel.anoyl & !sel.methyl %in% sel.oxo & !sel.methyl %in% sel.hydroxy) # 6
other_ACP_names[sel.methyl[subsel]]
# "12-methyl-tetra-decanoyl-ACP" "14-methyl-hexa-decanoyl-ACP"  "13-methyl-tetra-decanoyl-ACP" "15-methyl-hexa-decanoyl-ACP"  "12-methyl-tridecanoyl-ACP"    "14-methyl-pentadecanoyl-ACP" 
# 14-chain + methyl(12).         16-chain + methyl(14).         14-chain + methyl(13).          16-chain + methyl(15).         13-chain + methyl(12).        15-chain + methyl(14)

0 + 18 + 22 + 9 + 6 + 1 # 56

# ZERO monomethyl unsaturated (enoyl) bcfa + 18 oxo (carbonyl-containing) + 22 hydroxy (hydroxyl-containing) + 9 esters + 6 saturated (long, 13 to 16 chain length) + 1 x 2-methylbutyryl-ACP (saturated)






## T2D-SWE

# link compounds to functions
df.out <-readRDS( "df.out--tidy-compounds_indiv-cpp3d--forslund-t2d-swe-hostremoval-pawsey.RDS" )
str(df.out)
# 'data.frame':	507391 obs. of  9 variables:
# $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_2" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" ...
# $ rxn_id             : chr  "rxn25279" "rxn25279" "rxn25279" "rxn03060" ...
# $ cpd_id             : chr  "cpd24620" "cpd00001" "cpd25681" "cpd01501" ...
# $ cpd_name           : chr  "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" "2-methyl-trans-aconitate" "2-Methylcitrate" ...
# $ cpd_form           : chr  "C7H7O7" "H2O" "C7H5O6" "C7H7O7" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.3333 0.3333 0.3333 0.0556 0.0556 ...

phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval.RDS")
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 17912 taxa and 76 samples ]
# sample_data() Sample Data:       [ 76 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 17912 taxa by 4 taxonomic ranks ]
df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 17912     4
table(phy@sam_data$Status)
# ND CTRL T2D metformin- 
#   43             33 

length(row.names(phy@tax_table)) # 17912
head(row.names(phy@tax_table)) # "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"

dim(df.out) # 507391      9
sel <- which(df.out$superfocus_fxn %in% row.names(phy@tax_table)) # 507391


length(acp35) # 35

sel1 <- which(df.out$cpd_id %in% acp35) # 72

df.out[sel1, ]

# ACP-35 names
unique(df.out$cpd_name[sel1])
# [1] "(2E)-Decenoyl-[acp]"                   "But-2-enoyl-[acyl-carrier protein]"    "(2E)-Dodecenoyl-[acp]"                 "(2E)-Hexenoyl-[acp]"                   "4-methyl-trans-hex-2-enoyl-ACP"       
# [6] "4-methyl-hexanoyl-ACP"                 "6-methyl-trans-oct-2-enoyl-ACP"        "6-methyl-octanoyl-ACP"                 "8-methyl-trans-dec-2-enoyl-ACP"        "8-methyl-decanoyl-ACP"                
# [11] "10-methyl-trans-dodec-2-enoyl-ACP"     "10-methyl-dodecanoyl-ACP"              "12-methyl-trans-tetra-dec-2-enoyl-ACP" "14-methyl-trans-hexa-dec-2-enoyl-ACP"  "5-methyl-trans-hex-2-enoyl-ACP"       
# [16] "5-methyl-hexanoyl-ACP"                 "7-methyl-trans-oct-2-enoyl-ACP"        "7-methyl-octanoyl-ACP"                 "9-methyl-trans-dec-2-enoyl-ACP"        "9-methyl-decanoyl-ACP"                
# [21] "11-methyl-trans-dodec-2-enoyl-ACP"     "11-methyl-dodecanoyl-ACP"              "13-methyl-trans-tetra-dec-2-enoyl-ACP" "15-methyl-trans-hexa-dec-2-enoyl-ACP"  "4-methyl-trans-pent-2-enoyl-ACP"      
# [26] "4-methyl-pentanoyl-ACP"                "6-methyl-trans-hept-2-enoyl-ACP"       "6-methyl-heptanoyl-ACP"                "8-methyl-trans-non-2-enoyl-ACP"        "8-methyl-nonanoyl-ACP"                
# [31] "10-methyl-trans-undec-2-enoyl-ACP"     "10-methyl-undecanoyl-ACP"              "12-methyl-trans-tridec-2-enoyl-ACP"    "14-methyl-trans-pentadec-2-enoyl-ACP"  "trans-Octodec-2-enoyl-ACP"


## what about other ACPs?

#sel2 <- grep(pattern = "-[acyl-carrier protein]|-[acp]|-ACP", x = df.out$cpd_name)
sel2 <- grep(pattern = "acyl-carrier protein|acp|ACP", x = df.out$cpd_name) # qty 8357

subsel <- which( sel2 %in% sel1 ) # 72

length( df.out$cpd_name[sel2[-subsel]] ) # 8285
length( unique(df.out$cpd_name[sel2[-subsel]]) ) # 289

( other_ACP_names <- unique(df.out$cpd_name[sel2[-subsel]]) )
# [1] "Pimelyl-[acyl-carrier protein]"                                      "ACP"                                                                 "18-Carbamoyl-3,5,7,9,11,13,15,17-octaoxo-octadecanoyl-[acp]"        
# [4] "18-Carbamoyl-9-hydroxy-3,5,7,11,13,15,17-octaoxo-octadecanoyl-[acp]" "dihomo gamma linolenoyl-2-enoyl [acp]"                               "a dihomo gamma linolenoyl [acp]"                                    
# [7] "eicosatrienoyl-2-enoyl [acp]"                                        "EICOSATRIENOYL-ACP"                                                  "DOCOSAPENTAENOYL-2-ENOYL-ACP"                                       
# [10] "DOCOSAPENTAENOYL-ACP"                                                "Petroselinoyl-ACPs"                                                  "Petrosel-2-enoyl-ACPs"                                              
# [13] "(R)-3-Hydroxybutanoyl-[acyl-carrier protein]"                        "Acetoacetyl-ACP"                                                     "(R)-3-Hydroxydecanoyl-[acyl-carrier protein]"                       
# [16] "3-oxodecanoyl-acp"                                                   "(R)-3-Hydroxyoctanoyl-[acyl-carrier protein]"                        "3-oxooctanoyl-acp"                                                  
# [19] "3-oxohexadecanoyl-acp"                                               "Myristoyl-ACP"                                                       "3,5,7,9,11,13,15-Heptaoxo-hexadecanoyl-[acp]"                       
# [22] "3,5,7,9,11,13,15-Heptaoxo-octadecanoyl-[acp]"                        "3,5,7,9,11,13,15,17,19-Nonaoxo-eicosanoyl-[acp]"                     "3,5,7,9,11,13,15,17,19-Nonaoxo-henicosanoyl-[acp]"                  
# [25] "3,5,7,9,11,13,15,17,19-Nonaoxo-docosanoyl-[acp]"                     "Myristoyl-ACPs"                                                      "Lignoceroyl-ACPs"                                                   
# [28] "3-oxo-cerotoyl-ACPs"                                                 "R-3-hydroxycerotoyl-ACPs"                                            "Trans-D2-hexacos-2-enoyl-ACPs"                                      
# [31] "Cerotoyl-ACPs"                                                       "Palmitoyl-ACP"                                                       "Heptadecanoyl-ACPs"                                                 
# [34] "Long-Chain-Acyl-ACPs"                                                "Stearoyl-[acyl-carrier protein]"                                     "Polyketide-ACP-Proteins"                                            
# [37] "Hepta-oxo-hexadecanoyl-ACPs"                                         "R-3-hydroxyarachidoyl-ACPs"                                          "a trans-eicos-2-enoyl-[acp]"                                        
# [40] "3-oxo-arachidoyl-ACPs"                                               "Arachidoyl-ACPs"                                                     "3-oxo-behenoyl-ACPs"                                                
# [43] "a trans-docos-2-enoyl-[acp]"                                         "Behenoyl-ACPs"                                                       "3-oxo-lignoceroyl-ACPs"                                             
# [46] "R-3-hydroxylignoceroyl-ACPs"                                         "a trans-tetracos-2-enoyl-[acp]"                                      "(R)-3-Hydroxyacyl-[acyl-carrier protein]"                           
# [49] "Octanoyl-ACP"                                                        "Lipoyl-ACP"                                                          "OCTANOYL-ACP"                                                       
# [52] "3-Hydroxystearoyl-[acp]"                                             "(2E)-Octadecenoyl-[acp]"                                             "trans-2,3-Dehydroacyl-[acyl-carrier protein]"                       
# [55] "3-hydroxy-dihomo gamma linolenoyl [acp]"                             "3-hydroxy-eicosatrienoyl [acp]"                                      "3-HYDROXY-DOCOSAPENTAENOYL-ACP"                                     
# [58] "R-3-hydroxypetroselinoyl-ACPs"                                       "a cis,cis-delta19-37-C56:2-[acp]"                                    "a cis-delta19-37-hydroxy-38-methyl-C57:1-[acp]"                     
# [61] "a cis-methoxy-C59-meroacyl-[acp]"                                    "a trans-delta18-37-hydroxy-38-methyl-C58:1-[acp]"                    "a cis,cis-delta21,39-C58:2-[acp]"                                   
# [64] "a cis-delta21-39-hydroxy-40-methyl-C59:1-[acp]"                      "a cis-keto-C60-meroacyl-[acp]"                                       "a trans-delta20-39-hydroxy-40-methyl-C60:1-[acp]"                   
# [67] "a trans-keto-C61-meroacyl-[acp]"                                     "a cis,cis-delta19,31-C50:2-[acp]"                                    "a C52-alpha-meroacyl-[acp]"                                         
# [70] "Dodecanoyl-ACP"                                                      "Palmitoleoyl-ACPs"                                                   "Tetradecenoyl-ACP"                                                  
# [73] "Hexadecenoyl-ACP"                                                    "Octadecanoyl-ACP"                                                    "Octadecenoyl-ACP"                                                   
# [76] "Acyl-[acyl-carrier protein]"                                         "Saturated-Fatty-Acyl-ACPs"                                           "D-3-Hydroxyhexanoyl-[acp]"                                          
# [79] "3-Oxohexanoyl-[acp]"                                                 "D-3-Hydroxydodecanoyl-[acp]"                                         "3-oxododecanoyl-acp"                                                
# [82] "3-oxotetradecanoyl-acp"                                              "4-methyl-3-oxo-hexanoyl-ACP"                                         "4-methyl-3-hydroxy-hexanoyl-ACP"                                    
# [85] "6-methyl-3-oxo-octanoyl-ACP"                                         "6-methyl-3-hydroxy-octanoyl-ACP"                                     "8-methyl-3-oxo-decanoyl-ACP"                                        
# [88] "8-methyl-3-hydroxy-decanoyl-ACP"                                     "10-methyl-3-oxo-dodecanoyl-ACP"                                      "10-methyl-3-hydroxy-dodecanoyl-ACP"                                 
# [91] "12-methyl-3-oxo-tetra-decanoyl-ACP"                                  "12-methyl-3-hydroxy-tetra-decanoyl-ACP"                              "14-methyl-3-oxo-hexa-decanoyl-ACP"                                  
# [94] "14-methyl-3-hydroxy-hexa-decanoyl-ACP"                               "5-methyl-3-oxo-hexanoyl-ACP"                                         "5-methyl-3-hydroxy-hexanoyl-ACP"                                    
# [97] "7-methyl-3-oxo-octanoyl-ACP"                                         "7-methyl-3-hydroxy-octanoyl-ACP"                                     "9-methyl-3-oxo-decanoyl-ACP"                                        
# [100] "9-methyl-3-hydroxy-decanoyl-ACP"                                     "11-methyl-3-oxo-dodecanoyl-ACP"                                      "11-methyl-3-hydroxy-dodecanoyl-ACP"                                 
# [103] "13-methyl-3-oxo-tetra-decanoyl-ACP"                                  "13-methyl-3-hydroxy-tetra-decanoyl-ACP"                              "15-methyl-3-oxo-hexa-decanoyl-ACP"                                  
# [106] "15-methyl-3-hydroxy-hexa-decanoyl-ACP"                               "4-methyl-3-oxo-pentanoyl-ACP"                                        "4-methyl-3-hydroxy-pentanoyl-ACP"                                   
# [109] "6-methyl-3-oxo-heptanoyl-ACP"                                        "6-methyl-3-hydroxy-heptanoyl-ACP"                                    "8-methyl-3-oxo-nonanoyl-ACP"                                        
# [112] "8-methyl-3-hydroxy-nonanoyl-ACP"                                     "10-methyl-3-oxo-undecanoyl-ACP"                                      "10-methyl-3-hydroxy-undecanoyl-ACP"                                 
# [115] "12-methyl-3-oxo-tridecanoyl-ACP"                                     "12-methyl-3-hydroxy-tridecanoyl-ACP"                                 "14-methyl-3-oxo-pentadecanoyl-ACP"                                  
# [118] "14-methyl-3-hydroxy-pentadecanoyl-ACP"                               "3-Oxooctodecanoyl-ACP"                                               "3-Hydroxyoctodecanoyl-ACP"                                          
# [121] "apo-ACP"                                                             "3-Oxoacyl-[acyl-carrier protein]"                                    "3-Oxostearoyl-[acp]"                                                
# [124] "3-oxo-cis-dodec-5-enoyl-[acyl-carrier protein]"                      "(R)-3-hydroxy-cis-dodec-5-enoyl-[acyl-carrier protein]"              "3-oxo-cis-myristol-7-eoyl-[acyl-carrier protein]"                   
# [127] "(R)-3-hydroxy-cis-myristol-7-eoyl-[acyl-carrier protein]"            "3-oxo-cis-palm-9-eoyl-[acyl-carrier protein]"                        "(R)-3-hydroxy-cis-palm-9-eoyl-[acyl-carrier protein]"               
# [130] "3-Oxooctadecanoyl-[acyl-carrier protein]"                            "3-oxo-cis-vacc-11-enoyl-[acyl-carrier protein]"                      "(R)-3-hydroxy-cis-vacc-11-enoyl-[acyl-carrier protein]"             
# [133] "3-oxo-cis-D7-tetradecenoyl-ACPs"                                     "3-hydroxy-cis-D7-tetraecenoyl-ACPs"                                  "3-oxo-cis-D9-hexadecenoyl-ACPs"                                     
# [136] "3-hydroxy-cis-D9-hexaecenoyl-ACPs"                                   "3-Ketoglutaryl-ACP-methyl-ester"                                     "3-Hydroxyglutaryl-ACP-methyl-ester"                                 
# [139] "3-Ketopimeloyl-ACP-methyl-esters"                                    "3-hydroxypimeloyl-ACP-methyl-esters"                                 "3-oxo-dihomo gamma linolenoyl-[acp]"                                
# [142] "3-oxo-eicosatrienoyl [acp]"                                          "3-OXO-EICOSAPENTAENOYL-ACP"                                          "Beta-3-hydroxybutyryl-ACPs"                                         
# [145] "3-Hydroxy-octanoyl-ACPs"                                             "R-3-Hydroxypalmitoyl-ACPs"                                           "3-oxo-petroselinoyl-ACPs"                                           
# [148] "R-3-hydroxy-cis-vaccenoyl-ACPs"                                      "3-oxo-cis-vaccenoyl-ACPs"                                            "R-3-hydroxystearoyl-ACPs"                                           
# [151] "a 3-oxo-cis-delta5-dodecenoyl-[acp]"                                 "a 3-hydroxy cis delta5-dodecenoyl-[acp]"                             "a cis,cis-delta13,31-3-oxo-C50:2-[acp]"                             
# [154] "a cis,cis-delta13,31-3-hydroxyC50:2-[acp]"                           "a cis-delta19-3-oxo-C38:1-[acp]"                                     "a cis-delta19-3-hydroxyC38:1-[acp]"                                 
# [157] "a cis,cis-delta15,33-3-oxo-C52:2-[acp]"                              "a cis,cis-delta15,33-3-hydroxyC52:2-[acp]"                           "R-3-hydroxybehenoyl-ACPs"                                           
# [160] "a cis-delta11-3-oxo-C30:1-[acp]"                                     "a cis-delta11-3-hydroxyC30:1-[acp]"                                  "a cis,cis-delta5,23-3-oxo-C42:2-[acp]"                              
# [163] "a cis,cis-delta5,23-3-hydroxyC42:2-[acp]"                            "a cis,cis-delta17,35-3-oxo-C54:2-[acp]"                              "a cis,cis-delta17,35-3-hydroxyC54:2-[acp]"                          
# [166] "a cis,cis-delta13,25-3-oxo-C44:2-[acp]"                              "a cis,cis-delta13,25-3-hydroxyC44:2-[acp]"                           "a cis-delta9-3-oxo-C28:1-[acp]"                                     
# [169] "a cis-delta9-3-hydroxyC28:1-[acp]"                                   "a cis,cis-delta9,21-3-oxo-C40:2-[acp]"                               "a cis,cis-delta9,21-3-hydroxyC40:2-[acp]"                           
# [172] "a cis-delta13-3-oxo-C32:1-[acp]"                                     "a cis-delta13-3-hydroxyC32:1-[acp]"                                  "a cis,cis-delta7,19-3-oxo-C38:2-[acp]"                              
# [175] "a cis,cis-delta7,19-3-hydroxyC38:2-[acp]"                            "a cis-delta21-3-oxo-C40:1-[acp]"                                     "a cis-delta21-3-hydroxyC40:1-[acp]"                                 
# [178] "a cis-delta15-3-oxo-C34:1-[acp]"                                     "a cis-delta15-3-hydroxyC34:1-[acp]"                                  "a cis-delta7-3-oxo-C26:1-[acp]"                                     
# [181] "a cis-delta7-3-hydroxyC26:1-[acp]"                                   "a cis,cis-delta7,25-3-oxo-C44:2-[acp]"                               "a cis,cis-delta7,25-3-hydroxyC44:2-[acp]"                           
# [184] "a cis,cis-delta5,17-3-oxo-C36:2-[acp]"                               "a cis,cis-delta5,17-3-hydroxyC36:2[acp]"                             "a cis,cis-delta19,37-3-oxo-C56:2-[acp]"                             
# [187] "a cis,cis-delta19,37-3-hydroxy-C56:2-[acp]"                          "a cis,cis-delta15,27-3-oxo-C46:2-[acp]"                              "a cis,cis-delta15,27-3-hydroxyC46:2-[acp]"                          
# [190] "a cis,cis-delta9,27-3-oxo-C46:2-[acp]"                               "a cis,cis-delta9,27-3-hydroxyC46:2-[acp]"                            "a cis,cis-delta21,39-3-oxo-C58:2-[acp]"                             
# [193] "a cis,cis-delta21,39-3-hydroxyC58:2-[acp]"                           "a cis,cis-delta11,23-3-oxo-C42:2-[acp]"                              "a cis,cis-delta11,23-3-hydroxyC42:2-[acp]"                          
# [196] "a cis,cis-delta17,29-3-oxo-C48:2-[acp]"                              "a cis,cis-delta17,29-3-hydroxyC48:2-[acp]"                           "a cis-delta5-3-oxo-C24:1-[acp]"                                     
# [199] "a cis-delta5-3-hydroxyC24:1-[acp]"                                   "a cis,cis-delta11,29-3-oxo-C48:2-[acp]"                              "a cis,cis-delta11,29-3-hydroxyC48:2-[acp]"                          
# [202] "a cis,cis-delta19,31-3-oxo-C50:2-[acp]"                              "a cis,cis-delta19,31-3-hydroxyC50:2-[acp]"                           "Malonyl-acp-methyl-ester"                                           
# [205] "Pimelyl-[acyl-carrier protein] methyl ester"                         "trans-3-cis-5-dodecenoyl-[acyl-carrier protein]"                     "cis-dodec-5-enoyl-[acyl-carrier protein]"                           
# [208] "trans-3-cis-7-myristoleoyl-[acyl-carrier protein]"                   "(2E)-Hexadecenoyl-[acp]"                                             "trans-3-cis-9-palmitoleoyl-[acyl-carrier protein]"                  
# [211] "trans-octadec-2-enoyl-[acyl-carrier protein]"                        "trans-3-cis-11-vacceoyl-[acyl-carrier protein]"                      "Crotonyl-ACPs"                                                      
# [214] "Butyryl-ACP"                                                         "Decanoyl-ACP"                                                        "Trans-D2-decenoyl-ACPs"                                             
# [217] "(2E)-Tetradecenoyl-[acp]"                                            "Cis-vaccenoyl-ACPs"                                                  "a cis-vaccen-2-enoyl-[acp]"                                         
# [220] "Oleoyl-ACPs"                                                         "2E-9Z-octadeca-2-9-dienoyl-ACPs"                                     "2-Hexadecenoyl-[acyl-carrier protein]"                              
# [223] "cis-dec-3-enoyl-[acyl-carrier protein]"                              "Iso-C7:0 ACP"                                                        "Iso-C13:0 ACP"                                                      
# [226] "Iso-C6:0 ACP"                                                        "Iso-C14:0 ACP"                                                       "iso-C15:0 ACP"                                                      
# [229] "Fatty Acid (n-C5:0 ACP)"                                             "pentadecanoyl-ACP (n-C15:0ACP)"                                      "Iso-C16:0 ACP"                                                      
# [232] "Iso-C17:0 ACP"                                                       "heptadecanoyl-ACP (n-C17:0ACP)"                                      "heptadecenoyl ACP (C17:1ACP)"                                       
# [235] "Cis-Delta5-dodecenoyl-ACPs"                                          "Cis-Delta7-tetradecenoyl-ACPs"                                       "Glutaryl-ACP-methyl-esters"                                         
# [238] "Hexanoyl-ACP"                                                        "a cis-delta17-C36:1-[acp]"                                           "Holo-ACP-Synthases"                                                 
# [241] "a cis-delta19-C38:1-[acp]"                                           "a cis-delta9-C28:1-[acp]"                                            "a cis-delta7-C26:1-[acp]"                                           
# [244] "a cis-delta11-C30:1-[acp]"                                           "a cis-docos-3-enoyl-[acp]"                                           "a cis-delta5-C24:1-[acp]"                                           
# [247] "a cis-delta13-C32:1-[acp]"                                           "a cis-delta15-C34:1-[acp]"                                           "a cis-delta17-3-oxo-C36:1-[acp]"                                    
# [250] "Acetyl-ACP"                                                          "All-ACPs"                                                            "(2E)-Octenoyl-[acp]"                                                
# [253] "Trans-D3-cis-D7-tetradecenoyl-ACPs"                                  "Trans-D3-cis-D9-hexadecenoyl-ACPs"                                   "Enoylglutaryl-ACP-methyl-esters"                                    
# [256] "Enoylpimeloyl-ACP-methyl-esters"                                     "Trans-D3-cis-D5-dodecenoyl-ACPs"                                     "cis-3-Decenoyl-[acyl-carrier protein]"                              
# [259] "Cis-delta-3-decenoyl-ACPs"                                           "a cis,cis-delta11,29-C48:2-[acp]"                                    "a cis,cis-delta13,31-C50:2-[acp]"                                   
# [262] "a cis,cis-delta7,19-C38:2-[acp]"                                     "a cis,cis-delta3,15-C34:2-[acp]"                                     "a cis,cis-delta11,23-C42:2-[acp]"                                   
# [265] "a cis,cis-delta5,23-C42:2-[acp]"                                     "a cis,cis-delta17,35-C54:2-[acp]"                                    "a cis,cis-delta13,25-C44:2-[acp]"                                   
# [268] "a cis,cis-delta7,25-C44:2-[acp]"                                     "a cis,cis-delta9,21-C40:2-[acp]"                                     "a cis,cis-delta15,27-C46:2-[acp]"                                   
# [271] "a cis,cis-delta3,21-C40:2-[acp]"                                     "a cis,cis-delta9,27-C46:2-[acp]"                                     "a cis,cis-delta17,29-C48:2-[acp]"                                   
# [274] "a cis,cis-delta5,17-C36:2-[acp]"                                     "a cis,cis-delta15,33-C52:2-[acp]"                                    "(Z)-hexadec-11-enoyl-ACP"                                           
# [277] "(Z)-3-oxooctadec-13-enoyl-ACP"                                       "hexadecanoyl-acp"                                                    "12-methyl-tetra-decanoyl-ACP"                                       
# [280] "14-methyl-hexa-decanoyl-ACP"                                         "13-methyl-tetra-decanoyl-ACP"                                        "15-methyl-hexa-decanoyl-ACP"                                        
# [283] "12-methyl-tridecanoyl-ACP"                                           "14-methyl-pentadecanoyl-ACP"                                         "Octodecanoyl-ACP"                                                   
# [286] "2-methylbutyryl-ACP"                                                 "isovaleryl-ACP"                                                      "isobutyryl-ACP"                                                     
# [289] "a trans-methoxy-C60-meroacyl-[acp]" 


unique(df.out$f__in[sel2[-subsel]]) # 137 enzyme / sub-functions


length(other_ACP_names) # 289

# other methyl containing ACPs?
sel.other_mm_acps <- grep(pattern = "methyl", x = other_ACP_names) # qty 56
other_ACP_names[sel.other_mm_acps]
# [1] "a cis-delta19-37-hydroxy-38-methyl-C57:1-[acp]"   "a trans-delta18-37-hydroxy-38-methyl-C58:1-[acp]" "a cis-delta21-39-hydroxy-40-methyl-C59:1-[acp]"   "a trans-delta20-39-hydroxy-40-methyl-C60:1-[acp]"
# [5] "4-methyl-3-oxo-hexanoyl-ACP"                      "4-methyl-3-hydroxy-hexanoyl-ACP"                  "6-methyl-3-oxo-octanoyl-ACP"                      "6-methyl-3-hydroxy-octanoyl-ACP"                 
# [9] "8-methyl-3-oxo-decanoyl-ACP"                      "8-methyl-3-hydroxy-decanoyl-ACP"                  "10-methyl-3-oxo-dodecanoyl-ACP"                   "10-methyl-3-hydroxy-dodecanoyl-ACP"              
# [13] "12-methyl-3-oxo-tetra-decanoyl-ACP"               "12-methyl-3-hydroxy-tetra-decanoyl-ACP"           "14-methyl-3-oxo-hexa-decanoyl-ACP"                "14-methyl-3-hydroxy-hexa-decanoyl-ACP"           
# [17] "5-methyl-3-oxo-hexanoyl-ACP"                      "5-methyl-3-hydroxy-hexanoyl-ACP"                  "7-methyl-3-oxo-octanoyl-ACP"                      "7-methyl-3-hydroxy-octanoyl-ACP"                 
# [21] "9-methyl-3-oxo-decanoyl-ACP"                      "9-methyl-3-hydroxy-decanoyl-ACP"                  "11-methyl-3-oxo-dodecanoyl-ACP"                   "11-methyl-3-hydroxy-dodecanoyl-ACP"              
# [25] "13-methyl-3-oxo-tetra-decanoyl-ACP"               "13-methyl-3-hydroxy-tetra-decanoyl-ACP"           "15-methyl-3-oxo-hexa-decanoyl-ACP"                "15-methyl-3-hydroxy-hexa-decanoyl-ACP"           
# [29] "4-methyl-3-oxo-pentanoyl-ACP"                     "4-methyl-3-hydroxy-pentanoyl-ACP"                 "6-methyl-3-oxo-heptanoyl-ACP"                     "6-methyl-3-hydroxy-heptanoyl-ACP"                
# [33] "8-methyl-3-oxo-nonanoyl-ACP"                      "8-methyl-3-hydroxy-nonanoyl-ACP"                  "10-methyl-3-oxo-undecanoyl-ACP"                   "10-methyl-3-hydroxy-undecanoyl-ACP"              
# [37] "12-methyl-3-oxo-tridecanoyl-ACP"                  "12-methyl-3-hydroxy-tridecanoyl-ACP"              "14-methyl-3-oxo-pentadecanoyl-ACP"                "14-methyl-3-hydroxy-pentadecanoyl-ACP"           
# [41] "3-Ketoglutaryl-ACP-methyl-ester"                  "3-Hydroxyglutaryl-ACP-methyl-ester"               "3-Ketopimeloyl-ACP-methyl-esters"                 "3-hydroxypimeloyl-ACP-methyl-esters"             
# [45] "Malonyl-acp-methyl-ester"                         "Pimelyl-[acyl-carrier protein] methyl ester"      "Glutaryl-ACP-methyl-esters"                       "Enoylglutaryl-ACP-methyl-esters"                 
# [49] "Enoylpimeloyl-ACP-methyl-esters"                  "12-methyl-tetra-decanoyl-ACP"                     "14-methyl-hexa-decanoyl-ACP"                      "13-methyl-tetra-decanoyl-ACP"                    
# [53] "15-methyl-hexa-decanoyl-ACP"                      "12-methyl-tridecanoyl-ACP"                        "14-methyl-pentadecanoyl-ACP"                      "2-methylbutyryl-ACP" 

# contain additional hydroxyl groups, carbonyl (oxo), 

sel.methyl <- grep(pattern = "methyl", x = other_ACP_names) # 56
sel.enoyl <- grep(pattern = "enoyl", x = other_ACP_names) # 54
sel.match <- which(sel.methyl %in% sel.enoyl) # Empty <<< i.e., no other methyl unsaturated bcfa !

sel.oxo <- grep(pattern = "-oxo-", x = other_ACP_names) # 60
length(which(sel.methyl %in% sel.oxo)) # 18 <<<
sel.hydroxy <- grep(pattern = "-hydroxy-", x = other_ACP_names) # 33
length(which(sel.methyl %in% sel.hydroxy)) # 22 <<<

sel.esters <- grep(pattern = "ester", x = other_ACP_names) # 9
length(which(sel.methyl %in% sel.esters)) # 9 <<<

sel.anoyl <- grep(pattern = "anoyl", x = other_ACP_names) # 75
length(which(sel.methyl %in% sel.anoyl & !sel.methyl %in% sel.oxo & !sel.methyl %in% sel.hydroxy)) # 6 <<<
subsel <- which(sel.methyl %in% sel.anoyl & !sel.methyl %in% sel.oxo & !sel.methyl %in% sel.hydroxy) # 6
other_ACP_names[sel.methyl[subsel]]
# "12-methyl-tetra-decanoyl-ACP" "14-methyl-hexa-decanoyl-ACP"  "13-methyl-tetra-decanoyl-ACP" "15-methyl-hexa-decanoyl-ACP"  "12-methyl-tridecanoyl-ACP"    "14-methyl-pentadecanoyl-ACP" 
# 14-chain + methyl(12).         16-chain + methyl(14).         14-chain + methyl(13).          16-chain + methyl(15).         13-chain + methyl(12).        15-chain + methyl(14)

0 + 18 + 22 + 9 + 6 + 1  # 56

# ZERO monomethyl unsaturated (enoyl) bcfa + 18 oxo (carbonyl-containing) + 22 hydroxy (hydroxyl-containing) + 9 esters + 6 saturated (long, 13 to 16 chain length) + 1 x 2-methylbutyryl-ACP (saturated)







## T2D-CHN

# link compounds to functions
df.out <-readRDS( "df.out--tidy-compounds_indiv-cpp3d--forslund-t2d-chn-hostremoval-pawsey.RDS" )
str(df.out)
# 'data.frame':	554562 obs. of  9 variables:
# $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_3" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "4-hydroxyproline epimerase (EC 5.1.1.8)" "4-hydroxyproline epimerase (EC 5.1.1.8)" "Alanine racemase (EC 5.1.1.1)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id             : chr  "rxn02360" "rxn02360" "rxn00283" "rxn00283" ...
# $ cpd_id             : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ cpd_name           : chr  "trans-4-Hydroxy-L-proline" "cis-4-Hydroxy-D-proline" "L-Alanine" "D-Alanine" ...
# $ cpd_form           : chr  "C5H9NO3" "C5H9NO3" "C3H7NO2" "C3H7NO2" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.167 0.167 0.167 ...

phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal.RDS")
df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 19352     4

table(phy@sam_data$Diagnosis)
# ND CTRL T2D metformin- 
#   52             30 

min(taxa_sums(phy)) # 1.734989e-06


length(acp35) # 35

sel1 <- which(df.out$cpd_id %in% acp35) # 72

df.out[sel1, ]

# ACP-35 names
unique(df.out$cpd_name[sel1])
# [1] "(2E)-Decenoyl-[acp]"                   "But-2-enoyl-[acyl-carrier protein]"    "(2E)-Dodecenoyl-[acp]"                 "(2E)-Hexenoyl-[acp]"                   "4-methyl-trans-hex-2-enoyl-ACP"       
# [6] "4-methyl-hexanoyl-ACP"                 "6-methyl-trans-oct-2-enoyl-ACP"        "6-methyl-octanoyl-ACP"                 "8-methyl-trans-dec-2-enoyl-ACP"        "8-methyl-decanoyl-ACP"                
# [11] "10-methyl-trans-dodec-2-enoyl-ACP"     "10-methyl-dodecanoyl-ACP"              "12-methyl-trans-tetra-dec-2-enoyl-ACP" "14-methyl-trans-hexa-dec-2-enoyl-ACP"  "5-methyl-trans-hex-2-enoyl-ACP"       
# [16] "5-methyl-hexanoyl-ACP"                 "7-methyl-trans-oct-2-enoyl-ACP"        "7-methyl-octanoyl-ACP"                 "9-methyl-trans-dec-2-enoyl-ACP"        "9-methyl-decanoyl-ACP"                
# [21] "11-methyl-trans-dodec-2-enoyl-ACP"     "11-methyl-dodecanoyl-ACP"              "13-methyl-trans-tetra-dec-2-enoyl-ACP" "15-methyl-trans-hexa-dec-2-enoyl-ACP"  "4-methyl-trans-pent-2-enoyl-ACP"      
# [26] "4-methyl-pentanoyl-ACP"                "6-methyl-trans-hept-2-enoyl-ACP"       "6-methyl-heptanoyl-ACP"                "8-methyl-trans-non-2-enoyl-ACP"        "8-methyl-nonanoyl-ACP"                
# [31] "10-methyl-trans-undec-2-enoyl-ACP"     "10-methyl-undecanoyl-ACP"              "12-methyl-trans-tridec-2-enoyl-ACP"    "14-methyl-trans-pentadec-2-enoyl-ACP"  "trans-Octodec-2-enoyl-ACP" 


## what about other ACPs?

#sel2 <- grep(pattern = "-[acyl-carrier protein]|-[acp]|-ACP", x = df.out$cpd_name)
sel2 <- grep(pattern = "acyl-carrier protein|acp|ACP", x = df.out$cpd_name) # qty 9311

subsel <- which( sel2 %in% sel1 ) # 72

length( df.out$cpd_name[sel2[-subsel]] ) # 9239
length( unique(df.out$cpd_name[sel2[-subsel]]) ) # 293

( other_ACP_names <- unique(df.out$cpd_name[sel2[-subsel]]) )
# [1] "Pimelyl-[acyl-carrier protein]"                                      "ACP"                                                                 "18-Carbamoyl-3,5,7,9,11,13,15,17-octaoxo-octadecanoyl-[acp]"        
# [4] "18-Carbamoyl-9-hydroxy-3,5,7,11,13,15,17-octaoxo-octadecanoyl-[acp]" "dihomo gamma linolenoyl-2-enoyl [acp]"                               "a dihomo gamma linolenoyl [acp]"                                    
# [7] "eicosatrienoyl-2-enoyl [acp]"                                        "EICOSATRIENOYL-ACP"                                                  "DOCOSAPENTAENOYL-2-ENOYL-ACP"                                       
# [10] "DOCOSAPENTAENOYL-ACP"                                                "Petroselinoyl-ACPs"                                                  "Petrosel-2-enoyl-ACPs"                                              
# [13] "Myristoyl-ACP"                                                       "3,5,7,9,11,13,15-Heptaoxo-hexadecanoyl-[acp]"                        "3,5,7,9,11,13,15-Heptaoxo-octadecanoyl-[acp]"                       
# [16] "3,5,7,9,11,13,15,17,19-Nonaoxo-eicosanoyl-[acp]"                     "3,5,7,9,11,13,15,17,19-Nonaoxo-henicosanoyl-[acp]"                   "3,5,7,9,11,13,15,17,19-Nonaoxo-docosanoyl-[acp]"                    
# [19] "Myristoyl-ACPs"                                                      "Lignoceroyl-ACPs"                                                    "3-oxo-cerotoyl-ACPs"                                                
# [22] "R-3-hydroxycerotoyl-ACPs"                                            "Trans-D2-hexacos-2-enoyl-ACPs"                                       "Cerotoyl-ACPs"                                                      
# [25] "Acetoacetyl-ACP"                                                     "Palmitoyl-ACP"                                                       "Heptadecanoyl-ACPs"                                                 
# [28] "Long-Chain-Acyl-ACPs"                                                "Stearoyl-[acyl-carrier protein]"                                     "Polyketide-ACP-Proteins"                                            
# [31] "Hepta-oxo-hexadecanoyl-ACPs"                                         "R-3-hydroxyarachidoyl-ACPs"                                          "a trans-eicos-2-enoyl-[acp]"                                        
# [34] "3-oxo-arachidoyl-ACPs"                                               "Arachidoyl-ACPs"                                                     "3-oxo-behenoyl-ACPs"                                                
# [37] "a trans-docos-2-enoyl-[acp]"                                         "Behenoyl-ACPs"                                                       "3-oxo-lignoceroyl-ACPs"                                             
# [40] "R-3-hydroxylignoceroyl-ACPs"                                         "a trans-tetracos-2-enoyl-[acp]"                                      "(R)-3-Hydroxyacyl-[acyl-carrier protein]"                           
# [43] "(R)-3-Hydroxydecanoyl-[acyl-carrier protein]"                        "Octanoyl-ACP"                                                        "Lipoyl-ACP"                                                         
# [46] "OCTANOYL-ACP"                                                        "3-Hydroxystearoyl-[acp]"                                             "(2E)-Octadecenoyl-[acp]"                                            
# [49] "trans-2,3-Dehydroacyl-[acyl-carrier protein]"                        "3-hydroxy-dihomo gamma linolenoyl [acp]"                             "3-hydroxy-eicosatrienoyl [acp]"                                     
# [52] "3-HYDROXY-DOCOSAPENTAENOYL-ACP"                                      "R-3-hydroxypetroselinoyl-ACPs"                                       "a cis,cis-delta19-37-C56:2-[acp]"                                   
# [55] "a cis-delta19-37-hydroxy-38-methyl-C57:1-[acp]"                      "a cis-methoxy-C59-meroacyl-[acp]"                                    "a trans-delta18-37-hydroxy-38-methyl-C58:1-[acp]"                   
# [58] "a cis,cis-delta21,39-C58:2-[acp]"                                    "a cis-delta21-39-hydroxy-40-methyl-C59:1-[acp]"                      "a cis-keto-C60-meroacyl-[acp]"                                      
# [61] "a trans-delta20-39-hydroxy-40-methyl-C60:1-[acp]"                    "a trans-keto-C61-meroacyl-[acp]"                                     "a cis,cis-delta19,31-C50:2-[acp]"                                   
# [64] "a C52-alpha-meroacyl-[acp]"                                          "Dodecanoyl-ACP"                                                      "Palmitoleoyl-ACPs"                                                  
# [67] "Tetradecenoyl-ACP"                                                   "Hexadecenoyl-ACP"                                                    "Octadecanoyl-ACP"                                                   
# [70] "Octadecenoyl-ACP"                                                    "Acyl-[acyl-carrier protein]"                                         "Saturated-Fatty-Acyl-ACPs"                                          
# [73] "3-oxohexadecanoyl-acp"                                               "D-3-Hydroxyhexanoyl-[acp]"                                           "3-Oxohexanoyl-[acp]"                                                
# [76] "3-oxodecanoyl-acp"                                                   "(R)-3-Hydroxybutanoyl-[acyl-carrier protein]"                        "D-3-Hydroxydodecanoyl-[acp]"                                        
# [79] "3-oxododecanoyl-acp"                                                 "(R)-3-Hydroxyoctanoyl-[acyl-carrier protein]"                        "3-oxooctanoyl-acp"                                                  
# [82] "3-oxotetradecanoyl-acp"                                              "4-methyl-3-oxo-hexanoyl-ACP"                                         "4-methyl-3-hydroxy-hexanoyl-ACP"                                    
# [85] "6-methyl-3-oxo-octanoyl-ACP"                                         "6-methyl-3-hydroxy-octanoyl-ACP"                                     "8-methyl-3-oxo-decanoyl-ACP"                                        
# [88] "8-methyl-3-hydroxy-decanoyl-ACP"                                     "10-methyl-3-oxo-dodecanoyl-ACP"                                      "10-methyl-3-hydroxy-dodecanoyl-ACP"                                 
# [91] "12-methyl-3-oxo-tetra-decanoyl-ACP"                                  "12-methyl-3-hydroxy-tetra-decanoyl-ACP"                              "14-methyl-3-oxo-hexa-decanoyl-ACP"                                  
# [94] "14-methyl-3-hydroxy-hexa-decanoyl-ACP"                               "5-methyl-3-oxo-hexanoyl-ACP"                                         "5-methyl-3-hydroxy-hexanoyl-ACP"                                    
# [97] "7-methyl-3-oxo-octanoyl-ACP"                                         "7-methyl-3-hydroxy-octanoyl-ACP"                                     "9-methyl-3-oxo-decanoyl-ACP"                                        
# [100] "9-methyl-3-hydroxy-decanoyl-ACP"                                     "11-methyl-3-oxo-dodecanoyl-ACP"                                      "11-methyl-3-hydroxy-dodecanoyl-ACP"                                 
# [103] "13-methyl-3-oxo-tetra-decanoyl-ACP"                                  "13-methyl-3-hydroxy-tetra-decanoyl-ACP"                              "15-methyl-3-oxo-hexa-decanoyl-ACP"                                  
# [106] "15-methyl-3-hydroxy-hexa-decanoyl-ACP"                               "4-methyl-3-oxo-pentanoyl-ACP"                                        "4-methyl-3-hydroxy-pentanoyl-ACP"                                   
# [109] "6-methyl-3-oxo-heptanoyl-ACP"                                        "6-methyl-3-hydroxy-heptanoyl-ACP"                                    "8-methyl-3-oxo-nonanoyl-ACP"                                        
# [112] "8-methyl-3-hydroxy-nonanoyl-ACP"                                     "10-methyl-3-oxo-undecanoyl-ACP"                                      "10-methyl-3-hydroxy-undecanoyl-ACP"                                 
# [115] "12-methyl-3-oxo-tridecanoyl-ACP"                                     "12-methyl-3-hydroxy-tridecanoyl-ACP"                                 "14-methyl-3-oxo-pentadecanoyl-ACP"                                  
# [118] "14-methyl-3-hydroxy-pentadecanoyl-ACP"                               "3-Oxooctodecanoyl-ACP"                                               "3-Hydroxyoctodecanoyl-ACP"                                          
# [121] "apo-ACP"                                                             "3-Oxoacyl-[acyl-carrier protein]"                                    "3-Oxostearoyl-[acp]"                                                
# [124] "3-oxo-cis-dodec-5-enoyl-[acyl-carrier protein]"                      "(R)-3-hydroxy-cis-dodec-5-enoyl-[acyl-carrier protein]"              "3-oxo-cis-myristol-7-eoyl-[acyl-carrier protein]"                   
# [127] "(R)-3-hydroxy-cis-myristol-7-eoyl-[acyl-carrier protein]"            "3-oxo-cis-palm-9-eoyl-[acyl-carrier protein]"                        "(R)-3-hydroxy-cis-palm-9-eoyl-[acyl-carrier protein]"               
# [130] "3-Oxooctadecanoyl-[acyl-carrier protein]"                            "3-oxo-cis-vacc-11-enoyl-[acyl-carrier protein]"                      "(R)-3-hydroxy-cis-vacc-11-enoyl-[acyl-carrier protein]"             
# [133] "3-oxo-cis-D7-tetradecenoyl-ACPs"                                     "3-hydroxy-cis-D7-tetraecenoyl-ACPs"                                  "3-oxo-cis-D9-hexadecenoyl-ACPs"                                     
# [136] "3-hydroxy-cis-D9-hexaecenoyl-ACPs"                                   "3-Ketoglutaryl-ACP-methyl-ester"                                     "3-Hydroxyglutaryl-ACP-methyl-ester"                                 
# [139] "3-Ketopimeloyl-ACP-methyl-esters"                                    "3-hydroxypimeloyl-ACP-methyl-esters"                                 "3-oxo-dihomo gamma linolenoyl-[acp]"                                
# [142] "3-oxo-eicosatrienoyl [acp]"                                          "3-OXO-EICOSAPENTAENOYL-ACP"                                          "Beta-3-hydroxybutyryl-ACPs"                                         
# [145] "3-Hydroxy-octanoyl-ACPs"                                             "R-3-Hydroxypalmitoyl-ACPs"                                           "3-oxo-petroselinoyl-ACPs"                                           
# [148] "R-3-hydroxy-cis-vaccenoyl-ACPs"                                      "3-oxo-cis-vaccenoyl-ACPs"                                            "R-3-hydroxystearoyl-ACPs"                                           
# [151] "a 3-oxo-cis-delta5-dodecenoyl-[acp]"                                 "a 3-hydroxy cis delta5-dodecenoyl-[acp]"                             "a cis,cis-delta13,31-3-oxo-C50:2-[acp]"                             
# [154] "a cis,cis-delta13,31-3-hydroxyC50:2-[acp]"                           "a cis-delta19-3-oxo-C38:1-[acp]"                                     "a cis-delta19-3-hydroxyC38:1-[acp]"                                 
# [157] "a cis,cis-delta15,33-3-oxo-C52:2-[acp]"                              "a cis,cis-delta15,33-3-hydroxyC52:2-[acp]"                           "R-3-hydroxybehenoyl-ACPs"                                           
# [160] "a cis-delta11-3-oxo-C30:1-[acp]"                                     "a cis-delta11-3-hydroxyC30:1-[acp]"                                  "a cis,cis-delta5,23-3-oxo-C42:2-[acp]"                              
# [163] "a cis,cis-delta5,23-3-hydroxyC42:2-[acp]"                            "a cis,cis-delta17,35-3-oxo-C54:2-[acp]"                              "a cis,cis-delta17,35-3-hydroxyC54:2-[acp]"                          
# [166] "a cis,cis-delta13,25-3-oxo-C44:2-[acp]"                              "a cis,cis-delta13,25-3-hydroxyC44:2-[acp]"                           "a cis-delta9-3-oxo-C28:1-[acp]"                                     
# [169] "a cis-delta9-3-hydroxyC28:1-[acp]"                                   "a cis,cis-delta9,21-3-oxo-C40:2-[acp]"                               "a cis,cis-delta9,21-3-hydroxyC40:2-[acp]"                           
# [172] "a cis-delta13-3-oxo-C32:1-[acp]"                                     "a cis-delta13-3-hydroxyC32:1-[acp]"                                  "a cis,cis-delta7,19-3-oxo-C38:2-[acp]"                              
# [175] "a cis,cis-delta7,19-3-hydroxyC38:2-[acp]"                            "a cis-delta21-3-oxo-C40:1-[acp]"                                     "a cis-delta21-3-hydroxyC40:1-[acp]"                                 
# [178] "a cis-delta15-3-oxo-C34:1-[acp]"                                     "a cis-delta15-3-hydroxyC34:1-[acp]"                                  "a cis-delta7-3-oxo-C26:1-[acp]"                                     
# [181] "a cis-delta7-3-hydroxyC26:1-[acp]"                                   "a cis,cis-delta7,25-3-oxo-C44:2-[acp]"                               "a cis,cis-delta7,25-3-hydroxyC44:2-[acp]"                           
# [184] "a cis,cis-delta5,17-3-oxo-C36:2-[acp]"                               "a cis,cis-delta5,17-3-hydroxyC36:2[acp]"                             "a cis,cis-delta19,37-3-oxo-C56:2-[acp]"                             
# [187] "a cis,cis-delta19,37-3-hydroxy-C56:2-[acp]"                          "a cis,cis-delta15,27-3-oxo-C46:2-[acp]"                              "a cis,cis-delta15,27-3-hydroxyC46:2-[acp]"                          
# [190] "a cis,cis-delta9,27-3-oxo-C46:2-[acp]"                               "a cis,cis-delta9,27-3-hydroxyC46:2-[acp]"                            "a cis,cis-delta21,39-3-oxo-C58:2-[acp]"                             
# [193] "a cis,cis-delta21,39-3-hydroxyC58:2-[acp]"                           "a cis,cis-delta11,23-3-oxo-C42:2-[acp]"                              "a cis,cis-delta11,23-3-hydroxyC42:2-[acp]"                          
# [196] "a cis,cis-delta17,29-3-oxo-C48:2-[acp]"                              "a cis,cis-delta17,29-3-hydroxyC48:2-[acp]"                           "a cis-delta5-3-oxo-C24:1-[acp]"                                     
# [199] "a cis-delta5-3-hydroxyC24:1-[acp]"                                   "a cis,cis-delta11,29-3-oxo-C48:2-[acp]"                              "a cis,cis-delta11,29-3-hydroxyC48:2-[acp]"                          
# [202] "a cis,cis-delta19,31-3-oxo-C50:2-[acp]"                              "a cis,cis-delta19,31-3-hydroxyC50:2-[acp]"                           "Malonyl-acp-methyl-ester"                                           
# [205] "Pimelyl-[acyl-carrier protein] methyl ester"                         "Decanoyl-ACP"                                                        "Iso-C14:0 ACP"                                                      
# [208] "iso-C15:0 ACP"                                                       "Iso-C17:0 ACP"                                                       "pentadecenoyl-ACP (C15:1ACP)"                                       
# [211] "heptadecenoyl ACP (C17:1ACP)"                                        "pentadecanoyl-ACP (n-C15:0ACP)"                                      "heptadecanoyl-ACP (n-C17:0ACP)"                                     
# [214] "Iso-C13:0 ACP"                                                       "trans-3-cis-5-dodecenoyl-[acyl-carrier protein]"                     "cis-dodec-5-enoyl-[acyl-carrier protein]"                           
# [217] "trans-3-cis-7-myristoleoyl-[acyl-carrier protein]"                   "(2E)-Hexadecenoyl-[acp]"                                             "trans-3-cis-9-palmitoleoyl-[acyl-carrier protein]"                  
# [220] "trans-octadec-2-enoyl-[acyl-carrier protein]"                        "trans-3-cis-11-vacceoyl-[acyl-carrier protein]"                      "Crotonyl-ACPs"                                                      
# [223] "Butyryl-ACP"                                                         "Trans-D2-decenoyl-ACPs"                                              "(2E)-Tetradecenoyl-[acp]"                                           
# [226] "Cis-vaccenoyl-ACPs"                                                  "a cis-vaccen-2-enoyl-[acp]"                                          "Oleoyl-ACPs"                                                        
# [229] "2E-9Z-octadeca-2-9-dienoyl-ACPs"                                     "2-Hexadecenoyl-[acyl-carrier protein]"                               "cis-dec-3-enoyl-[acyl-carrier protein]"                             
# [232] "Iso-C7:0 ACP"                                                        "Iso-C6:0 ACP"                                                        "Fatty Acid (n-C5:0 ACP)"                                            
# [235] "Iso-C16:0 ACP"                                                       "Cis-Delta5-dodecenoyl-ACPs"                                          "Cis-Delta7-tetradecenoyl-ACPs"                                      
# [238] "Glutaryl-ACP-methyl-esters"                                          "Hexanoyl-ACP"                                                        "a cis-delta17-C36:1-[acp]"                                          
# [241] "Holo-ACP-Synthases"                                                  "a cis-delta19-C38:1-[acp]"                                           "a cis-delta9-C28:1-[acp]"                                           
# [244] "a cis-delta7-C26:1-[acp]"                                            "a cis-delta11-C30:1-[acp]"                                           "a cis-docos-3-enoyl-[acp]"                                          
# [247] "a cis-delta5-C24:1-[acp]"                                            "a cis-delta13-C32:1-[acp]"                                           "a cis-delta15-C34:1-[acp]"                                          
# [250] "a cis-delta17-3-oxo-C36:1-[acp]"                                     "Acetyl-ACP"                                                          "All-ACPs"                                                           
# [253] "(2E)-Octenoyl-[acp]"                                                 "Trans-D3-cis-D7-tetradecenoyl-ACPs"                                  "Trans-D3-cis-D9-hexadecenoyl-ACPs"                                  
# [256] "Enoylglutaryl-ACP-methyl-esters"                                     "Enoylpimeloyl-ACP-methyl-esters"                                     "Trans-D3-cis-D5-dodecenoyl-ACPs"                                    
# [259] "cis-3-Decenoyl-[acyl-carrier protein]"                               "Cis-delta-3-decenoyl-ACPs"                                           "a cis,cis-delta11,29-C48:2-[acp]"                                   
# [262] "a cis,cis-delta13,31-C50:2-[acp]"                                    "a cis,cis-delta7,19-C38:2-[acp]"                                     "a cis,cis-delta3,15-C34:2-[acp]"                                    
# [265] "a cis,cis-delta11,23-C42:2-[acp]"                                    "a cis,cis-delta5,23-C42:2-[acp]"                                     "a cis,cis-delta17,35-C54:2-[acp]"                                   
# [268] "a cis,cis-delta13,25-C44:2-[acp]"                                    "a cis,cis-delta7,25-C44:2-[acp]"                                     "a cis,cis-delta9,21-C40:2-[acp]"                                    
# [271] "a cis,cis-delta15,27-C46:2-[acp]"                                    "a cis,cis-delta3,21-C40:2-[acp]"                                     "a cis,cis-delta9,27-C46:2-[acp]"                                    
# [274] "a cis,cis-delta17,29-C48:2-[acp]"                                    "a cis,cis-delta5,17-C36:2-[acp]"                                     "a cis,cis-delta15,33-C52:2-[acp]"                                   
# [277] "(Z)-hexadec-11-enoyl-ACP"                                            "(Z)-3-oxooctadec-13-enoyl-ACP"                                       "hexadecanoyl-acp"                                                   
# [280] "12-methyl-tetra-decanoyl-ACP"                                        "14-methyl-hexa-decanoyl-ACP"                                         "13-methyl-tetra-decanoyl-ACP"                                       
# [283] "15-methyl-hexa-decanoyl-ACP"                                         "12-methyl-tridecanoyl-ACP"                                           "14-methyl-pentadecanoyl-ACP"                                        
# [286] "Octodecanoyl-ACP"                                                    "2-methylbutyryl-ACP"                                                 "isovaleryl-ACP"                                                     
# [289] "isobutyryl-ACP"                                                      "Octadecynoyl-ACP"                                                    "a trans-methoxy-C60-meroacyl-[acp]"                                 
# [292] "Delta6-hexadecenoyl-ACPs"                                            "Delta4-hexadecenoyl-ACPs" 


unique(df.out$f__in[sel2[-subsel]]) # 146 enzyme / sub-functions


length(other_ACP_names) # 293

# other methyl ACPs?
sel.other_mm_acps <- grep(pattern = "methyl", x = other_ACP_names)
other_ACP_names[sel.other_mm_acps]
# [1] "a cis-delta19-37-hydroxy-38-methyl-C57:1-[acp]"   "a trans-delta18-37-hydroxy-38-methyl-C58:1-[acp]" "a cis-delta21-39-hydroxy-40-methyl-C59:1-[acp]"   "a trans-delta20-39-hydroxy-40-methyl-C60:1-[acp]"
# [5] "4-methyl-3-oxo-hexanoyl-ACP"                      "4-methyl-3-hydroxy-hexanoyl-ACP"                  "6-methyl-3-oxo-octanoyl-ACP"                      "6-methyl-3-hydroxy-octanoyl-ACP"                 
# [9] "8-methyl-3-oxo-decanoyl-ACP"                      "8-methyl-3-hydroxy-decanoyl-ACP"                  "10-methyl-3-oxo-dodecanoyl-ACP"                   "10-methyl-3-hydroxy-dodecanoyl-ACP"              
# [13] "12-methyl-3-oxo-tetra-decanoyl-ACP"               "12-methyl-3-hydroxy-tetra-decanoyl-ACP"           "14-methyl-3-oxo-hexa-decanoyl-ACP"                "14-methyl-3-hydroxy-hexa-decanoyl-ACP"           
# [17] "5-methyl-3-oxo-hexanoyl-ACP"                      "5-methyl-3-hydroxy-hexanoyl-ACP"                  "7-methyl-3-oxo-octanoyl-ACP"                      "7-methyl-3-hydroxy-octanoyl-ACP"                 
# [21] "9-methyl-3-oxo-decanoyl-ACP"                      "9-methyl-3-hydroxy-decanoyl-ACP"                  "11-methyl-3-oxo-dodecanoyl-ACP"                   "11-methyl-3-hydroxy-dodecanoyl-ACP"              
# [25] "13-methyl-3-oxo-tetra-decanoyl-ACP"               "13-methyl-3-hydroxy-tetra-decanoyl-ACP"           "15-methyl-3-oxo-hexa-decanoyl-ACP"                "15-methyl-3-hydroxy-hexa-decanoyl-ACP"           
# [29] "4-methyl-3-oxo-pentanoyl-ACP"                     "4-methyl-3-hydroxy-pentanoyl-ACP"                 "6-methyl-3-oxo-heptanoyl-ACP"                     "6-methyl-3-hydroxy-heptanoyl-ACP"                
# [33] "8-methyl-3-oxo-nonanoyl-ACP"                      "8-methyl-3-hydroxy-nonanoyl-ACP"                  "10-methyl-3-oxo-undecanoyl-ACP"                   "10-methyl-3-hydroxy-undecanoyl-ACP"              
# [37] "12-methyl-3-oxo-tridecanoyl-ACP"                  "12-methyl-3-hydroxy-tridecanoyl-ACP"              "14-methyl-3-oxo-pentadecanoyl-ACP"                "14-methyl-3-hydroxy-pentadecanoyl-ACP"           
# [41] "3-Ketoglutaryl-ACP-methyl-ester"                  "3-Hydroxyglutaryl-ACP-methyl-ester"               "3-Ketopimeloyl-ACP-methyl-esters"                 "3-hydroxypimeloyl-ACP-methyl-esters"             
# [45] "Malonyl-acp-methyl-ester"                         "Pimelyl-[acyl-carrier protein] methyl ester"      "Glutaryl-ACP-methyl-esters"                       "Enoylglutaryl-ACP-methyl-esters"                 
# [49] "Enoylpimeloyl-ACP-methyl-esters"                  "12-methyl-tetra-decanoyl-ACP"                     "14-methyl-hexa-decanoyl-ACP"                      "13-methyl-tetra-decanoyl-ACP"                    
# [53] "15-methyl-hexa-decanoyl-ACP"                      "12-methyl-tridecanoyl-ACP"                        "14-methyl-pentadecanoyl-ACP"                      "2-methylbutyryl-ACP" 

# contain additional hydroxyl groups, carbonyl (oxo), 

sel.methyl <- grep(pattern = "methyl", x = other_ACP_names) # 56
sel.enoyl <- grep(pattern = "enoyl", x = other_ACP_names) # 57
sel.match <- which(sel.methyl %in% sel.enoyl) # Empty <<< i.e., no other methyl unsaturated bcfa !

sel.oxo <- grep(pattern = "-oxo-", x = other_ACP_names) # 60
length(which(sel.methyl %in% sel.oxo)) # 18 <<<
sel.hydroxy <- grep(pattern = "-hydroxy-", x = other_ACP_names) # 33
length(which(sel.methyl %in% sel.hydroxy)) # 22 <<<

sel.esters <- grep(pattern = "ester", x = other_ACP_names) # 9
length(which(sel.methyl %in% sel.esters)) # 9 <<<

sel.anoyl <- grep(pattern = "anoyl", x = other_ACP_names) # 75
length(which(sel.methyl %in% sel.anoyl & !sel.methyl %in% sel.oxo & !sel.methyl %in% sel.hydroxy)) # 6 <<<
subsel <- which(sel.methyl %in% sel.anoyl & !sel.methyl %in% sel.oxo & !sel.methyl %in% sel.hydroxy) # 6
other_ACP_names[sel.methyl[subsel]]
# "12-methyl-tetra-decanoyl-ACP" "14-methyl-hexa-decanoyl-ACP"  "13-methyl-tetra-decanoyl-ACP" "15-methyl-hexa-decanoyl-ACP"  "12-methyl-tridecanoyl-ACP"    "14-methyl-pentadecanoyl-ACP" 
# 14-chain + methyl(12).         16-chain + methyl(14).         14-chain + methyl(13).          16-chain + methyl(15).         13-chain + methyl(12).        15-chain + methyl(14)

0 + 18 + 22 + 9 + 6 + 1 # 56

# ZERO monomethyl unsaturated (enoyl) bcfa + 18 oxo (carbonyl-containing) + 22 hydroxy (hydroxyl-containing) + 9 esters + 6 saturated (long, 13 to 16 chain length) + 1 x 2-methylbutyryl-ACP (saturated)


#-------------------------


##########################
##########################
##########################
##########################

#### Australian Microbiome Initiative (AMI) disturbed vs natural soils
#    additional validation dataset
#    using same phyloseq data object as per previous CPP method in Liddicoat et al 2024

#### AMI - disturbed vs natural - build reaction search in parallel - get_reactions & compounds
#-------------------------

# use same phyloseq data object for SUPER_FOCUS functions as used in Liddicoat et al 2024
phy <- readRDS("phy-phyloseq-object-AMI.RDS")

phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 37335 taxa and 84 samples ]
# sample_data() Sample Data:       [ 84 samples by 92 sample variables ]
# tax_table()   Taxonomy Table:    [ 37335 taxa by 4 taxonomic ranks ]

table(phy@sam_data$nat_vs_dist)
# Disturbed   Natural 
# 29        55 

## convert each row in functional tax_table to "mean van Krevelen distance to health-associated compounds"

df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 37335     4


get_rxns_and_compounds_indiv <- function( df.tax, subsys.lut, rxns.lut, rxn_pathways.lut ) {
  
  rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
  rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
  
  ###testing
  # #for (i in 1:dim(df.tax)[1]) {
  # for (i in 1:1000) {
  #i<-1
  #i<-22889
  #i<-22894
  #i<-31768
  #i<-9406
  #i<-9422
  #i<-8
  #1<-13
  #i<-7
  
  sub1 <- df.tax$subsys_L1[i]
  sub2 <- df.tax$subsys_L2[i]
  sub3 <- df.tax$subsys_L3[i]
  
  fxn.temp <- df.tax$fxn[i]
  fxn.superfocus.rowlabel <- row.names(df.tax)[i]
  
  # store results corresponding to each Superfocus row
  fxn.list <- list()
  fxn.list[[ fxn.superfocus.rowlabel  ]] <- list()
  
  # check for multiple functions/reactions?
  flag1 <- grepl(pattern = "_/_|/", x = fxn.temp)
  flag2 <- grepl(pattern = "_@_", x = fxn.temp)
  if (!any(flag1,flag2)==TRUE) {
    # no multiples
    fxns <- fxn.temp
  } else if (flag1==TRUE) {
    fxns <- unlist( strsplit(fxn.temp, split = "_/_") )  ###### WHAT ABOUT SPLIT FOR "/" WITHOUT UNDERSCORES ??
  } else {
    fxns <- unlist( strsplit(fxn.temp, split = "_@_") )
  }
  # remove underscores
  ( fxns <- gsub(pattern = "_", replacement = " ", x = fxns) )
  
  # process each fxn & store attributes
  #df.fxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=1:length(fxns),`f__in`=fxns, matching_method=NA, min_adist_modelSEED=NA, min_amatch_modelSEED=NA, rxns=NA, tot_mean_OC_x=NA, tot_mean_HC_y=NA , tot_mean_NC_z=NA )
  
  df.fxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=1:length(fxns),`f__in`=fxns, matching_method=NA, rxns=NA) #, compound_name=NA, compound_form=NA, compound_molar_prop=NA, OC_x=NA, HC_y=NA , NC_z=NA )
  
  # # do round brackets interfere with search? - YES
  # lookfor <- "option 4 (this one)"
  # lookuplist <- c("option 1", "option 2", "option 3 (this one)", "option 4 (this one)")
  # grep(pattern = lookfor, x = lookuplist)
  
  # Identify '/' separators with no '_'  ??
  
  for (f in 1:length(fxns)) {  # this accounts for multiple functions/reactions reported in Superfocus outputs
    #f<-1
    #f<-2
    f.in <- fxns[f]
    
    # these concatenated expressions will be used to look for exact match using hierarchy in ModelSEED Subsystem table
    full_hier_target <- paste0(sub1,"__",sub2,"__",sub3,"__",f.in)
    full_hier_list <- paste0(subsys.lut$Class,"__",subsys.lut$Subclass,"__",gsub("_"," ",subsys.lut$Name),"__",subsys.lut$Role)
    
    ## data cleaning
    
    # trim off '_#' and '_##' tags
    trim_nchar <- str_locate(string = f.in, pattern = " # | ## ")[1]
    if (!is.na(trim_nchar) & length(trim_nchar)==1) {
      f.in <- substring(text = f.in , first = 1, last = trim_nchar-1)
    }
    
    # Eliminate unwanted parsing of regular expressions: '[', ']','***', '(', ')'
    f.in <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\} ", replacement ="." , x = f.in) # used later
    
    #rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
    #rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
    
    full_hier_target <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_target)
    full_hier_list <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_list)
    
    sel.rx <- grep(pattern = full_hier_target, x = full_hier_list)
    
    ## ALTERNATIVE #1 == FULL HIERACHICAL MATCH
    if (length(sel.rx)>=1) {
      df.fxns$matching_method[f] <- "Exact hierachy match"
      df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
      
    } else if (str_detect(string = fxns[f], pattern = " \\(EC ")) {  ## ALTERNATIVE #2 == MATCHING ECs
      # search by EC id if present
      
      f.in <- fxns[f] # this goes back to string with brackets for EC
      ## LOOK FOR MULTIPLE ECs ??????????
      # 22889
      # 22894
      # 31768
      
      how_many_ECs <- str_count(string = f.in, pattern = "\\(EC.*?\\)")
      
      ECs <- as.character( str_extract_all(string = f.in, pattern = "\\(EC.*?\\)", simplify = TRUE) )
      #class(ECs)
      ECs <- gsub(pattern = "\\(EC |\\)", replacement = "", x = ECs)
      ECs.collapse <- paste0(ECs, collapse = "|")
      
      sel.rx <- which(rxns.lut$ec_numbers == ECs.collapse)
      
      if (length(how_many_ECs)==0 | length(ECs)==0) {
        # there was a glitch, database typo, or some error in identifying the EC number
        df.fxns$matching_method[f] <- "No match found"
        df.fxns$rxns[f] <- NA
        
      } else if (length(sel.rx)>=1) {
        # combined EC hits identified
        df.fxns$matching_method[f] <- "EC number"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if (length(which(rxns.lut$ec_numbers %in% ECs)) >=1) {
        # treat EC hits individually
        sel.rx <- which(rxns.lut$ec_numbers %in% ECs) # look 1st where ECs are exact matches for EC numbers in Reactions lookup table
        
        df.fxns$matching_method[f] <- "EC number"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if (length(grep(pattern = ECs, x = rxns.lut$ec_numbers)) >=1) {
        # this allows EC to be part of a combination of EC numbers that are listed in Reactions lookup table
        sel.rx <- grep(pattern = ECs, x = rxns.lut$ec_numbers)
        
        df.fxns$matching_method[f] <- "EC number"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else {
        # it had an EC number but couldn't find a match in the EC numbers listed in Reaction lookup table
        df.fxns$matching_method[f] <- "No match found"
        df.fxns$rxns[f] <- NA
        
      }
      # END EC matching
      
      
    } else {  ## ALTERNATIVE 3 == FXN NAME MATCHING
      ## otherwise attempt to match function name - a) first look for exact matches   ########## then b) closest match above a threshold
      # 1. 'reactions' table by name: rxns.lut$name
      # 2. 'reactions' table by aliases: rxns.lut$aliases
      # 3. 'Model SEED Subsystems' table by Role: subsys.lut$Role
      # 4. 'Unique_ModelSEED_Reaction_Pathways' table by External ID: rxn_pathways.lut$External_rxn_name
      
      if ( length( grep(pattern = f.in, x = rxns.lut$name) )>=1 ) {
        # 1a - exact match - rxns.lut$name
        sel.rx <- grep(pattern = f.in, x = rxns.lut$name)
        #rxns.lut$name[sel.rx]
        df.fxns$matching_method[f] <- "Matched Reactions name"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if ( length( grep(pattern = f.in, x = rxns.lut$aliases) )>=1 ) {
        # 2a - exact match - rxns.lut$aliases
        sel.rx <- grep(pattern = f.in, x = rxns.lut$aliases)
        #rxns.lut$aliases[sel.rx]
        #rxns.lut$name[sel.rx]
        
        df.fxns$matching_method[f] <- "Matched Reactions aliases"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if ( length( grep(pattern = f.in, x = subsys.lut$Role) )>=1 ) {
        # 3a - exact match - subsys.lut$Role
        sel.rx <- grep(pattern = f.in, x = subsys.lut$Role)
        #subsys.lut$Role[sel.rx]
        #subsys.lut$Reaction[sel.rx]
        
        df.fxns$matching_method[f] <- "Matched Subsytem role"
        df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
        
      } else if ( length( grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name) )>=1 ) {
        # 4a - exact match - rxn_pathways.lut$External_rxn_name
        sel.rx <- grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name)
        
        df.fxns$matching_method[f] <- "Matched ModelSEED Reaction pathways"
        df.fxns$rxns[f] <- paste0( unique(rxn_pathways.lut$rxn_id[sel.rx]), collapse = ";")
        
        
      } else {
        df.fxns$matching_method[f] <- "No match found"
        df.fxns$rxns[f] <- NA
        
      }
      
      ## DON'T RUN PARTIAL MATCHING AT THIS STAGE
      
      
    } # END function - reaction search
    
    #fxn.list[[ fxn.superfocus.rowlabel  ]][[ f ]][[ "fxns" ]] <- df.fxns
    
    print(paste0("completed fxn ", f))
    
    
    ## now investigate these reactions ...
    # Reactions lookup table: 
    # - "equation": Definition of reaction expressed using compound IDs and after protonation
    # Compounds lookup table:
    # - "formula": Standard chemical format (using Hill system) in protonated form to match reported charge
    #df.fxns
    
    
    #if (df.fxns$matching_method == "No match found") {
    if (df.fxns$rxns[f] == "" | is.na(df.fxns$rxns[f])) {
      
      df.Rxns <- NA
      df.Compounds <- NA
      
    } else { # reaction(s) were identified
      
      # consider reactions for this f.in only (possibly > 1 f.in per Superfocus row)
      f.in.rxns <- unique(unlist(str_split(string = df.fxns$rxns[f], pattern = ";")))
      
      df.Rxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel, f=f, f__in=f.in,rxn_id= f.in.rxns,
                            rxn_name=NA, rxn_eqn=NA, rxn_defn=NA,compds=NA,compd_coef=NA, chem_formx=NA ) #, OC_ratios=NA, HC_ratios=NA, NC_ratios=NA, coefwtmean_OC_x=NA, coefwtmean_HC_y=NA, coefwtmean_NC_z=NA)
      
      #df.Compounds <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=f, f__in=f.in,rxn_id= f.in.rxns, compound_name=NA, compound_form=NA, compound_molar_prop=NA, OC_x=NA, HC_y=NA , NC_z=NA )
      
      for (r in 1:dim(df.Rxns)[1]) {
        #r<-1
        #this_rxn <- "rxn00004"
        this_rxn <- df.Rxns$rxn_id[r]
        sel <- which(rxns.lut$id == this_rxn)
        ( df.Rxns$rxn_name[r] <- rxns.lut$name[sel] )
        ( df.Rxns$rxn_eqn[r] <- rxns.lut$equation[sel] )
        ( df.Rxns$rxn_defn[r] <- rxns.lut$definition[sel] )
        
        # extract compound info
        
        #df.Rxns$rxn_eqn[r]
        #[1] "(1) cpd00010[0] + (1) cpd29672[0] <=> (1) cpd00045[0] + (1) cpd11493[0]"
        #[1] "(45) cpd00144[0] + (45) cpd00175[0] <=> (45) cpd00014[0] + (45) cpd00091[0] + (1) cpd15634[0]"
        
        ( compds.idx <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd")[[1]][,"start"] )
        # 5 23 43 61
        # 6 25 46 65 83
        
        ( compds <- as.character( str_extract_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd.....", simplify = TRUE) ) )
        # "cpd00010" "cpd29672" "cpd00045" "cpd11493"
        
        if (length(compds)>=1) {
          
          df.Rxns$compds[r] <- paste0(compds, collapse = ";")
          
          ## get compound coefficients?
          start_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\(")[[1]][,"start"]
          end_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\)")[[1]][,"start"]
          ( compd.coeff <- as.numeric( substring(text = df.Rxns$rxn_eqn[r], first = start_brackets+1, last = end_brackets-1)) )
          
          df.Rxns$compd_coef[r] <- paste0(compd.coeff, collapse = ";")
          
          # get formulas of compounds
          
          formx <-filter(compounds.lut, id %in% compds )
          row.names(formx) <- formx$id
          ( formx.char <- formx[compds, ]$formula )
          # "C21H32N7O16P3S" "HOR"            "C10H11N5O10P2"  "C11H22N2O7PRS" 
          # "C15H19N2O18P2"      "C17H25N3O17P2"      "C9H12N2O12P2"       "C9H11N2O9P"         "C630H945N45O630P45"
          # "C7H7O7" "H2O"    "C7H5O6"
          df.Rxns$chem_formx[r] <- paste0(formx.char, collapse = ";")
          
          ( compd.names <- formx[compds, ]$name )
          # "2-methyl-trans-aconitate" "cis-2-Methylaconitate"
          
          
          # # O:C ratio (replace Inf with NA when C not present)
          # OC_ratio <- atomic_no(form = formx.char, atom = "O")/atomic_no(form = formx.char, atom = "C")
          # OC_ratio[is.infinite(OC_ratio)] <- NA
          # 
          # # H:C ratio (replace Inf with NA when C not present)
          # HC_ratio <- atomic_no(form = formx.char, atom = "H")/atomic_no(form = formx.char, atom = "C")
          # HC_ratio[is.infinite(HC_ratio)] <- NA
          # 
          # # N:C ratio (replace Inf with NA when C not present)
          # NC_ratio <- atomic_no(form = formx.char, atom = "N")/atomic_no(form = formx.char, atom = "C")
          # NC_ratio[is.infinite(NC_ratio)] <- NA
          
          temp.df.Compounds <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=f, f__in=f.in,rxn_id= f.in.rxns[r], 
                                          cpd_id=compds, cpd_name=compd.names, cpd_form=formx.char, cpd_molar_prop=compd.coeff #, 
                                          #OC_x=OC_ratio, HC_y=HC_ratio , NC_z=NC_ratio 
          )
          
        } else {
          # No specified reaction equation or chemical formula info
          df.Rxns$compds[r] <- NA
          df.Rxns$compd_coef[r] <- NA
          df.Rxns$chem_formx[r] <- NA
          
          temp.df.Compounds <- NA
          
        }
        
        if (r==1) { df.Compounds <- temp.df.Compounds }
        
        if (r>1 & is.data.frame(df.Compounds) & is.data.frame(temp.df.Compounds)) { df.Compounds <- rbind(df.Compounds, temp.df.Compounds) }
        
        # clean up - if there are additional reactions?
        temp.df.Compounds <- NA
        
      } # END loop for r - rxn_id's per f/f.in
      
    } # END else loop when reactions identified
    
    # store results corresponding to each sub-reaction of each Superfocus row
    fxn.list[[ fxn.superfocus.rowlabel  ]][[ "fxns" ]] <- df.fxns
    
    if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]] <- list() } # set this only once
    fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]][[ f ]] <- df.Rxns
    
    if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]] <- list() } # set this only once
    fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]][[ f ]] <- df.Compounds
    
    
  } # END loop - f in 1:length(fxns)) - to account for multiple functions/reactions reported in each row of Superfocus outputs
  
  
  #return(fxn.list)
  
  saveRDS(object = fxn.list, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/AMI-R-working-files-indiv/fxn-list-",fxn.superfocus.rowlabel,".rds") ) # use readRDS()
  
  #print(paste0("COMPLETED ROW ",i," OF SUPERFOCUS FUNCTIONAL TAXA  # # # # # # # # # # # # # # # # # # # # #"))
  
} # END function to be run in parallel for each superfocus row


# # # # # # # # # # # # # # # # # #


time.start <- Sys.time()
cl<-makeCluster( detectCores()-1 )
registerDoParallel(cl)

#foreach(i=1:100 , .packages=c('stringr', 'dplyr')) %dopar%
foreach(i=1:dim(df.tax)[1] , .packages=c('stringr', 'dplyr')) %dopar%  #
  get_rxns_and_compounds_indiv( df.tax=df.tax, subsys.lut=subsys.lut, rxns.lut=rxns.lut, rxn_pathways.lut=rxn_pathways.lut )

stopCluster(cl)
time.finish <- Sys.time()



time.start #  
#[1] "2025-10-11 14:20:48 ACDT"
time.finish # 
#[1] "2025-10-11 17:43:17 ACDT"




## assemble results

# results saved here
#saveRDS(object = fxn.list, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/AMI-R-working-files-indiv/fxn-list-",fxn.superfocus.rowlabel,".rds") ) # use readRDS()


modelSEED_rxn_result_dir <- "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/AMI-R-working-files-indiv"


dim(df.tax)
# 37335     4





# Check for secondary fxns ?? in list objects ??????????


# read first output
i<-1
#temp <- readRDS(paste0(modelSEED_rxn_result_dir,"/fxn-list-fxn_",i,".rds"))
fxn.superfocus.rowlabel <- row.names(df.tax)[i]
temp <- readRDS(paste0(modelSEED_rxn_result_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))

print( length(temp) )
print( names(temp) )
# "fxn_1"
print( paste0( class(temp[[1]])," ",names(temp)," of length ", length(temp[[1]])," named ",paste0(names(temp[[1]]), collapse = " & ") ))
# "list fxn_1 of length 3 named fxns & rxns & compounds"
class ( temp[[1]][["fxns"]] ) # data.frame
class ( temp[[1]][["rxns"]] ) # list
class ( temp[[1]][["compounds"]] ) # list

length( temp[[1]][["compounds"]] ) # 1
names( temp[[1]][["compounds"]] ) # NULL
class( temp[[1]][["compounds"]][[1]] ) # "logical"
is.na( temp[[1]][["compounds"]][[1]] )

i<-2
fxn.superfocus.rowlabel <- row.names(df.tax)[i]
temp <- readRDS(paste0(modelSEED_rxn_result_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))

length(temp) # 1
names(temp) # "fxn_2"
print( paste0( class(temp[[1]])," ",names(temp)," of length ", length(temp[[1]])," named ",paste0(names(temp[[1]]), collapse = " & ") ))
# "list fxn_2 of length 3 named fxns & rxns & compounds"
class ( temp[[1]][["fxns"]] ) # data.frame
class ( temp[[1]][["rxns"]] ) # list
class ( temp[[1]][["compounds"]] ) # list

length( temp[[1]][["compounds"]] ) # 1
names( temp[[1]][["compounds"]] ) # NULL
class( temp[[1]][["compounds"]][[1]] ) # "data.frame"

df.out <- temp[[1]][["compounds"]][[1]]

names(df.out) #
# [1] "superfocus_fxn" "f"              "f__in"          "rxn_id"         "cpd_id"         "cpd_name"       "cpd_form"       "cpd_molar_prop" 

# i<-13
# fxn.superfocus.rowlabel <- row.names(df.tax)[i]
# temp <- readRDS(paste0(modelSEED_rxn_result_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))
# 
# length(temp) # 1
# names(temp) # "fxn_13"
# print( paste0( class(temp[[1]])," ",names(temp)," of length ", length(temp[[1]])," named ",paste0(names(temp[[1]]), collapse = " & ") ))
# # "list fxn_13 of length 3 named fxns & rxns & compounds"
# class ( temp[[1]][["fxns"]] ) # data.frame
# class ( temp[[1]][["rxns"]] ) # list
# class ( temp[[1]][["compounds"]] ) # list
# 
# length( temp[[1]][["compounds"]] ) # 3
# names( temp[[1]][["compounds"]] ) # NULL
# class( temp[[1]][["compounds"]][[1]] ) # "data.frame"
# class( temp[[1]][["compounds"]][[2]] ) # "data.frame"
# class( temp[[1]][["compounds"]][[3]] ) # "data.frame"



# total results written to disk - 30125
dim(df.tax) # 37335     4
num_results_files <- dim(df.tax)[1]
#num_results_files <- 100


# assemble all compound data outputs
# start with blank row

df.out <- data.frame(superfocus_fxn=NA, f=NA, f__in=NA, rxn_id=NA, cpd_id=NA, cpd_name=NA, cpd_form=NA, cpd_molar_prop=NA #, 
                     #OC_x=NA, HC_y=NA, NC_z=NA
)

for (i in 1:num_results_files) {
  #for (i in 25760:num_results_files) {
  #i<-1
  #i<-2
  #i<-13
  fxn.superfocus.rowlabel <- row.names(df.tax)[i]
  temp <- readRDS(paste0(modelSEED_rxn_result_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))
  
  f_no <- length( temp[[1]][["compounds"]] )
  
  for (f in 1:f_no) {
    #f<-2
    # only add non-NA results
    if (is.data.frame( temp[[1]][["compounds"]][[f]] )) {
      
      df.temp <- temp[[1]][["compounds"]][[f]]
      ok <- complete.cases(df.temp)
      df.temp <- df.temp[ which(ok==TRUE), ] # updated version will include some compounds with vK coordinates that are NA. vK coordinates are considered later
      df.out <- rbind(df.out,df.temp)
    }
  }
  
  print(paste0("added df ",i," of ",num_results_files ))
  
}


str(df.out)
# 'data.frame':	1340067 obs. of  8 variables:

saveRDS(object = df.out, file = "df.out--get_rxns_and_compounds_indiv--ami.RDS")
df.out <- readRDS(file = "df.out--get_rxns_and_compounds_indiv--ami.RDS")

# remove NA first row
head(df.out)
#   superfocus_fxn  f                                                                         f__in   rxn_id   cpd_id
# 1            <NA> NA                                                                          <NA>     <NA>     <NA>
# 2           fxn_2  1                                                   2-methylaconitate isomerase rxn25278 cpd25681
# 3           fxn_2  1                                                   2-methylaconitate isomerase rxn25278 cpd02597
# 11          fxn_3  1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd24620
# 31          fxn_3  1 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117) rxn25279 cpd25681
# 14          fxn_4  1                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79) rxn03060 cpd01501
#                                        cpd_name cpd_form cpd_molar_prop      OC_x      HC_y NC_z
# 1                                          <NA>     <NA>             NA        NA        NA   NA
# 2                      2-methyl-trans-aconitate   C7H5O6              1 0.8571429 0.7142857    0
# 3                         cis-2-Methylaconitate   C7H5O6              1 0.8571429 0.7142857    0
# 11 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7              1 1.0000000 1.0000000    0
# 31                     2-methyl-trans-aconitate   C7H5O6              1 0.8571429 0.7142857    0
# 14                              2-Methylcitrate   C7H7O7              1 1.0000000 1.0000000    0

df.out <- df.out[-1, ]


# check for different cpd_molar_prop ??
hist(df.out$cpd_molar_prop)

dim(df.out) # 1340066       8


# normalise molar_prop to cpd_relabun so total of 1 per superfocus function !!

df.out$cpd_molar_prop_norm <- NA

length(unique(df.out$superfocus_fxn)) # 19717

phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 37335 taxa and 84 samples ]
# sample_data() Sample Data:       [ 84 samples by 92 sample variables ]
# tax_table()   Taxonomy Table:    [ 37335 taxa by 4 taxonomic ranks ]

100*(length(unique(df.out$superfocus_fxn)) / ntaxa(phy)) # 52.81104 % of functions represented
100*(19717/37335) # 52.81104

fxns_found <- unique(df.out$superfocus_fxn)

for (k in 1:length(fxns_found)) {
  #k<-1
  this_fxn <- fxns_found[k]
  sel <- which(df.out$superfocus_fxn == this_fxn)
  
  sum_molar_prop <- sum( df.out$cpd_molar_prop[sel], na.rm = TRUE)
  # calculate 
  
  df.out$cpd_molar_prop_norm[sel] <- df.out$cpd_molar_prop[sel]/sum_molar_prop
  
  print(paste0("completed ",k))
  
}

sum(df.out$cpd_molar_prop_norm) # 19717


sample_sums(phy)
# 12465 12469 12471 12473 12475 12477 12479 12491 12493 12495 12497 12501 12503 12505 12507 12509 12511 12513 12515 12517 12519 12521 12523 12525 12527 
# 100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100 
# 12560 12562 12564 12566 12568 12570 12572 12574 12576 12578 12580 12582 12590 12614 12616 12618 12620 12622 12624 13260 13274 13278 13729 13731 13733 
# 100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100 
# 13735 13906 19471  7845  7861  7863  7888  7890  7900  7902  8076  8078  8084  8088  8126  8138  8140  8142  8144  8146  8150  8152  9430  9434  9436 
# 100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100   100 
# 9438  9440  9442  9446  9458  9462  9464  9466  9468 
# 100   100   100   100   100   100   100   100   100 

dim(df.out) # 1340066       9




getwd() # "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R"

saveRDS(object = df.out, file = "df.out--tidy-compounds_indiv--cpp3d-ami.RDS")

df.out <-readRDS( "df.out--tidy-compounds_indiv--cpp3d-ami.RDS" )


#-------------------------

#### AMI
#    Collate CPP (compound rel abun %)
#    Test disturbed vs natural for compounds with trending metabolism: sugars, lignin, BCFA-ACPs
#-------------------------

this_study <- "-ami-"
header <- "cpp3d-cpdtrend"
phy <- readRDS("phy-phyloseq-object-AMI.RDS")
df.out <-readRDS( "df.out--tidy-compounds_indiv--cpp3d-ami.RDS" )
dim(df.out) # 1340066       9

df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 37335     4


phy@sam_data


sample_names(phy)
identical( sample_names(phy), colnames( as.matrix( phy@otu_table)) ) # TRUE

df.OTU <- as.data.frame( phy@otu_table ) # this is Superfocus functional relative abundance data represented in phyloseq OTU abundance table
dim(df.OTU) # 37335    84
df.OTU[1:5, 1:8]
#              12465        12469        12471        12473        12475        12477        12479        12491
# fxn_1 4.085207e-04 4.205200e-04 2.964560e-04 3.198079e-04 3.852306e-04 3.076140e-04 4.778209e-04 2.779646e-04
# fxn_2 2.288676e-05 2.821954e-05 0.000000e+00 1.035635e-05 4.866925e-05 3.717924e-05 2.275338e-05 1.166878e-05
# fxn_3 5.276470e-05 5.603135e-05 1.824345e-05 2.051317e-05 1.165885e-05 2.124528e-05 6.929437e-05 0.000000e+00
# fxn_4 3.152034e-05 2.634775e-05 3.420646e-05 2.616284e-05 1.391540e-04 7.170283e-05 1.965064e-04 0.000000e+00
# fxn_5 1.636470e-03 1.417673e-03 1.334052e-03 1.783856e-03 2.032855e-03 1.633231e-03 1.566880e-03 1.736613e-03

mean( df.OTU[ , "12465"] ) # 0.002678452
sum( df.OTU[ , "12465"] ) # 100

sample_sums(phy) # all values of 100


# loop through each sample to get potential metabolism of compounds


get_cpd_relabun_per_sample <- function(phy_in, dat.cpd) {
  #i<-1
  #phy_in = phy
  #dat.cpd = df.out
  
  this_samp <- sample_names(phy_in)[i]
  df.OTU <- as.data.frame( phy_in@otu_table[ ,this_samp] )
  
  dat.cpd$sample <- this_samp
  
  dat.cpd$cpd_rel_abun_norm <- NA
  
  fxns_all <- row.names(df.OTU)
  
  for (k in 1:length(fxns_all)) {
    #k<-1
    this_fxn <- fxns_all[k]
    sel <- which(dat.cpd$superfocus_fxn == this_fxn)
    
    if (length(sel)>=1) {
      dat.cpd$cpd_rel_abun_norm[sel] <- df.OTU[this_fxn, ]*dat.cpd$cpd_molar_prop_norm[sel]
      
    }
  } # END rel abun values for all relevant functions added
  
  saveRDS(object = dat.cpd, file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/AMI-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
  
} # END


time.start <- Sys.time()
cl<-makeCluster( detectCores()-1 )
registerDoParallel(cl)

foreach(i=1: length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
  get_cpd_relabun_per_sample( phy_in = phy, dat.cpd = df.out)

stopCluster(cl)
time.finish <- Sys.time()

# output 1
i<-1
this_samp <- sample_names(phy)[i]
dat <- readRDS ( file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/AMI-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
head(dat)

for ( i in 2:length(sample_names(phy)) ) {
  #i<-1
  this_samp <- sample_names(phy)[i]
  temp <- readRDS ( file = paste0("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/R-working-files/AMI-R-working-files-indiv/dat.cpd-",this_samp,".rds") ) # use readRDS()
  
  dat <- rbind(dat, temp)
  
  print(paste0("completed ",i))
}


saveRDS(object = dat, file = "dat.cpd-long-all-samps-cpp3d-ami.rds" )
dat <- readRDS("dat.cpd-long-all-samps-cpp3d-ami.rds")

rm(temp)

str(dat)
# 'data.frame':	112565544 obs. of  11 variables:
#   $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_3" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylaconitate isomerase" "2-methylaconitate isomerase" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" ...
# $ rxn_id             : chr  "rxn25278" "rxn25278" "rxn25279" "rxn25279" ...
# $ cpd_id             : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ cpd_name           : chr  "2-methyl-trans-aconitate" "cis-2-Methylaconitate" "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" ...
# $ cpd_form           : chr  "C7H5O6" "C7H5O6" "C7H7O7" "H2O" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.333 0.333 0.333 ...
# $ sample             : chr  "12465" "12465" "12465" "12465" ...
# $ cpd_rel_abun_norm  : num  1.14e-05 1.14e-05 1.76e-05 1.76e-05 1.76e-05 ...

#dat$combos <- paste0(dat$OC_x,"__",dat$HC_y,"__",dat$NC_z)

sum(dat$cpd_rel_abun_norm) # 5520.761
sum(dat$cpd_rel_abun_norm)/nsamples(phy) # 65.72334 = average 65.7% functional relative abundance per sample

length(which(is.na(dat$cpd_rel_abun_norm))==TRUE) # 0
length(which( dat$cpd_rel_abun_norm > 0) == TRUE) #  75050608
length(which( dat$cpd_rel_abun_norm == 0) == TRUE) # 37514936

# so this step does collect some zero relative abundances from the 'otu-table' 

names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"               "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm" "sample"              "cpd_rel_abun_norm" 


length(unique(dat$cpd_id)) # 8515


## Collate compounds within each sample 


# ## Steps run on Pawsey HPC
# 
# # # # # # # # # # # # #
# #
# # collate_compounds() - R script for cpp3d
# #
# # AMI (rerun with new method, using AMI phyloseq data object from 2024 bioenergetic mapping study)
# # Craig Liddicoat - Flinders University
# # Running on Pawsey Setonix
# # # # # # # # # # # # #
# 
# # Add a new path
# .libPaths(c("/software/projects/pawsey1216/cliddicoat/setonix/2024.05/r/4.4.1",
#             "/software/projects/pawsey1216/cliddicoat/setonix/2024.05/r/4.3", .libPaths()))
# 
# R.Version()
# 
# # load packages
# #library(readxl); packageVersion("readxl")
# library(parallel); packageVersion("parallel")
# library(doParallel); packageVersion("doParallel")
# library(dplyr); packageVersion("dplyr")
# library(stringr); packageVersion("stringr")
# library(phyloseq); packageVersion("phyloseq") # '1.44.0'
# 
# no_forks <- 8
# 
# workdir <- "/scratch/pawsey1216/cliddicoat/ami"
# setwd(workdir)
# temp_dir <- "/scratch/pawsey1216/cliddicoat/ami/working"
# 
# this_study <- "-AMI"
# phy <- readRDS("phy-phyloseq-object-AMI.RDS") # phyloseq object for SUPER-FOCUS functions
# dat <- readRDS("dat.cpd-long-all-samps-cpp3d-ami.rds")
# 
# str(dat)
# 
# sum(dat$cpd_rel_abun_norm)
# sum(dat$cpd_rel_abun_norm)/nsamples(phy) # average functional relative abundance per sample
# 
# length(which(is.na(dat$cpd_rel_abun_norm))==TRUE)
# length(which( dat$cpd_rel_abun_norm > 0) == TRUE)
# length(which( dat$cpd_rel_abun_norm == 0) == TRUE)
# 
# # so this step does collect some zero relative abundances from the 'otu-table' 
# 
# names(dat)
# length(unique(dat$cpd_id)) # 8549
# 
# 
# ## Collate compounds within each sample 
# 
# unique_cpd <- unique(dat$cpd_id)
# samp_names <- sample_names(phy)
# 
# collate_compounds <- function(dat.cpd, unique_cpd, samp) {
#   #i<-1
#   #samp = samp_names[i]
#   #dat.cpd = dat[which(dat$sample == samp_names[i]), ]
#   
#   this_samp <- samp
#   
#   cpd_data <- data.frame(cpd_id = unique_cpd, sample=this_samp, cpd_rel_abun=NA)
#   
#   for (c in 1:length(unique_cpd)) {
#     #c<-1
#     this_cpd <- unique_cpd[c]
#     sel.cpd <- which(dat.cpd$cpd_id == this_cpd)
#     
#     if (length(sel.cpd) >=1) {
#       cpd_data$cpd_rel_abun[c] <- sum(dat.cpd$cpd_rel_abun_norm[sel.cpd])
#     }
#     
#   } # END all compounds
#   
#   saveRDS(object = cpd_data, file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") ) # use readRDS()
#   
# } # END
# 
# 
# #time.start <- Sys.time()
# #cl<-makeCluster( detectCores()-2 )
# #cl<-makeCluster( no_cores )
# 
# # this makes clusters on Unix-like system (may need to use other alternative for Windows)
# cl<-makeForkCluster(nnodes = no_forks)      # no of nodes will depend on your HPC facility
# registerDoParallel(cl)
# 
# foreach(i=1:length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
#   collate_compounds(dat.cpd = dat[which(dat$sample == samp_names[i]), ], unique_cpd = unique_cpd, samp = samp_names[i])
# 
# stopCluster(cl)
# #time.finish <- Sys.time()
# 
# 
# # output 1
# i<-1
# this_samp <- sample_names(phy)[i]
# dat.cpd.collate <- readRDS ( file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") ) # use readRDS()
# head(dat.cpd.collate)
# 
# for ( i in 2:length(sample_names(phy)) ) {
#   #i<-1
#   this_samp <- sample_names(phy)[i]
#   temp <- readRDS ( file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") ) # use readRDS()
#   
#   dat.cpd.collate <- rbind(dat.cpd.collate, temp)
#   
#   print(paste0("completed ",i))
# }
# 
# 
# str(dat.cpd.collate)
# 
# sum(dat.cpd.collate$cpd_rel_abun)
# 
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# 
# saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-ami.rds" )
# 
# # END

str(dat.cpd.collate)
# 'data.frame':	715260 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "12465" "12465" "12465" "12465" ...
# $ cpd_rel_abun: num  0.000155 0.017153 0.000263 4.789375 0.012866 ...
sum(dat.cpd.collate$cpd_rel_abun)
# [1] 5520.761
sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# [1] 65.72334




# data object retrieved from Pawsey HPC

dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ami.rds")


str(dat.cpd.collate)
# 'data.frame':	715260 obs. of  3 variables:
# $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "12465" "12465" "12465" "12465" ...
# $ cpd_rel_abun: num  0.000155 0.017153 0.000263 4.789375 0.012866 ...

sum(dat.cpd.collate$cpd_rel_abun) # 5520.761
sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample)) # 65.72334

#saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-ami.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ami.rds")

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000012 0.000162 0.007719 0.001458 8.223545

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf -4.9338 -3.7894    -Inf -2.8363  0.9151 


# log10 abun
dat.cpd.collate$log10_abun <- dat.cpd.collate$cpd_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(dat.cpd.collate$log10_abun == 0) # qty 43765
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(dat.cpd.collate$log10_abun[ -subsel.zero ])
  dat.cpd.collate$log10_abun[ subsel.zero ] <- zero_replace
}
dat.cpd.collate$log10_abun <- log10(dat.cpd.collate$log10_abun)

hist(dat.cpd.collate$log10_abun); summary( dat.cpd.collate$log10_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -8.6204 -4.9338 -3.7894 -4.1318 -2.8363  0.9151



# dat.sel <- dat.cpd.collate[ which(dat.cpd.collate$cpd_rel_abun > 0), ]
# n <- length(unique(dat.sel$sample))
# sum(dat.sel$abun) # 982.155
# sum(dat.sel$abun)/n # 65.477
# dat.sel$abun <- dat.sel$abun/n


# make group variable from sample name

dat.cpd.collate$group <- NA

class(dat.cpd.collate$sample) # "character"

for (i in 1:length(sample_names(phy))) {
  #i<-1
  this_samp <- as.character(sample_names(phy)[i])
  sel <- which(dat.cpd.collate$sample == this_samp)
  dat.cpd.collate$group[sel] <- phy@sam_data$nat_vs_dist[i]
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group) # "Natural"   "Disturbed"
dat.cpd.collate$group <- factor(dat.cpd.collate$group, levels = c("Disturbed", "Natural"), ordered = TRUE)


saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-ami.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-ami.rds")



str(dat.cpd.collate)
# 'data.frame':	715260 obs. of  5 variables:
# $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "12465" "12465" "12465" "12465" ...
# $ cpd_rel_abun: num  0.000155 0.017153 0.000263 4.789375 0.012866 ...
# $ log10_abun  : num  -3.81 -1.77 -3.58 0.68 -1.89 ...
# $ group       : Ord.factor w/ 2 levels "Disturbed"<"Natural": 2 2 2 2 2 2 2 2 2 2 ...


length( unique(dat.cpd.collate$cpd_id) ) # 8515
8515*84 # 715260



### Test trends in disturbed vs natural groups for key compounds: Sugars, Lignin, BCFA-ACPs
# L-arabinose, melibiose, D-fructose (, 6-Phosphosucrose) ; Lignin ; 35 x BCFA-ACPs

dat.cpd.ami <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-ami.rds")
str(dat.cpd.ami)
# same as above

unique(dat.cpd.ami$group)
# [1] Natural   Disturbed
# Levels: Disturbed < Natural

# compound info is in here
names(df.comp)
# [1] "id"        "abbrev"    "name"      "form"      "OC_ratio"  "HC_ratio"  "NC_ratio"  "PC_ratio"  "NP_ratio"  "O_count"   "N_count"   "P_count"   "S_count"   "mass"     
# [15] "SC_ratio"  "MgC_ratio" "ZnC_ratio" "KC_ratio"  "CaC_ratio" "MnC_ratio" "FeC_ratio" "CoC_ratio" "CuC_ratio" "MoC_ratio"
head(df.comp$abbrev) # "h2o"   "atp"   "nad"   "nadh"  "nadph" "nadp"
head(df.comp$name) # "H2O"   "ATP"   "NAD"   "NADH"  "NADPH" "NADP" 
head(df.comp$form) # "H2O"           "C10H13N5O13P3" "C21H26N7O14P2" "C21H27N7O14P2" "C21H26N7O17P3" "C21H25N7O17P3"


## "arabinose"

sel.cpd <- which(df.comp$name == "L-Arabinose") # 
df.comp[sel.cpd, ]
#           id abbrev        name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio
# 223 cpd00224 arab-L L-Arabinose C5H10O5
this_var <- "L-Arabinose"

# AMI
sel <- which(dat.cpd.ami$cpd_id == "cpd00224") # 84
subsel.dist <- which(dat.cpd.ami$group[sel] == "Disturbed")
subsel.nat <- which(dat.cpd.ami$group[sel] == "Natural")
dat.cpd.ami[sel ,]
dat.cpd.ami[sel[subsel.dist] ,]
dat.cpd.ami[sel[subsel.nat] ,]
x <- dat.cpd.ami$cpd_rel_abun[sel[subsel.dist]] # n = 29
y <- dat.cpd.ami$cpd_rel_abun[sel[subsel.nat]] # n = 55
#x <- dat.cpd.ami$log10_abun[sel[subsel.dist]] # n =
#y <- dat.cpd.ami$log10_abun[sel[subsel.nat]] # n =

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Disturbed vs Natural\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.ami[sel, ] , aes(x=group, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Group")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(a)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-AMI-v5.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## D-fructose 
## Fructose
sel.cpd <- which(df.comp$name == "D-Fructose") # 
df.comp[sel.cpd, ]
#          id abbrev       name    form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio
# 82 cpd00082    fru D-Fructose C6H12O6
this_var <- "D-Fructose"

# AMI
sel <- which(dat.cpd.ami$cpd_id == "cpd00082") # 84
subsel.dist <- which(dat.cpd.ami$group[sel] == "Disturbed")
subsel.nat <- which(dat.cpd.ami$group[sel] == "Natural")
dat.cpd.ami[sel ,]
dat.cpd.ami[sel[subsel.dist] ,]
dat.cpd.ami[sel[subsel.nat] ,]
x <- dat.cpd.ami$cpd_rel_abun[sel[subsel.dist]] # n = 29
y <- dat.cpd.ami$cpd_rel_abun[sel[subsel.nat]] # n = 55
#x <- dat.cpd.ami$log10_abun[sel[subsel.dist]] # n =
#y <- dat.cpd.ami$log10_abun[sel[subsel.nat]] # n =

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Disturbed vs Natural\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.ami[sel, ] , aes(x=group, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Group")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(b)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-AMI-v5.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## "Melibiose" = cpd03198 
sel.cpd <- which(df.comp$name == "Melibiose") # 
df.comp[sel.cpd, ]
# id    abbrev      name      form  OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count     mass SC_ratio MgC_ratio
# 3139  cpd03198     melib Melibiose C12H22O11 0.9166667 1.833333        0        0      NaN      11       0       0       0      342        0         0
# 27105 cpd29840 Melibiose Melibiose      null       NaN      NaN
this_var <- "Melibiose"

# AMI
sel <- which(dat.cpd.ami$cpd_id == "cpd03198") # 84
subsel.dist <- which(dat.cpd.ami$group[sel] == "Disturbed")
subsel.nat <- which(dat.cpd.ami$group[sel] == "Natural")
dat.cpd.ami[sel ,]
dat.cpd.ami[sel[subsel.dist] ,]
dat.cpd.ami[sel[subsel.nat] ,]
x <- dat.cpd.ami$cpd_rel_abun[sel[subsel.dist]] # n = 29
y <- dat.cpd.ami$cpd_rel_abun[sel[subsel.nat]] # n = 55
#x <- dat.cpd.ami$log10_abun[sel[subsel.dist]] # n =
#y <- dat.cpd.ami$log10_abun[sel[subsel.nat]] # n =

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Disturbed vs Natural\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.ami[sel, ] , aes(x=group, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Group")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(c)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-AMI-v5.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")





## "6-Phosphosucrose"
sel.cpd <- which(df.comp$name == "6-Phosphosucrose") # 
df.comp[sel.cpd, ]
# id abbrev             name       form OC_ratio HC_ratio NC_ratio   PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio
# 1658 cpd01693  suc6p 6-Phosphosucrose C12H21O14P 1.166667     1.75        0 0.08333333
this_var <- "6-Phosphosucrose"

# AMI
sel <- which(dat.cpd.ami$cpd_id == "cpd01693") # 84
subsel.dist <- which(dat.cpd.ami$group[sel] == "Disturbed")
subsel.nat <- which(dat.cpd.ami$group[sel] == "Natural")
dat.cpd.ami[sel ,]
dat.cpd.ami[sel[subsel.dist] ,]
dat.cpd.ami[sel[subsel.nat] ,]
x <- dat.cpd.ami$cpd_rel_abun[sel[subsel.dist]] # n = 29
y <- dat.cpd.ami$cpd_rel_abun[sel[subsel.nat]] # n = 55
#x <- dat.cpd.ami$log10_abun[sel[subsel.dist]] # n =
#y <- dat.cpd.ami$log10_abun[sel[subsel.nat]] # n =

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "greater", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Disturbed vs Natural\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.ami[sel, ] , aes(x=group, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Group")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(d)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-AMI-v5.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## lignin
sel.cpd <- which(df.comp$name == "Lignin") # 
df.comp[sel.cpd, ]
#             id abbrev   name form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio CaC_ratio MnC_ratio
# 12548 cpd12745 Lignin Lignin null
this_var <- "Lignin"

# AMI
sel <- which(dat.cpd.ami$cpd_id == "cpd12745") # 84
subsel.dist <- which(dat.cpd.ami$group[sel] == "Disturbed")
subsel.nat <- which(dat.cpd.ami$group[sel] == "Natural")
dat.cpd.ami[sel ,]
dat.cpd.ami[sel[subsel.dist] ,]
dat.cpd.ami[sel[subsel.nat] ,]
x <- dat.cpd.ami$cpd_rel_abun[sel[subsel.dist]] # n = 29
y <- dat.cpd.ami$cpd_rel_abun[sel[subsel.nat]] # n = 55
#x <- dat.cpd.ami$log10_abun[sel[subsel.dist]] # n =
#y <- dat.cpd.ami$log10_abun[sel[subsel.nat]] # n =

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Disturbed vs Natural\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot(data = dat.cpd.ami[sel, ] , aes(x=group, y = cpd_rel_abun))+ # log10_abun cpd_rel_abun
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Group")+ ylab("CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(e)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-AMI-v5.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## BCFA-ACPs - Sum rel abun of 35 x  BCFA-ACPs

this_var <- "BCFA-ACPs"

length(t2d.zoom1.decreasing.proteins)
# [1] 35

t2d.zoom1.decreasing.proteins
# [1] "cpd11499" "cpd11503" "cpd11507" "cpd11511" "cpd11524" "cpd11528" "cpd11532" "cpd11536" "cpd11549" "cpd11553" "cpd11557" "cpd11561" "cpd11475"
# [14] "cpd11465" "cpd11469" "cpd11473" "cpd11498" "cpd11502" "cpd11506" "cpd11510" "cpd11514" "cpd11518" "cpd11523" "cpd11527" "cpd11531" "cpd11535"
# [27] "cpd11539" "cpd11543" "cpd11548" "cpd11552" "cpd11556" "cpd11560" "cpd11564" "cpd11568" "cpd11572"

df <- dat.cpd.ami
length(unique(df$sample)) # 84
length(unique(df$cpd_id)) # 8515

sel <- which(df$cpd_id %in% t2d.zoom1.decreasing.proteins)
df <- df[sel, ]
length(unique(df$cpd_id)) # 35
35*84 # 2940

str(df)
# 'data.frame':	2940 obs. of  5 variables:
# $ cpd_id      : chr  "cpd11465" "cpd11469" "cpd11473" "cpd11475" ...
# $ sample      : chr  "12465" "12465" "12465" "12465" ...
# $ cpd_rel_abun: num  0.000831 0.000831 0.000831 0.000932 0.000831 ...
# $ log10_abun  : num  -3.08 -3.08 -3.08 -3.03 -3.08 ...
# $ group       : Ord.factor w/ 2 levels "Disturbed"<"Natural": 2 2 2 2 2 2 2 2 2 2 ...

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group[i] <- as.character(unique(df$group[subsel]))
  
  print(paste0("completed ",i))
}

str(res)
# 'data.frame':	82 obs. of  3 variables:

unique(res$group) # "Natural"   "Disturbed"
res$group <- factor(res$group, levels = c("Disturbed", "Natural"), ordered = TRUE)


x <- res$sum_rel_abun[ res$group == "Disturbed" ] # n = 29
y <- res$sum_rel_abun[ res$group == "Natural" ] # n = 55

# Wilcoxon-Mann-Whitney Test
wmw.test <- wilcox.test(x, y, alternative = "less", paired = FALSE) # based on cumulative rel abun %
pval <- ifelse(test = wmw.test$p.value < 0.001, yes = paste0("P-value < 0.001"), no = paste0("P-value = ",round(wmw.test$p.value,3)) )
test_result <- paste0("Disturbed vs Natural\nWilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,3),"; ",pval)

p <- ggplot( data = res, aes(x = group, y = sum_rel_abun) )+ 
  ggtitle( this_var)+
  
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, alpha = 0.3)+
  xlab("Group")+ ylab("Sum CPP rel abun (%)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = rel(1.1)))+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 3 , lineheight = 0.85)+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p
grid.text(label = "(f)", x = unit(0.035, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Lookfor-Compounds-rel-abun-Wilcoxon-Mann-Whitney-",gsub(pattern="/",replacement="-", x = this_var),"-AMI-v5.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



#-------------------------

#### Check overlap in trends as found matching in Soil-resto, T2D-SWE, T2D-CHN?
#-------------------------

# this object captures compounds identified with matching potential metabolism trends across Degraded soil, T2D-SWE, T2D-CHN
temp.match.degrad <- readRDS("temp.match.degrad-Sig-results-EcoDegradation-consistent-T2D-SWE-CHN-Hostremoval-data-linked-fxns-subsyL3-v6.RDS")


table(temp.match.degrad$trend_with_eco_degradation)
# Decreasing Increasing 
# 60         16 

check_cpds <- unique(temp.match.degrad$cpd) # qty 76
identical(temp.match.degrad$cpd, check_cpds) # TRUE

# check trends for these in AMI disturbed vs natural soils?
dat.cpd.ami <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-ami.rds")
str(dat.cpd.ami)
# 'data.frame':	715260 obs. of  5 variables:
# $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "12465" "12465" "12465" "12465" ...
# $ cpd_rel_abun: num  0.000155 0.017153 0.000263 4.789375 0.012866 ...
# $ log10_abun  : num  -3.81 -1.77 -3.58 0.68 -1.89 ...
# $ group       : Ord.factor w/ 2 levels "Disturbed"<"Natural": 2 2 2 2 2 2 2 2 2 2 ...

length(unique(dat.cpd.ami$sample)) # 84

#data_in <- dat.cpd.collate
data_in <- dat.cpd.ami
sel <- which(data_in$cpd_id %in% check_cpds) # qty 6384
76*84 # 6384
data_in <- data_in[sel, ]
dim(data_in) # 6384  5

unique(data_in$group)
# [1] Natural   Disturbed
# Levels: Disturbed < Natural

dat.test <- data.frame(cpd = unique(data_in$cpd_id), data_for_this_cpd=NA , 
                       #median_t2d = NA, median_normal = NA,
                       median_dist = NA, median_nat = NA,  
                       alt = NA,
                       p_val = NA, W_statistic = NA, hl_effect_wilcox = NA,
                       trend_with_dist = NA
)

for (i in 1:dim(dat.test)[1]) {
  #i<-1
  this_cpd <- dat.test$cpd[i]
  sel <- which(data_in$cpd_id == this_cpd)
  # prepare data
  df = data.frame(group = as.character( data_in$group[sel]), 
                  value = data_in$log10_abun[sel], 
                  value_perc = data_in$cpd_rel_abun[sel] )
  
  x <- df$value[df$group == "Disturbed"]
  y <- df$value[df$group == "Natural"]
  
  if ( length(which(x == min(x) )) > 0.5*length(x) & length(which(y == min(y) )) > 0.5*length(y) ) {
    # low data: if at least one dataset does not have at least 50% non-zero cases
    dat.test$data_for_this_cpd[i] <- "low data"
    
  } else {
    
    # # test for homogeneity of variances
    # var.test(x,y) # e.g., j<-100: F = 2.3477, num df = 32, denom df = 42, p-value = 0.009872
    
    dat.test$median_dist[i] <- median( x ) # log10 abun  
    dat.test$median_nat[i] <- median( y ) # log10 abun
    
    alt <- NA
    if (median( x ) > median( y ) ) {
      alt <- "greater"
    } else {
      alt <- "less"
    }
    dat.test$alt[i] <- alt
    
    # Wilcoxon-Mann-Whitney Test
    wmw.test <- wilcox.test(x, y, alternative = alt, paired = FALSE, conf.int = TRUE) # based on log10 abun
    
    dat.test$p_val[i] <- wmw.test$p.value
    dat.test$W_statistic[i] <- wmw.test$statistic
    
    # https://search.r-project.org/CRAN/refmans/DescTools/html/HodgesLehmann.html
    # https://aakinshin.net/posts/r-hodges-lehmann-problems
    
    dat.test$hl_effect_wilcox[i] <- wmw.test$estimate
    #dat.test$hl_effect_hlfxn[i] <- hl(x, y)
    
    if (!(is.na(dat.test$p_val[i])|is.na(dat.test$W_statistic[i]))) {
      if (dat.test$p_val[i] <= 0.05 & alt == "greater") { dat.test$trend_with_dist[i] <- "Increasing" }
      else if (dat.test$p_val[i] <= 0.05 & alt == "less") { dat.test$trend_with_dist[i] <- "Decreasing" }
      else { dat.test$trend_with_dist[i] <- "No trend" }
    }
  }
  print(paste0("Completed ",i))
}


sel.low <- which(dat.test$data_for_this_cpd == "low data") # empty
# length(unique(dat.cpd.collate$cpd_id)) #
# length(unique(dat.cpd.collate$cpd_id)) - length(sel.low) #
# sel.nonNA <- which(!is.na(dat.test$p_val)) # xx applicable tests

head(dat.test)

#dat.test %>% arrange(., p_val) %>% slice(1:50)

identical(dat.test$cpd, check_cpds) # FALSE
# rearrange order to match temp.match.degrad$cpd
length(dat.test$cpd) # 76
length(check_cpds) #76
row.names(dat.test) <- dat.test$cpd
dat.test <- dat.test[check_cpds, ]
identical(dat.test$cpd, check_cpds) # TRUE
identical(temp.match.degrad$cpd, dat.test$cpd)
names(dat.test)
# [1] "cpd"               "data_for_this_cpd" "median_dist"       "median_nat"        "alt"               "p_val"            
# [7] "W_statistic"       "hl_effect_wilcox"  "trend_with_dist"    
names(temp.match.degrad)
# [1] "cpd"                                   "median_log10_cpp_relabun_t2d_SWE"      "median_log10_cpp_relabun_normal_SWE"  
# [4] "alt_wilcox_SWE"                        "p_val_wilcox_SWE"                      "W_statistic_SWE"                      
# [7] "hl_effect_wilcox_SWE"                  "trend_with_t2d_SWE"                    "minuslog10_p_val_wilcox_SWE"          
# [10] "cpd_names"                             "cpd_forms"                             "OC_x"                                 
# [13] "HC_y"                                  "NC_z"                                  "z_layer"                              
# [16] "sub_fxns_SWE"                          "subsys_L3_SWE"                         "no_linked_superfocus_fxns_SWE"        
# [19] "no_linked_subsys_L3_SWE"               "p_val_kendall_Restoration"             "kendall_tau_Restoration"              
# [22] "trend_with_age_Restoration"            "minuslog10_p_val_kendall_Restoration"  "sub_fxns_Restoration"                 
# [25] "subsys_L3_Restoration"                 "no_linked_superfocus_fxns_Restoration" "no_linked_subsys_L3_Restoration"      
# [28] "median_log10_cpp_relabun_t2d_CHN"      "median_log10_cpp_relabun_normal_CHN"   "alt_wilcox_CHN"                       
# [31] "p_val_wilcox_CHN"                      "W_statistic_CHN"                       "hl_effect_wilcox_CHN"                 
# [34] "trend_with_t2d_CHN"                    "minuslog10_p_val_wilcox_CHN"           "sub_fxns_CHN"                         
# [37] "subsys_L3_CHN"                         "no_linked_superfocus_fxns_CHN"         "no_linked_subsys_L3_CHN"              
# [40] "trend_with_eco_degradation"            "minuslog10_p_val_wilcox_T2D_mean"      "no_linked_subsys_L3_T2D_mean"
temp.match.degrad$trend_with_eco_degradation
temp.match.degrad$trend_with_t2d_SWE
temp.match.degrad$trend_with_t2d_CHN
dat.test$trend_with_dist

identical(dat.test$trend_with_dist, temp.match.degrad$trend_with_eco_degradation) # FALSE
# i.e. not a complete match
identical(temp.match.degrad$trend_with_t2d_SWE, temp.match.degrad$trend_with_eco_degradation) # TRUE
identical(temp.match.degrad$trend_with_t2d_CHN, temp.match.degrad$trend_with_eco_degradation) # TRUE

sel <- which(dat.test$cpd %in% t2d.zoom1.decreasing.proteins) # 35
identical(dat.test$trend_with_dist[sel], temp.match.degrad$trend_with_eco_degradation[sel]) # TRUE
# i.e. matching trends for all 35 BCFA-ACPs

sel <- which(dat.test$trend_with_dist == temp.match.degrad$trend_with_eco_degradation) # qty 65
100*length(sel)/length(dat.test$trend_with_dist) 
# i.e., 65 out of 76 or 85.5 % of compounds had matching potential metabolism trends in AMI disturbed soils
# which compounds had matching trends in AMI disturbed soils?
temp.match.degrad$cpd_names[sel]
# [1] "L-Arabinose"                                               "3-Octaprenyl-4-hydroxybenzoate"                            "2-Octaprenylphenol"                                       
# [4] "Pyrroloquinoline-quinone"                                  "PQQH2"                                                     "R-CN"                                                     
# [7] "3-chlorobenzyl alcohol"                                    "Melibiose"                                                 "Protocatechuate"                                          
# [10] "D-Fructose"                                                "6-Phosphosucrose"                                          "Phenolic-Donors"                                          
# [13] "Phenoxyl-rad-of-phenolic-donors"                           "betanidin quinone"                                         "Lignin"                                                   
# [16] "Sinapyl alcohol"                                           "p-Coumaryl alcohol"                                        "Baicalein"                                                
# [19] "6,7-dehydrobaicalein"                                      "Molybdoenzyme molybdenum cofactor"                         "Chenodiol"                                                
# [22] "Molybdate"                                                 "Adenylated molybdopterin"                                  "(2E)-Decenoyl-[acp]"                                      
# [25] "But-2-enoyl-[acyl-carrier protein]"                        "(2E)-Dodecenoyl-[acp]"                                     "(2E)-Hexenoyl-[acp]"                                      
# [28] "4-methyl-trans-hex-2-enoyl-ACP"                            "4-methyl-hexanoyl-ACP"                                     "6-methyl-trans-oct-2-enoyl-ACP"                           
# [31] "6-methyl-octanoyl-ACP"                                     "8-methyl-trans-dec-2-enoyl-ACP"                            "8-methyl-decanoyl-ACP"                                    
# [34] "10-methyl-trans-dodec-2-enoyl-ACP"                         "10-methyl-dodecanoyl-ACP"                                  "12-methyl-trans-tetra-dec-2-enoyl-ACP"                    
# [37] "14-methyl-trans-hexa-dec-2-enoyl-ACP"                      "5-methyl-trans-hex-2-enoyl-ACP"                            "5-methyl-hexanoyl-ACP"                                    
# [40] "7-methyl-trans-oct-2-enoyl-ACP"                            "7-methyl-octanoyl-ACP"                                     "9-methyl-trans-dec-2-enoyl-ACP"                           
# [43] "9-methyl-decanoyl-ACP"                                     "11-methyl-trans-dodec-2-enoyl-ACP"                         "11-methyl-dodecanoyl-ACP"                                 
# [46] "13-methyl-trans-tetra-dec-2-enoyl-ACP"                     "15-methyl-trans-hexa-dec-2-enoyl-ACP"                      "4-methyl-trans-pent-2-enoyl-ACP"                          
# [49] "4-methyl-pentanoyl-ACP"                                    "6-methyl-trans-hept-2-enoyl-ACP"                           "6-methyl-heptanoyl-ACP"                                   
# [52] "8-methyl-trans-non-2-enoyl-ACP"                            "8-methyl-nonanoyl-ACP"                                     "10-methyl-trans-undec-2-enoyl-ACP"                        
# [55] "10-methyl-undecanoyl-ACP"                                  "12-methyl-trans-tridec-2-enoyl-ACP"                        "14-methyl-trans-pentadec-2-enoyl-ACP"                     
# [58] "trans-Octodec-2-enoyl-ACP"                                 "Phenylcarbinol"                                            "Acetone"                                                  
# [61] "PRPn"                                                      "5-phospho-alpha-D-ribosyl 1,2-cyclic phosphate"            "alpha-D-ribose-1-[N-(phosphonomethyl)glycine] 5-phosphate"
# [64] "5-Ureido-4-imidazole carboxylate"                          "ARSENOBETAINE"                

# sig results ... (no p-adjustment)
sel.sig <- which(dat.test$p_val <= 0.05) # 65

write.csv(x = dat.test, file = "dat.test-cpp3d-indiv-Compounds-AMI-check-matching-trends-Dist-vs-Nat-Wilcox-v6.csv")

saveRDS(dat.test, file = "dat.test-cpp3d-indiv-Compounds-AMI-check-matching-trends-Dist-vs-Nat-Wilcox-v6.RDS")





# link compounds to functions
df.out <-readRDS( "df.out--tidy-compounds_indiv--cpp3d-AMI.RDS" )
str(df.out)
# 'data.frame':	1340066 obs. of  9 variables:
# $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_3" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylaconitate isomerase" "2-methylaconitate isomerase" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" ...
# $ rxn_id             : chr  "rxn25278" "rxn25278" "rxn25279" "rxn25279" ...
# $ cpd_id             : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ cpd_name           : chr  "2-methyl-trans-aconitate" "cis-2-Methylaconitate" "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" ...
# $ cpd_form           : chr  "C7H5O6" "C7H5O6" "C7H7O7" "H2O" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.333 0.333 0.333 ...

phy <- readRDS("phy-phyloseq-object-AMI.RDS")
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 37335 taxa and 84 samples ]
# sample_data() Sample Data:       [ 84 samples by 92 sample variables ]
# tax_table()   Taxonomy Table:    [ 37335 taxa by 4 taxonomic ranks ]
df.tax <- as.data.frame(phy@tax_table)
head(row.names(df.tax))
dim(df.tax) # 37335     4

table(phy@sam_data$nat_vs_dist)
# Disturbed   Natural 
# 29        55 

min(taxa_sums(phy)) # 1.712637e-06

## tidy sig test results table and add function info
#temp <- dat.test.sig
temp <- dat.test # in this case consider all 76 compounds at the intersection of restoration gradient, T2D (SWE), T2D (CHN)
head(temp)
names(temp)
# [1] "cpd"               "data_for_this_cpd" "median_dist"       "median_nat"        "alt"               "p_val"             "W_statistic"      
# [8] "hl_effect_wilcox"  "trend_with_dist"  
temp <- temp[ ,c(
  "cpd", "median_dist", "median_nat", "alt", "p_val", "W_statistic", "hl_effect_wilcox", "trend_with_dist")]
names(temp) <- c(
  "cpd", "median_dist_AMI", "median_nat_AMI", "alt_AMI", "p_val_AMI", "W_statistic_AMI", "hl_effect_wilcox_AMI", "trend_with_dist_AMI")

dim(temp) #  76 8


# add functions
temp$fxns <- NA
#temp$superfocus_fxns <- NA
temp$sub_fxns <- NA
temp$subsys_L3 <- NA
temp$no_linked_superfocus_fxns <- NA
temp$no_linked_subsys_L3 <- NA

for (i in 1:dim(temp)[1]) {
  #i<-1
  sel <- which(df.out$cpd_id == temp$cpd[i]) # 2
  ( list_fxns <- unique( df.out$superfocus_fxn[sel] ) ) # "fxn_9950"  "fxn_9951"
  temp$fxns[i] <- paste0(list_fxns, collapse = "  |  ")
  temp$sub_fxns[i] <- paste0( unique(df.out$f__in[sel]), collapse = "  |  ")
  
  sel <- which(row.names(df.tax) %in% list_fxns) # 2
  
  #temp$superfocus_fxns[i] <- paste0( unique(df.tax$fxn[sel]), collapse = "  |  ")
  temp$subsys_L3[i] <- paste0( unique(df.tax$subsys_L3[sel]), collapse = "  |  ")
  
  temp$no_linked_superfocus_fxns[i] <- length( unique(df.tax$fxn[sel]) )
  temp$no_linked_subsys_L3[i] <- length( unique(df.tax$subsys_L3[sel]) )
  
  print(paste0("completed ",i))
}
#temp$superfocus_fxns <- gsub(pattern = ",", replacement = ";", x = temp$superfocus_fxns)
temp$sub_fxns <- gsub(pattern = ",", replacement = ";", x = temp$sub_fxns)

head(temp)


saveRDS(object = temp, file = "temp-results-Confirm-Check-AMI-Disturbed-vs-Natural-data-linked-fxns-subsyL3-v6.RDS")
temp <- readRDS("temp-results-Confirm-Check-AMI-Disturbed-vs-Natural-data-linked-fxns-subsyL3-v6.RDS")

dim(temp) # 76 13
row.names(temp)
length(which(row.names(temp) %in% check_cpds)) # 76

colnames(phy@tax_table) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn" 
length(unique(phy@tax_table[ , "subsys_L1" ])) #
length(unique(phy@tax_table[ , "subsys_L2" ])) #
length(unique(phy@tax_table[ , "subsys_L3" ])) #
length(unique(phy@tax_table[ , "fxn" ])) #

write.table(x = temp, file = "temp-results-Confirm-Check-AMI-Disturbed-vs-Natural-data-linked-fxns-subsyL3-v6.tsv", sep = "\t", quote = FALSE, row.names = FALSE )
write.table(x = temp[ order(temp$p_val_AMI, decreasing = FALSE) , ], file = "temp-results-Confirm-Check-AMI-Disturbed-vs-Natural-data-linked-fxns-subsyL3-ORDERED-v6.tsv", sep = "\t", quote = FALSE, row.names = FALSE )

# write table in the same order of compounds as this previous output:
# Sig-Matching-SUMMARY-EcoDegradation-T2D-SWE-CHN-Hostremoval-compounds-ordered-by-kendall-tau-restoration-v6.tsv
tab <- read.table(file = "Sig-Matching-SUMMARY-EcoDegradation-T2D-SWE-CHN-Hostremoval-compounds-ordered-by-kendall-tau-restoration-v6.tsv",sep = "\t", header = TRUE)
head(tab)
order_cpds <- rownames(tab)
order_cpds # this matches Supplementary table
# [1] "cpd00224" "cpd00082" "cpd11597" "cpd03198" "cpd01693" "cpd12672" "cpd12673" "cpd12674" "cpd12697" "cpd27089" "cpd27073" "cpd27088"
# [13] "cpd29537" "cpd29536" "cpd29314" "cpd29315" "cpd26519" "cpd01663" "cpd22764" "cpd00197" "cpd00435" "cpd21090" "cpd11475" "cpd26508"
# [25] "cpd32092" "cpd19503" "cpd11465" "cpd11469" "cpd11473" "cpd11498" "cpd11499" "cpd11502" "cpd11503" "cpd11506" "cpd11507" "cpd11510"
# [37] "cpd11511" "cpd11514" "cpd11518" "cpd11523" "cpd11524" "cpd11527" "cpd11528" "cpd11531" "cpd11532" "cpd11535" "cpd11536" "cpd11539"
# [49] "cpd11543" "cpd11548" "cpd11549" "cpd11552" "cpd11553" "cpd11556" "cpd11557" "cpd11560" "cpd11561" "cpd11564" "cpd11568" "cpd11572"
# [61] "cpd03280" "cpd03444" "cpd00178" "cpd00542" "cpd11574" "cpd27851" "cpd27852" "cpd25482" "cpd12745" "cpd01554" "cpd01722" "cpd06913"
# [73] "cpd34699" "cpd03443" "cpd00097" "cpd00986"

head(temp)
write.table(x = temp[ order_cpds , ], file = "temp-results-Confirm-Check-AMI-Disturbed-vs-Natural-data-linked-fxns-subsyL3-ORDER-BY-SUPP-TABLE-v6.tsv", sep = "\t", quote = FALSE, row.names = FALSE )


#-------------------------

#### sample summary stats - AMI
#-------------------------

## AMI
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-indiv-ExtraData-ami.rds")

str(dat.cpd.collate)
# 'data.frame':	715260 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "12465" "12465" "12465" "12465" ...
# $ cpd_rel_abun: num  0.000155 0.017153 0.000263 4.789375 0.012866 ...
# $ log10_abun  : num  -3.81 -1.77 -3.58 0.68 -1.89 ...
# $ group       : Ord.factor w/ 2 levels "Disturbed"<"Natural": 2 2 2 2 2 2 2 2 2 2 ...

length( unique(dat.cpd.collate$cpd_id) ) # 8515
length( unique(dat.cpd.collate$cpd_id[ dat.cpd.collate$cpd_rel_abun > 0] ) ) # 8515
length( unique(dat.cpd.collate$sample) ) # 84

data_in <- dat.cpd.collate

head(data_in)
# cpd_id sample cpd_rel_abun log10_abun   group
# 1 cpd25681  12465 0.0001549552 -3.8097940 Natural
# 2 cpd02597  12465 0.0171531030 -1.7656573 Natural
# 3 cpd24620  12465 0.0002631835 -3.5797413 Natural
# 4 cpd00001  12465 4.7893746559  0.6802788 Natural
# 5 cpd01501  12465 0.0128658550 -1.8905613 Natural
# 6 cpd00851  12465 0.0095119192 -2.0217318 Natural

dim(data_in) # 715260      5
8515*84 # 715260

unique_samps <- unique(data_in$sample)

no_compounds <- numeric(length = length(unique_samps))
sample_sum_relabun <- numeric(length = length(unique_samps))

for (i in 1:length(unique_samps)) {
  #i<-1
  this_samp <- unique_samps[i]
  sel <- which(data_in$sample == this_samp)
  
  values <- data_in$cpd_rel_abun[sel]
  values <- values[values > 0]
  
  no_compounds[i] <- length( values )
  sample_sum_relabun[i] <- sum(values)
  print(paste0("completed ",i))
}

mean(no_compounds) # 7993.988
sd(no_compounds) #  127.2065

mean(sample_sum_relabun) # 65.72334
sd(sample_sum_relabun) # 1.239307



phy <- readRDS("phy-phyloseq-object-AMI.RDS")
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 37335 taxa and 84 samples ]
# sample_data() Sample Data:       [ 84 samples by 92 sample variables ]
# tax_table()   Taxonomy Table:    [ 37335 taxa by 4 taxonomic ranks ]

head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 84
NonZeroFxns

mean(NonZeroFxns) # 22958.61
sd(NonZeroFxns) # 1171.752
#rm(NonZeroFxns)

# check
unique_samps <- sample_names(phy)

unique_samps[1] # "12465"
a <- prune_samples(unique_samps[1], phy)
min(taxa_sums(a)) # 0
# prune taxa that have zero sequence reads
a <- prune_taxa(taxa = taxa_sums(a) > 0, x = a)
a

unique_samps[15] # "12507"
a <- prune_samples(unique_samps[15], phy)
min(taxa_sums(a)) # 0
# prune taxa that have zero sequence reads
a <- prune_taxa(taxa = taxa_sums(a) > 0, x = a)
a



#-------------------------

## END OF FILE
