# Written by Chul Lee (c) Seoul National Univ. email: chul.bioinfo@gmail.com
# Library
import sys
import os
import glob

# Global variables
AA_SET        = ["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","U","O","-",'*',"X"] 
CODONTABLE_dic  = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L','ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P','ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A','TAT':'Y','TAC':'Y','TAA':'Z','TAG':'Z','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E','TGT':'C','TGC':'C','TGA':'Z','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',"---":"-",'CTN':'L','GTN':'V','TCN':'S','CCN':'P','ACN':'T','GCN':'A','CGN':'R','GGN':'G'} # Stop codon = 'Z', NNN = 'X'


# Functions
def read_seq_file(fNAME_seq):
  sID_nSeq_dic = {}
  fpin = open(fNAME_seq,'r')
  for line in fpin:
    if line[0]==">":
      sID = line[1:].strip('\n')
      sID_nSeq_dic.setdefault(sID,'')
    else:
      sID_nSeq_dic[sID]+=line.strip()  
  fpin.close()
  return(sID_nSeq_dic)


def read_tree_file(fNAME_tree):
  fpin = open(fNAME_tree,'r')
  nTree = ''
  for line in fpin:
    nTree += line.strip()
  fpin.close()
  return(nTree)

    
def read_tre_removing_branchlength(nTree):
  line = nTree.strip(";")
  if ":" in line:
    tmp_node_winfo_list = []
    node_winfo_list = line.split(")")
    for node_winfo in node_winfo_list:
      if ',' in node_winfo:
        nodeID_list = []
        for eachnode_winfo in node_winfo.split(","):
          nodeID = eachnode_winfo.split(":")[0]
          nodeID_list.append(nodeID)
        tmp_node_winfo_list.append(','.join(nodeID_list))
      else:
        nodeID = node_winfo.split(":")[0]
        tmp_node_winfo_list.append(nodeID)
    line = ')'.join(tmp_node_winfo_list)+';'
  return line


def terminal_node_from_tree(nTree):
  ## list of terminal nodes
  terNode_list = []
  for terInfo in nTree.split(","):
    if "(" in terInfo:
      sID = terInfo.split("(")[-1]
    elif ")" in terInfo:
      sID = terInfo.split(")")[0]
    else:
      pass
    terNode_list.append(sID)
  return(terNode_list)


def ancestral_node_from_tree(nTree):
  ## list of ancestral nodes
  ancNode_list = []
  for ancInfo in nTree.split(")"):
    if "(" in ancInfo:
      pass
    else:
      if ',' in ancInfo:
        sID = ancInfo.split(",")[0]
      elif ';' in ancInfo:
        sID = ancInfo.split(";")[0]
      else:
        sID = ancInfo.strip()
      ancNode_list.append(sID)
  return(ancNode_list)


def reordering_targetID_list(terNode_list, targetID_list, MonophyleticPair_list):
  tmp_targetID_list = []
  mono_list = []
  poly_list = []
  for targetID in terNode_list:
    if targetID in targetID_list:
      if targetID in MonophyleticPair_list:
        mono_list.append(targetID)
      else:
        poly_list.append(targetID)
  for targetID in mono_list:
    tmp_targetID_list.append(targetID)
  for targetID in poly_list:
    tmp_targetID_list.append(targetID)
  return(tmp_targetID_list)


def make_mID_sID_dic(partial_nTree, mID_sID_dic): ## dictionary maternal and sister nodes
  # mother ID
  mID = partial_nTree.split(")")[-1]
  mID = mID.replace(";",'')
  
  # set partial_tree excluding root
  tmp_partial_nTree = partial_nTree[1:].split(")"+mID)[0]

  # set sister ID 1 and 2
  if not ")" in tmp_partial_nTree:            # case: "s1 , s2"    (mID = a1)
    sID_1 = tmp_partial_nTree.split(",")[0]
    sID_2 = tmp_partial_nTree.split(",")[1]
    mID_sID_dic.setdefault(mID,[sID_1,sID_2])
    return(mID_sID_dic)
  else: #"( and )" is present in partial tree
    sID_info_list = tmp_partial_nTree.split(",")
    FirstID_info  = sID_info_list[0]
    LastID_info   = sID_info_list[-1]
    if FirstID_info.count("(") > 0:
      if LastID_info.count(")") > 0:          # case "(s1,s2)a1 , (s3,s4)a2"
        p1_tree = FirstID_info
        for i in range(1,len(sID_info_list),1):
          p1_tree += ','+sID_info_list[i]
          if p1_tree.count("(") == p1_tree.count(")"):
            p2_tree = tmp_partial_nTree.split(p1_tree+',')[1]
            if not p2_tree.count("(") == p2_tree.count(")"):
              print("# Errorneous_Topology_Of_2nd_sister_lineage:",tmp_partial_nTree)
              sys.exit()
            sID_1 = p1_tree.split(")")[-1]
            sID_2 = p2_tree.split(")")[-1]
            break
        mID_sID_dic.setdefault(mID,[sID_1,sID_2])
        mID_sID_dic = make_mID_sID_dic(p1_tree, mID_sID_dic)
        mID_sID_dic = make_mID_sID_dic(p2_tree, mID_sID_dic)
        return(mID_sID_dic)
      else:                                   # case "(s1,s2)a1 , s3"
        sID_2 = sID_info_list[-1]
        p1_tree = tmp_partial_nTree.split(","+sID_2)[0]
        sID_1 = p1_tree.split(")")[-1]
        mID_sID_dic.setdefault(mID,[sID_1,sID_2])
        mID_sID_dic = make_mID_sID_dic(p1_tree, mID_sID_dic)
        return(mID_sID_dic)
    else:
      if LastID_info.count(")") > 0:          # case "s1 , (s2,s3)a1"
        sID_1 = sID_info_list[0]
        p2_tree = tmp_partial_nTree.split(sID_1+",")[1]
        sID_2 = p2_tree.split(")")[-1]
        mID_sID_dic.setdefault(mID,[sID_1,sID_2])
        mID_sID_dic = make_mID_sID_dic(p2_tree, mID_sID_dic)
        return(mID_sID_dic)
      else:
        print("Error-parsing_tree:",tmp_partial_nTree)
        sys.exit()


def make_sID_mID_dic(mID_sID_dic):
  sID_mID_dic = {}
  for mID in mID_sID_dic.keys():
    for sID in mID_sID_dic[mID]:
      sID_mID_dic.setdefault(sID,mID)
  return(sID_mID_dic)


def remove_ID_from_list(tmp_targetID_list,tmp_targetID):
  iPos_LastElement  = len(tmp_targetID_list)-1
  iPos_tmp_targetID = tmp_targetID_list.index(tmp_targetID) 
  tmp_targetID_list[iPos_tmp_targetID]  = tmp_targetID_list[iPos_LastElement]
  tmp_targetID_list[iPos_LastElement]   = tmp_targetID
  tmp_targetID_list = tmp_targetID_list[:-1]
  return(tmp_targetID_list)
  
  
def find_monophyletic_clade(sID_mID_dic, ref_targetID, targetID_list, targetID_mID_dic):
  ## First time of the reference target
  iFlag = 0
  if not ref_targetID in targetID_mID_dic.keys():
    for tmp_targetID in targetID_list:
      if not ref_targetID == tmp_targetID:
        ### case finding monophyletic pairs within target list
        if sID_mID_dic[ref_targetID] == sID_mID_dic[tmp_targetID]:
          iFlag= 1
          if ref_targetID in targetID_list:
            targetID_list = remove_ID_from_list(targetID_list, ref_targetID)
          targetID_list = remove_ID_from_list(targetID_list, tmp_targetID)
          mID = sID_mID_dic[ref_targetID]
          targetID_mID_dic.setdefault(ref_targetID,mID)
          targetID_mID_dic.setdefault(tmp_targetID,mID)
          return(targetID_list, targetID_mID_dic)
    ### case without monophyletic pairs within target list
    if iFlag == 0:
      mID = sID_mID_dic[ref_targetID]
      targetID_mID_dic.setdefault(ref_targetID,mID)
      targetID_list = remove_ID_from_list(targetID_list, ref_targetID)
      return(targetID_list, targetID_mID_dic)
  ## Next time of the reference target
  else:
    ref_mID = targetID_mID_dic[ref_targetID]
    for tmp_targetID in targetID_list:
      ### case finding additional monophyletic species in target list
      if sID_mID_dic[ref_mID] == sID_mID_dic[tmp_targetID]:
        iFlag = 1
        targetID_list = remove_ID_from_list(targetID_list, tmp_targetID)
        mID = sID_mID_dic[ref_mID]
        for targetID in targetID_mID_dic.keys():
          targetID_mID_dic[targetID] = mID
        targetID_mID_dic.setdefault(tmp_targetID,mID)
        return(targetID_list, targetID_mID_dic)
    ### case without additional monophyletic species within target list
    if iFlag == 0:
      return(targetID_list, targetID_mID_dic)

  
def make_MRCA_targetID_dic(sID_mID_dic, ref_targetID, targetID_list, targetID_mID_dic, mID_targetID_dic, terNode_list, MonophyleticPair_list):
  iLen_targetID_list = len(targetID_list)
  targetID_list, targetID_mID_dic = find_monophyletic_clade(sID_mID_dic, ref_targetID, targetID_list, targetID_mID_dic)
  tmp_iLen_targetID_list = len(targetID_list)

  if not iLen_targetID_list == tmp_iLen_targetID_list:
    ### Recursive function
    targetID_list, targetID_mID_dic, mID_targetID_dic = make_MRCA_targetID_dic(sID_mID_dic, ref_targetID, targetID_list, targetID_mID_dic, mID_targetID_dic, terNode_list, MonophyleticPair_list)
    targetID_list = reordering_targetID_list(terNode_list, targetID_list,MonophyleticPair_list)
    return(targetID_list, targetID_mID_dic, mID_targetID_dic)
  else:
    for targetID in targetID_mID_dic.keys():
      mID = targetID_mID_dic[targetID]
      mID_targetID_dic.setdefault(mID,[])
      mID_targetID_dic[mID].append(targetID)
    targetID_mID_dic = {}
    return(targetID_list, targetID_mID_dic, mID_targetID_dic)
    

def scan_all_targetID(sID_mID_dic, targetID_list, targetID_mID_dic , mID_targetID_dic, terNode_list, MonophyleticPair_list):
  ref_targetID = targetID_list[0]
  targetID_list, targetID_mID_dic, mID_targetID_dic = make_MRCA_targetID_dic(sID_mID_dic, ref_targetID, targetID_list, targetID_mID_dic, mID_targetID_dic, terNode_list, MonophyleticPair_list)
  if len(targetID_list) > 0:
    ### Recursive function
    targetID_list, targetID_mID_dic , mID_targetID_dic = scan_all_targetID(sID_mID_dic, targetID_list, targetID_mID_dic , mID_targetID_dic, terNode_list, MonophyleticPair_list)
  return(targetID_list, targetID_mID_dic , mID_targetID_dic)


def parsing_MonophyleticPair_list(mID_sID_dic,terNode_list):
  MonophyleticPair_list = []
  for mID in mID_sID_dic.keys():
    iCNT_terNode = 0
    for sID in mID_sID_dic[mID]:
      if sID in terNode_list:
        iCNT_terNode +=1
    if iCNT_terNode == 2:
      for sID in mID_sID_dic[mID]:
        MonophyleticPair_list.append(sID)
  return(MonophyleticPair_list)
    

def Anc_Finder(nTree, nTargets):
  targetID_list = nTargets.split(",")
  nTree = read_tre_removing_branchlength(nTree)
  terNode_list = terminal_node_from_tree(nTree)
  ancNode_list = ancestral_node_from_tree(nTree)
  mID_sID_dic = {}
  mID_sID_dic = make_mID_sID_dic(nTree, mID_sID_dic)
  sID_mID_dic = make_sID_mID_dic(mID_sID_dic)
  MonophyleticPair_list = parsing_MonophyleticPair_list(mID_sID_dic,terNode_list)
  targetID_mID_dic = {}
  mID_targetID_dic = {}
  targetID_list, targetID_mID_dic , mID_targetID_dic = scan_all_targetID(sID_mID_dic, targetID_list, targetID_mID_dic , mID_targetID_dic, terNode_list, MonophyleticPair_list)
  mID_mmID_dic = {}
  for mID in mID_targetID_dic.keys():
    if len(mID_targetID_dic[mID]) > 1:
      mmID = sID_mID_dic[mID]
      mID_mmID_dic.setdefault(mID,mmID)
  return(mID_targetID_dic, mID_mmID_dic)


def parse_target_from_tree(nTargets,nTree):
  target_list = nTargets.split(",")
  terNode_list = terminal_node_from_tree(nTree)
  tmp_target_list = []
  for sID in terNode_list:
    if sID in target_list:
      tmp_target_list.append(sID)
  return(tmp_target_list)

def parse_target_others_from_list(nTargets,sID_list):
  target_list = nTargets.split(",")
  tmp_target_list = []
  tmp_others_list = []
  for sID in sID_list:
    if sID in target_list:
      tmp_target_list.append(sID)
    else:
      tmp_others_list.append(sID)
  return(tmp_target_list, tmp_others_list)

  
def set_options():
  # default options
  nTargets = "TAEGU,GEOFO,CORBR,MELUN,NESNO,CALAN"
  templete_nTree = "(((((((((((((TAEGU,GEOFO)60,CORBR)59,MANVI)58,(MELUN,NESNO)61)57,FALPE)56,CARCR)55,(((((((MERNU,PICPU)68,BUCRH)67,APAVI)66,LEPDI)65,COLST)64,TYTAL)63,((HALLE,HALAL)70,CATAU)69)62)54,((((((PELCR,EGRGA)76,NIPNI)75,PHACA)74,(FULGL,(PYGAD,APTFO)78)77)73,GAVST)72,(PHALE,EURHE)79)71)53,((CHAVO,BALRE)81,OPHHO)80)52,(((CALAN,CHAPE)84,CAPCA)83,((CHLUN,TAUER)86,CUCCA)85)82)51,(((MESUN,PTEGU)89,COLLI)88,(PHORU,PODCR)90)87)50,((GALGA,MELGA)92,ANAPL)91)49,(TINMA,STRCA)48)ROOT;"
  iPATH = "../example/"
  fNAME_seq_format = ".fasta"
  fNAME_tree_format = ".tre"
  oPATH = "./"
  # end of default options

  # imported options
  option_list = ["-l=",    # input: List of target species
                 "-ttf=",  # input: file name of templete tree
                 "-ip=",   # input: input path
                 "-tfmt=", # input: file format of binary tree with labelled ancestral and terminal nodes (output of RAxML: )
                 "-sfmt=", # input: file format of reconstruced sequences with ancestral and terminal nodes (output of RAxML: )
                 "-op=",   # output: output path
                 "-h"]     # help

  description_list = ["START OF LOG FILE",
                      "USAGE:  ./ConVarFinder.py [-options]",
                      "|------------------------------ HELP: -----------------------------+",
                      "|-h    help                                                        |",
                      "|-l=   strings of list of target species as comma seperated (csv)  |",
                      "|-ttf=  input file of templete tree                                |",
                      "|-tfmt=  file format of binary tree with ancestral nodes           |",
                      "|-sfmt=  file format of input with ancestral reconstructions       |",
                      "|-ip=  input path                                                  |",
                      "|-op=  output path                                                 |",
                      "+------------------------------------------------------------------+"]

  # update options
  if len(sys.argv) > 1:
    for i in len(1,len(sys.argv),1):
      nOpt = sys.argv[i]
      if not nOpt.split("=")[0]+'=' in option_list:
        print("# Argument error:", nOpt)
        sys.exit()
      if "-l=" in nOpt:
        nTargets = nOpt.split("-l=")[1]
      if "-ttf=" in nOpt:
        templete_fNAME_tree = nOpt.split("-ttf=")[1]
        templete_nTree = read_tree_file(templete_fNAME_tree)
        templete_nTree = read_tre_removing_branchlength(templete_nTree)
      if "-sfmt=" in nOpt:
        fNAME_seq_format = nOpt.split("-sfmt=")[1]
      if "-tfmt=" in nOpt:
        fNAME_seq_format = nOpt.split("-tfmt=")[1]
      if "-ip=" in nOpt:
        iPATH = nOpt.split("-ip=")[1]
      if "-op=" in nOpt:
        oPATH = nOpt.split("-op=")[1]
      if "-h" in nOpt:
        for nDescription in description_list:
          print(nDescrtiption)
        sys.exit()

  flist_seq = glob.glob(iPATH+"*"+fNAME_seq_format)
  flist_tree = glob.glob(iPATH+"*"+fNAME_tree_format)
  templete_mID_targetID_dic, templete_mID_mmID_dic = Anc_Finder(templete_nTree, nTargets)
  templete_sID_list = terminal_node_from_tree(templete_nTree)
  target_list, others_list = parse_target_others_from_list(nTargets, templete_sID_list)
  templete_sID_list = target_list + others_list
  
  # check input files
  iFlag_err = 0
  for fNAME_seq in flist_seq:
    nGene = os.path.basename(fNAME_seq).split(fNAME_seq_format)[0]
    fNAME_tree = os.path.dirname(fNAME_seq) +"/"+ nGene + fNAME_tree_format
    if os.path.isfile(fNAME_tree) == False:
      print("Error-Absent tree:",fNAME_tree)
      iFlag_err = 1
  for fNAME_tree in flist_tree:
    nGene = os.path.basename(fNAME_tree).split(fNAME_tree_format)[0]
    fNAME_seq = os.path.dirname(fNAME_tree) +"/"+ nGene + fNAME_seq_format
    if os.path.isfile(fNAME_seq) == False:
      print("Error-Absent sequence:",fNAME_seq)
      iFlag_err = 1
  if iFlag_err == 1:
    sys.exit()
  option_list = [nTargets,
                 templete_mID_targetID_dic,
                 iPATH,
                 fNAME_seq_format,
                 fNAME_tree_format,
                 oPATH,
                 templete_sID_list]
  return(option_list)

def ConVarFinder(fNAME_seq, terNode_list, mID_targetID_dic, mID_mmID_dic, local_target_list, templete_sID_list, fpout):
  # terNode_list      : terminal nodes of gene tree of partial family tree
  # templete_sID_list : terminal nodes of raw family tree
  sID_nSeq_dic = read_seq_file(fNAME_seq)
  for sID in sID_nSeq_dic.keys():
    iLen_nSeq = len(sID_nSeq_dic[sID])
    break
  if not iLen_nSeq%3 == 0:
    print("Errorneous sequence - not codon-wise:",fNAME_seq)
    sys.exit()
  iPos_Codon = -2
  iPos_AA = 0
  for i in range(0,iLen_nSeq,3):
    iPos_Codon  += 3
    iPos_AA     += 1
    # set codon and aa lists of target species and the others
    tmp_target_codon_list = []
    tmp_others_codon_list = []
    tmp_target_aa_list = []
    tmp_others_aa_list = []
    for sID in terNode_list:
      if sID in local_target_list:
        nCodon_target = sID_nSeq_dic[sID][i:i+3]
        tmp_target_codon_list.append(nCodon_target)
        nAA_target = CODONTABLE_dic[nCodon_target]
        tmp_target_aa_list.append(nAA_target)
      else:
        nCodon_others = sID_nSeq_dic[sID][i:i+3]
        tmp_others_codon_list.append(nCodon_others)
        nAA_others = CODONTABLE_dic[nCodon_others]
        tmp_others_aa_list.append(nAA_others)
        
    # codon substitutions
    tmp_target_codon_list = list(set(tmp_target_codon_list))
    tmp_others_codon_list = list(set(tmp_others_codon_list))
    iFlag_Codon_Sub = 0
    for nCodon_target in tmp_target_codon_list:
      if nCodon_target in tmp_others_codon_list:
        iFlag_Codon_Sub = 1

    # amino acid substitutions
    tmp_target_aa_list = list(set(tmp_target_aa_list))
    tmp_others_aa_list = list(set(tmp_others_aa_list))
    iFlag_AA_Sub = 0
    for nAA_target in tmp_target_aa_list:
      if nAA_target in tmp_others_aa_list:
        iFlag_AA_Sub = 1
        
    if iFlag_Codon_Sub == 0: # if Target-specific codon substitutions
      # nGene                gene name
      # nPhyRel              phylogenetic relationship ex) "mmID4>mID2>target1,target2,target3/..."
      # nPos_codon           start position of codon substitutions (bp)
      # nType_EvoDir_Codon   evolutionary direction: Convergent, Parallel, and Divergent
      # nType_CodonSub       Type_subsitutions: identical, and different
      # nList_EvoDir_Codon   subsitutions ex) "ATG>ATC>ATC,ATC,ATC/..."
      # nCodon_targets       codon set of targets
      # nCodon_others        codon set of the others
      # nPos_AA              position of aa substitutions (pep)
      # nType_EvoDir_AA      evolutionary direction: Convergent, Parallel, and Divergent
      # nType_AASub          Type_subsitutions: identical, and different
      # nList_EvoDir_AA      subsitutions ex) "M>T>T,T,T/..."
      # nAA_targets          AA set of targets
      # nAA_others           AA set of others
      mCodon_list     = []
      mAA_list        = []
      mm_m_Codon_list = []
      mm_m_AA_list    = []
      mm_m_target_ID_list     = []
      mm_m_target_Codon_list  = []
      mm_m_target_AA_list     = []
      for mID in mID_targetID_dic.keys():
        mCodon = sID_nSeq_dic[mID][i:i+3]
        mAA    = CODONTABLE_dic[mCodon]
        tmp_targets_id_list     = mID_targetID_dic[mID]
        if len(tmp_targets_id_list) > 1:
          tmp_targets_codon_list  = []
          tmp_targets_aa_list     = []
          for targetID in tmp_targets_id_list:
            targetCodon = sID_nSeq_dic[targetID][i:i+3]
            tmp_targets_codon_list.append(targetCodon)
            targetAA    = CODONTABLE_dic[targetCodon]
            tmp_targets_aa_list.append(targetAA)
          mCodon_list.append(mCodon)
          mAA_list.append(mAA)
          mmID          = mID_mmID_dic[mID]
          mmCodon = sID_nSeq_dic[mmID][i:i+3]
          mmAA    = CODONTABLE_dic[mmCodon]
          mm_m_Codon_list.append(mmCodon+'>'+mCodon)
          mm_m_AA_list.append(mmAA+'>'+mAA)
        
          mm_m_target_ID_list.append(      mmID+'>'+mID+'>'+','.join(tmp_targets_id_list)           )
          mm_m_target_Codon_list.append(   mmCodon+'>'+mCodon+'>'+','.join(tmp_targets_codon_list)  )
          mm_m_target_AA_list.append(      mmAA+'>'+mAA+'>'+','.join(tmp_targets_aa_list)           )

        else:
          sID = mID_targetID_dic[mID][0]
          sCodon = sID_nSeq_dic[sID][i:i+3]
          sAA = CODONTABLE_dic[sCodon]
          mCodon_list.append(sCodon)
          mAA_list.append(sAA)
          mm_m_Codon_list.append(mCodon+'>'+sCodon)
          mm_m_AA_list.append(mAA+'>'+sAA)
          mm_m_target_ID_list.append(    mID+'>'+sID       )
          mm_m_target_Codon_list.append( mCodon+'>'+sCodon )
          mm_m_target_AA_list.append(    mAA+'>'+sAA       )
          
      mCodon_list     = list(set(mCodon_list))
      mAA_list        = list(set(mAA_list))
      mm_m_Codon_list = list(set(mm_m_Codon_list))
      mm_m_AA_list    = list(set(mm_m_AA_list))

      # set local variables
      ## for gene and species information
      nGene      = os.path.basename(fNAME_seq).split(".")[0]
      nPhyRel    = '/'.join(mm_m_target_ID_list)
      ## for codon information
      nPos_codon = str(iPos_Codon)
      #nType_EvoDir_Codon 
      if len(mCodon_list) == 1:
        if len(mm_m_Codon_list) ==1:
          nType_EvoDir_Codon  = "P" # Parallel:    AAA > TTT / AAA > TTT
        else:
          nType_EvoDir_Codon  = "C" # Convergent:  AAA > TTT / GGG > TTT
      else:
        nType_EvoDir_Codon    = "D" # Divergent:   AAA > TTT / AAA > GGG
      #nType_CodonSub
      if len(tmp_target_codon_list) == 1:
        nType_CodonSub = "i" # Identical
      else:
        nType_CodonSub = "d" # Different
      nList_EvoDir_Codon = '/'.join(mm_m_target_Codon_list)
      nCodon_targets     = ','.join(tmp_targets_codon_list)
      nCodon_others      = ','.join(tmp_others_codon_list)
      
      ## for amino acid information
      if iFlag_AA_Sub == 0: # if Target-specific amino acid substitutions
        nPos_AA            = str(iPos_AA)
        #nType_EvoDir_AA 
        if len(mAA_list) == 1:
          if len(mm_m_AA_list) ==1:
            nType_EvoDir_AA  = "P" # Parallel:    AAA > TTT / AAA > TTT
          else:
            nType_EvoDir_AA  = "C" # Convergent:  AAA > TTT / GGG > TTT
        else:
          nType_EvoDir_AA    = "D" # Divergent:   AAA > TTT / AAA > GGG
        #nType_AASub
        if len(tmp_target_aa_list) == 1:
          nType_AASub = "i" # Identical
        else:
          nType_AASub = "d" # Different
        nList_EvoDir_AA = '/'.join(mm_m_target_AA_list)
        nAA_targets     = ','.join(tmp_target_aa_list)
        nAA_others      = ','.join(tmp_others_aa_list)

      else: # if Target-specific amino acid substitutions
        nPos_AA         = ''
        nType_EvoDir_AA = ''
        nType_AASub     = ''
        nList_EvoDir_AA = ''
        nAA_targets     = ''
        nAA_others      = ''
      # site-wise codon and amino acid sequence information of each species
      nCodon_line = ''
      nAA_line    = ''
      for sID in templete_sID_list:
        if sID in sID_nSeq_dic.keys():
          nCodon = sID_nSeq_dic[sID][i:i+3]
          nAA = CODONTABLE_dic[nCodon]
        else:
          nCodon = ''
          nAA    = ''
        nCodon_line += '\t' + nCodon
        nAA_line    += '\t' + nAA
      tmpline = nGene +'\t'+ nPhyRel +'\t'+ nPos_codon +'\t'+ nType_EvoDir_Codon +'\t'+ nType_CodonSub +'\t'+ nList_EvoDir_Codon +'\t'+ nCodon_targets +'\t'+ nCodon_others
      tmpline +=                      '\t'+ nPos_AA    +'\t'+ nType_EvoDir_AA    +'\t'+ nType_AASub    +'\t'+ nList_EvoDir_AA    +'\t'+ nAA_targets    +'\t'+ nAA_others
      tmpline += nCodon_line + nAA_line +'\n'
      fpout.write(tmpline)
    else: # if not Target-specific codon substitutions
      pass
  return(fpout)
      
      
      

def main(option_list, nGene, fpout):
  # local vaiables for each gene
  nGene                     = nGene
  nTargets                  = option_list[0]
  templete_mID_targetID_dic = option_list[1]
  iPATH                     = option_list[2]
  fNAME_seq_format          = option_list[3]
  fNAME_tree_format         = option_list[4]
  templete_sID_list         = option_list[6]

  fNAME_seq   = iPATH + nGene + fNAME_seq_format
  fNAME_tree  = iPATH + nGene + fNAME_tree_format

  #targetID_list       = nTargets.split(',') 
  nGene_nTree         = read_tree_file(fNAME_tree)
  nGene_nTree         = read_tre_removing_branchlength(nGene_nTree)
  nGene_terNode_list  = terminal_node_from_tree(nGene_nTree)
  nGene_target_list   = parse_target_from_tree(nTargets,nGene_nTree)

  # check_clade
  iCNT_clade = 0
  for mID in templete_mID_targetID_dic.keys():
    iFlag_clade = 0
    for targetID in templete_mID_targetID_dic[mID]:
      if targetID in nGene_target_list:
        iFlag_clade += 1
    if iFlag_clade > 0:
      iCNT_clade += 1
  if not iCNT_clade == len(templete_mID_targetID_dic.keys()):
    return(fpout)
  else:
    # Phylogenetic relationships
    nGene_mID_targetID_dic, nGene_mID_mmID_dic = Anc_Finder(nGene_nTree, ','.join(nGene_target_list))

    # Convergent variant finder
    fpout = ConVarFinder(fNAME_seq, nGene_terNode_list, nGene_mID_targetID_dic, nGene_mID_mmID_dic, nGene_target_list, templete_sID_list, fpout)
  return(fpout)
  
  


if __name__ == "__main__":
  # Input
  option_list       = set_options()
  nTargets          = option_list[0]
  oPATH             = option_list[5]
  iPATH             = option_list[2]
  fNAME_seq_format  = option_list[3]
  templete_sID_list = option_list[6]
  oNAME             = oPATH + "ConVarFinder_"+nTargets + '.txt'


  # analysis for each gene
  flist_seq = glob.glob(iPATH+"*"+fNAME_seq_format)
  fpout = open(oNAME,'w')
  # "Gene"\n"               gene name
  # Phy.relation"\n"        phylogenetic relationship ex) "mmID4>mID2>target1,target2,target3/..."
  # Pos_codon"\n"           start position of codon substitutions (bp)
  # Type_EvoDir_Codon"\n'   evolutionary direction: Convergent, Parallel, and Divergent
  # Type_CodonSub'\n"       Type_subsitutions: identical, and different
  # List_EvoDir_Codon"\n"   subsitutions ex) "ATG>ATC>ATC,ATC,ATC/..."
  # Codon_targets"\n"       codon set of targets
  # Codon_others"\n"        codon set of the others
  # Pos_AA"\n"              position of aa substitutions (pep)
  # Type_EvoDir_AA"\n'      evolutionary direction: Convergent, Parallel, and Divergent
  # Type_AASub'\n"          Type_subsitutions: identical, and different
  # List_EvoDir_AA"\n"      subsitutions ex) "M>T>T,T,T/..."
  # AA_targets"\n"          AA set of targets
  # AA_others"\n            AA set of others
  # "\t".join(templete_sID_list)\n"   list of codon
  # '\t".join(templete_sID_list)+'\n' list of aa
  
  nHeader = "Gene"+'\t'+"Phy.relation"+'\t'+"Pos_codon"+'\t'+"Type_EvoDir_Codon"+'\t'+'Type_CodonSub'+'\t'+"List_EvoDir_Codon"+'\t'+"Codon_targets"+'\t'+"Codon_others"+'\t'+"Pos_AA"+'\t'+"Type_EvoDir_AA"+'\t'+'Type_AASub'+'\t'+"List_EvoDir_AA"+'\t'+"AA_targets"+'\t'+"AA_others"+'\t'+"\t".join(templete_sID_list)+'\t'+"\t".join(templete_sID_list)+'\n'
  fpout.write(nHeader)
  for fNAME_seq in flist_seq:
    nGene = os.path.basename(fNAME_seq).split(fNAME_seq_format)[0]
    fpout = main(option_list, nGene, fpout)
  fpout.close()
