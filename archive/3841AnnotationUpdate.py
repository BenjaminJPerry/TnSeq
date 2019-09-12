# This script reads in an old 3841 annotation file
#     Updates the metadata string in column 7 to have the appropriate
#     information in the appropriate fields (ie. ID=RL0000;Name=hemA;Product=Heme-related protein;)
#     As well, This script will identify any
#     It creates 3 lines for each entry: gene, mRNA, CDS
#
#     Finally, it makes a prot_table used for transit INSeq analysis to have the
#     gene names, and appropriate gene IDs (locus tags)
import os


def read_in_tsv_file(infile):  # Make a function to parse the input files ###


    """" (infile -> list[][])
        Reads in a tab seperated file, removes whitespace and " characters.
        Lines in infile -> rows in list[][]
        Tab delimited elements in infile -> columns in list[][]
    """
    input = infile
    data = []
    for line in input:
        items = line.strip().split('\t')  # Split the line into a string using tabs. Stripping.
        items = list([s.strip() for s in items])  # Stripping whitespace from each element
        items = list([s.strip('"') for s in items])
        data.append(items)
    return data

def get_file():
    """"() -> string
    Returns a user input file
    """
    the_file = input("Designate file: ")
    return the_file


### gff3 file format: seqname   source  feature start   end score   strand  frame   *attributes*
#    *attribute* format:
#    gene    :ID=RL0001;Name=genA;
#    mRNA    :ID=RL0001.t01;Name=RL0001.t01;Parent=RL0001
#    CDS     :"ID=RL0001.p01;Name=RL0001.p01;Parent=RL0001.t01;Dbxref=GI:115254415,EnsemblGenomes-Gn:RL0001,
#               EnsemblGenomes-Tr:CAK05489,GOA:Q1MNF6,InterPro:IPR005177,InterPro:IPR026565;Name=RL0001;
#               Note=conserved hypothetical protein;codon_start=1;product=Hypothetical UPF0085 protein R00001.;
#               protein_id=CAK05489.1;transl_table=11"

### Read in Old NCBI gff3 annotation file
print ("Old GFF3 Annotation File")
# gffInput = get_file()

gffInput = 'RLV3841.Old.txt'

print ("GFF3 File: " + gffInput + '\n')
gffIn = open(str(gffInput), "r")
gff_old = read_in_tsv_file(gffIn)
### Building a library of the old gff3 annotation. Using the protein table.
#   We want the key to be the locus, and the return to be the entire gff3 line for that locus
gff_dict = {}
gff_errors=[]
for entry1 in gff_old: # make a locus_tag:entry1 dictionary for gff file
    if "#" in entry1[0]:
        continue
    elif entry1[2] == "gene":
        attributes = dict(item.split("=") for item in entry1[8].split(";"))
        if 'locus_tag' in attributes:
            locus = attributes['locus_tag'].strip()
            gff_dict.update({locus:entry1})
        else:
            gff_errors.append(entry1)# Some of the new loci are not in the old annotation... So we have to deal with them separately
#            locus = attributes['locus_tag'].strip()
#            gff_new_loci_dict.update({locus: entry1})
    else:
        continue
### gff_dict now contains a dictionary with 'old_locus_tag as' the key and the gff entry for that locus as the value


### Protein Table Format: 
#   Name    Accession   Start   Stop    Strand  GeneID  Locus   Locus_Tag   Protein_Product Length  Protein_Name
### Read in NCBI protein table file
print ("NCBI Exported Protein Table")
#protInput = get_file()

protInput = 'RLV3841ProteinTable.txt'

print ("Protien Table: " + protInput + '\n')
protIn = open(str(protInput), "rU")
prot_table = read_in_tsv_file(protIn)
### Build a prot_table dictionary using locus_tag:prot_table entry2
prot_dict={}
for entry2 in prot_table: # make a locus_tag:entry2 dictionary for gff file
    if "#" in entry2[0]:
        continue
    else:
         locus=entry2[7].strip()
         prot_dict.update({locus:entry2})
### prot_dict now contains a dictionary of the protein table with locus as the key


### pfam table format:
#   #Pfam ID	Pfam    Name	RNB MEDIAN  NC MEDIAN	Transcript induction in nodule		LPD	% total seqs with SignalP	Antismash Pfam (Seconday metabolite)	PID	PCPfam	SignifEnrich
# Read in pfam table
print("Pfam Annotation Table")
#pfamInput = get_file()

pfamInput = 'RLV3841PfamSummary.txt'

print("Pfam Annotation Table: " + pfamInput + '\n')
pfamIn = open(str(pfamInput), 'rU')
pfam_table = read_in_tsv_file(pfamIn)
# Build a dictionary of pfam file with old_locus as key
pfam_dict={}
for entry3 in pfam_table:
    if "#" in entry3[0]:
        continue
    else:
        locus = entry3[1].strip()
        pfam_dict.update({locus:entry3})
### pfam_dict now contains a dictionary of the pfam file with locus as the tag


### Plant Interaction Domain (PID) table format:
#   #Pfam ID	Pfam Name	RNB MEDIAN	NC MEDIAN	Transcript induction in nodule		LPD	% total seqs with SignalP	Antismash Pfam (Seconday metabolite)	PID	PCPfam	SignifEnrich
#read in PID summary table
print("PID Summary Table")
#pidInput = get_file()

pidInput = 'PIDOverRepresentedPfams.txt'

print ("PID Table: " + pidInput + '\n')
pidIn = open(pidInput, 'rU')
pid_table = read_in_tsv_file(pidIn)
pid_dict={}
for entry4 in pid_table:
    if "#" in entry4[0]:
        continue
    else:
        pfam = entry4[0].strip()
        pid_dict.update({pfam:entry4})
### pid_dict now contains a dictionary of the PID file with pfam as the key


###Summary: Summary of data files in memory
# Old Annotation Gff3 file: gff_old (key: ID)
# Dictionary of Gff3 file: gff_dict (key: ID)
# Dictionary of NCBI protein table file: prot_dict (key: locus_tag)
# Dictionary of IMG pfam annotation file: pfam_dict (key: old_locus_tag)
# Dictionary of proposed PID pfams: pid_dict (key: pfam)

### Use the dictionaries to update the gff3 file for VF39
gff_new=[]
# Loop over the unedited gff3 file
for entry5 in gff_old:
    if "#" in entry5[0]:
        gff_new.append(entry5)
    else:
        entry5_attributes = dict(item.split("=") for item in entry5[8].split(";"))  # dict(att1)=1, dict(att2)=2, dict(att3)=3 ...
        if entry5[2] == 'gene' and entry5_attributes["gene_biotype"] == 'rRNA':
            gff_new.append(entry5)
        elif entry5[2] == 'gene' and entry5_attributes["gene_biotype"] == 'tRNA':
            gff_new.append(entry5)
        elif entry5[2] == 'gene' and entry5_attributes["gene_biotype"] == 'RNase_P_RNA':
            gff_new.append(entry5)
        elif entry5[2] == "gene" and entry5_attributes["gene_biotype"] == "protein_coding":
            gff_entry = entry5  # seqname;source;feature;start;end;score;strand;frame;attributes
            # attributes for 'gene' entry
            ID = str(entry5_attributes["locus_tag"])
            Name = ""
            Gbkey = "Gene"
            Gene_biotype = "protein_coding"
            Old_locus_tag = ''

            if 'old_locus_tag' in entry5_attributes:
                Old_locus_tag = str(entry5_attributes["old_locus_tag"])
            else:
                Old_locus_tag = 'none annotated check NCBI for conflict'

            COG = ""
            EnzyComm = ""

            # Find entry information in protein and pfam dictionaries
            entry_prot_table = []
            entry_pfam_table = []

            if ID in prot_dict:  # Protien dictionary is keyed on locus_tag(Id=locus_tag now)
                entry_prot_table = prot_dict[ID]
            else:
                entry_prot_table = "-"

            if Old_locus_tag in pfam_dict:  # pfam dictionar is keyed on old RL locus tags (IMG's fault)
                entry_pfam_table = pfam_dict[Old_locus_tag]
            else:
                entry_pfam_table = "-"

            # Updating the 'Name' attribute
            if entry_prot_table == "-":  # Find the gene Name (gapD) for ID; if none, use ID.
                Name = ID
            elif entry_prot_table[6] == '-':
                Name = ID
            else:
                Name = entry_prot_table[6].strip()

            # Updating the new 'COG' attribute
            COG = entry_pfam_table
            COG_tmp = []
            if COG == "-":
                COG = "none annotated"
            elif COG[4] == "-":
                COG = "none annotated"
            elif "<<>>" in COG[4]:
                COG = COG[4]
                COG = list([s.strip() for s in COG.split("<<>>")])
                for i in COG:
                    item = i.strip().split("===")
                    item = "::".join(item)
                    COG_tmp.append(item)
                COG = "<<>>".join(COG_tmp)
            elif "===" in COG[4]:
                COG = COG[4]
                COG_tmp = list([s.strip() for s in COG.split("===")])
                COG_tmp = "::".join(COG_tmp)
                COG = COG_tmp
            else:
                COG = "irregular parse signal"

            # Updating the new 'EnzyComm' attribute
            EnzyComm = entry_pfam_table
            EnzyComm_tmp = []
            if EnzyComm == "-":
                EnzyComm = "none annotated"
            elif EnzyComm[6] == "-":
                EnzyComm = "none annotated"
            elif "<<>>" in EnzyComm[6]:
                EnzyComm = EnzyComm[6]
                EnzyComm = list([s.strip() for s in EnzyComm.split("<<>>")])
                for i in EnzyComm:
                    if "-===" in i:
                        item = i.strip().split("-===")
                        item = "::".join(item)
                        EnzyComm_tmp.append(item)
                    elif "===" in i:
                        item = i.strip().split("===")
                        item = "::".join(item)
                        EnzyComm_tmp.append(item)
                    else:
                        EnzyComm_tmp.append("irregular parse signal")
                EnzyComm = "<<>>".join(EnzyComm_tmp)
            elif "-===" in EnzyComm[6]:
                EnzyComm = EnzyComm[6]
                EnzyComm_tmp = list([s.strip() for s in EnzyComm.split("-===")])
                EnzyComm_tmp = "::".join(EnzyComm_tmp)
                EnzyComm = EnzyComm_tmp
            else:
                EnzyComm = EnzyComm[6]
                EnzyComm_tmp = list([s.strip() for s in EnzyComm.split("===")])
                EnzyComm_tmp = "::".join(EnzyComm_tmp)
                EnzyComm = EnzyComm_tmp

            #  All fields needed for 'gene' entry attributes have been updated
            #   ID = locus_tag
            #   Name = 'xzyA' else locus_tag
            #   Gbkey = 'gene'
            #   Gene_biotype = 'protein_coding'
            #   old_locus_tag = 'RL####'
            #   COG = 'COG#####::Some Description(<<>>COG#####::Some Description)'
            #   EnzyComm = 'EC:#.#.#.::Something(<<>>EC:#.#.#.::Something Else)'

            # Attributes for new 'mRNA' entry
            mRNA_ID = ID + ".t01"
            mRNA_Parent = ID
            mRNA_gbkey = "mRNA"

            # Attributes for new 'CDS' entry
            CDS_ID = ID + ".p01"
            CDS_Parent = mRNA_ID

            Dbxref = ""
            if entry_prot_table == "-":  # Updating the Dbxref attribute
                Dbxref = "none"
            elif entry_prot_table[8] == '-':
                Dbxref = "none"
            else:
                Dbxref = "Genbank:" + entry_prot_table[8]

            CDS_Name = Name  # 'Name' attribute defined previously
            CDS_gbkey = "CDS"

            Product = ""
            if entry_prot_table == "-":  # Updating the 'product' attribute with protein table product
                Product = "none annotated"
            elif entry_prot_table[10] == '-':
                Product = "none annotated"
            else:
                Product = entry_prot_table[10]

            SigRNBD = []
            PID = []
            Pfams = []
            # Updating the new 'pfam' attribute and the SigRNBD and PID attributes
            if entry_pfam_table == "-":
                Pfams = str("none annotated")
                SigRNBD = str("none annotated")
                PID = str("none annotated")
            elif entry_pfam_table[5] == '-':
                Pfams = str("none annotated")
                SigRNBD = str("none annotated")
                PID = str("none annotated")
            elif "<<>>" in entry_pfam_table[5]:
                pfams = entry_pfam_table[5]
                pfams = list([s.strip() for s in pfams.split("<<>>")])
                # looping over each 'pfam00000===something'
                for each in pfams:
                    item = each.strip().split("===")  # item = ['pfam00000', 'something']
                    if item[0] in pid_dict:
                        PID_item = pid_dict[item[0]]
                        SigRNBD.append(PID_item[0])
                        if PID_item[8] == "PID":
                            PID.append(PID_item[0])
                    # Make a sanity check outside the loop on the length of SigRNBD and PID; if len() == 0: PID = 'none...'
                    item = "::".join(item)
                    Pfams.append(item)
                Pfams = "<<>>".join(Pfams)
            else:
                pfam = entry_pfam_table[5]
                item = pfam.strip().split("===")  # item = ['pfam00000', 'something']
                if item[0] in pid_dict:
                    PID_item = pid_dict[str(item[0])]
                    SigRNBD = PID_item[0]
                    if PID_item[8] == "PID":
                        PID = PID_item[0]
                Pfams = "::".join(item)

            # Finalize the SigRNBD
            if len(SigRNBD) == 0 and type(SigRNBD) == list:
                SigRNBD = "none annotated"
            elif len(SigRNBD) >= 2 and type(SigRNBD) == list:
                SigRNBD = "<<>>".join(SigRNBD)
            else:
                pass

            # and PID attributes
            if len(PID) == 0 and type(PID) == list:
                PID = "none annotated"
            elif len(PID) >= 2 and type(PID) == list:
                PID = "<<>>".join(PID)
            else:
                pass
            # All items should be updated now, ready for merger and appending to the updated gff file
            #   gff_gene
            # Attributes
            #   ID
            #   Name
            #   Gbkey = 'gene'
            #   Gene_biotype = 'protein_coding'
            #   Old_locus_tag = 'RL####'
            #   COG = 'COG#####::Some Description(<<>>COG#####::Some Description)'
            #   EnzyComm = 'EC:#.#.#.::Something(<<>>EC:#.#.#.::Something Else)'
            gff_gene = entry5[0:8]
            gff_gene[7] = '.'
            gene_attributes = []
            ID = "ID=" + ID
            Name  = "Name=" + Name
            Gbkey = "gbkey=" + Gbkey
            Gene_biotype = "gene_biotype=" + Gene_biotype
            Old_locus_tag = "old_locus_tag=" + Old_locus_tag
            COG = "COG=" + COG
            EnzyComm = "EnzyComm=" + EnzyComm

            gene_attributes.append(ID)
            gene_attributes.append(Name)
            gene_attributes.append(Gbkey)
            gene_attributes.append(Gene_biotype)
            gene_attributes.append(Old_locus_tag)
            gene_attributes.append(COG)
            gene_attributes.append(EnzyComm)
            gene_attributes = ';'.join(gene_attributes)
            gff_gene.append(gene_attributes)

            gff_new.append(gff_gene)

            # Attributes for new 'CDS' entry
            #   gff_CDS = []  # list for updated information
            #   CDS_ID = ID + '.p01'
            #   CDS_Parent = mRNA_ID
            #   Dbxref = ''
            #   Product
            #   Pfams
            #   SigRNBD
            #   PID
            gff_CDS = entry5[0:8]
            gff_CDS[7] = '0'
            CDS_attributes = []

            CDS_ID = "ID=" + CDS_ID
            CDS_Parent = "Parent=" + CDS_Parent
            CDS_Name = "Name=" + CDS_Name
            Dbxref = "dbxref=" + Dbxref
            Product = "product=" + Product
            Pfams = "pfam=" + str(Pfams)
            SigRNBD = "SigRNBD=" + str(SigRNBD)
            PID = "PID=" + str(PID)

            CDS_attributes.append(CDS_ID)
            CDS_attributes.append(CDS_Parent)
            CDS_attributes.append(CDS_Name)
            CDS_attributes.append(Dbxref)
            CDS_attributes.append(Product)
            CDS_attributes.append(Pfams)
            CDS_attributes.append(SigRNBD)
            CDS_attributes.append(PID)

            CDS_attributes = ";".join(CDS_attributes)

            gff_CDS[2] = "CDS"
            gff_CDS.append(CDS_attributes)

            gff_new.append(gff_CDS)

            # Attributes for new 'mRNA' entry
            #   gff_mRNA = []  # list for updated information
            #   mRNA_ID = ID + '.t01'
            #   mRNA_Parent = ID
            #   mRNA_gbkey = 'mRNA'
            gff_mRNA = entry5[0:8]
            gff_mRNA[7] = '0'
            mRNA_attributes = []
            mRNA_ID = "ID=" + mRNA_ID
            mRNA_Parent = "Parent=" + mRNA_Parent
            mRNA_gbkey = "gbkey=" + mRNA_gbkey

            mRNA_attributes.append(mRNA_ID)
            mRNA_attributes.append(mRNA_Parent)
            mRNA_attributes.append(mRNA_gbkey)

            mRNA_attributes = ";".join(mRNA_attributes)

            gff_mRNA.append(mRNA_attributes)
            gff_mRNA[2] = "mRNA"
            gff_new.append(gff_mRNA)
        elif entry5[2] == "mRNA":  # Pass over old mRNA entries
            continue
        elif entry5[2] == "CDS":  # Pass over old CDS entries
            continue
        else:  # Push anything else to the new gff3
            gff_new.append(entry5)

# make an output file

wkdir = os.getcwd()
outfile = str(wkdir + "RLV3841.updated.gff3")
gffout = open(outfile, "w")
print('Printing to ' + str(gffout) + '\n')

# gff_new contains the updated gff file in a list[][]
for line3 in gff_new:
    print_line = str("\t".join(line3))
    gffout.write(print_line + '\n')

gffout.close()

exit()

