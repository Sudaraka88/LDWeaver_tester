if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY

### GFF3 RUN with SNP_ONLY_ALIGNMENT
rm(list=ls())
dset = "ec_test_gff"
if(file.exists(dset)) unlink(dset, recursive = T)

gff_path = "ref.gff3"
ref_fasta_path = "ref.fa"

aln_path = "out.fa.gz"
pos = as.numeric(readLines("out.fa.pos"))

LDWeaver::LDWeaver(dset = dset,
                   aln_path = aln_path, aln_has_all_bases = F,
                   pos = pos, 
                   gff3_path = gff_path,
                   ref_fasta_path = ref_fasta_path,
                   perform_SR_analysis_only = T, SnpEff_Annotate = T, save_additional_outputs = F)

### GBK RUN with SNP_ONLY_ALIGNMENT
rm(list=ls())
dset = "ec_test_gbk"
if(file.exists(dset)) unlink(dset, recursive = T)

gbk_path = "ref.gbk"
aln_path = "out.fa.gz"
pos = as.numeric(readLines("out.fa.pos"))

LDWeaver::LDWeaver(dset = dset,
                   aln_path = aln_path, aln_has_all_bases = F,
                   pos = pos,
                   gbk_path = gbk_path,
                   perform_SR_analysis_only = T, SnpEff_Annotate = F, save_additional_outputs = F)

### CHECK FOR FILE SIMILARITY
s1 = read.csv("ec_test_gbk/Temp/sr_links.tsv", sep = "\t", header = F)
s2 = read.csv("ec_test_gff/Temp/sr_links.tsv", sep = "\t", header = F)
all.equal(s1, s2)

# s2 = read.csv("ec_test_gff/Temp/sr_links.tsv", sep = "\t", header = F)
sd1 = readRDS("ec_test_gbk/Additional_Outputs/snp_ACGTN.rds")
sd2 = readRDS("ec_test_gff/Additional_Outputs/snp_ACGTN.rds")
all.equal(sd1, sd2)

cd1 = readRDS("ec_test_gbk/Additional_Outputs/cds_var.rds")
cd2 = readRDS("ec_test_gff/Additional_Outputs/cds_var.rds")
all.equal(cd1, cd2)

h1 = readRDS("ec_test_gbk/Additional_Outputs/hdw.rds")
h2 = readRDS("ec_test_gff/Additional_Outputs/hdw.rds")
all.equal(cd1, cd2)

p1 = readRDS("ec_test_gbk/Additional_Outputs/parsed_gbk.rds")
p2 = readRDS("ec_test_gff/Additional_Outputs/parsed_gff3.rds")

t1 = LDWeaver::read_TopHits("ec_test_gbk/Tophits/sr_tophits.tsv")
t2 = LDWeaver::read_TopHits("ec_test_gff/Tophits/sr_tophits.tsv")
all.equal(t1, t2)

sra1 = LDWeaver::read_AnnotatedLinks("ec_test_gbk/Annotated_links/sr_links_annotated.tsv")
sra2 = LDWeaver::read_AnnotatedLinks("ec_test_gff/Annotated_links/sr_links_annotated.tsv")
all.equal(sra1, sra2)
