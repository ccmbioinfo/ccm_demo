import re
import os
import argparse

parser = argparse.ArgumentParser(description='process mmseqs .tab output into a3m')
parser.add_argument('-i', '--input', help="input .tab file", type=str, action="store")
parser.add_argument('-o', '--output', help="ouptut .a3m file", type=str, action="store")
args=parser.parse_args()

mmseq_tab=open(args.input)
mmseq_data=mmseq_tab.readlines()
mmseq_out=open(args.output, "w")
mmseq_count=0
for mmseq_line in mmseq_data:
  mmseq_line=mmseq_line.replace("\r", "").replace("\n", "").split("\t")
  mmseq_count+=1
  #Parse cigar
  mmseq_cigar=re.findall(r'(\d+)([MDI])?', mmseq_line[7])
  #Pad start of alignment
  alignment_seq="-"*(int(mmseq_line[2])-1)
  alignment_index=int(mmseq_line[4])-1;
  #Loop through cigar re-writing sequence
  for cigar_entry in mmseq_cigar:
    for sequence_index in range(0, int(cigar_entry[0])):
      if cigar_entry[1]=="M":
        alignment_seq+=mmseq_line[6][alignment_index:alignment_index+1].upper()
        alignment_index+=1
      elif cigar_entry[1]=="D":
        alignment_seq+=mmseq_line[6][alignment_index:alignment_index+1].lower()
        alignment_index+=1
      elif cigar_entry[1]=="I":
        alignment_seq+="-"
  #Pad end of alignment
  alignment_seq=alignment_seq+("-"*(int(mmseq_line[1])-int(mmseq_line[3])))
  #Print result
  mmseq_out.write(">%s\n" % mmseq_line[0])
  mmseq_out.write("%s\n" % alignment_seq)
mmseq_out.close()