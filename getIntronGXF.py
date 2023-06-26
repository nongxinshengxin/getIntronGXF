import argparse
import sys
import re


def check_format(gxf,format):
    with open(gxf) as fr:
        if "exon" not in fr.read():
            print("Your GFF/GTF file lost the type of exon, please check !")
            sys.exit()
        if format=="gff" and "Parent" not in fr.read():
            print("Your GFF file format is wrong, please check !")
            sys.exit()
        if format=="gtf" and "transcript_id" not in fr.read():
            print("Your GTF file format is wrong, please check !")
            sys.exit()

def read_gxf(gxf):
    with open(gxf) as fr:
        for line in fr:
            yield line
                

def find_intron(gxf,format,newgxf):
    exon_dic={}
    info_dic={}
    for line in read_gxf(gxf):
        if not line.startswith("#"):
            lin=line.strip().split("\t")
            type=lin[2]
            start=int(lin[3])
            end=int(lin[4])
            attribution=lin[8]
            if type=="exon" and format=="gff":
                attri_list=attribution.split(";")
                for i in attri_list:
                    if re.search(r"Parent=",i):
                        trans_id=re.search(r"(?<=Parent=).+",i).group()
                        exon_dic.setdefault(trans_id,[]).append(start-1)
                        exon_dic.setdefault(trans_id,[]).append(end+1)
                        info_dic.setdefault(trans_id,[]).append(lin)
            if type=="exon" and format=="gtf":
                attri_list=attribution.split(";")
                for i in attri_list:
                    if re.search(r"transcript_id",i):
                        trans_id=re.search(r"(?<=transcript_id ).+",i).group()
                        exon_dic.setdefault(trans_id,[]).append(start-1)
                        exon_dic.setdefault(trans_id,[]).append(end+1)
                        info_dic.setdefault(trans_id,[]).append(lin)
    with open(newgxf,"w") as w:
    
        for key in exon_dic.keys():
            if len(exon_dic[key])>2:
                start_list=sorted(exon_dic[key])[1:-1][::2]
                end_list=sorted(exon_dic[key])[1:-1][1::2]
                for i in range(len(start_list)):
                    intron_start=start_list[i]
                    intron_end=end_list[i]
                    newlin=info_dic[key][0]
                    newlin[2]="intron"
                    newlin[3]=str(intron_start)
                    newlin[4]=str(intron_end)
                    newlin="\t".join(newlin)
                    w.write(newlin+"\n")



def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('gxfFile',type=str,help='input GFF/GTF file')
    parser.add_argument('introGXFfile',type=str,help='output GFF/GTF file')
    parser.add_argument('--format',type=str,default="gff", help='the format of inputfile. format takes the string "gtf" or " gff" (default).')

    args=parser.parse_args()

    check_format(args.gxfFile,args.introGXFfile)
    find_intron(args.gxfFile,args.format,args.introGXFfile)

if __name__== '__main__':
    main()

