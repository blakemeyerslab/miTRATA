#!/usr/local/bin/python3
# coding=utf-8

#Created by Parth Patel, DBI @ University of Delaware, Newark, Delaware 19717
#Date created: 10/08/2014 Modified on 11/2/2015 Line number 350,297


import os
import sys
if sys.version < '2.6':
    print ('You are using a version of Python that this program does not support. Please update to the latest version!')
    sys.exit(1)
import subprocess, multiprocessing
from multiprocessing import Process, Queue, Pool
import logging
import re #pattern search https://docs.python.org/2/howto/regex.html
#import Bio #importing biopython
#from Bio import SeqIO #import Bio-python
#plot specific packages
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from matplotlib import rc
#from pylab import *
import datetime,time,timeit #import Date and Time
import argparse
from PyPDF2 import PdfFileReader, PdfFileMerger #working with PDF in python
import random

###################################################################################################################
# STEP-1                                                                                                          #
#Given a tag, finds the longest subsequence in it which aligns to a reference genome (corresponding BOWTIE INDEX) #
#Given a file containing short sequences, finds the longest aligned subsequence from the 5' end                   #
###################################################################################################################

nproc = 'Y'
nthread = 6
global DIR_NAME
Date=str(datetime.date.today())#get Today's date
Time=str(datetime.datetime.now().time())#get Current time  )#get Current time
DIR_NAME="5pLCS_results"+"_"+Date+"_"+Time

global OUTPUT_DIR
#OUTPUT_DIR_NAME="Tailing_output_"+Date+"_"+Time

def RunBowtie(InputFile, BOWTIE_INDEX, Aligned, Unaligned, AlignmentOutput, pas):
    # command = "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/glibc-2.14/lib; bowtie  -v 0 "
    command = "bowtie  -v 0 "			+\
              "--threads %s " % nthread         +\
                        "%s " % BOWTIE_INDEX    +\
                     "-r %s " % InputFile       +\
                    "%s " % AlignmentOutput     +\
               "--al %s " % Aligned             +\
               "--un %s " % Unaligned         
    
    if os.system(command) == 0:
        print ("Bowtie Completed Successfully - Pass %d " % pas)
       
def pLCS_finder(inputfile):

    print ("Working on "+inputfile+"\n")
    InputFile=inputfile
    Date=str(datetime.date.today())#get Today's date
    OutputFile=os.path.splitext(inputfile)[0]+"_"+Date+".txt" #get the input filename without the extension and append Today's date
    bowtie_index=bowtieindex
#    print (bowtie_index)
   
    #open output file for writing
    IN= open(InputFile,'r')
    # Input validation
    if not os.path.isfile(InputFile):
        print ("Cannot find specified input file")
        sys.exit(1)
 
    if not os.path.isfile(bowtie_index+'.1.ebwt'):
        print ("Cannot find specified bowtie index")
        sys.exit(1)
     
     
    abundances={}
    tail_length={}
 
    for index, line in enumerate(IN.readlines()):
        fields = line.strip().split('\t')
        if index == 0:
            # if the line is the first line make sure its a header line before we skip (bug 1482)
            if "".join(sorted(set(fields[0]))) == 'ACGT' and fields[1].isnumeric():
                # first line is also a data line
                seq, value = fields[0], int(fields[1])
                abundances[seq] = value
        else:
            seq, value = fields[0], int(fields[1])
            abundances[seq] = value
     
    LogFile = os.path.splitext(OutputFile)[0] + '.size_distribution'
    tempInput = 'tempInput' +InputFile+ '.txt'
    AlignmentOutput = 'Bowtie' +InputFile+ '.out'
    aligned = 'aligned_tags' +InputFile+ '.out'
    unaligned = 'unaligned_tags' +InputFile+ '.out'
    pas = 0
 
    logging.basicConfig(filename=LogFile, level=logging.NOTSET, format="%(message)s")
 
    next_run = open(tempInput, 'w')
    next_run.write("\n".join(tag for tag in abundances.keys())+'\n')
    next_run.close()
 
    print ("Abundance is %d long and tail length is %d long" %(len(abundances), len(tail_length)))
#     
#     ###############################################################################################
    while len(abundances) != len(tail_length):
        
#        print ("Abundances")
       RunBowtie(tempInput, bowtie_index, aligned,unaligned, AlignmentOutput, pas)            
#           count = 0
       list_of_aligned = {}
       try:
           for tag in open(aligned):
               list_of_aligned[tag[:-1]] = 1 # Create an entry for each tag that is aligned. The value bears no significance 
       except IOError:
           pass
 
#        print ("\n\n\n\n\n\n\n\n\n\nDone Listing")
#        exit(0)

 
       for tag in abundances.keys():
            if tag not in tail_length:
                check_tag = tag                
                if pas < 0:
                    check_tag = tag[:pas]    # For the special case of pas = 0, the intended slicing to trim would not work. For every other pass, this statement trims.
                if check_tag in list_of_aligned:
                    tail_length[tag] = pas
 
#        print ("\n\n\n\n\n\n\n\n\nfound %d tails until this pass" % len(tail_length))
#     
       pas = int(pas)-1
       next_run = open(tempInput, 'w')
       unal = open(unaligned)
       next_run.write("\n".join(tag[:-2] for tag in unal.readlines())+'\n')
       unal.close()
       next_run.close()
#        print("\n\n\n\n\n\n\n\n\nClosed")
    #############################################################################################
    
 
    command = "rm %s %s %s %s" %(tempInput, AlignmentOutput, aligned, unaligned)
    if os.system(command) == 0:
        print ("Deleting temporary Files.........")
 
    taglength_count = {} 
    headlength_count = {}
    taglength_sum = {}
    headlength_sum = {}    
    root_of_head_length_sum = {}
     
#    Date=str(datetime.date.today())#get Today's date
#    Time=str(datetime.datetime.now().time())#get Current time  )#get Current time
#    dirName="5pLCS_results"+"_"+Date+"_"+Time
    dirName=DIR_NAME
    if not os.path.exists(dirName): #create a directory for storing the results of STEP1
        os.makedirs(dirName)
    current_path= os.getcwd() #get the path of current directory
    new_path= os.path.join(current_path,dirName) #add newly created results' folder to an exisiting path 
    os.chdir(new_path) #change the directory    
     
    Output=open(OutputFile, 'w')
    Output.write("tag\tabundance\thead\thead_abundance\ttail\n")
    for tag in abundances.keys():
        Output.write("%s\t" % tag)
        Output.write("%d\t" % abundances[tag])
        if tail_length[tag] != 0:        
            Output.write("%s\t" % tag[:tail_length[tag]] )        
            head_abundance = abundances.get(tag[:tail_length[tag]], 0 )        
            if head_abundance == 0:
                Output.write("0\t")
            else:
                Output.write("%d\t" % head_abundance)
            Output.write("%s\n" % tag[len(tag)+tail_length[tag]:] )        
        else:
            Output.write("%s\t" % tag)
            Output.write("%d\t" % abundances[tag])    
            #Output.write("None\n")
            Output.write("\n")
         
         
     
            # The three try/except clauses have to be maintained separately and as they are:
 
        try:
            taglength_count[ len(tag) ] = taglength_count[len(tag) ] + 1
            taglength_sum[ len(tag) ] = taglength_sum[ len(tag) ] + abundances[tag]
        except KeyError:
            taglength_count[ len(tag)] = 1
            taglength_sum[ len(tag) ] = abundances[tag]
             
        try:
            headlength_count[ len(tag) + tail_length[tag] ] = headlength_count[ len(tag) + tail_length[tag] ] + 1
        except:
            headlength_count[ len(tag) + tail_length[tag] ] = 1
 
        try:
            headlength_sum[ len(tag) + tail_length[tag] ] = headlength_sum[ len(tag) + tail_length[tag] ] + abundances.get(tag[:tail_length[tag]], 0 )    
        except:            
            headlength_sum[ len(tag) + tail_length[tag] ] = abundances.get(tag[:tail_length[tag]], 0 )
 
        try:
            root_of_head_length_sum[ len(tag) ] = root_of_head_length_sum[ len(tag) ] + abundances.get(tag[:tail_length[tag]], 0 )
        except KeyError:
            root_of_head_length_sum[ len(tag) ] = abundances.get(tag[:tail_length[tag]], 0 )
 
    Output.close()
     
    logging.info( "(Size-Group)\t(Number of tags)\t(Sum of tags)\t(Number of heads)\t(Sum of heads)\t(Sum of heads with tag in size)\n" )
    for key in range( min( headlength_count.keys() ), max( taglength_count.keys() )+1 ):
        logging.info( "%d\t%d\t%d\t%d\t%d\t%d" %( key, taglength_count.get(key,0), taglength_sum.get(key,0), headlength_count.get(key,0), headlength_sum.get(key,0), root_of_head_length_sum.get(key,0) ) )
 
    logging.shutdown()
    os.chdir(current_path)#go back to main directory
    print("\n\npLCS_finder Done")
    return dirName


####################################################################################
# STEP-2 For multiple libraries from the 1st step, merge 5GMC/Tail from above into one   #
####################################################################################

# USAGE = " python2.6 Merge_5pLCS.py  -[Name of the results Directory from first step]  -[output_file_name] "

def merge_5plCS(DIR,outfile):

    seqdir = DIR #"results"#str(sys.argv[1]) #seqdir directory_with_chrom_seq
    outputfile = outfile#"step2_output_HEN1.txt"#str(sys.argv[2]) #output_file_name
    
    #open output file for writing
    out= open(outputfile,'w')
    
    #library id starting from 1
    lib_id=1
    
    
    for filename in os.listdir(seqdir): #list all the files in the current directory
    #    print os.listdir(seqdir)
    #for filename in os.listdir(os.getcwd()): #list all the files in the current directory
        lib_name= os.path.splitext(filename)[0] #get the filename without the extension
    #    print lib_name + "\n" + filename
        in_file = open(os.path.join(seqdir,filename),'r') # open file for reading     
        
        
        f = in_file.readlines()
        firstLine = f.pop(0) #removes the first line or the header row
        for line in f:
            
            out.write(str(lib_id))
            out.write("\t")         
            out.write(lib_name)
            out.write("\t") 
            out.write(line)
            out.write("\n")  
            
        in_file.close() #closing the inputfile to fetch the next one
        lib_id=lib_id+1 #update the library id
#        format_Output(input_file)
        
       
    print ("Second Step is DONE..Files Merged.")
    out.close() #closing the output file 
    return outputfile

####################################################################################
# STEP-3 For multiple libraries from the 1st step, merge 5GMC/Tail from above into one   #
####################################################################################

# USAGE = " python2.6 Tailing_analysis_DB_v6.py  -[Name of the miRNA/siRNA file] -[Merged file from the second step] "

def tailling_Analysis(miRNA_file,mergred_file):
    miRNA_file_name = miRNA_file#"new_miRNA_Test1.fa"#str(sys.argv[1]) #seqdir directory_with_chrom_seq
    merged_file_name = mergred_file#"step2_output_HEN1.txt"#str(sys.argv[2]) #output_file_name
    
    miRNA_NAME=[]
    miRNA_TAG=[]

    fh_in = open(miRNA_file_name,'r')
    miRNAData = fh_in.read().split('>')
    miRList = [] ## Store miRNA NAME and TAG as tuple
    for i in miRNAData[1:]:
        block = i.split('\n')
        ID = block[0].split(' ')[0]
        seq = block[1]
        miRList.append((ID,seq))
#        print (ID,seq)
    
    miRListS = sorted(miRList) #Sorted by miRNA NAME
    for i in miRListS:    
        ID,SEQ=i
        miRNA_NAME.append(ID)
        miRNA_TAG.append(SEQ)
   
       ##convert all U'S into T's for miRNA_TAG  
    for index in range(0,len(miRNA_TAG)):
       #    print miRNA_TAG[index]+"->>"+ re.sub("U","T",str(miRNA_TAG[index])) 
       miRNA_TAG[index]=re.sub("U","T",str(miRNA_TAG[index]))
       miRNA_TAG[index]=miRNA_TAG[index].upper() # create an upper case if miRNAs are not already in uppercase, Added on 11/2/2015
#       print (miRNA_TAG[index],miRNA_NAME[index])

    print ("miRNA laoding DONE!") 
        
    
    
#    output_dirName=OUTPUT_DIR_NAME
#    if not os.path.exists(output_dirName): #create a directory for storing the results of STEP1
#        os.makedirs(output_dirName)
#    current_path= os.getcwd() #get the path of current directory   
    
    merged_input=open(merged_file_name, 'r') # READ merged file from step 2
    
#    new_path= os.path.join(current_path,output_dirName) #add newly created results' folder to an exisiting path 
#    os.chdir(new_path) #change the directory    
    os.chdir(OUTPUT_DIR)
    out= open("Tail_truncated_"+miRNA_file_name+merged_file_name,'w')    
    out.write("Lib_ID \t Lib_NAME \t Tag \t Tag_abundance  \t Sub_tag \t Sub_tag_hits \t Sub_tag_abun \t Tail \t Tail length\n") # tail length added by Parth after the 1st review of this paper. 8/5/2015  
    out2=open("Tail_summary_truncated_"+miRNA_file_name+merged_file_name,'w')
    #out2.write("miRNA_NAME \t miRNA_TAG \t miR_abundance \t Lib_ID \t Lib_NAME \t sum_tailed \t percent \n") didn't use percent because it's giving value more than 100% e.g. 466.66
    out2.write("miRNA \t miRNA sequence \t Abundance \t library # \t library \t Sum of abundances \t tailing ratio = Sum of abundances/Abundance \n")
    
    hash_miR= {}
    miR_abun= {}
    hash_count= {}
    hash_lib_id= {}
    hash_miR_rows=[] #Created this list which acts as a Arrays of Hashes
    
    
    for lines in merged_input: #reading each line from the file
    
        
        if not lines.strip(): # To remove empty line in your file (usually the last line)
            continue
        #print (lines)
        lib_id, lib_name, tag ,tag_abun, sub_tag, sub_tag_abun,tail = lines.split('\t',21)
    #   1 test1 CCCAGGTCCAGACATAGTAAGGATTGACAGACTGAGATCACTTTCTTGATTC 1 CCCAGGTCCAGAC 0 ATAGTAAGGATTGACAGACTGAGATCACTTTCTTGATTC
       
    #    print lib_id, lib_name, tag ,tag_abun, sub_tag, sub_tag_abun,tail     
       
        hash_lib_id [lib_id]=lib_name    
        
        if (len(tag)>30):
            continue
        if (len(sub_tag)<10):
            continue      
        
    
        #search for truncated miRNA
        for index in range(0,len(miRNA_NAME)):
            ID=lib_id+"_"+miRNA_NAME[index]            
       
           
            if((re.match(str(miRNA_TAG[index]),sub_tag,re.IGNORECASE)) or (re.match(sub_tag,str(miRNA_TAG[index]),re.IGNORECASE))): # the pattern at the start of the string, returning a match object, or None if no match was found. re.IGNORECASE Added on 11/2/2015 to handle lower case tag count file.
    
            
                if (ID in hash_miR): #Append matched pattern to an exisiting miRNA entry in the hash(dictionary)
                    
                    hash_miR[ID].append([lib_id, lib_name, tag ,tag_abun, sub_tag, sub_tag_abun,tail])
                   
                else: #if miRNA entry doesn't exist, create an empty entry in the hash and append the pattern
                    
                    hash_miR[ID] = []                      
                    hash_miR[ID].append([lib_id, lib_name, tag ,tag_abun, sub_tag, sub_tag_abun,tail])            
              
    
                hash_count[sub_tag]=+1           
                
    
            if(miRNA_TAG[index]==tag):
                miR_abun[ID]=tag_abun
    
      
    print ("Tail matching Done!\n")
    
    
    for index in range(0,len(miRNA_NAME)):
        
        sortedKeys=sorted(hash_lib_id.keys())# sort the hash_lib_id on thier keys
    
        
        for key in sortedKeys:
            lib_id=key
            ID=lib_id+"_"+miRNA_NAME[index]
    #        print sorted(miR_abun.keys())
            if(ID in miR_abun):
                out.write (">"+miRNA_NAME[index]+"_"+miRNA_TAG[index]+"_"+lib_id+"_"+ hash_lib_id[lib_id]+"_"+ miR_abun[ID]+"\n")
            else:
                out.write (">"+miRNA_NAME[index]+"_"+miRNA_TAG[index]+"_"+lib_id+"_"+ hash_lib_id[lib_id]+"_"+ "_" +"\n")
    #            continue
            
            sum_tailed=0
            tail_ratio=0
            
            if(ID in hash_miR): #Check if we have the miRNA entry in the hash or not        
                VALUES = hash_miR[ID]
    #            print VALUES
            
    
                
                
                for hash_miR_rows in VALUES:
    #                print hash_miR_rows   
                    lib_id, lib_name, tag ,tag_abun, sub_tag,sub_tag_abun,tail=hash_miR_rows  
                
                    sub_tag_hits=hash_count[sub_tag]   # compute sub_tag_hits
                    tail=tail.strip()# remove "\n" char from string
                    hash_miR_rows= lib_id, lib_name, tag ,tag_abun, sub_tag,sub_tag_hits,sub_tag_abun,tail,len(tail) # tail length added by Parth after the 1st review of this paper. 8/5/2015  
                
                    #print (hash_miR_rows)
    
                    delimiter='\t';
                    output = delimiter.join(str(x) for x in hash_miR_rows)
                    out.write("%s \n" % output) # tail length added by Parth after the 1st review of this paper. 8/5/2015           
    
                    sum_tailed=sum_tailed+int(tag_abun)
         
                if(ID in miR_abun):
                    if (int(miR_abun[ID])<1):
#                        percent = "N\A"
                        tail_ratio="N\A"
                    else:
                        tail_ratio= (float(sum_tailed)/int(miR_abun[ID]))#Calculate the ratio- CORRECT WAY
#                        tail_ratio= 100*(float(sum_tailed)/int(miR_abun[ID]))#Calculate the ratio -WRONG WAY
#                    percent=round(tail_ratio,2) #Round number to 2 digits after decimal point
                    tail_ratio=round(tail_ratio,2) #Round number to 2 digits after decimal point
#                    out2.write("%s \t %s \t %d \t %s \t %s \t %d \t %f \n" % (miRNA_NAME[index] ,miRNA_TAG[index],int(miR_abun[ID]),lib_id,hash_lib_id[lib_id],sum_tailed,percent))
                    out2.write("%s \t %s \t %d \t %s \t %s \t %d \t %f \n" % (miRNA_NAME[index] ,miRNA_TAG[index],int(miR_abun[ID]),lib_id,hash_lib_id[lib_id],sum_tailed,tail_ratio))
                else:
                    out2.write("%s \t %s \t %s \t %s \t %s \t %d \t %s \n" % (miRNA_NAME[index] ,miRNA_TAG[index],"0",lib_id,hash_lib_id[lib_id],sum_tailed,"N\A"))
    
    # closing out put file   
    out.close()
    out2.close()   
    
    print ("Step 3 is DONE!\n")
    return ("Tail_truncated_"+miRNA_file_name+merged_file_name)

####################################################################################
# STEP- 4 Format output file from the 3rd step                                            #
####################################################################################

# USAGE = " python2.6 Format_Taling_v2.py  -[Name of the outputfile from step3] "
#>ath-miR156a_TGACAGAAGAGAGTGAGCAC_10_hen1_8_660
#10 hen1_8  TGACAGAAGAGAGTGACCACA   1   TGACAGAAGAGAGTGA    216 0   CCACA
#10    hen1_8    TGACAGAAGAGAGTGAGCACTTTT    16    TGACAGAAGAGAGTGAGCAC    726    660    TTTT

def format_Output(input_file):

    inputfile = input_file #"Tail_truncated_new_miRNA_Test1.fastep2_output_HEN1.txt"#str(sys.argv[1]) #input_file_name
    outputfile="FORTMATTED_STEP4_OUTPUT.txt"
    out= open(outputfile,'w')
    
    #open output file for writing
    IN = open(inputfile,'r')
    
    #defining hash variable for miRNA
    hash_miRNA={}
    hash_lib={}
    hash_abun={}
    
    miRNA_name={}
    
    lib_id=[]
    
    
    LINES = IN.readlines()
    firstLine = LINES.pop(0) #removes the first line or the header row 
    print (firstLine)
    for lines in LINES: #reading each line from the file 
    
        if not lines.strip():# To remove empty line in your file (usually the last line)
            continue  
       
        if(re.search('\>',lines)):                 
            
            lines=re.sub('\>',"",lines) # remove '>' from the line
            splitLine = re.split('_',lines,maxsplit=3) #split first 3 item with "_" delimiter        
            name= splitLine[0] #storing miRNA name
            miRNA_seq=splitLine[1] #storing miRNA sequence
            lib_id = splitLine[2] #storing libary id   
    #        print name,miRNA_seq,lib_id
            hash_miRNA[name]=miRNA_seq
            continue
        
    #    lines="10  hen1_8    TGACAGAAGAGAGTGAGCACTTTT    16    TGACAGAAGAGAGTGAGCAC    726    660    TTTT"
    #    lines="1    hen1-1-rep1    TCGGACCAGGCTTCACTTTTTT    4    TCGGACCAGGCTTCA    1    46    CTTTTTT"
     
        lib_id, lib_name,sRNA_seq,sRNA_abun,sub_tag,sub_hit,sub_abun,tail,tail_len = lines.split('\t',9)  # tail length added by Parth after the 1st review of this paper. 8/5/2015        
    #  OR  lib_id,lib_name,sRNA_seq,sRNA_abun,sub_tag,sub_hit,sub_abun,tail = re.split('\t',lines,maxsplit=8)
    #    print lib_id, lib_name,sRNA_seq,sRNA_abun,sub_tag,sub_hit,sub_abun,tail
        
    #    miRNA_name[name]=sRNA_seq,sub_tag,miRNA_seq,tail
        if (name in miRNA_name): 
            miRNA_name[name].append([sRNA_seq,sub_tag,miRNA_seq,tail,tail_len])#Append matched pattern to an exisiting miRNA entry in the hash(dictionary) # tail length added by Parth after the 1st review of this paper. 8/5/2015  
        else:
            miRNA_name[name]=[]
            miRNA_name[name].append([sRNA_seq,sub_tag,miRNA_seq,tail,tail_len])#if miRNA entry doesn't exist, create an empty entry in the hash and append the pattern   # tail length added by Parth after the 1st review of this paper. 8/5/2015  
    #    print miRNA_name[name]
    #    print lib_name,sRNA_seq,miRNA_name
        hash_lib[lib_name]=lib_name
    #    print hash_lib [lib_name]
    #    print hash_lib[lib_name]   
        hash_abun[lib_name,sRNA_seq]=sRNA_abun
    #    print hash_abun[lib_name,sRNA_seq]
        
        
     
    sortedKeys=sorted(hash_miRNA.keys()) # sort the hash_miRNA on thier keys
    for key in sortedKeys: 
        
        hash_sRNA={}
        
    #    print ">"+key+"_"+hash_miRNA[key]+"\n"
        out.write(">"+key+"_"+hash_miRNA[key]+"\n")
    #    print "Complete_Sequence \t Sub_tag \t miRNA_seq \t Tail"
        out.write("Complete_Sequence \t"+" Sub_tag \t"+" miRNA_seq \t"+"Tail \t"+"Tail length")
        sortedKeys_1=sorted(hash_lib.keys())# sort the hash_lib on thier keys
        for lib in sortedKeys_1:
            out.write("\t"+lib)
        out.write("\n")
        
     
        
            
        if(key in miRNA_name):
                
            VALUES=miRNA_name[key]
            for miRNA_rows in VALUES:
                sRNA_seq,sub_tag,miRNA_seq,tail,tail_len=miRNA_rows # tail length added by Parth after the 1st review of this paper. 8/5/2015  
            
                delimiter='\t';
                output = delimiter.join(str(x) for x in miRNA_rows)
    #            print output
                count=0
                for lib_entry in sortedKeys_1:
    #                count=+1
    #                print "count ="+str(count)+"\n"
                    if((lib_entry,sRNA_seq) in hash_abun):
                        if(int(hash_abun[lib_entry,sRNA_seq])<1):
                            hash_abun[lib_entry,sRNA_seq]=0
                        output=output.replace("\n","")
                        output=output+"\t"+hash_abun[lib_entry,sRNA_seq] #Appending sRNA_seq aundance based on lib_ID
                        print (output)
                    else:
                        output=output.replace("\n","")
                        output=output+"\t"+"0" #Appending sRNA_seq aundance = 0  based on IF[lib_entry,sRNA_seq] DOES NOT EXIST in hash_sRNA
                                     
                hash_sRNA[sRNA_seq]=output      
            
        
        for sRNA_keys in (hash_sRNA.keys()):
            out.write(hash_sRNA[sRNA_keys]+"\n")
    #    output=""
    
    out.write(">EOF_NNNNNNNNNN\n")
    #closing files
    IN.close()
    out.close()
    print  ("Step-4 is DONE! \n")
    return (outputfile)

####################################################################################
# STEP-5  Generate plots from the 4th step                                            #
####################################################################################

def genrate_Plots(inputfile):

    ######################### User input #############################
    #input file
    input_file = inputfile
    #only U tail? (1=yes, 0= no)
    OnlyU = 1
    #exclude canonical miRNA? (1=yes, 0= no)
    exclude_miR = 0
    # Define number of Columns and Rows for each figure
    Num_Column = 4
    Num_Row = 7
    Total_lib = 0
    # Calculation Range, default 10
    Cal_range = 10
    # plotting range for truncation and tailing (e.g. miR163 need larger range because it's 24nt)
    Plot_range = 9
    #use n to initiate drawfigure
    n = 1
    # Amplification Factor
    AF = 1000
    # figure size
    Figure_size = 5
    
    # Draw figure
    def drawfigure(miRname, miRsize, lib_name, miR, miRsum):
    
        fig = plt.figure(figsize=(Figure_size * Num_Column, Figure_size * Num_Row))
#        fig.suptitle(r'%s, size %snt'%(miRname,miRsize),color='k',fontsize=18, ha='center')
        fig.suptitle(r'%s, size %s-nt'%(miRname,miRsize),x=0.30,y=0.925,color='k',fontsize=18, ha='right') #Added x=0.35 and y=0.925 by Parth to center the suptitle- Link http://stackoverflow.com/questions/8248467/matplotlib-tight-layout-doesnt-take-into-account-figure-suptitle
        fig.subplots_adjust(left=0.125, right=1, bottom=0, top=0.9, wspace=0.4, hspace=.5)# add white space between subplots
        
        def set_properties(ax):
            for i in range(Plot_range + 1):
                ax.plot((Plot_range,-1+i),(Plot_range-i,-1),'#AFC7C7',linewidth=1,zorder=2)#, linestyle= 'dashdot'
                ax.plot((Plot_range-i,-1),(Plot_range,-1+i),'#AFC7C7',linewidth=1,zorder=2)#, linestyle= 'dashdot'
            ax.set_yscale('linear')
            ax.set_xscale('linear')
            plt.xlabel('Length of Truncation', fontsize=16)  #added by Parth after the 1st review of this paper. 8/5/2015 
            plt.ylabel('Length of Tailing', fontsize=16)  #added by Parth after the 1st review of this paper. 8/5/2015 
            ax.set_xlim(-1, Plot_range) 
            ax.set_ylim(-1, Plot_range)
            ax.set_xlim(ax.get_xlim()[::-1])
            ax.set_xticks(range(-1,Plot_range))
            ax.set_yticks(range(-1,Plot_range))
            ax.grid(True)
    # calculate the propotion of each position in matrix
        for i in range (Cal_range):
            for j in range (Cal_range):
                for lib in range (Num_libs):
                    if miRsum[lib] >0:
                        miR[lib][i][j] = float(miR[lib][i][j])/miRsum[lib]*AF
                    else:
                        miR[lib][i][j] = 0
    # plotting with miR name and library name
        ax = [ 0 for x in range(Num_libs)]
        for lib in range (Num_libs):
            random_color="#"+("%06x"%random.randint(0,16777215)) # Added by Parth to plot each library with different random colors
#            random_color = "#%06x" % random.randint(0,0xFFFFFF)
            ax[lib] = fig.add_subplot(Num_Row,Num_Column,lib+1)
            ax[lib].set_title(r'%s (%d reads)' %(lib_name[lib], miRsum[lib]),color='k',fontsize=14, ha='center')
            for i in range (Cal_range):
                for j in range (Cal_range):
#                    ax[lib].scatter(i, j, c='red', s=miR[lib][i][j],linewidth=1,zorder=7)
                    ax[lib].scatter(i, j, c=random_color, s=miR[lib][i][j],linewidth=1,zorder=7)
            set_properties(ax[lib])
    
    #    plt.show()
        plt.autoscale()
        plt.savefig('%s-%s.png' % (miRname,miRsize),bbox_inches='tight') #bbox_inches='tight' is used to reduce left and right margins in matplotlib plot
        plt.savefig('%s-%s.pdf' % (miRname,miRsize),bbox_inches='tight')
        plt.close()
    
    
    
    
     
    #save multiple pdf into single pdf.
    def mergePDF():
        
        files_dir = os.getcwd() #get the path of current directory
        pdf_files = [f for f in os.listdir(files_dir) if f.endswith("pdf")]
        merger = PdfFileMerger()
        for filename in pdf_files:
            merger.append(PdfFileReader(os.path.join(files_dir, filename), "rb"))
        merger.write(os.path.join(files_dir, "merged_full.pdf"))
        
        
    ###################### Main Script ################################
    
    # Open file
    f = open (input_file, 'r')
    for line in f:
        if line.startswith( '>' ):
            if n >1:
                print (miRsum)
                drawfigure (miRname, miRsize, lib_name, miR, miRsum)
            head = [(x) for x in line.split('_')]
            miRname = head[0][1:]
            miRseq = head[1].rstrip('\n')
            miRsize = len(miRseq)
            print (miRname, miRseq, miRsize,"\n")
            n = n + 1
            continue
    # Find number of libraries. The related data is in the second row for each miRNA section, start with "Complete"
        elif line.startswith( 'Complete' ):
            data = [(x) for x in line.split('\t')]
            #number of libraries
            Num_libs = len(data) -5
            lib_name = data[5:]
    # Define a 10x10 two-dimention array for each miRNA 
            miR = [[[0 for x in range(Cal_range)] for x in range(Cal_range)] for x in range(Num_libs)]
            miRsum = [ 0 for x in range(Num_libs)]
            for lib in range (Num_libs):
                miR[lib] = [[0 for col in range(Cal_range)] for row in range(Cal_range)]
                miRsum[lib] = 0
            continue
    # Recording truncation and tailing data
        else:
            data = [(x) for x in line.split('\t')]
            Seq = data[0]
            Com = data[1]
            miRNA = data[2]
            Tail_seq = data[3]
    # decide if only include U tail
            if OnlyU == 1:
                if Tail_seq.find( 'A' ) > -1:
                    continue
                if Tail_seq.find( 'G' ) > -1:
                    continue
                if Tail_seq.find( 'C' ) > -1:
                    continue
    # calculate truncation and tailing length
            if Tail_seq.startswith( 'None' ):
                if len(Seq) >= len(miRNA):
                    Tail = len(Seq) - len(miRNA)
                    Truncation = 0
                elif len(Seq) < len(miRNA):
                    Tail = 0;
                    Truncation = len(miRNA) - len(Seq)
            else:
                if len(Com) <= len(miRNA):
                    Tail = len(Tail_seq)
                    Truncation = len(miRNA) - len (Com)
                elif len(Com) > len(miRNA):
                    Tail = len(Tail_seq) + (len(Com) - len(miRNA))
                    Truncation = 0
    # remove data points that are out of calculation range.
            if (Tail >(Cal_range-1) or Truncation > (Cal_range-1)):
                continue
    # decide if to include cononical miRNA
            if exclude_miR == 1:
                if (Tail ==0 and Truncation ==0):
                    continue
    # calculate abundnace for each position on truncation and tailing matrix
            abun = {}
            for lib in range (Num_libs):
                abun[lib] = int(data[5+lib])
                miRsum[lib] += abun[lib]
                miR[lib][Truncation][Tail] += abun[lib]
    f.close()
    mergePDF()


def PP(module,alist):
    print('***********Parallel instance of %s is being executed*********' % (module))
    
    start = time.time()
    ##PP is being used for Bowtie mappings - This will avoid overflooding of processes to server
    nprocPP = round((nproc/int(nthread))+1) ## 1 added so as to avoid 0 processor being allocated in serial mode
    print('\nnprocPP:%s\n' % (nprocPP))
    npool = Pool(int(nprocPP))
    npool.map(module, alist)

def main(inputFiles,Bowtie_index,miRNA_List):
    
    #inputfiles = ["MP_hen1-1-rep1.txt","MP_hen1-8-rep1.txt"] #str(sys.argv[1:3])#input files
#    outputfile = ["hen1_1_rep1_TEST_10_13_2014.txt","hen1_8_rep1_TEST_10_13_2014.txt"]#str(sys.argv[3:5]) #output_file_name
#    global bowtieindex
#    bowtieindex= "/alldata/Genomic/Arabidopsis/TAIR9/BowtieGenomicIndexes/AT_TAIR9_genome" #str(sys.argv[3])
    
#     for i in inputfiles:
#         dirname = pLCS_finder(i)   
    
    inputfiles=inputFiles
    global bowtieindex
    bowtieindex=Bowtie_index
    PP(pLCS_finder,inputfiles) ## Added by Atul
    #for i in inputfiles:
        #pLCS_finder(i)
    print ("STEP1 DONE!!")
    outfile="step2_output.txt"
    DIRNAME=DIR_NAME#"5pLCS_results_2015-08-05_15:21:19.830100"#DIR_NAME #directory name
    Second_step_output=merge_5plCS(DIRNAME,outfile) 
#    miRNA_file="new_miRNA_Test.fa"#str(sys.argv[4])
    miRNA_file=miRNA_List
    Third_step_output= tailling_Analysis(miRNA_file,Second_step_output)
    Fourth_step_output=format_Output(Third_step_output)
    genrate_Plots(Fourth_step_output)
    print("Tailing Pipeline is SUCCESSFULLY DONE!")
    
if __name__ == '__main__':   #calls the main function
    if nproc == 'Y':
        nproc = int(multiprocessing.cpu_count()*0.80)
    else:
        nproc == int(nproc)
        
    #create empty variables    
    inputFiles=[]
    Bowtie_index=""
    miRNA_List=""
    
    #adding positional parser for command line argument
    parser = argparse.ArgumentParser(description='Getting input for Tailing Analysis')
#    parser.add_argument('tag_countfiles',nargs='+',help='load tag_count files')
    parser.add_argument('tag_count_PATH',nargs=1, help='Provide path to the tag_count files')
    parser.add_argument('bowtie_PATH',nargs=1,help='Provide path to the bowtie index')
    parser.add_argument('mirna_FILE',nargs=1,help='Load miRNA file')
    parser.add_argument('output_PATH',nargs=1, help='Provide path to the output of Tailing Pipeline')
    args=parser.parse_args()    
#    inputFiles=args.tag_countfiles 
    
    INPUT_DIR=str(args.tag_count_PATH[0])
    for fn in os.listdir(INPUT_DIR):
        if fn.endswith(".txt"):
            inputFiles.append(fn)        
    Bowtie_index=str(args.bowtie_PATH[0])    
    miRNA_List=str(args.mirna_FILE[0])              
    # print (inputFiles)
    # print (Bowtie_index+"\n")
    # print (miRNA_List+"\n")  
    OUTPUT_DIR=str(args.output_PATH[0])   
    start_time = time.clock() # note start time of the exceution added by Parth after review of the paper- 8/6/2015
    main(inputFiles,Bowtie_index,miRNA_List)
    print (time.clock() - start_time, "seconds") # note start time of the exceution added by Parth after review of the paper- 8/6/2015    
    sys.exit()




