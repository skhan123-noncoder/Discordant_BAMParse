import pandas as pd
import pysam
import argparse
import logging
import logging.config
import os
import time
import numpy as np

start_time = time.time()

logger=logging.getLogger(__name__)
logging.basicConfig(filename='run.log', filemode='w', level=logging.DEBUG)

def main():
    usr_input=parse_args()
    try:
        discordant_reads=extract_discordant_reads(usr_input.sorted_bam_file,usr_input.outdir)
        analysis=analyse_discordant_reads(discordant_reads)
        
        write_output=writeout(discordant_reads, analysis[0], analysis[1], analysis[2], usr_input.outdir) 

    except Exception as e:
        print("Something wrong! Please check log file")
        logger.exception("Unexpected exception! %s",e)
        

def parse_args():
    
    parser=argparse.ArgumentParser(prog='BamDiscordantParser', description="Provides information about paired-end reads which are mapped to different chromosomes", epilog='none')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-b', '--bam_file', action='store', dest='sorted_bam_file', type=str, required=True, help='Sorted and indexed BAM file to be analysed')
    parser.add_argument('-o', '--output_folder', action='store', dest='outdir', type=str, help="name of the output directory", default="out_dir")
    args=parser.parse_args()

    logger.debug('Imported arguments from user')
    
    return(args)

def extract_discordant_reads(BAM_file, outdir):

    #Generate the output folder for storing results
    os.mkdir(os.path.join(outdir)+os.path.dirname(outdir))

    #Read the bam file using Pysam
    bamfile= pysam.AlignmentFile(BAM_file, "rb")

    print('Bam file read successfully')

    #Generate lists to store essential values
    query_name=[]
    flag=[]
    reference_name=[]
    reference_start=[]
    mapping_quality=[]
    cigarstring=[]
    next_reference_name=[]
    seq=[]

    #FLAG information for reads which are discordant
    discordant_flags = [113, 177, 81, 161, 97, 145, 65, 129]

    #Append BAM info if reads are discordant
    for read in bamfile.fetch():
        if read.flag in discordant_flags:
            query_name.append(read.query_name)
            flag.append(read.flag)
            reference_name.append(read.reference_name)
            reference_start.append(read.reference_start)
            mapping_quality.append(read.mapping_quality)
            cigarstring.append(read.cigarstring)
            next_reference_name.append(read.next_reference_name)
            seq.append(read.seq)

    #Store the information from BAM reads as a dataframe
    bam_df= pd.DataFrame({"query_name":query_name,
                      "chr_name":reference_name,
                      "start_pos":reference_start,
                      "MapQ":mapping_quality,
                      "chr_of_mate":next_reference_name,
                      "seq":seq,
                      "CigarString":cigarstring})

    #Sort the reads
    bam_df=bam_df.sort_values('query_name', ignore_index=True)

    print('Discordant reads formatted as dataframe')

    return(bam_df)

def analyse_discordant_reads(bam_df):

    print('Analysing discordant reads now')
    #Find the reads which have a mate pair and store them as a list
    duplicated_reads=bam_df[bam_df.duplicated(subset=['query_name'],
                                             keep=False)]
    dup_np=list(duplicated_reads['query_name'])

    #Convert the bam_df to a np_array for easier handling
    bam_np = np.array(bam_df)

    #Generate empty lists to store the reads
    same_chr_mate1=[]
    same_chr_mate2=[]
    diff_chr_mate1=[]
    diff_chr_mate2=[]
    singles=[]

    #Define a variable to access the numpy array
    i=0
    while i<len(bam_np):

        #This block will process only those reads which have a mate pair
        if bam_np[i][0] in dup_np:

            #The bam_np is sorted by read name. Hence checking for only the next read
            if bam_np[i][1]==bam_np[i+1][1]:

                #Store the values of the read and the mate which are on the same chr
                same_chr_mate1.append(bam_np[i])
                same_chr_mate2.append(bam_np[i+1])
            else:

                #Store the values of the read and the mate which are on the different chr
                diff_chr_mate1.append(bam_np[i])
                diff_chr_mate2.append(bam_np[i+1])
            
            #Increase the value of i by 2. Will skip the mate pair read
            i=i+2
        
        #This block will take care of reads lacking a mate pair
        else:
            singles.append(bam_np[i])

            #Increase the value of i by 1 here. Makes sure the next unique read is analysed
            i=i+1
    
    #Concatenate the list with reads belonging to same chr into a dataframe and label columns
    same_chr_df = pd.concat([pd.DataFrame(same_chr_mate1),pd.DataFrame(same_chr_mate2)], axis=1)

    if len(same_chr_df)>0:
        same_chr_df.columns=["Read1_query_name", "Read1_chr_name", "Read1_start_pos", "Read1_MapQ", "Read1_mate_chr",
                      "Read1_sequence", "Read1_cigarstring", "Read2_query_name", "Read2_chr_name", "Read2_start_pos",
                     "Read2_MapQ", "Read2_mate_chr",
                      "Read2_sequence", "Read2_cigarstring"]

        same_chr_df = same_chr_df.drop(["Read1_mate_chr", "Read2_mate_chr"], axis=1)

    #Concatenate the list with reads belonging to diff chr into a dataframe and label columns
    diff_chr_df = pd.concat([pd.DataFrame(diff_chr_mate1),pd.DataFrame(diff_chr_mate2)], axis=1)

    if len(diff_chr_df)>0:
        diff_chr_df.columns=["Read1_query_name", "Read1_chr_name", "Read1_start_pos", "Read1_MapQ", "Read1_mate_chr",
                      "Read1_sequence", "Read1_cigarstring", "Read2_query_name", "Read2_chr_name", "Read2_start_pos",
                     "Read2_MapQ", "Read2_mate_chr",
                      "Read2_sequence", "Read2_cigarstring"]

        diff_chr_df = diff_chr_df.drop(["Read1_mate_chr", "Read2_mate_chr"], axis=1)

    #Make a dataframe for single reads
    singles_df = pd.DataFrame(singles)
    
    return(same_chr_df, diff_chr_df, singles_df)


def writeout(bam_df, same_df, diff_df, singles_df, out_folder):

    print('Writing output files')

    #Write all Discordant read information
    bam_df.to_csv(f"{os.path.join(out_folder)}/All_discordant.txt", sep='\t', index=False)

    #Write Discordant reads aligning to the same chromosome
    same_df.to_csv(f"{os.path.join(out_folder)}/Discordant_reads_aligning_to_same_chr.txt", sep='\t', index=False)

    #Write Discordant reads aligning to the different chromosomes
    diff_df.to_csv(f"{os.path.join(out_folder)}/Discordant_reads_aligning_to_diff_chr.txt", sep='\t', index=False)

    #Write Singles read info
    if len(singles_df)>0:
        singles_df.to_csv(f"{os.path.join(out_folder)}/Single_reads.txt", sep='\t', index=False)

    logger.debug('Total time for the entire process in seconds is %s' %(time.time()-start_time))


if __name__=='__main__':
    main()

