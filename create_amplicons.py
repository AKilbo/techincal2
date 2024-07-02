import argparse
import requests
import pandas as pd
from Bio.Seq import Seq
import primer3


#blast was taking too long, and crashing my docker images, so I looked around for something that was faster and found uscs blat
#this function takes a sequence and returns the coordinates of the sequence in the human genome
def blat_coords(seq):
    url = f"https://genome.ucsc.edu/cgi-bin/hgBlat?userSeq={seq}&type=DNA&db=hs1&output=json"
    response = requests.get(url)
    json_response = response.json()

    #convert the json response to a dataframe
    cords_df = pd.DataFrame(json_response['blat'], columns=json_response['fields'])
    cords_df = cords_df[cords_df['matches'] == 20]
    #change the tstart column that is tstart + 1 because the UCSC genome browser has some weird indexing
    cords_df['tStart'] = cords_df['tStart'] + 1
    return cords_df

#this function gets the sequence around the coordinates of the guide RNA that was found by the blat function
def get_sequence_around_coords(cords_df, flanking_length=250):
    chrom = cords_df['tName'].iloc[0]
    start = int(cords_df['tStart'].iloc[0]) - flanking_length
    end = int(cords_df['tEnd'].iloc[0]) + flanking_length
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hs1;chrom={chrom};start={start};end={end}"
    response = requests.get(url)
    dna= response.json()['dna']
    if cords_df['strand'].iloc[0] == '-':
        dna = str(Seq(dna).reverse_complement())
    return dna.upper() 
    
#this function takes the sequence around the guide RNA and the guide RNA sequence and returns the primers    
def get_primers(dna,guide_seq, target_pos=250, target_len=20):
    # Set the sequence arguments
    seq_args = {
        'SEQUENCE_ID': 'example',
        'SEQUENCE_TEMPLATE': dna,
        'SEQUENCE_INCLUDED_REGION': [0,len(dna)],
        # the target position is the position of the guide RNA in the sequence
        'SEQUENCE_TARGET': [target_pos,target_len]
        
    }
    #set size of the final product
    global_args = {'PRIMER_PRODUCT_SIZE_RANGE':[150,250]}

    # Call the function
    results = primer3.bindings.design_primers(seq_args, global_args)


    df = pd.DataFrame()
    df['guide_seq'] = [guide_seq]
    df['primer_left'] = results['PRIMER_LEFT_0_SEQUENCE']
    df['primer_right'] = results['PRIMER_RIGHT_0_SEQUENCE']
    df['primer_left_gc'] = results['PRIMER_LEFT_0_GC_PERCENT']
    df['primer_right_gc'] = results['PRIMER_RIGHT_0_GC_PERCENT']
    df['primer_left_start'] = results['PRIMER_LEFT_0'][0]
    df['primer_left_end'] = results['PRIMER_LEFT_0'][0] + len(results['PRIMER_LEFT_0_SEQUENCE'])
    df['primer_right_start'] = results['PRIMER_RIGHT_0'][0]
    df['primer_right_end'] = results['PRIMER_RIGHT_0'][0] + len(results['PRIMER_RIGHT_0_SEQUENCE'])
    df['primer_left_length'] = results['PRIMER_LEFT_0'][1]
    df['primer_right_length'] = results['PRIMER_RIGHT_0'][1]
    #return first row of df
    return df.iloc[0]
    # return sequence

#this function takes the primer dataframe and the coordinates dataframe and returns the coordinates of the primers
# theres some funky UCSC genome browser coordinate stuff, so I played around with adding 1 at various places 
def get_primer_coords(primer_df,cords_df, dna, guide_rna_seq):
    #add code to get the correct primer coordinates for ones on the reverse strand, by flipping the left and right primers coordinates
    # I don't think this will get the correct coordinates for the reverse strand currently

    #get the coordinates of the left primer based on the coordinates of the guide RNA
    left_start_diff = 250 - primer_df['primer_left_start']
    left_start_coord = cords_df['tStart'].iloc[0] - left_start_diff
    left_end_coord = left_start_coord + len(primer_df['primer_left'])

    #get the coordinates of the right primer based on the coordinates of the guide RNA
    right_start_spacing = primer_df['primer_right_start'] - primer_df['primer_left_start']
    right_start_coord = left_start_coord + right_start_spacing + 1
    right_end_coord = right_start_coord - primer_df['primer_right_length'] +1 
   
    #format the coordinates as a string
    left_primer_coords = str(cords_df['tName'].iloc[0]) + ':' + str(left_start_coord) + '-' + str(left_end_coord)
    right_primer_coords = str(cords_df['tName'].iloc[0]) + ':' + str(right_start_coord) + '-' + str(right_end_coord)


    return left_primer_coords, right_primer_coords

#this is the main function that reads the tsv file and calls the other functions
def main():
    #parse the arguments
    parser = argparse.ArgumentParser(description='Create guides from a TSV file')
    parser.add_argument('tsv_path', type=str, help='The path to the TSV file')
    args = parser.parse_args()

    grna = pd.read_csv(args.tsv_path, sep='\t')
    output_df = pd.DataFrame()

    #iterate over the rows of the tsv file
    for index, row in grna.iterrows():


        guide_seq = row['guide_seq']

        #get th genomic coordinates of the guide RNA
        cords_df= (blat_coords(guide_seq))
        #get the sequence around the guide RNA        
        full_seq = get_sequence_around_coords(cords_df)

        #get the primers using primer3
        primer_df = get_primers(full_seq,guide_seq)

        #ready the data for the output
        guide_rna_name = row['guide_name']
        guide_rna_seq = row['guide_seq']
        guide_rna_coords = str(cords_df['tName'].iloc[0]) +":" +  str(cords_df['tStart'].iloc[0]) + "-" + str(cords_df['tEnd'].iloc[0])
        primer_left = primer_df['primer_left']
        primer_right = primer_df['primer_right']

        #get the coordinates of the primers in the genome using the coordinates of the guide RNA
        primer_left_coords, primer_right_coords = get_primer_coords(primer_df,cords_df, full_seq, guide_rna_seq)

        primer_left_gc = primer_df['primer_left_gc']
        primer_right_gc = primer_df['primer_right_gc']
        
        amplicon_seq = full_seq[primer_df['primer_left_start']:primer_df['primer_right_start']+1]
        
        data = {
                'guide_rna_name': [guide_rna_name],
                'guide_rna_seq': [guide_rna_seq],
                'guide_rna_coords': [guide_rna_coords],
                'strand' : [cords_df['strand'].iloc[0]],
                'primer_left_seq': [primer_left],
                'primer_left_coords': [primer_left_coords],
                'primer_left_gc': [primer_left_gc],
                'primer_right_seq': [primer_right],
                'primer_right_coords': [primer_right_coords],
                'primer_right_gc': [primer_right_gc],
                'amplicon_seq': [amplicon_seq]
            }

        df = pd.DataFrame(data)

        output_df = pd.concat([output_df, df], ignore_index=True)

    print('see output.csv for the results')
    print(output_df)
    output_df.to_csv('output.csv', index=False)


if __name__ == "__main__":
    main()

    