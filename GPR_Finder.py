import os
import requests
import argparse
import json
import cobra
import pandas as pd
from Bio import Entrez, SeqIO
from glob import glob

def main():
    ## Define arguments for running the script
    parser = argparse.ArgumentParser()
    parser.add_argument('-template_model',type=str,default='iML1515',
                        help='BiGG ID for the model you want to use as a template')
    parser.add_argument('-target_genome',type=str,
                        help='Path to file of the genome (GenBank format e.g .gb) for the genome of your model (target)')
    parser.add_argument('-target_model',type=str,
                        help='Path to file of your model (target) in .json, .mat or .xml(SBML) format')
    parser.add_argument('-pid',type=int,default=60,
                        help='PID Threshold for filtering the blast results')
    parser.add_argument('-num_threads',type=int,default=1,
                        help='Number of threads for the BLAST algorithm')
    args = parser.parse_args()

    ## Make the necessary directories for the pipeline
    if os.path.isdir('genomes') == False:
        os.system('mkdir genomes')
    if os.path.isdir('prots') == False:
        os.system('mkdir prots')
    if os.path.isdir('bbh') == False:
        os.system('mkdir bbh')
    if os.path.isdir('results') == False:
        os.system(f'mkdir results')

    print(f'Starting GPRFinder. Template: {args.template_model} vs Target: {args.target_model}')

    ## Execute the pipeline
    target_genome_id = args.target_genome.split('.gb')[0].split('/')[1]
    # Donwload the genome of the requested BiGG model
    template_genome_id = Get_Model_Genome(args.template_model)
    # Parse genomes into fasta format
    parse_genome(template_genome_id)
    parse_genome(target_genome_id)
    # Make protein BLAST DB
    make_blast_db(target_genome_id, folder='prots', db_type='prot')
    make_blast_db(template_genome_id, folder='prots', db_type='prot')
    # Execute BLASTP to get BBHs
    get_bbh(template_genome_id,target_genome_id, in_folder='bbh', threads=args.num_threads)
    # Get missing genes
    missingGenes, GPRMissing = get_missing_genes(args.target_model, target_genome_id, 
                                                 template_genome_id, args.template_model, 
                                                 args.pid) 
    file_prefix = f'{args.template_model}_vs_{args.target_model.split(".")[0]}'
    missingGenes = missingGenes.rename(columns = {'gene':'template_gene', 'subject':'target_gene'})
    missingGenes.to_csv(f'results/{file_prefix}_missingGenes.csv', index=False)
    GPRMissing.to_csv(f'results/{file_prefix}_GPRsMissing.csv', index=False)
    print('--------------------------------------------------')
    print(f'GPRFinder done! Found {GPRMissing.shape[0]} GPR association rules')
    print(f'Results files located at results/ directory with the following prefix: {file_prefix}')
    print('--------------------------------------------------')

## Define function to get genome of reference model with BiGG and NCBI API
def Get_Model_Genome(referenceModelID, out_folder='genomes'):
    # Download template model from BiGG database in JSON format
    out_file = f'{referenceModelID}.json'
    if os.path.isfile(out_file) == False:
        os.system(f'curl -O http://bigg.ucsd.edu/static/models/{referenceModelID}.json')

    # Download model information from BiGG database and save data into a dict
    r = requests.get(f'http://bigg.ucsd.edu/api/v2/models/{referenceModelID}')
    dct = r.json()

    if 'ncbi_accession' in dct['genome_ref_string']:
        try:
            model_genome_id = dct['genome_name']
        except:
            print(f'Oops.. No genome was found for model {referenceModelID}')
    else:
        print(f'Oops.. No genome was found for model {referenceModelID}')

    # Download genome file (GenBank format) from NCBI's nucleotide database
    out_file = f"{out_folder}/{model_genome_id}.gb"
    files = glob(f'{out_folder}/*.gb')
    if out_file in files:
        print(f'Genome {model_genome_id}.gb from {referenceModelID} already downloaded')
    
    else:
        print(f'Downloading genome for model {referenceModelID}')
        Entrez.email = ""
        handle = Entrez.efetch(db="nuccore", id = model_genome_id, rettype="gbwithparts", retmode="text")
        fout = open(out_file,'w')
        fout.write(handle.read())
        fout.close()
        print(f'Downloaded genome: {model_genome_id}.gb')
    return(model_genome_id)

## Define function to parse genbank files into protein and nucleotide format
def parse_genome(id, type='prot', in_folder = 'genomes',out_folder='prots'):
    # Define input and output file
    in_file = f'{in_folder}/{id}.gb'
    out_file = f'{out_folder}/{id}.fa'
    files = glob(f'{out_folder}/*.fa')

    if out_file in files:
        print(f'Genome {id}.gb already parsed')
    else:
        # Parse the genome
        print(f'Parsing {id}.gb...')
        handle = open(in_file)
        fout = open(out_file,'w')
        x = 0
        records = SeqIO.parse(handle, "genbank")
        for record in records:
            for f in record.features:
                if f.type=='CDS':
                    seq=f.extract(record.seq)

                    if type=='nucl':
                        seq=str(seq)
                    else:
                        seq=str(seq.translate())

                    if 'locus_tag' in f.qualifiers.keys():
                        locus = f.qualifiers['locus_tag'][0]
                    elif 'gene' in f.qualifiers.keys():
                        locus = f.qualifiers['gene'][0]
                    else:
                        locus = 'gene_%i'%x
                        x+=1
                    fout.write('>%s\n%s\n'%(locus, seq))
        fout.close()
        print(f'Parsing done, results located in {out_folder}/{id}.fa')

## Define function for making a blast db
def make_blast_db(id,folder='prots',db_type='prot'):
    out_file = f'{folder}/{id}.fa.pin'
    files = glob(f'{folder}/*.fa.pin')

    if out_file in files:
        print (id, 'already has a blast db')
        return
    if db_type=='nucl':
        ext='fna'
    else:
        ext='fa'

    cmd_line = f'makeblastdb -in {folder}/{id}.{ext} -dbtype {db_type}'
    print ('making blast db with following command line...')
    print (cmd_line)
    os.system(cmd_line)

def run_blastp(seq,db,in_folder='prots', out_folder='bbh', out=None,outfmt=6,evalue=0.001,threads=1):
    if out==None:
        out = f'{out_folder}/{seq}_vs_{db}.txt'
        print(out)

    files =glob('%s/*.txt'%out_folder)
    if out in files:
        print (seq, 'already blasted')
        return

    print (f'blasting {seq} vs {db}')
    db = f'{in_folder}/{db}.fa'
    seq = f'{in_folder}/{seq}.fa'
    cmd_line = f'blastp -db {db} -query {seq} -out {out} -evalue {evalue} -outfmt {outfmt} -num_threads {threads}'
    print ('running blastp with following command line...')
    print (cmd_line)
    os.system(cmd_line)
    return out

## Define a function to get sequence length
def get_gene_lens(query, in_folder='prots'):
    file = '%s/%s.fa'%(in_folder, query)
    handle = open(file)
    records = SeqIO.parse(handle, "fasta")
    out = []

    for record in records:
        out.append({'gene':record.name, 'gene_length':len(record.seq)})
    out = pd.DataFrame(out)
    return out

## Define a function to get Bi-Directional BLASTp Best Hits
def get_bbh(query, subject, in_folder='bbh', threads=1):
    #Utilize the defined protein BLAST function
    run_blastp(query, subject, threads=threads)
    run_blastp(subject, query, threads=threads)
    query_lengths = get_gene_lens(query, in_folder='prots')
    subject_lengths = get_gene_lens(subject, in_folder='prots')

    #Define the output file of this BLAST
    out_file = f'{in_folder}/{query}_vs_{subject}_parsed.csv'
    files=glob(f'{in_folder}/*_parsed.csv')
    if out_file in files:
        print(f'{out_file} already parsed')
        return

    #Combine the results of the protein BLAST into a dataframe
    print ('parsing BBHs for', query, subject)
    cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 
            'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 
            'subjectEnd', 'eVal', 'bitScore']
    bbh=pd.read_csv(f'{in_folder}/{query}_vs_{subject}.txt', sep = '\t', names=cols)
    bbh = pd.merge(bbh, query_lengths)
    bbh['COV'] = bbh['alnLength']/bbh['gene_length']

    bbh2=pd.read_csv(f'{in_folder}/{subject}_vs_{query}.txt', sep='\t', names=cols)
    bbh2 = pd.merge(bbh2, subject_lengths)
    bbh2['COV'] = bbh2['alnLength']/bbh2['gene_length']
    out = pd.DataFrame()

    # Filter the genes based on coverage
    bbh = bbh[bbh.COV>=0.25]
    bbh2 = bbh2[bbh2.COV>=0.25]

    #Delineate the best hits from the BLAST
    for g in bbh.gene.unique():
        res = bbh[bbh.gene==g]
        if len(res)==0:
            continue
        best_hit = res.loc[res.PID.idxmax()]
        best_gene = best_hit.subject
        res2 = bbh2[bbh2.gene==best_gene]
        if len(res2)==0:
            continue
        best_hit2 = res2.loc[res2.PID.idxmax()]
        best_gene2 = best_hit2.subject
        if g==best_gene2:
            best_hit['BBH'] = '<=>'
        else:
            best_hit['BBH'] = '->'
        out=pd.concat([out, pd.DataFrame(best_hit).transpose()])

    #Save the final file to a designated CSV file
    out.to_csv(out_file, index=False)

def get_missing_genes(targetModelFile, targetGenomeID, templateGenomeID, templateModelID, PIDThreshold):
    # Load target model
    ext = targetModelFile.split('.')[1]
    if ext == 'json':
        model = cobra.io.load_json_model(targetModelFile)
    elif ext == 'mat':
        model = cobra.io.load_matlab_model(targetModelFile)
    elif ext == 'xml':
        model = cobra.io.read_sbml_model(targetModelFile)
    
    else:
        print('Error loading model, please provide your model in json, matlab or SBML (xml) format')
    
    # Load template model
    templateModel = cobra.io.load_json_model(f'{templateModelID}.json')
    # Get genes of target model
    targetGenes = []
    for g in model.genes:
        targetGenes.append(g.id)
    target_genes_df = pd.DataFrame(targetGenes)
    target_genes_df.columns = ['subject']
    
    # Get genes of template model
    templateGenes = []
    for g in templateModel.genes:
        templateGenes.append(g.id)
    template_genes_df = pd.DataFrame(templateGenes)
    template_genes_df.columns = ['gene']

    # Import blast results and merge with target model genes
    blast = pd.read_csv(f'bbh/{templateGenomeID}_vs_{targetGenomeID}_parsed.csv')
    merged_df = pd.merge(blast, target_genes_df, on = 'subject')
    not_merged = blast[~blast.subject.isin(merged_df.subject)].dropna()
    filtered_target_genes = not_merged[['gene', 'subject', 'PID', 'eVal', 'COV']]

    # Filter results with genes present in the template model and according to the PID threshold specified
    missingGenes = pd.merge(filtered_target_genes, template_genes_df, on = 'gene')
    missingGenes = missingGenes[missingGenes.PID >= PIDThreshold]

    # Generate dataframe with GPR information
    GeneReactionsMissing = pd.DataFrame(columns = ['template_gene','target_gene',
                                      'reactionID','reaction_name','subsystem', 
                                      'reaction_string'])
    
    for index,row in missingGenes.iterrows():
        for r in templateModel.genes.get_by_id(row['gene']).reactions:
            info_list = []
            info_list.append(row['gene'])
            info_list.append(row['subject']) 
            info_list.append(r.id)
            info_list.append(r.name)
            info_list.append(r.subsystem)
            info_list.append(r.reaction)
            GeneReactionsMissing.loc[len(GeneReactionsMissing)] = info_list
    return(missingGenes, GeneReactionsMissing)

if __name__ == '__main__':
    main()