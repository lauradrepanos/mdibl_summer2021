#Laura Drepanos MDIBL 7/1/21. Experiment jcoffman008
#Usage: Run in directory with countmatrix with command python3 countMatrixCondense.py -i countmatrixfilename
#WHAT IT DOES: reads count matrix and writes a new .txt file that collapses all rows with the same gene name, keeping the count sum under the smaller geneID 
#PURPOSE: addressing the issue of one gene (gene name) being split up into different versions (GeneIDs) that causes misleading differential expression results. etc
#updated version of collapsingGeneIDs.py from 6/16/21 to enable use on user-specified count matrix 

import argparse

def matrixFindDuplicates(infile):
    #Takes in a count matrix and defines the rows of a new count matrix that sums rows with the same gene name and collapses them under the first Gene ID
    #going through the count matrix to find duplicates
    f = open(infile, 'r') 
    header=f.readline() #pass header
    rowskeep=[] #store the rows that will go into the new version of the count matrix here
    genenames=[] #store gene names here to identify when we come across a repeat gene name (check after to see how it handles case)
    num_duplicates=0 #tracking how many rows are not added b/c repeat gene name
    cols= header.split('\t')
    num_cols= len(cols)+1 #the number of columns a row will have after ID and gene name are split *IF the ID has an associated gene name* 
    while True: #go through each line until the code breaks (the row is empty)
        contents=f.readline()
        if len(contents)<1:
            break #break out of while loop if the row is empty
        terms = contents.split('\t')
        ID_and_gene = terms[0].split('_',1)
        terms=ID_and_gene+terms[1:]
        genename=terms[1]
        genename=genename.lower() #this line makes it so that a gene name registers as a duplicate if one is all caps and one is lowercase
        if genename in genenames and len(terms)==num_cols: #if the gene listed in this row is already associated with another geneID...
            num_duplicates +=1
            index=genenames.index(genename) #finding the row to add onto 
            for c in range(2,25): #add these counts on to the counts from the first occurence of this gene in the matrix 
                rowskeep[index][c]=float(rowskeep[index][c])+float(terms[c])
            #print(rowskeep[index][0]+"\t"+ terms[0]+"\t"+genename) --this prints the list of gene name duplicates
        else: #if this is the first occurance of this gene name
            genenames.append(genename)
            rowskeep.append(terms)
    f.close()
    return (rowskeep,header,num_cols)

def matrixWrite(rowskeep,infile,header,num_cols):
    #Builds the new count matrix with rows defined in matrixFindDuplicates
    outfile=infile[:-16]+"_nogeneduplicates.gene.counts.txt"
    newcountmatrix= open(outfile, 'wt')
    newcountmatrix.write(header)
    for row in rowskeep:
        #first rejoining GeneID and genename so that new count matrix is identical to old ones (just without the gene name duplicates)
        if len(row)==(num_cols): #this is here because some gene IDs have no associated gene name after a _ , the line below would cause these two have the first count as the gene name
            row=[str(row[0])+"_"+str(row[1])]+row[2:]
        #newcountmatrix.write("\n")
        newcountmatrix.write(str(row[0]))
        for cell in row[1:]:
            newcountmatrix.write("\t"+str(cell))
    newcountmatrix.close()


def main():
    parser = argparse.ArgumentParser(
        description='Input file name of count matrix to remove gene duplicates from (will generate a new .txt file).')

    parser.add_argument("-i", "--infile", required=True, dest="infile",
                        help="Enter counts matrix, should end in .gene.counts.txt")

    parsed_args= parser.parse_args()
    print("Reading count matrix for gene name duplicates...")
    (rowskeep,header,num_cols)=matrixFindDuplicates(parsed_args.infile)
    print("Writing new count matrix...")
    matrixWrite(rowskeep,parsed_args.infile,header,num_cols)
    print("New count matrix generation successful")

if __name__ == '__main__':
    main()
