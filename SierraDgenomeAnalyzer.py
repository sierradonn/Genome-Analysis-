from sequenceAnalysis import NucParams
import sys 
from sequenceAnalysis import FastAreader
def main(filename= None): # change to stdin
    myReader = FastAreader(filename)
    myNuc = NucParams()

    for head, seq in myReader.readFasta():
        myNuc.addSequence(seq)

    # Calculate sequence length and GC content
    seqLength = myNuc.nucCount()  
    gcContent = (myNuc.nucComposition()['G'] + myNuc.nucComposition()['C']) / seqLength * 100
    # using the nucleotide composition find percentage of G and C nucleotides in the sequence 

    print(f'Sequence length = {seqLength / 1e6:.2f} Mb')
    print(f'GC content = {gcContent:.1f}%')
    print('')

    # Sort codons in alphabetical order, by Amino Acid
    sortedCodons = sorted(myNuc.codonComposition().items(), key=lambda x: myNuc.rnaCodonTable[x[0]])

    # Iterate over sorted codons and print relative codon usage for each codon
    for codon, count in sortedCodons:
        aa = myNuc.rnaCodonTable[codon]
        val = (count / myNuc.aaComposition().get(aa) ) * 100
       
        print('{:s} : {:s} {:5.1f} ({:6d})'.format(codon, aa, val, count))
    
if __name__ == "__main__":
    main('tass2.fa')  # Change this if you want to use stdin or a different filename