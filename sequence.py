import sys
import getopt
from Bio import SeqIO

def usage():
    message = """
    reads a fasta file and builds a dictionary 
    """
    print(message)

def displayRecords(filename):
    file = SeqIO.parse(filename, "fasta")
    for record in file:
        print(record.seq)
    
def lengthOfSequences(filename):
    file = SeqIO.parse(filename, "fasta")
    length = len([record for record in file])
    print(f"Length of sequences in a file is {length} ") 

def longestSequences(filename):
    file = SeqIO.parse(filename, "fasta")
    max = 0
    for record in file:
        if len(record.seq) > max:
            max = len(record.seq)
    
    print(f"Max Length of sequences in a file is {max} ") 


o, a = getopt.getopt(sys.argv[1:], 'l:h:g', ['length', 'longest', 'identifier', 'display'])
opts = {}
seqlen = 0

for k,v in o:
    opts[k] = v

# print(f"O : ${o} \na : ${a}\nOpts : ${opts} ")

if '-h' in opts.keys():
    usage()
    sys.exit()

if len(a) < 1:
    usage()
    sys.exit("Input fasta file is missing")

if '-g' in opts.keys():
    print("Run")
    lengthOfSequences(a[0])
    sys.exit()

if '-l' in opts.keys():
    if int(opts['l']) < 0:
        print("Length of sequence should be positive")
        sys.exit()
    seqlen = opts['l']
    print(seqlen)

if "--length" in opts.keys():
    lengthOfSequences(a[0])
    sys.exit()

if "--longest" in opts.keys():
    longestSequences(a[0])
    print("Test pass longest")

if "--identifier" in opts.keys():
    print("Test pass identifier")

if "--display" in opts.keys():
    displayRecords(a[0])
