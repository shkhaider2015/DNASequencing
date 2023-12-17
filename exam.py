import sys
import getopt
from Bio import SeqIO
from collections import defaultdict

def usage():
    message = """
    reads a fasta file and builds a dictionary 
    """
    print(message)

def displayRecords(filename):
    file = SeqIO.parse(filename, "fasta")
    for record in file:
        print(record)

def recordLength(filename):
    file = SeqIO.parse(filename, "fasta")
    length = len([record for record in file])
    print(f"Length of sequences in a file is {length} ") 

def lengthOfSequences(filename):
    file = SeqIO.parse(filename, "fasta")
    lengths = [len(record) for record in file]
    print(lengths)

def longestSequence(filename):
    file = SeqIO.parse(filename, "fasta")
    max = 0
    for record in file:
        if len(record.seq) > max:
            max = len(record.seq)

def shortestSequence(filename):
    file = SeqIO.parse(filename, "fasta")
    min = 5000
    for record in file:
        if len(record.seq) < min:
            min = len(record.seq)

def find_orfs(sequence):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []

    for i in range(len(sequence)):
        if sequence[i:i+3] == start_codon:
            for j in range(i+3, len(sequence), 3):
                codon = sequence[j:j+3]
                if codon in stop_codons:
                    orf = sequence[i:j+3]
                    orfs.append((i+1, j+3, orf))
                    break

    return orfs

def find_longest_orf_in_sequence(sequences):
    longest_orf = ""
    longest_sequence_id = ""
    longest_orf_length = 0
    longest_orf_start = 0

    for seq_id, sequence in sequences.items():
        for frame in range(1, 4):
            current_orfs = find_orfs(sequence[frame-1:])
            for start, end, orf in current_orfs:
                if len(orf) > longest_orf_length:
                    longest_orf = orf
                    longest_sequence_id = seq_id
                    longest_orf_length = len(orf)
                    longest_orf_start = start

    return longest_orf_length, longest_sequence_id, longest_orf_start, longest_orf

def find_all_orfs(filename):
    file = SeqIO.parse(filename, "fasta")
    orfs = [ find_orfs(record.seq) for record in file  ]
    print(orfs)

def find_longest_orf_details(sequence):
    length, seq_id, start, orf = find_longest_orf_in_sequence(sequence)
    print(f"Length of the longest ORF: {length}")
    print(f"Identifier of the sequence containing the longest ORF: {seq_id}")
    print(f"For sequence {seq_id}, the longest ORF is: {orf}")
    print(f"Starting position of the longest ORF: {start}")


def find_repeats(sequence, n):
    repeats = defaultdict(int)

    for i in range(len(sequence) - n + 1):
        repeat = sequence[i:i+n]
        repeats[repeat] += 1

    return repeats

def find_most_frequent_repeat(sequences, n):
    all_repeats = defaultdict(int)

    for seq_id, sequence in sequences.items():
        repeats = find_repeats(sequence, n)
        for repeat, count in repeats.items():
            all_repeats[repeat] += count

    most_frequent_repeat = max(all_repeats, key=all_repeats.get)
    frequency = all_repeats[most_frequent_repeat]

    return most_frequent_repeat, frequency

def most_frequent_repeat(sequenc, repeat_length):
    most_frequent_repeat, frequency = find_most_frequent_repeat(sequenc, repeat_length)

    print(f"Most frequent repeat of length {repeat_length}: {most_frequent_repeat}")
    print(f"Frequency of the most frequent repeat: {frequency}")

o, a = getopt.getopt(sys.argv[1:], 'l:h', ['list-orf'])
opts = {}
seqlen = 0

for k,v in o:
    opts[k] = v


if '-h' in opts.keys():
    usage()
    sys.exit()

if len(a) < 1:
    usage()
    sys.exit("Input fasta file is missing")

if '-l' in opts.keys():
    if int(opts['l']) < 0:
        print("Length of sequence should be positive")
        sys.exit()
    seqlen = opts['l']
    print(seqlen)

if '--list-orf' in opts.keys():
    print("Listing ...")
    find_all_orfs(a[0])