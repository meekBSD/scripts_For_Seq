import sys
import gzip

filename = sys.argv[1]

def get_reads(file_obj):
    items = [""] * 5
    for number , line in enumerate(file_obj):
        slot = number % 4
        items[0] = number
        items[slot + 1] = line.rstrip()
        if slot == 3:
            yield items

f = gzip.open(filename, mode='rb', compresslevel = 9, encoding = None, errors = None, newline = None)

for r in get_reads(f):
    seq = str(r[2])

    for n,b in enumerate(seq):
	
        if b == 'N': 
            ## print(str(r[1]))
            print(str(r[4])[n])

f.close()
    
