import random


b_list = ['A', 'C', 'T',  'G']
low_Q = ['(', ')', '*']

High_Q = [ 'D', 'E', 'B']

outFastq =open('Test_A_R1.fastq', 'w')

for i in range(3000):
    R_id = '@simulated_{0:06d}'.format(i)

    seq = ""
    qual = ""
    for j in range(100):
        base = random.randint(0,3)
        seq += b_list[base]

        q_ind = random.randint(0,2)
        qual  += High_Q[q_ind]
    for z in range(50):
        base = random.randint(0,3)
        seq += b_list[base]
        q_ind  = random.randint(0,2)

        qual  += low_Q[q_ind]

    outFastq.write('{0}\n{1}\n+\n{2}\n'.format(R_id, seq, qual))


Nbase_L = ['A', 'C', 'N', 'G', 'T']
for i in range(8000,10000):
    R_id = '@simulated_{0:06d}'.format(i)

    seq = ""
    qual = ""
    for j in range(150):
        base = random.randint(0,4)
        seq += Nbase_L[base]

        q_ind = random.randint(0,2)
        qual  += High_Q[q_ind]

    outFastq.write('{0}\n{1}\n+\n{2}\n'.format(R_id, seq, qual))

for i in range(12000,14000):
    R_id = '@simulated_{0:06d}'.format(i)

    seq = ""
    qual = ""
    for j in range(100):
        base = random.randint(0,4)
        seq += Nbase_L[base]

        q_ind = random.randint(0,2)
        qual  += High_Q[q_ind]
    for z in range(50):
        base = random.randint(0,4)
        seq += Nbase_L[base]
        
        q_ind  = random.randint(0,2)
        qual  += low_Q[q_ind]
    outFastq.write('{0}\n{1}\n+\n{2}\n'.format(R_id, seq, qual))

outFastq.close()
