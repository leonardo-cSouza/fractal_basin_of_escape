import sys
import math

# Configuration Parameters
N = int(sys.argv[1])
da = int(sys.argv[2])
db = int(sys.argv[3])

# Read data and calculate i1 and i2 counts for each box
M = {}
i = 1
j = 1

for line in sys.stdin:
    values = line.strip().split()
    M[(i, j)] = int(values[0])
    j += 1
    if j > N:
        j = 1
        i += 1

i0t = 0
i1t = 0
i2t = 0
output_data = []

for a in range(1, N+1, da):
    for b in range(1, N+1, db):
        i0 = 0
        i1 = 0
        i2 = 0
        for k in range(a, a+da):
            for l in range(b, b+db):
                if M[k, l] == 0:
                    i0 += 1
                elif M[k, l] == 1:
                    i1 += 1
                elif M[k, l] == 2:
                    i2 += 1

        if i0 == 0:
            output_data.append((i1, i2))
            i1t += i1
            i2t += i2

# Calculate entropy
S = 0
for i in range(len(output_data)):
    n = output_data[i][0] + output_data[i][1]
    if n == 0:
        print("DIVISÃO POR ZERO!")
        exit(1)
    
    P1 = output_data[i][0] / n
    P2 = output_data[i][1] / n
    S1 = -math.log(P1**P1)
    S2 = -math.log(P2**P2)
    S += S1 + S2

S /= len(output_data)
print(f"{S:.5f}")
