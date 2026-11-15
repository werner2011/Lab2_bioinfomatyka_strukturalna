import argparse 
import numpy as np 
from Bio import SeqIO 
from pathlib import Path 

def needleman_wunsch(seq1, seq2): 
    match = 1 
    mismatch = -1 
    gap = -2 

    n = len(seq2) 
    m = len(seq1) 

    M = np.zeros((n+1, m+1), dtype=int) 

    M[:,0] = np.arange(n+1) * gap 
    M[0, :] = np.arange(m+1) * gap 

    for i in range(1, n+1):
        for j in range(1, m+1): 
            if seq2[i-1] == seq1[j-1]: 
                s = match 
            else:
                s = mismatch 

            diagonal = M[i-1, j-1] + s 
            up = M[i-1,j] + gap 
            left = M[i, j-1] + gap 
            M[i,j] = max(diagonal, up, left) 

    dopasowania = [] 

    def traceback(i, j, a, b): 
        if i==0 and j ==0 :
            dopasowania.append((a[::-1], b[::-1])) 
            return 

        if i>0 and j>0:
            s = match if seq2[i-1] == seq1[j-1] else mismatch 
            if M[i,j] == M[i-1, j-1] + s: 
                traceback(i-1, j-1, a+seq1[j-1], b+seq2[i-1]) 

        if i > 0 and M[i,j] == M[i-1, j] + gap: 
            traceback(i-1, j, a+ "-", b + seq2[i-1]) 

        if j>0 and M[i,j] == M[i,j-1] + gap: 
            traceback(i, j-1, a+seq1[j-1], b + '-') 

    traceback(n,m,"","") 
    return M, dopasowania 

def main(): 
    parser = argparse.ArgumentParser() 
    parser.add_argument("fasta") 
    parser.add_argument("--output", default="output.txt") 
    parser.add_argument("--matrix", action="store_true") 
    args = parser.parse_args() 
    path = Path(args.fasta) 

    if not path.is_file(): 
        parser.error("Nie znaleziono pliku") 

    rekordy = list(SeqIO.parse(path, "fasta")) 
    if len(rekordy) !=2: 
        parser.error("Plik musi zawierac dwie sekwencje!") 

    query = str(rekordy[0].seq).upper()
    target = str(rekordy[1].seq).upper() 

    M, dopasowania = needleman_wunsch(query, target) 
    if args.matrix:
        print("Macierz punktacji: ")
        print(M) 

    with open(args.output, "w") as file: 
        file.write(f"Query: {query}\n") 
        file.write(f"Target: {target}\n")
        file.write(f"Score: {M[-1, -1]}\n")
        file.write(f"Liczba dopasowan: {len(dopasowania)}\n\n")
        
        for k, (a, b) in enumerate(dopasowania, 1):
            file.write(f"Dopasowanie {k}:\n")
            file.write(f"{a}\n{b}\n\n")

    print(f"\nScore: {M[-1,-1]}") 
    print(f"Liczba dopasowan: {len(dopasowania)}") 
    print(f"Wyniki zapisano do {args.output}") 

if __name__=="__main__":
    main() 