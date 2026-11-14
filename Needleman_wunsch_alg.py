import argparse # do obsługi argumentów 
import numpy as np # tworzenie macierzy 
from Bio import SeqIO # odczytywanie sekwencji biologicznych z FASTA
from pathlib import Path # obsługa ścieżek, katalogów i plików 
#Algorytm Needlemana-Wunscha dla globalnego dopasowywania sekwencji 

def needleman_wunsch(seq1, seq2): #funkcja przyjmuje dwie sekwencje
    #seq1 - umieszczona w kolumnach macierzy 
    #seq2 - umieszczona w wierszach macierzy 

    #parametry punktacji 
    match = 1 #gdy nukleotydy są identyczne 
    mismatch = -1 #gdy są różne nukleotydy 
    gap = -2 # wstawienie przerwy 

    n = len(seq2) #wiersze 
    m = len(seq1) #kolumny 

    #utworzenie macierzy NumPy na początku dwuwymiarowa tablica wypełniona samymi zerami 
    M = np.zeros((n+1, m+1), dtype=int) 
    #dodajemy +1, bo potrzebujemy dodatkowego wiersza i kolumny 

    #M[:, 0] - wybór wszystkich wierszy z kolumnny[0] 
    #M[0,:] - wszystkie kolumny z wiersza[0] 
    #funckja np.arange(n+1) - generuje liczby całkowite od 0 do n
    #służy do utworzenia ciągu 0, -2, -4, -6, -8 
    M[:, 0] = np.arange(n+1) * gap #utworzenie dla pustych prefiksów 0, -2, -4, -6, -8 .. 
    M[0, :] = np.arange(m+1) * gap 

    #zagnieżdżone pętle to implementacja wzoru ze slajdu 
    for i in range(1, n+1): #przechodzi po wierszach 
        for j in range(1, m+1): #przechodzi po kolumnach 

            #porównywanie czy nukleotydy są identyczne? 
            if seq2[i-1] == seq1[j-1]: 
                s = match 
            else: 
                s = mismatch 

            #obliczenie możliwych wartości z 3 kierunków: 
            # przekątna 
            # góra 
            # dół 
            diagonal = M[i-1, j-1] + s 
            up = M[i-1,j] + gap 
            left = M[i, j-1] + gap 

            M[i,j] = max(diagonal, up, left) #wybór maksymalnej wartości z 3 kierunków 

    dopasowania = [] #lista, która przechowuje znalezione dopasowania 

    #Funkcja rekurencyjna (odwołuje sie sama do siebie), do cofania się w macierzy 
    #Przyjmuje 4 argumenty: 
    # i indeks wiersza 
    # j indeks kolumny 
    # a odtworzona sekwencja1
    # b odtworzona sekwencja2 
    def traceback(i, j, a, b): 

        if i == 0 and j==0: 
            dopasowania.append((a[::-1], b[::-1])) #przechodzi przez znaki w odwrotnej kolejności 
            #bo przechodzimy od ostatnich nukleotydów do pierwszych 
            return 


        #Możliwości przejścia podczas cofania się: 
        if i>0 and j>0: #przejscie po przekatnej 
            s = match if seq2[i-1] == seq1[j-1] else mismatch 
            if M[i, j] == M[i - 1, j - 1] + s:
                traceback(i - 1, j - 1, a + seq1[j - 1], b + seq2[i - 1])

        if i > 0 and M[i, j] == M[i - 1, j] + gap: #przejscie do góry 
            #czy komórka mogła zostać uzyskana z komórki znajdującej się nad nią? 
            traceback(i - 1, j, a + "-", b + seq2[i - 1])

        if j > 0 and M[i, j] == M[i, j - 1] + gap: #przejście w lewo 
            traceback(i, j - 1, a + seq1[j - 1], b + "-")

    traceback(n, m, "", "") #puste znaki na dopasowaną sekwencje 

    return M, dopasowania #zwrocenie macierzy i dopasowań 

def main(): 
    #Argumenty programu 
    parser = argparse.ArgumentParser() 
    parser.add_argument("fasta") 
    parser.add_argument("--output", default="output.txt")
    parser.add_argument("--matrix", action="store_true") 
    args = parser.parse_args() 

    path = Path(args.fasta) 

    if not path.is_file(): 
        parser.error("Nie znaleziono pliku!") 

    rekordy = list(SeqIO.parse(path, "fasta")) 

    if len(rekordy) !=2: 
        parser.error("Plik musi zawierac dwie sekwencje") 

    query = str(rekordy[0].seq).upper() 
    target = str(rekordy[1].seq).upper() 

    M, dopasowania = needleman_wunsch(query, target) 
    if args.matrix:
        print("Macierz punktacji:") 
        print(M)

    with open(args.output, "w") as file:
        file.write(f"Query: {query}\n")
        file.write(f"Target: {target}\n")
        file.write(f"Score: {M[-1, -1]}\n")
        file.write(f"Liczba dopasowań: {len(dopasowania)}\n\n")

        for k, (a, b) in enumerate(dopasowania, 1):
            file.write(f"Dopasowanie {k}:\n{a}\n{b}\n\n")

    print(f"\nScore: {M[-1, -1]}")
    print(f"Liczba dopasowań: {len(dopasowania)}")
    print(f"Wyniki zapisano do {args.output}")

if __name__=="__main__": 
    main() 
