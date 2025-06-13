# Pour rouler le programme il suffit de rentrer la commande python3 (python si vous utilisez windows) plast.py -i AGCGGGGTAGAGGAATTGGTTTACTCATCAGGCTCATGACCTGAAGACTGCAGGTTCGAATCCTGT CCCCGCCT -db tRNAs.fasta -E 4 -ss 0.001 -seed 11111111111
#Changer la séquence pour chaque unknown 
import argparse
import math
from Bio import SeqIO #(besoin de Biopython package)

# Toutes les valeurs constantes utiles pour les calculs 
MATCH_SCORE = 5
MISMATCH_PENALTY = -4
DEFAULT_E = 4
DEFAULT_SS = 0.001
LAMBDA = 0.192
K = 0.176

#1.1 Banque de données
def parse_fasta(file):
    return {record.id: str(record.seq) for record in SeqIO.parse(file, "fasta")} 


#1.2 Graines (seeds) et Kmer
#Genère des kmers et leurs positions en utilisant la graine standard 11111111111
def generate_kmers(sequence, seed):
    
    k = len(seed)
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        if all(seed[j] == '1' for j in range(k)):  
            kmers.append((kmer, i))
    return kmers

#1.3 Recherche exacte des mots de tailles k dans la banque de donnée tRNAs.fasta 
#Trouve tous les matchs exactes pour les kmers et chaque match représentent les HSP initiaux 
def find_exact_matches(kmers, db_sequences):
    
    matches = []
    for db_id, db_seq in db_sequences.items():
        for kmer, query_pos in kmers:
            start = 0
            while (index := db_seq.find(kmer, start)) != -1:
                matches.append((db_id, kmer, query_pos, index))
                start = index + 1  
                
    return matches 


#1.4 Heuristique gloutonne pour l'extension des HSP
# Étendre dans les 2 directions, arrêter directement lorsqu'on atteint le seuil -E 4 

def extend_hsp(query_seq, db_seq, q_start, d_start, e_thresh):
    
    q_len, d_len = len(query_seq), len(db_seq)
    q_left, d_left = q_start, d_start  # Position de départ des HSP
    q_right, d_right = q_start + 1, d_start + 1  # Match initial

    max_score, current_score = 0, 0
    max_region = (q_start,  d_start)  

    # Étendre vers la gauche 
    while q_left > 0 and d_left > 0:
        q_left -= 1
        d_left -= 1
        new_score = current_score + (MATCH_SCORE if query_seq[q_left] == db_seq[d_left] else MISMATCH_PENALTY) # Les matchs et mismatchs ont les valeurs fixées de l'énoncé 
        if new_score > max_score:
            max_score = new_score
            max_region = (q_left, q_right , d_left, d_right )  
        current_score = new_score
        
        if max_score - current_score >= e_thresh:
            q_left += 1  
            d_left += 1
            break  
   
    current_score = max_score

    # Etendre vers la droite 
    while q_right < q_len and d_right < d_len:
        new_score = current_score + (MATCH_SCORE if query_seq[q_right] == db_seq[d_right] else MISMATCH_PENALTY)
        if new_score > max_score:
            max_score = new_score
            max_region = (q_left, q_right + 1, d_left, d_right + 1)  
        current_score = new_score
        q_right += 1
        d_right += 1
        if max_score - current_score >= e_thresh:
            break  
        
    
    return max_score, max_region


#1.5 Fusion des HSP chevauchants 
def merge_overlapping_hsps(hsp_list):
    if not hsp_list:
        return []
    
    hsp_list.sort(key=lambda x: (x[0], x[2]))   
    merged_hsps = []
    current_hsp = hsp_list[0]
    
    for hsp in hsp_list[1:]:
        q_start, q_end, d_start, d_end, score = current_hsp
        next_q_start, next_q_end, next_d_start, next_d_end, next_score = hsp
        
        if q_end >= next_q_start and d_end >= next_d_start:
            # Fusionner les HSP chevauchant 
            current_hsp = (
                min(q_start, next_q_start), 
                max(q_end, next_q_end), 
                min(d_start, next_d_start), 
                max(d_end, next_d_end), 
                max(score, next_score)  
            )
        else:
            merged_hsps.append(current_hsp)
            current_hsp = hsp
    
    merged_hsps.append(current_hsp)
    return merged_hsps

#1.6 Statistiqu sur les HSPs et cutoff 
#Calculer le bitscore avec le score brut en utilisant la forumle et les valeurs données 
def calculate_bitscore(raw_score):
    return round((LAMBDA * raw_score - math.log(K)) / math.log(2))
    
#Calculer le evalue pour chaque HSP aec un seuil de signification de 1e-3 ou 0.001 
def calculate_evalue(bitscore, total_db_length, query_length):
    return total_db_length * query_length * 2 ** -bitscore

#Fonction main PLAST 
def plast(query_seq, db_sequences, seed, e_thresh, ss_thresh):
    results = []
    kmers = generate_kmers(query_seq, seed)
    matches = find_exact_matches(kmers, db_sequences)
    
    total_db_length = sum(len(seq) for seq in db_sequences.values())
    query_length = len(query_seq)
    
    hsp_by_sequence = {}

    for db_id, _, q_start, d_start in matches:
        #Extension des HSPs
        score, region = extend_hsp(query_seq, db_sequences[db_id], q_start, d_start, e_thresh)
        q_start, q_end, d_start, d_end = region
        hsp = (q_start, q_end, d_start, d_end, score)

        if db_id not in hsp_by_sequence:
            hsp_by_sequence[db_id] = []
        hsp_by_sequence[db_id].append(hsp)

    # Fusionner les hsp et calculer leurs métriques 
    for db_id, hsps in hsp_by_sequence.items():
        merged_hsps = merge_overlapping_hsps(hsps)
        for hsp in merged_hsps:
            q_start, q_end, d_start, d_end, score = hsp
            bitscore = calculate_bitscore(score)
            evalue = calculate_evalue(bitscore, total_db_length, query_length)

            if evalue < ss_thresh:
                results.append((db_id, score, bitscore, evalue, (q_start, q_end, d_start, d_end)))

    # DOnner les résultats en ordre de signification 
    results.sort(key=lambda x: x[2], reverse=True)
    return results


#1.7 Output du programme PLAST en utilisant la structure demandée 
def PLAST_output(results, query_seq, db_sequences):
    
    total = len(results)
    output = []
    for db_id, score, bitscore, evalue, region in results:
        q_start, q_end, d_start, d_end = region
        query_fragment = query_seq[q_start:q_end]
        db_fragment = db_sequences[db_id][d_start:d_end]

        output.append(f">{db_id}")
        output.append(f"# Best HSP score:{score:.2f}, bitscore:{bitscore:.2f}, evalue: {evalue:.2e}")
        output.append(f"{q_start} {query_fragment} {q_end-1}")
        output.append(f"{d_start} {db_fragment} {d_end-1}")
        output.append("-" * 40)
    output.append(f"Total : {total}")
    return "\n".join(output)

# Format des inputs et messages d'erreurs si une des valeurs est manquante 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PLAST: Primitive Local Alignment Search Tool")
    parser.add_argument("-i", "--input", required=True, help="Input la séquence nucléotide")
    parser.add_argument("-db", "--database", required=True, help=" Fichier de base de données FASTA")
    parser.add_argument("-E", type=float, default=DEFAULT_E, help="Seuil pour la descente")
    parser.add_argument("-ss", type=float, default=DEFAULT_SS, help="Seuil de signification")
    parser.add_argument("-seed", required=True, help="Format seed, e.x., '11111111111'")

    args = parser.parse_args()
    query_seq = args.input
    db_sequences = parse_fasta(args.database)
    results = plast(query_seq, db_sequences, args.seed, args.E, args.ss)
    output = PLAST_output(results, query_seq, db_sequences)
    print(output)
